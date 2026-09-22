"""Versioned spinup and configurable process experiments; v1 code is unchanged."""
from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Callable
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

from . import temperature_precipitation as v1
from .tp_processes import (ProcessOptions, visit_temperature, visit_gpp,
                           visit_quadratic_supply, visit_c3_soil_scalar)


@dataclass(frozen=True)
class Model:
    parameters: v1.Parameters = v1.Parameters()
    control: v1.Controls = v1.Controls()
    processes: ProcessOptions = ProcessOptions()

    def __post_init__(self):
        if self.processes.photo_temperature == "visit_tem" and not (
                self.processes.photo_tmin < self.parameters.t_opt < self.processes.photo_tmax):
            raise ValueError("VISIT temperature option needs Tmin < Topt < Tmax")

    def fluxes(self, state, climate):
        p, c, o = self.parameters, self.control, self.processes
        f = v1.fluxes(state, climate, p, c)
        x = np.maximum(np.asarray(state), 0.)  # v1 trial-state extension only
        r = v1.relative_water(x, p)
        t = climate.temperature if c.photo_temperature is None else c.photo_temperature
        ft = np.exp(-.5*((t-p.t_opt)/p.t_width)**2)
        if o.photo_temperature == "visit_tem":
            ft = visit_temperature(t, o.photo_tmin, p.t_opt, o.photo_tmax)
        if not c.temperature_response:
            ft = 1.
        rg = np.sqrt(r[0]*r[2])
        fw = rg**p.photo_power/(rg**p.photo_power+p.photo_half**p.photo_power)
        if o.photo_water == "visit_soil_with_plant_gate":
            fw *= visit_c3_soil_scalar(x[7], p.soil_fc, o.soil_half)
        if not c.water_response:
            fw = 1.
        g = p.gpp_max*ft*fw*(-np.expm1(-p.k_leaf*x[0]))
        if o.photo_canopy == "visit_monsi":
            ref = visit_gpp(1., 1., 1., 1., o.canopy_light_load, 1000.)
            g = p.gpp_max*visit_gpp(ft*fw, 1., 1., 1., o.canopy_light_load, p.k_leaf*x[0])/ref
        el, eo = f.transpiration, f.evaporation
        if o.evap_supply == "visit_quadratic":
            te = climate.temperature if c.evap_temperature is None else c.evap_temperature
            qe = p.q10_evap**((te-p.t_ref)/10) if c.temperature_response else 1.
            potential_leaf = p.transpiration*qe*(-np.expm1(-p.k_area*x[0]))
            potential_soil = p.evaporation*qe
            el = visit_quadratic_supply(x[4]/o.supply_time_days, potential_leaf, o.supply_curvature)
            eo = visit_quadratic_supply(max(0., x[7]-p.soil_dry)/o.supply_time_days, potential_soil, o.supply_curvature)
        drainage = f.drainage
        if o.drainage == "visit_baseflow_plus_excess":
            drainage += o.baseflow_rate*x[7]
        return replace(f, gpp=float(g), carbon_input=np.r_[np.asarray(p.allocation)*g, 0.],
                       transpiration=float(el), evaporation=float(eo), drainage=float(drainage))

    def rhs(self, time, state, climate):
        f = self.fluxes(state, climate)
        dw = v1.WATER_INCIDENCE@f.transport
        dw[:3] -= f.water_release
        dw[0] -= f.transpiration
        dw[3] += climate.precipitation+np.asarray(self.parameters.zeta)@f.water_release-f.evaporation-f.drainage-f.runoff
        return np.r_[f.carbon_input+f.carbon_matrix@np.maximum(np.asarray(state)[:4], 0.), dw]

    def budgets(self, state, climate):
        f = self.fluxes(state, climate)
        dx = self.rhs(0., state, climate)
        return np.array([sum(dx[:4])-f.nep+f.carbon_export,
            sum(dx[4:])-climate.precipitation+f.transpiration+f.evaporation+f.drainage+f.runoff+f.water_export])

    def jacobian(self, state, climate, step=1e-5):
        x = np.asarray(state)
        j = np.empty((8, 8))
        for k in range(8):
            d = np.zeros(8)
            d[k] = step*max(1., abs(x[k]))
            if x[k] >= d[k]:
                j[:, k] = (self.rhs(0, x+d, climate)-self.rhs(0, x-d, climate))/(2*d[k])
            else:
                j[:, k] = (self.rhs(0, x+d, climate)-self.rhs(0, x, climate))/d[k]
        return j

    def equilibrium(self, climate, guess):
        fit = least_squares(lambda y: self.rhs(0, y*v1.STATE_SCALE, climate)/v1.STATE_SCALE,
                            np.asarray(guess)/v1.STATE_SCALE, bounds=(0, np.inf),
                            gtol=None, ftol=1e-13, xtol=1e-13, max_nfev=2000)
        x = fit.x*v1.STATE_SCALE
        residual = np.max(np.abs(self.rhs(0, x, climate)))
        if not fit.success or residual > 1e-8 or v1.domain_violation(x, self.parameters) > 1e-8:
            raise RuntimeError("unverified coupled equilibrium")
        j = self.jacobian(x, climate)
        return v1.Equilibrium(x, float(residual), np.linalg.eigvals(j), float(np.linalg.cond(j)), fit.nfev)


@dataclass(frozen=True)
class Solver:
    method: str = "Radau"
    max_step: float = 10.
    rtol: float = 1e-8
    atol: float = 1e-10
    bound_tolerance: float = 1e-7
    retries: int = 3

    def __post_init__(self):
        if self.method not in ("Radau", "BDF", "DOP853"):
            raise ValueError("unsupported solver")
        if not all(np.isfinite(v) and v > 0 for v in (self.max_step, self.rtol, self.atol, self.bound_tolerance)):
            raise ValueError("positive finite tolerances and max_step required")
        if self.max_step >= 28:
            raise ValueError("internal timestep must be shorter than a month")
        if not isinstance(self.retries, int) or self.retries < 0:
            raise ValueError("invalid retries")


@dataclass(frozen=True)
class Trajectory:
    times: np.ndarray
    states: np.ndarray
    audit: dict


def integrate(model, initial, times, climate: v1.Climate | Callable, solver=Solver()):
    """Reject negative accepted states; retry smaller steps, never add mass.

    Audits every accepted step, requested output, and three interior dense
    samples per step. These sampled checks are not an interval proof. The
    differential equations' invariant-domain tests supply the structural check.
    Implicit Newton trial states may be negative; they are not accepted states.
    """
    initial, times = np.asarray(initial, dtype=float), np.asarray(times, dtype=float)
    if initial.shape != (8,) or not np.all(np.isfinite(initial)) or initial.min() < 0:
        raise ValueError("initial pools must be eight finite nonnegative values")
    if v1.domain_violation(initial, model.parameters) > solver.bound_tolerance:
        raise ValueError("initial water exceeds capacity")
    if times.ndim != 1 or len(times) < 2 or not np.all(np.isfinite(times)) or np.any(np.diff(times) <= 0):
        raise ValueError("strictly increasing finite times required")
    forcing = climate if callable(climate) else lambda t: climate
    attempts = []
    for attempt in range(solver.retries+1):
        max_step = solver.max_step/2**attempt
        sol = solve_ivp(lambda t, x: model.rhs(t, x, forcing(t)), (times[0], times[-1]), initial,
                        method=solver.method, dense_output=True, max_step=max_step,
                        rtol=solver.rtol, atol=solver.atol)
        if not sol.success:
            attempts.append({"reason": sol.message, "max_step": max_step})
            continue
        inner = np.concatenate([sol.t[:-1]+q*np.diff(sol.t) for q in (.25, .5, .75)])
        audit_times = np.unique(np.r_[sol.t, times, inner])
        sampled = sol.sol(audit_times).T
        minima = sampled.min(axis=0)
        cap = sampled[:, :3]*np.asarray(model.parameters.omega)+model.parameters.epsilon
        excess = max(0., float((sampled[:, 4:7]-cap).max()), float(sampled[:, 7].max()-model.parameters.soil_sat))
        if not np.all(np.isfinite(sampled)) or minima.min() < 0 or excess > solver.bound_tolerance:
            attempts.append({"reason": "negative pool or capacity violation", "minima": minima.tolist(),
                             "max_capacity_excess": excess, "max_step": max_step})
            continue
        return Trajectory(times, sol.sol(times).T, {"minima": dict(zip(v1.STATE_NAMES, minima.tolist())),
            "accepted_steps": len(sol.t)-1, "inspected_points": len(audit_times),
            "max_capacity_excess": excess, "rejected_attempts": attempts,
            "max_step": max_step, "mass_added_by_clipping": 0.,
            "scope": "all accepted steps, requested outputs and quarter-step dense samples; not implicit Newton trial states"})
    raise RuntimeError(f"integration failed physical-domain guard: {attempts}")


@dataclass(frozen=True)
class SpinupSettings:
    max_years: int = 800
    chunk_years: int = 25
    consecutive_years: int = 3
    annual_scaled_change: float = 1e-7
    rhs_absolute: float = 1e-8

    def __post_init__(self):
        for v in (self.max_years, self.chunk_years, self.consecutive_years):
            if not isinstance(v, int) or v < 1:
                raise ValueError("spinup lengths must be positive integer years")
        for v in (self.annual_scaled_change, self.rhs_absolute):
            if not np.isfinite(v) or v <= 0:
                raise ValueError("invalid spinup tolerance")


def spinup(model, initial, climate, settings=SpinupSettings(), solver=Solver()):
    """Actual integration from a user seed; never substitute a root solution.

    Annual checks use 365.25-day diagnostic years under constant climate.
    Stop after a chunk whose last consecutive years pass both criteria.
    """
    states, times, checks, audits = [np.asarray(initial)], [0.], [], []
    streak = 0
    for start in range(0, settings.max_years, settings.chunk_years):
        end = min(start+settings.chunk_years, settings.max_years)
        t = np.arange(start, end+1)*365.25
        result = integrate(model, states[-1], t, climate, solver)
        audits.append(result.audit)
        for k in range(1, len(t)):
            x, prev = result.states[k], result.states[k-1]
            drift = np.max(np.abs(x-prev)/np.maximum(1., abs(x)))
            residual = np.max(abs(model.rhs(t[k], x, climate)))
            passed = drift < settings.annual_scaled_change and residual < settings.rhs_absolute
            streak = streak+1 if passed else 0
            checks.append({"year": t[k]/365.25, "annual_scaled_change": float(drift),
                           "rhs_absolute": float(residual), "passes": bool(passed)})
        times.extend(t[1:])
        states.extend(result.states[1:])
        if streak >= settings.consecutive_years:
            return {"converged": True, "times": np.array(times), "states": np.array(states),
                    "checks": checks, "audits": audits, "years": end}
    return {"converged": False, "times": np.array(times), "states": np.array(states),
            "checks": checks, "audits": audits, "years": settings.max_years}
