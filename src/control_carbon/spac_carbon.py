"""One-way Manzoni + four carbon pools, independent synthetic theory.

State: [psi_x, psi_leaf, leaf, wood, root, soil]; MPa and kg C/m2 ground.
Time: day; D: mol/mol. Assimilation is GPP; live losses are explicitly
partitioned into litter transfer and respiration. Not VISIT or Leuning.
"""

from dataclasses import dataclass
from typing import Callable

import numpy as np
from scipy.optimize import brentq
from scipy.integrate import solve_ivp

from .manzoni_hydraulics import HydraulicParameters, conductances, hydraulic_rhs, hydraulic_jacobian
from .manzoni_extensions import (
    SigmoidShape, sigmoid_conductances, sigmoid_rhs, sigmoid_jacobian,
)

YEAR_DAYS = 365.25


@dataclass(frozen=True)
class CarbonParameters:
    assimilation_max: float
    allocation: tuple[float, float, float]
    losses: tuple[float, float, float, float]
    litter_fractions: tuple[float, float, float]
    specific_leaf_area: float
    carbon_dry_fraction: float
    extinction: float

    def __post_init__(self):
        for values, size in [
            (self.allocation, 3), (self.losses, 4), (self.litter_fractions, 3),
        ]:
            array = np.asarray(values, dtype=float)
            if array.shape != (size,) or not np.all(np.isfinite(array)):
                raise ValueError("carbon coefficient arrays must have correct shape and finite values")
        scalars = [self.assimilation_max, self.specific_leaf_area,
                   self.carbon_dry_fraction, self.extinction]
        if not np.all(np.isfinite(scalars)) or min(scalars) <= 0:
            raise ValueError("assimilation, SLA, carbon fraction and extinction must be finite positive")
        if self.carbon_dry_fraction > 1:
            raise ValueError("carbon dry-mass fraction must not exceed one")
        if min(self.allocation) <= 0 or not np.isclose(sum(self.allocation), 1., atol=1e-12, rtol=0):
            raise ValueError("positive live allocation fractions must sum to one")
        if min(self.losses) <= 0 or min(self.litter_fractions) < 0 or max(self.litter_fractions) > 1:
            raise ValueError("loss rates must be positive and litter fractions in [0,1]")

    @property
    def canopy_coefficient(self):
        """q=extinction*SLA/fC in (kg C/m2)^-1; SLA in m2/kg dry."""
        return self.extinction*self.specific_leaf_area/self.carbon_dry_fraction


@dataclass(frozen=True)
class TemperatureResponse:
    reference: float
    optimum: float
    width: float
    loss_q10: float

    def __post_init__(self):
        if not np.all(np.isfinite([self.reference, self.optimum, self.width, self.loss_q10])):
            raise ValueError("temperature response parameters must be finite")
        if self.width <= 0 or self.loss_q10 <= 0:
            raise ValueError("temperature width and Q10 must be positive")

    def factors(self, temperature: float) -> tuple[float, float]:
        """Synthetic Gaussian GPP and Q10 losses, normalized at reference."""
        if not np.isfinite(temperature):
            raise ValueError("temperature must be finite (Celsius)")
        with np.errstate(over="raise", under="raise", invalid="raise"):
            assimilation = np.exp(
                ((self.reference-self.optimum)**2-(temperature-self.optimum)**2)/self.width**2,
            )
            loss = self.loss_q10**((temperature-self.reference)/10.)
        if not np.all(np.isfinite([assimilation, loss])) or min(assimilation, loss) <= 0:
            raise ValueError("temperature coefficients are outside representable positive range")
        return float(assimilation), float(loss)


def carbon_matrix(p: CarbonParameters, loss_factor: float = 1.) -> np.ndarray:
    if not np.isfinite(loss_factor) or loss_factor <= 0:
        raise ValueError("loss factor must be finite positive")
    k = np.asarray(p.losses)*loss_factor
    matrix = -np.diag(k)
    matrix[3, :3] = k[:3]*p.litter_fractions
    return matrix


def _carbon(state):
    state = np.asarray(state, dtype=float)
    if state.shape != (4,) or not np.all(np.isfinite(state)) or np.any(state < 0):
        raise ValueError("carbon state must contain four finite nonnegative stocks")
    return state


def _environment(h, assimilation_factor):
    if not np.all(np.isfinite([h, assimilation_factor])) or not 0 <= h <= 1 or assimilation_factor <= 0:
        raise ValueError("hydraulic factor must lie in [0,1] and assimilation factor be positive")


def assimilation(leaf, h, p: CarbonParameters, assimilation_factor=1.):
    _environment(h, assimilation_factor)
    if not np.isfinite(leaf) or leaf < 0:
        raise ValueError("leaf carbon must be finite nonnegative")
    return p.assimilation_max*assimilation_factor*h*(-np.expm1(-p.canopy_coefficient*leaf))


def carbon_rhs(state, h, p: CarbonParameters, assimilation_factor=1., loss_factor=1.):
    state = _carbon(state)
    gpp = assimilation(state[0], h, p, assimilation_factor)
    return carbon_matrix(p, loss_factor)@state + np.r_[p.allocation, 0.]*gpp


def carbon_budget(state, h, p: CarbonParameters, assimilation_factor=1., loss_factor=1.):
    state = _carbon(state)
    carbon_matrix(p, loss_factor)
    gpp = assimilation(state[0], h, p, assimilation_factor)
    live_respiration = np.dot((1-np.array(p.litter_fractions))*np.array(p.losses[:3])*loss_factor, state[:3])
    heterotrophic = p.losses[3]*loss_factor*state[3]
    npp = gpp-live_respiration
    nep = npp-heterotrophic
    return {
        "gpp": float(gpp), "autotrophic_respiration": float(live_respiration),
        "npp": float(npp), "heterotrophic_respiration": float(heterotrophic),
        "nep": float(nep),
        "residual": float(np.sum(carbon_rhs(state, h, p, assimilation_factor, loss_factor))-nep),
    }


def frozen_coefficient_capacity(state, h, p: CarbonParameters, assimilation_factor=1., loss_factor=1.):
    """Freeze current nonlinear input: a map, not generally an equilibrium."""
    state = _carbon(state)
    u = np.r_[p.allocation, 0.]*assimilation(state[0], h, p, assimilation_factor)
    return np.linalg.solve(-carbon_matrix(p, loss_factor), u)


def reproduction_number(h, p: CarbonParameters, assimilation_factor=1., loss_factor=1.):
    _environment(h, assimilation_factor)
    carbon_matrix(p, loss_factor)
    return p.allocation[0]*p.assimilation_max*assimilation_factor*h*p.canopy_coefficient/(p.losses[0]*loss_factor)


def positive_carbon_equilibrium(h, p: CarbonParameters, assimilation_factor=1., loss_factor=1.):
    """Unique positive equilibrium if R0>1; otherwise no positive root."""
    r0 = reproduction_number(h, p, assimilation_factor, loss_factor)
    if r0 <= 1:
        return None
    # z=q*L; divide by z so the zero root cannot be mistaken for the positive one.
    def ratio(z):
        return r0*(-np.expm1(-z))/z-1 if z > 0 else r0-1
    z = brentq(ratio, 0., r0, xtol=1e-14, rtol=1e-14)
    state = np.zeros(4)
    state[0] = z/p.canopy_coefficient
    gpp = assimilation(state[0], h, p, assimilation_factor)
    state[:3] = np.array(p.allocation)*gpp/(np.array(p.losses[:3])*loss_factor)
    state[3] = np.dot(np.array(p.litter_fractions)*np.array(p.losses[:3]), state[:3])/p.losses[3]
    return state


def carbon_jacobian(state, h, p: CarbonParameters, assimilation_factor=1., loss_factor=1.):
    state = _carbon(state)
    _environment(h, assimilation_factor)
    matrix = carbon_matrix(p, loss_factor)
    slope = p.assimilation_max*assimilation_factor*h*p.canopy_coefficient*np.exp(-p.canopy_coefficient*state[0])
    matrix[:, 0] += np.r_[p.allocation, 0.]*slope
    return matrix


def hydraulic_factor(leaf_potential, hydraulic: HydraulicParameters, shape: SigmoidShape | None = None):
    gx, gs = (conductances(leaf_potential, hydraulic) if shape is None
              else sigmoid_conductances(leaf_potential, hydraulic, shape))
    return float(gs/hydraulic.gs_max)


def _combined_state(state):
    state = np.asarray(state, dtype=float)
    if state.shape != (6,):
        raise ValueError("combined state must have six entries")
    _carbon(state[2:])
    return state


def combined_rhs(state, root, vpd, hydraulic, carbon, assimilation_factor=1., loss_factor=1., shape=None):
    state = _combined_state(state)
    water = (hydraulic_rhs(state[:2], root, vpd, hydraulic) if shape is None
             else sigmoid_rhs(state[:2], root, vpd, hydraulic, shape))
    h = hydraulic_factor(state[1], hydraulic, shape)
    return np.r_[water, carbon_rhs(state[2:], h, carbon, assimilation_factor, loss_factor)]


def combined_jacobian(state, root, vpd, hydraulic, carbon, assimilation_factor=1., loss_factor=1., shape=None):
    state = _combined_state(state)
    water = (hydraulic_jacobian(state[:2], root, vpd, hydraulic) if shape is None
             else sigmoid_jacobian(state[:2], root, vpd, hydraulic, shape))
    h = hydraulic_factor(state[1], hydraulic, shape)
    if shape is None:
        dh = -1/(2*hydraulic.psi50_s) if h > 0 else 0.
    else:
        if state[1] == 0:
            if shape.alpha_s < 1:
                raise ValueError("hydraulic stress derivative is singular at zero potential")
            dh = 1/(-hydraulic.psi50_s) if shape.alpha_s == 1 else 0.
        else:
            dh = shape.alpha_s*h*(1-h)/(-state[1])
    jac = np.zeros((6, 6))
    jac[:2, :2] = water
    jac[2:, 2:] = carbon_jacobian(state[2:], h, carbon, assimilation_factor, loss_factor)
    jac[2:, 1] = np.r_[carbon.allocation, 0.]*carbon.assimilation_max*assimilation_factor*dh*(-np.expm1(-carbon.canopy_coefficient*state[2]))
    return jac


def simulate(
    initial, times, forcing: Callable[[float], tuple[float, float, float, float]],
    hydraulic: HydraulicParameters, carbon: CarbonParameters,
    *, shape: SigmoidShape | None = None, rtol=1e-9,
):
    """Integrate positive carbon in log coordinates, with no state clipping.

    forcing(day) -> (root MPa, D mol/mol, assimilation multiplier, loss multiplier).
    This integrator requires strictly positive initial carbon; zero equilibria
    are analyzed directly, not initialized with an artificial numerical floor.
    """
    initial = _combined_state(initial)
    times = np.asarray(times, dtype=float)
    if times.ndim != 1 or len(times) < 2 or not np.all(np.isfinite(times)) or np.any(np.diff(times) <= 0):
        raise ValueError("sample times must be finite strictly increasing")
    if not np.isfinite(rtol) or not 0 < rtol < 1 or np.any(initial[2:] <= 0):
        raise ValueError("rtol must lie in (0,1) and initial carbon be strictly positive")
    def rhs(t, z):
        stocks = np.exp(z[2:])
        if not np.all(np.isfinite(stocks)) or np.any(stocks == 0):
            raise FloatingPointError("log carbon trajectory exceeds representable range")
        root, vpd, assimilation_factor, loss_factor = forcing(t)
        f = combined_rhs(np.r_[z[:2], stocks], root, vpd, hydraulic, carbon,
                         assimilation_factor, loss_factor, shape=shape)
        return np.r_[f[:2], f[2:]/stocks]
    result = solve_ivp(
        rhs, (times[0], times[-1]), np.r_[initial[:2], np.log(initial[2:])],
        method="Radau", rtol=rtol, atol=rtol*.01, dense_output=True,
    )
    if not result.success:
        raise RuntimeError(result.message)
    sampled = result.sol(times)
    states = np.vstack([sampled[:2], np.exp(sampled[2:])]).T
    if np.any(result.y[:2] > 0) or np.any(states[:, :2] > 0) or not np.all(np.isfinite(states)) or np.any(states[:, 2:] <= 0):
        raise RuntimeError("accepted or sampled trajectory left the physical/representable domain")
    return {"time_days": times, "state": states, "nfev": result.nfev}
