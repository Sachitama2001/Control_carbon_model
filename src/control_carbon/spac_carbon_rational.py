"""Independent rational-canopy E0--E4 variant; existing exponential ODEs unchanged."""

from dataclasses import dataclass

import numpy as np
from numpy.polynomial import Polynomial
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from .manzoni_hydraulics import _state as hydraulic_state, conductances, hydraulic_jacobian
from .spac_carbon import _combined_state, carbon_matrix
from .spac_carbon_equilibria import stock_ratios
from .spac_carbon_mpa import PressureSPAC
from .spac_carbon_stability import _Interval


@dataclass(frozen=True)
class RationalSPAC:
    base: PressureSPAC

    def __post_init__(self):
        h, f, c = self.base.hydraulic, self.base.feedback, self.base.carbon
        strengths = (f.leaf_supply, f.leaf_demand, f.wood_supply)
        if strengths not in [(0., 0., 0.), (1., 0., 0.), (0., 1., 0.),
                             (1., 1., 0.), (1., 1., 1.)]:
            raise ValueError("Rational analysis requires a full E0--E4 feedback variant")
        if (self.base.root_potential_mpa != 0 or h.psi50_x != -1.5 or h.psi50_s != -2
                or f.wood_exponent != 1 or c.canopy_coefficient != 5):
            raise ValueError("Algebraic classification requires root0, thresholds -3/-4, q5, beta1")
        if not .25 < self.opening(0.) < .5:
            raise ValueError("Algebraic classification requires 1/4 < zero-leaf opening < 1/2")

    @property
    def label(self):
        f = self.base.feedback
        return {
            (0., 0., 0.): "E0", (1., 0., 0.): "E1", (0., 1., 0.): "E2",
            (1., 1., 0.): "E3", (1., 1., 1.): "E4",
        }[(f.leaf_supply, f.leaf_demand, f.wood_supply)]

    @property
    def wet_leaf(self):
        c = self.base.carbon
        return c.allocation[0]*c.assimilation_max/c.losses[0] - 1/c.canopy_coefficient

    def opening(self, leaf):
        if not np.isfinite(leaf) or leaf < 0:
            raise ValueError("Leaf carbon must be finite nonnegative")
        c = self.base.carbon
        return c.losses[0]*(1+c.canopy_coefficient*leaf)/(
            c.allocation[0]*c.assimilation_max*c.canopy_coefficient
        )

    @staticmethod
    def _state(state):
        state = _combined_state(state)
        hydraulic_state(state[:2])
        return state

    def assimilation(self, state):
        state = self._state(state)
        _, gs = conductances(state[1], self.base.hydraulic)
        c = self.base.carbon
        return c.assimilation_max*(gs/self.base.hydraulic.gs_max)*(
            c.canopy_coefficient*state[2]/(1+c.canopy_coefficient*state[2])
        )

    def factors(self, state):
        f = self.base.feedback
        return (
            (state[2]+f.leaf_floor)/(f.leaf_reference+f.leaf_floor) if f.leaf_supply else 1.,
            state[2]/f.leaf_reference if f.leaf_demand else 1.,
            state[3]/f.wood_reference if f.wood_supply else 1.,
        )

    def fluxes(self, state, vpd_mpa):
        state = self._state(state)
        gx, gs = conductances(state[1], self.base.hydraulic)
        fl, fe, fw = self.factors(state)
        return np.array([2*gx*fw*(self.base.root_potential_mpa-state[0]),
                         2*gx*fl*(state[0]-state[1]), gs*fe*self.base.mole_fraction(vpd_mpa)])

    def capacity(self, state):
        c = self.base.carbon
        return np.linalg.solve(-carbon_matrix(c), np.r_[c.allocation, 0.]*self.assimilation(state))

    def rhs(self, state, vpd_mpa):
        erx, exl, evaporation = self.fluxes(state, vpd_mpa)
        h, c = self.base.hydraulic, self.base.carbon
        return np.r_[(erx-exl)/h.cx, (exl-evaporation)/h.cl,
                     carbon_matrix(c)@np.asarray(state)[2:]
                     +np.r_[c.allocation, 0.]*self.assimilation(state)]

    def jacobian(self, state, vpd_mpa):
        state = self._state(state)
        h, c, f = self.base.hydraulic, self.base.carbon, self.base.feedback
        d = self.base.mole_fraction(vpd_mpa)
        hydraulic_jacobian(state[:2], self.base.root_potential_mpa, d, h)
        jac = np.zeros((6, 6))
        gx, gs = conductances(state[1], h)
        dgx = -h.gx_max/(2*h.psi50_x) if gx > 0 else 0.
        dgs = -h.gs_max/(2*h.psi50_s) if gs > 0 else 0.
        fl, fe, fw = self.factors(state)
        dfl = f.leaf_supply/(f.leaf_reference+f.leaf_floor)
        dfe, dfw = f.leaf_demand/f.leaf_reference, f.wood_supply/f.wood_reference
        x, leaf = state[:2]
        jac[0, 0] = -2*gx*(fw+fl)/h.cx
        jac[0, 1] = (2*dgx*(fw*(self.base.root_potential_mpa-x)-fl*(x-leaf))+2*gx*fl)/h.cx
        jac[0, 2] = -2*gx*dfl*(x-leaf)/h.cx
        jac[0, 3] = 2*gx*dfw*(self.base.root_potential_mpa-x)/h.cx
        jac[1, 0] = 2*gx*fl/h.cl
        jac[1, 1] = (2*dgx*fl*(x-leaf)-2*gx*fl-dgs*fe*d)/h.cl
        jac[1, 2] = (2*gx*dfl*(x-leaf)-gs*dfe*d)/h.cl
        ql = c.canopy_coefficient*state[2]
        jac[2:, 2:] = carbon_matrix(c)
        jac[2:, 1] = np.r_[c.allocation, 0.]*c.assimilation_max*dgs/h.gs_max*ql/(1+ql)
        jac[2:, 2] += (np.r_[c.allocation, 0.]*c.assimilation_max*gs/h.gs_max
                       *c.canopy_coefficient/(1+ql)**2)
        return jac

    def curve_polynomials(self):
        """Return N(L), Q(L) with D=N/Q; no denominator roots are admissible."""
        c, h, f = self.base.carbon, self.base.hydraulic, self.base.feedback
        leaf = Polynomial([0., 1.])
        opening = Polynomial([self.opening(0.),
                              c.losses[0]/(c.allocation[0]*c.assimilation_max)])
        numerator = 4*h.gx_max/(3*h.gs_max)*(1-opening)*(4*opening-1)
        denominator = opening
        if self.label == "E1":
            numerator *= 2*(leaf+f.leaf_floor)
            denominator *= leaf+f.leaf_reference+2*f.leaf_floor
        elif self.label == "E2":
            numerator *= f.leaf_reference
            denominator *= leaf
        elif self.label == "E3":
            numerator *= 2*f.leaf_reference*(leaf+f.leaf_floor)
            denominator *= leaf*(leaf+f.leaf_reference+2*f.leaf_floor)
        elif self.label == "E4":
            ratio = stock_ratios(self.base)[1]
            numerator *= 2*ratio*f.leaf_reference*(leaf+f.leaf_floor)
            denominator *= (ratio*leaf*(f.leaf_reference+f.leaf_floor)
                            +f.wood_reference*(leaf+f.leaf_floor))
        return numerator, denominator

    def equilibrium_polynomial(self, vpd_mpa):
        numerator, denominator = self.curve_polynomials()
        return self.base.mole_fraction(vpd_mpa)*denominator-numerator

    def curve(self, leaf):
        if not np.isfinite(leaf) or not 0 < leaf <= self.wet_leaf:
            raise ValueError("Equilibrium curve requires 0 < C_L <= wet endpoint")
        psi = -4*(1-self.opening(leaf)) if leaf != self.wet_leaf else 0.
        stocks = stock_ratios(self.base)*leaf
        fl, _, fw = self.factors(np.r_[0., psi, stocks])
        state = np.r_[fl*psi/(fw+fl), psi, stocks]
        numerator, denominator = self.curve_polynomials()
        value = numerator(leaf)/denominator(leaf) if leaf != self.wet_leaf else 0.
        derivative = ((numerator.deriv()(leaf)*denominator(leaf)
                       -numerator(leaf)*denominator.deriv()(leaf))/denominator(leaf)**2)
        return state, self.base.pressure(value), derivative*self.base.atmospheric_pressure_mpa

    def shape_bounds(self):
        """Outward constant margins completing the baseline E2/E3/E4 proofs."""
        c, f = self.base.carbon, self.base.feedback
        zero = (_Interval(c.losses[0])
                /(_Interval(c.allocation[0])*_Interval(c.assimilation_max)
                  *_Interval(c.canopy_coefficient)))
        slope = _Interval(c.losses[0])/(_Interval(c.allocation[0])*_Interval(c.assimilation_max))
        cut = _Interval(".03")
        upper = slope*(1/(zero*zero)-4)/(5-4*zero-1/zero)
        ratio = (_Interval(c.allocation[1])*_Interval(c.losses[0])
                 /(_Interval(c.allocation[0])*_Interval(c.losses[1])))
        lf, lr, wr = _Interval(f.leaf_floor), _Interval(f.leaf_reference), _Interval(f.wood_reference)
        lowers = {
            "E2": 1/cut, "E3": lf/(cut*(cut+lf)),
            "E4": ratio*(lr+lf)*lf/(
                (cut+lf)*(ratio*cut*(lr+lf)+wr*(cut+lf))),
        }
        if (zero+slope*cut).lo <= .5 or any((lower-upper).lo <= 0 for lower in lowers.values()):
            raise ValueError("Declared global monotonicity margins did not certify")
        return {"cut_leaf": ".03", "positive_log_base_upper": upper.strings(),
                "negative_log_factor_lower": {label: bound.strings() for label, bound in lowers.items()},
                "strict_gaps": {label: (bound-upper).strings() for label, bound in lowers.items()}}

    def landmarks(self):
        self.shape_bounds()
        numerator, denominator = self.curve_polynomials()
        limit = (self.base.pressure(numerator(0)/denominator(0))
                 if denominator(0) else None)
        fold = None
        if self.label in ["E0", "E1"]:
            c = self.base.carbon
            half = (.5-self.opening(0.))*c.allocation[0]*c.assimilation_max/c.losses[0]
            fold = half if self.label == "E0" else brentq(
                lambda leaf: self.curve(leaf)[2], half, self.wet_leaf, xtol=1e-14)
        return {"zero_leaf_vpd_limit_mpa": limit, "wet_leaf": self.wet_leaf,
                "fold_leaf": fold, "fold_vpd_mpa": None if fold is None else self.curve(fold)[1]}

    def equilibria(self, vpd_mpa):
        self.base.mole_fraction(vpd_mpa)
        points = self.landmarks()
        fold, limit = points["fold_leaf"], points["zero_leaf_vpd_limit_mpa"]
        if vpd_mpa == 0:
            return [self.curve(self.wet_leaf)[0]]
        if fold is not None and vpd_mpa == points["fold_vpd_mpa"]:
            return [self.curve(fold)[0]]
        if limit is None:
            low = self.wet_leaf/2
            while self.curve(low)[1] < vpd_mpa:
                low /= 2
                if low < np.finfo(float).tiny:
                    raise FloatingPointError("Positive equilibrium below representable range")
            nodes = [low, self.wet_leaf]
        else:
            nodes = [0., self.wet_leaf] if fold is None else [0., fold, self.wet_leaf]
        def value(leaf):
            return (limit if leaf == 0 else self.curve(leaf)[1])-vpd_mpa
        roots = []
        for left, right in zip(nodes[:-1], nodes[1:]):
            if value(left)*value(right) < 0:
                leaf = brentq(value, left, right, xtol=np.finfo(float).tiny)
                roots.append(self.curve(leaf)[0])
        return sorted(roots, key=lambda state: -state[2])

    def simulate(self, initial, times, vpd_mpa, *, rtol=1e-9):
        initial, times = np.asarray(initial), np.asarray(times)
        self.rhs(initial, vpd_mpa)
        if (np.any(initial[2:] <= 0) or initial[2] <= self.base.feedback.minimum_leaf
                or times.ndim != 1 or len(times) < 2 or np.any(~np.isfinite(times))
                or np.any(np.diff(times) <= 0) or not np.isfinite(rtol) or not 0 < rtol < 1):
            raise ValueError("Simulation needs positive stocks, increasing times, and 0 < rtol < 1")
        def rhs(t, z):
            stocks = np.exp(z[2:])
            if not np.all(np.isfinite(stocks)) or np.any(stocks == 0):
                raise FloatingPointError("Log-carbon state exceeds representable range")
            values = self.rhs(np.r_[z[:2], stocks], vpd_mpa)
            return np.r_[values[:2], values[2:]/stocks]
        def boundary(t, z):
            return z[2]-np.log(self.base.feedback.minimum_leaf)
        boundary.terminal, boundary.direction = True, -1
        run = solve_ivp(rhs, (times[0], times[-1]), np.r_[initial[:2], np.log(initial[2:])],
                        method="Radau", rtol=rtol, atol=rtol*.01, events=boundary, dense_output=True)
        if not run.success:
            raise RuntimeError(run.message)
        sampled = np.unique(np.r_[times[times <= run.t[-1]], run.t[-1]])
        z = run.sol(sampled)
        states = np.vstack([z[:2], np.exp(z[2:])]).T
        if (np.any(~np.isfinite(states)) or np.any(states[:, :2] > 0)
                or np.any(states[:, 2:] <= 0)
                or np.any(states[:, 2] < self.base.feedback.minimum_leaf*(1-1e-8))):
            raise RuntimeError("Accepted trajectory left the declared physical/domain bounds")
        return {"status": "leaf_domain_exit" if run.status else "completed",
                "day": sampled, "state": states}
