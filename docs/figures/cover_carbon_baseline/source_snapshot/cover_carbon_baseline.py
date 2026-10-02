"""Area-exclusive forest--grass cover coupled to ten carbon pools.

This is the analytical baseline specified in
``docs/forest_grass_cover_climate_research_plan.md``. The cover vector field
is reused from the published-model analysis; carbon parameters are explicit
and are not calibrated defaults.
"""

from dataclasses import dataclass

import numpy as np

from .forest_grass_tipping import (
    Parameters as CoverParameters,
    drift as cover_drift,
    drift_prime as cover_drift_prime,
    interior_roots,
    phi,
    phi_prime,
)


GROUPS = ("forest", "grass")
LIVING_POOLS = ("leaf", "stem", "root")
POOL_NAMES = ("leaf", "stem", "root", "dead", "humus")
STATE_SIZE = 11


def _validated_array(name, values, shape):
    array = np.asarray(values, dtype=float)
    if array.shape != shape or not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be a finite array with shape {shape}.")
    array = array.copy()
    array.setflags(write=False)
    return array


@dataclass(frozen=True)
class CarbonParameters:
    """Fixed carbon coefficients; rows of 2-D arrays are forest then grass."""

    production: np.ndarray
    allocation: np.ndarray
    turnover: np.ndarray
    transfer_emission: np.ndarray
    decomposition: np.ndarray
    humification: np.ndarray
    humus_loss: np.ndarray

    def __post_init__(self):
        shapes = {
            "production": (2,),
            "allocation": (2, 3),
            "turnover": (2, 3),
            "transfer_emission": (2, 3),
            "decomposition": (2,),
            "humification": (2,),
            "humus_loss": (2,),
        }
        for name, shape in shapes.items():
            object.__setattr__(
                self, name, _validated_array(name, getattr(self, name), shape)
            )

        if np.any(self.production <= 0):
            raise ValueError("Production rates must be positive.")
        if np.any(self.allocation <= 0) or not np.allclose(
            self.allocation.sum(axis=1), 1.0, rtol=0.0, atol=1e-12
        ):
            raise ValueError("Each group must have positive allocations summing to one.")
        if np.any(self.turnover <= 0):
            raise ValueError("Living-pool turnover rates must be positive.")
        if np.any((self.transfer_emission < 0) | (self.transfer_emission > 1)):
            raise ValueError("Transfer emission fractions must lie in [0, 1].")
        if np.any(self.transfer_emission[1] != 0):
            raise ValueError("Grass-to-forest conversion has zero emission in this baseline.")
        if np.any(self.decomposition <= 0) or np.any(self.humus_loss <= 0):
            raise ValueError("Decomposition and humus-loss rates must be positive.")
        if np.any((self.humification <= 0) | (self.humification > 1)):
            raise ValueError("Humification fractions must lie in (0, 1].")


def _check_cover(g):
    g = float(g)
    if not np.isfinite(g) or not 0.0 <= g <= 1.0:
        raise ValueError("Grass cover must be finite and lie in [0, 1].")
    return g


def _check_carbon_state(carbon):
    carbon = np.asarray(carbon, dtype=float)
    if carbon.shape != (10,) or not np.all(np.isfinite(carbon)):
        raise ValueError("Carbon state must be a finite vector with ten entries.")
    return carbon


def transition_rates(g, cover=CoverParameters()):
    """Return forest-to-grass and grass-to-forest rates, in that order."""
    g = _check_cover(g)
    return np.array([phi(g, cover), cover.alpha * (1.0 - g)], dtype=float)


def carbon_forcing(g, carbon_parameters, cover=CoverParameters()):
    """Return the ten-pool NPP input at the current area fractions."""
    g = _check_cover(g)
    group_npp = carbon_parameters.production * np.array([1.0 - g, g])
    forcing = np.zeros(10)
    forcing[:3] = carbon_parameters.allocation[0] * group_npp[0]
    forcing[5:8] = carbon_parameters.allocation[1] * group_npp[1]
    return forcing


def carbon_matrix(g, carbon_parameters, cover=CoverParameters()):
    """Return the block-diagonal 10x10 transfer matrix at fixed cover."""
    rates = transition_rates(g, cover)
    matrix = np.zeros((10, 10))
    for group in range(2):
        start = 5 * group
        living = slice(start, start + 3)
        dead = start + 3
        humus = start + 4
        turnover = carbon_parameters.turnover[group]
        emission = carbon_parameters.transfer_emission[group]
        rate = rates[group]

        matrix[living, living] = np.diag(-(turnover + rate))
        matrix[dead, start : start + 3] = (
            turnover + (1.0 - emission) * rate
        )
        matrix[dead, dead] = -carbon_parameters.decomposition[group]
        matrix[humus, dead] = (
            carbon_parameters.humification[group]
            * carbon_parameters.decomposition[group]
        )
        matrix[humus, humus] = -carbon_parameters.humus_loss[group]
    return matrix


def carbon_rhs(g, carbon, carbon_parameters, cover=CoverParameters()):
    """Evaluate carbon dynamics, preserving the specified area-conversion flows."""
    carbon = _check_carbon_state(carbon)
    return (
        carbon_matrix(g, carbon_parameters, cover) @ carbon
        + carbon_forcing(g, carbon_parameters, cover)
    )


def carbon_capacity(g, carbon_parameters, cover=CoverParameters()):
    """Frozen-cover equilibrium of the linear carbon subsystem."""
    matrix = carbon_matrix(g, carbon_parameters, cover)
    forcing = carbon_forcing(g, carbon_parameters, cover)
    return np.linalg.solve(-matrix, forcing)


def carbon_potential(g, carbon, carbon_parameters, cover=CoverParameters()):
    """Frozen-cover capacity minus the current carbon state."""
    return carbon_capacity(g, carbon_parameters, cover) - _check_carbon_state(carbon)


def lifted_equilibrium(g, carbon_parameters, cover=CoverParameters()):
    """Return the full 11-state equilibrium over a specified cover equilibrium."""
    g = _check_cover(g)
    state = np.empty(STATE_SIZE)
    state[0] = g
    state[1:] = carbon_capacity(g, carbon_parameters, cover)
    return state


def interior_lifted_equilibria(carbon_parameters, cover=CoverParameters()):
    """Lift all interior equilibria found by the established cover analysis."""
    return [lifted_equilibrium(g, carbon_parameters, cover)
            for g in interior_roots(cover)]


def rhs(state, carbon_parameters, cover=CoverParameters()):
    """Evaluate the coupled 11-state autonomous baseline."""
    state = np.asarray(state, dtype=float)
    if state.shape != (STATE_SIZE,) or not np.all(np.isfinite(state)):
        raise ValueError("State must be a finite vector with eleven entries.")
    g = _check_cover(state[0])
    derivative = np.empty(STATE_SIZE)
    derivative[0] = cover_drift(g, cover)
    derivative[1:] = carbon_rhs(g, state[1:], carbon_parameters, cover)
    return derivative


def _carbon_matrix_derivative(g, carbon_parameters, cover):
    rates_derivative = np.array([phi_prime(g, cover), -cover.alpha])
    matrix_derivative = np.zeros((10, 10))
    for group in range(2):
        start = 5 * group
        living = slice(start, start + 3)
        dead = start + 3
        rate_derivative = rates_derivative[group]
        matrix_derivative[living, living] = np.diag(
            np.full(3, -rate_derivative)
        )
        matrix_derivative[dead, start : start + 3] = (
            (1.0 - carbon_parameters.transfer_emission[group])
            * rate_derivative
        )
    return matrix_derivative


def jacobian(state, carbon_parameters, cover=CoverParameters()):
    """Return the analytic Jacobian, exposing its block-triangular structure."""
    state = np.asarray(state, dtype=float)
    if state.shape != (STATE_SIZE,) or not np.all(np.isfinite(state)):
        raise ValueError("State must be a finite vector with eleven entries.")
    g = _check_cover(state[0])
    matrix = carbon_matrix(g, carbon_parameters, cover)
    matrix_derivative = _carbon_matrix_derivative(g, carbon_parameters, cover)
    forcing_derivative = np.zeros(10)
    forcing_derivative[:3] = (
        -carbon_parameters.allocation[0] * carbon_parameters.production[0]
    )
    forcing_derivative[5:8] = (
        carbon_parameters.allocation[1] * carbon_parameters.production[1]
    )

    result = np.zeros((STATE_SIZE, STATE_SIZE))
    result[0, 0] = cover_drift_prime(g, cover)
    result[1:, 0] = matrix_derivative @ state[1:] + forcing_derivative
    result[1:, 1:] = matrix
    return result


def carbon_budget(g, carbon, carbon_parameters, cover=CoverParameters()):
    """Return total-stock budget terms and their numerical closure residual."""
    carbon = _check_carbon_state(carbon)
    derivative = carbon_rhs(g, carbon, carbon_parameters, cover)
    by_group = carbon.reshape(2, 5)
    rates = transition_rates(g, cover)
    npp = float(np.sum(carbon_forcing(g, carbon_parameters, cover)))
    heterotrophic_respiration = float(np.sum(
        (1.0 - carbon_parameters.humification)
        * carbon_parameters.decomposition
        * by_group[:, 3]
        + carbon_parameters.humus_loss * by_group[:, 4]
    ))
    transition_emission = float(np.sum(
        carbon_parameters.transfer_emission
        * rates[:, None]
        * by_group[:, :3]
    ))
    stock_change = float(np.sum(derivative))
    residual = stock_change - (
        npp - heterotrophic_respiration - transition_emission
    )
    return {
        "npp": npp,
        "heterotrophic_respiration": heterotrophic_respiration,
        "transition_emission": transition_emission,
        "stock_change": stock_change,
        "residual": residual,
    }
