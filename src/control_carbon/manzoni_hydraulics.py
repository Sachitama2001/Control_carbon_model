"""Independent reproduction of Manzoni et al. (2014), equations (1)--(4).

Potentials: MPa; time: day; fluxes: m water/day per ground area.
VPD is a mole fraction (mol/mol), NOT kPa. Constant capacitances are
in m water/MPa per ground area. This is neither VISIT nor a carbon model.
"""

from dataclasses import dataclass

import numpy as np


PAPER_DOI = "10.1002/2013WR015236"


@dataclass(frozen=True)
class HydraulicParameters:
    gx_max: float
    gs_max: float
    psi50_x: float
    psi50_s: float
    cx: float
    cl: float

    def __post_init__(self) -> None:
        values = np.array([
            self.gx_max, self.gs_max, self.psi50_x, self.psi50_s,
            self.cx, self.cl,
        ])
        if not np.all(np.isfinite(values)):
            raise ValueError("hydraulic parameters must be finite")
        if min(self.gx_max, self.gs_max, self.cx, self.cl) <= 0:
            raise ValueError("conductances and capacitances must be positive")
        if max(self.psi50_x, self.psi50_s) >= 0:
            raise ValueError("half-loss potentials must be negative")


def _potential(value: float) -> None:
    if not np.isfinite(value) or value > 0:
        raise ValueError("water potential must be finite and nonpositive (MPa)")


def _forcing(psi_root: float, vpd: float) -> None:
    _potential(psi_root)
    if not np.isfinite(vpd) or vpd < 0:
        raise ValueError("VPD must be finite and nonnegative (mol/mol)")


def vpd_mole_fraction(vpd_kpa: float, pressure_kpa: float) -> float:
    """Explicit gas-pressure conversion D = VPD/atmospheric pressure."""
    if not np.isfinite(vpd_kpa) or vpd_kpa < 0:
        raise ValueError("VPD must be finite and nonnegative (kPa)")
    if not np.isfinite(pressure_kpa) or pressure_kpa <= 0:
        raise ValueError("atmospheric pressure must be finite and positive (kPa)")
    return vpd_kpa / pressure_kpa


def conductances(psi_leaf: float, p: HydraulicParameters) -> tuple[float, float]:
    """Equations (3),(4), with the paper's zero lower bounds."""
    _potential(psi_leaf)
    return (
        p.gx_max * max(0.0, 1.0 - psi_leaf / (2.0 * p.psi50_x)),
        p.gs_max * max(0.0, 1.0 - psi_leaf / (2.0 * p.psi50_s)),
    )


def _state(state: np.ndarray) -> np.ndarray:
    state = np.asarray(state, dtype=float)
    if state.shape != (2,):
        raise ValueError("state must contain [psi_xylem, psi_leaf]")
    _potential(state[0])
    _potential(state[1])
    return state


def _fluxes(
    state: np.ndarray, psi_root: float, vpd: float, gx: float, gs: float,
) -> tuple[float, float, float]:
    return (
        2.0 * gx * (psi_root - state[0]),
        2.0 * gx * (state[0] - state[1]),
        gs * vpd,
    )


def _storage_rhs(fluxes: tuple[float, float, float], p: HydraulicParameters) -> np.ndarray:
    erx, exl, evaporation = fluxes
    return np.array([(erx - exl) / p.cx, (exl - evaporation) / p.cl])


def _state_jacobian(
    state: np.ndarray, psi_root: float, vpd: float, p: HydraulicParameters,
    gx: float, dgx: float, dgs: float,
) -> np.ndarray:
    x, leaf = state
    return np.array([
        [-4.0 * gx / p.cx,
         2.0 * (dgx * (psi_root - 2.0 * x + leaf) + gx) / p.cx],
        [2.0 * gx / p.cl,
         (2.0 * dgx * (x - leaf) - 2.0 * gx - dgs * vpd) / p.cl],
    ])


def water_fluxes(
    state: np.ndarray, psi_root: float, vpd: float, p: HydraulicParameters,
) -> tuple[float, float, float]:
    """Return root-to-xylem, xylem-to-leaf and atmospheric fluxes."""
    state = _state(state)
    _forcing(psi_root, vpd)
    gx, gs = conductances(state[1], p)
    return _fluxes(state, psi_root, vpd, gx, gs)


def hydraulic_rhs(
    state: np.ndarray, psi_root: float, vpd: float, p: HydraulicParameters,
) -> np.ndarray:
    """Full two-storage ODE, equations (1),(2); no state clipping."""
    return _storage_rhs(water_fluxes(state, psi_root, vpd, p), p)


def hydraulic_jacobian(
    state: np.ndarray, psi_root: float, vpd: float, p: HydraulicParameters,
) -> np.ndarray:
    """Full state Jacobian away from conductance-curve kinks."""
    water_fluxes(state, psi_root, vpd, p)
    leaf = state[1]
    if leaf == 2.0 * p.psi50_x or leaf == 2.0 * p.psi50_s:
        raise ValueError("Jacobian is not defined at a piecewise-linear kink")
    gx, gs = conductances(leaf, p)
    dgx = -p.gx_max / (2.0 * p.psi50_x) if gx > 0 else 0.0
    dgs = -p.gs_max / (2.0 * p.psi50_s) if gs > 0 else 0.0
    return _state_jacobian(state, psi_root, vpd, p, gx, dgx, dgs)


def steady_fluxes(
    psi_leaf: float, psi_root: float, vpd: float, p: HydraulicParameters,
) -> tuple[float, float]:
    """Supply and demand with conducting xylem at steady state (Fig. 2b)."""
    _forcing(psi_root, vpd)
    gx, gs = conductances(psi_leaf, p)
    return gx * (psi_root - psi_leaf), gs * vpd


def transpiring_equilibria(
    psi_root: float, vpd: float, p: HydraulicParameters,
) -> tuple[float, ...]:
    """Leaf potentials on the active branches of equation (5), wettest first.

    Reject algebraic roots outside BOTH positive-conductance intervals.
    Nontranspiring equilibria/continua are deliberately not enumerated.
    """
    _forcing(psi_root, vpd)
    a, b, r = -2.0 * p.psi50_x, -2.0 * p.psi50_s, -psi_root
    t = p.gs_max * vpd / p.gx_max
    coefficient = a + r + a * t / b
    constant = a * (r + t)
    discriminant = coefficient**2 - 4.0 * constant
    tolerance = 64.0 * np.finfo(float).eps * max(
        coefficient**2, 4.0 * abs(constant), 1.0,
    )
    if discriminant < -tolerance:
        return ()
    if abs(discriminant) <= tolerance:
        roots = (coefficient / 2.0,)
    else:
        large = (coefficient + np.sqrt(discriminant)) / 2.0
        roots = (constant / large, large)
    return tuple(
        -float(u) for u in roots
        if r <= u < min(a, b)
    )


def hydraulic_fold(
    psi_root: float, p: HydraulicParameters,
) -> tuple[float, float] | None:
    """Return (critical D in mol/mol, leaf MPa) for an interior active fold.

    Safe trait ordering and roots beyond xylem closure have no such fold.
    This is not a claim that arbitrary soil drying is harmless.
    """
    _potential(psi_root)
    a, b, r = -2.0 * p.psi50_x, -2.0 * p.psi50_s, -psi_root
    if b <= a or r >= a:
        return None
    u = b - np.sqrt((b - a) * (b - r))
    critical_t = b / a * (np.sqrt(b - r) - np.sqrt(b - a))**2
    return p.gx_max / p.gs_max * critical_t, -float(u)
