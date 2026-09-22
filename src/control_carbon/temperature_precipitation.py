"""Eight-state T/P theoretical model; not a native VISIT or calibrated site model.

Equations: docs/temperature_precipitation_tipping_research_plan.md §§6–8.
Explicit closure changes: docs/temperature_precipitation_implementation.md.
Time: days; carbon: Mg C ha-1; water: mm per ground area. All defaults are
uncalibrated numerical hypotheses. No native model coefficients are imported.
"""

from __future__ import annotations

from calendar import monthrange
from collections.abc import Callable
from dataclasses import asdict, dataclass, field, fields

import numpy as np
from numpy.typing import ArrayLike
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares
from scipy.special import expit

STATE_NAMES = ("C_L", "C_S", "C_R", "C_O", "W_L", "W_S", "W_R", "W_O")
STATE_UNITS = ("Mg C ha-1",) * 4 + ("mm",) * 4
STATE_SCALE = np.array([10., 200., 40., 200., 5., 30., 10., 300.])
WATER_INCIDENCE = np.array([[0., 0., 1.], [0., 1., -1.],
                            [1., -1., 0.], [-1., 0., 0.]])


def assumed(value, unit: str, section: str, meaning: str):
    return field(default=value, metadata={"unit": unit, "source": section,
                 "meaning": meaning, "status": "uncalibrated assumption",
                 "ecological_range": None})


@dataclass(frozen=True)
class Parameters:
    allocation: tuple = assumed((.3, .4, .3), "1", "plan §6.2", "葉・幹・根への配分")
    turnover: tuple = assumed((.002, .00008, .0006), "day-1", "plan §6.2", "通常ターンオーバー")
    respiration: tuple = assumed((.0006, .00003, .00015), "day-1", "plan §6.2", "維持呼吸の基準率")
    soil_respiration: float = assumed(.0002, "day-1", "plan §6.2", "土壌呼吸の基準率")
    eta: tuple = assumed((1., 1., 1.), "1", "plan §6.2", "枯死Cの土壌移行割合")
    gpp_max: float = assumed(.08, "Mg C ha-1 day-1", "plan §6.1", "最大群落GPP")
    k_leaf: float = assumed(.4, "ha (Mg C)-1", "plan §6.1", "葉量に対する受光飽和")
    t_ref: float = assumed(27., "degC", "plan §6.2", "温度応答の基準")
    t_opt: float = assumed(27., "degC", "plan §6.1", "GPP最適温度")
    t_width: float = assumed(8., "degC", "plan §6.1", "GPP温度応答の幅")
    q10_resp: float = assumed(2., "1", "plan §6.2", "植物呼吸Q10")
    q10_soil: float = assumed(2., "1", "plan §6.2", "土壌呼吸Q10")
    photo_half: float = assumed(.3, "1", "plan §6.1", "GPP相対水分の半飽和値")
    photo_power: float = assumed(2., "1", "plan §6.1", "GPP水分応答の指数")
    soil_half: float = assumed(.3, "1", "implementation closure", "分解水分応答の半飽和値")
    omega: tuple = assumed((.4, .1, .2), "mm ha (Mg C)-1", "plan §7.1", "器官C当たり水容量")
    epsilon: tuple = assumed((.02, .02, .02), "mm", "plan §7.1", "容量の正則化")
    soil_dry: float = assumed(0., "mm", "plan §7.1", "土壌の利用可能水下限")
    soil_fc: float = assumed(240., "mm", "plan §7.1", "土壌圃場容水量")
    soil_sat: float = assumed(360., "mm", "implementation closure", "土壌飽和貯留量")
    conductance: tuple = assumed((25., 20., 15.), "mm day-1", "plan §7.2", "OR・RS・SLの相対水分輸送係数")
    flow_smoothing: float = assumed(.01, "1", "implementation closure", "順方向C1輸送の平滑幅")
    transpiration: float = assumed(3., "mm day-1", "plan §7.2", "基準蒸散強度")
    evaporation: float = assumed(.6, "mm day-1", "plan §7.2", "基準土壌蒸発強度")
    evap_half: float = assumed(.2, "1", "implementation closure", "蒸散・蒸発の半飽和水分")
    k_area: float = assumed(.4, "ha (Mg C)-1", "plan §7.2", "蒸散面積飽和係数")
    q10_evap: float = assumed(1.5, "1", "implementation closure", "温度による蒸発散需要の集約Q10")
    drainage: float = assumed(.025, "day-1", "implementation closure", "圃場容水量超過水の排水率")
    runoff_power: float = assumed(8., "1", "implementation closure", "供給に対する流出率の貯留指数")
    zeta: tuple = assumed((1., 1., 1.), "1", "plan §7.2", "植物放出水の土壌還流割合")
    heat_max: float = assumed(.001, "day-1", "plan §8.1", "慢性的熱枯死の最大率")
    heat_half: float = assumed(35., "degC", "plan §8.1", "月代表温度に対する枯死半値温度")
    heat_slope: float = assumed(.7, "degC-1", "plan §8.1", "慢性的熱枯死の傾き")
    drought_max: float = assumed(.001, "day-1", "plan §8.2", "水分枯死の最大率")
    drought_half: float = assumed(.35, "1", "plan §8.2", "枯死半値相対水分")
    drought_slope: float = assumed(15., "1", "plan §8.2", "水分枯死の傾き")

    def __post_init__(self):
        for f in fields(self):
            v = np.asarray(getattr(self, f.name), dtype=float)
            if not np.all(np.isfinite(v)):
                raise ValueError(f"non-finite parameter: {f.name}")
            if isinstance(f.default, tuple) and v.shape != (3,):
                raise ValueError(f"three organ/edge values required: {f.name}")
            if f.name not in ("t_ref", "t_opt", "heat_half") and np.any(v < 0):
                raise ValueError(f"negative parameter: {f.name}")
        if not np.isclose(sum(self.allocation), 1., atol=1e-12, rtol=0):
            raise ValueError("allocation must sum to one")
        for key in ("eta", "zeta"):
            if np.any(np.asarray(getattr(self, key)) > 1):
                raise ValueError(f"{key} must lie in [0,1]")
        for key in ("t_width", "photo_half", "photo_power", "soil_half",
                    "epsilon", "flow_smoothing", "evap_half", "q10_resp",
                    "q10_soil", "q10_evap", "runoff_power"):
            if np.any(np.asarray(getattr(self, key)) <= 0):
                raise ValueError(f"positive parameter required: {key}")
        if not 0 <= self.soil_dry < self.soil_fc < self.soil_sat:
            raise ValueError("require dry < field capacity < saturation")


@dataclass(frozen=True)
class Climate:
    temperature: float = 27.
    precipitation: float = 6.

    def __post_init__(self):
        if not np.isfinite(self.temperature) or not np.isfinite(self.precipitation):
            raise ValueError("climate must be finite")
        if self.precipitation < 0:
            raise ValueError("negative precipitation")

    @classmethod
    def from_month(cls, temperature, precipitation_mm, year, month):
        """Convert a monthly total to a daily flux using Gregorian month length."""
        return cls(temperature, precipitation_mm / monthrange(year, month)[1])


@dataclass(frozen=True)
class Controls:
    mortality: str = "none"  # none -> heat -> water -> additive
    photo_temperature: float | None = None
    respiration_temperature: float | None = None
    evap_temperature: float | None = None
    heat_temperature: float | None = None
    drought_rate: float | None = None
    temperature_response: bool = True
    water_response: bool = True

    def __post_init__(self):
        if self.mortality not in ("none", "heat", "water", "additive"):
            raise ValueError("unknown mortality hierarchy")
        for key in ("photo_temperature", "respiration_temperature",
                    "evap_temperature", "heat_temperature", "drought_rate"):
            value = getattr(self, key)
            if value is not None and not np.isfinite(value):
                raise ValueError(f"non-finite clamp: {key}")
        if self.drought_rate is not None and self.drought_rate < 0:
            raise ValueError("negative clamped mortality")


def capacity(c, p: Parameters):
    return np.asarray(p.omega) * np.asarray(c)[:3] + p.epsilon


def relative_water(x, p: Parameters):
    return np.r_[x[4:7] / capacity(x[:4], p),
                 np.clip((x[7]-p.soil_dry)/(p.soil_fc-p.soil_dry), 0., 1.)]


def positive_flow(z, width):
    pos = np.maximum(z, 0.)
    return pos * pos / (pos + width)


def _temperature(clamp, actual):
    return actual if clamp is None else clamp


@dataclass(frozen=True)
class Fluxes:
    gpp: float
    respiration: np.ndarray
    soil_respiration: float
    mortality_rates: np.ndarray
    mortality: np.ndarray
    turnover: np.ndarray
    carbon_export: float
    transport: np.ndarray
    transpiration: float
    evaporation: float
    drainage: float
    runoff: float
    water_release: np.ndarray
    water_export: float
    heat_hazard: float
    drought_hazard: float
    carbon_input: np.ndarray
    carbon_matrix: np.ndarray

    @property
    def nep(self):
        return self.gpp - self.respiration.sum() - self.soil_respiration


def fluxes(state: ArrayLike, climate: Climate, p=Parameters(), control=Controls()):
    """Pure flux function. Positive-part extension is only for solver trial states.

    Accepted trajectories are never clipped and are independently domain-checked.
    All budget identities hold also for this explicitly defined extension.
    """
    raw = np.asarray(state, dtype=float)
    if raw.shape != (8,) or not np.all(np.isfinite(raw)):
        raise ValueError("eight finite states required")
    x = np.maximum(raw, 0.)
    c, w = x[:4], x[4:]
    r = relative_water(x, p)
    tc = _temperature(control.photo_temperature, climate.temperature)
    tr = _temperature(control.respiration_temperature, climate.temperature)
    te = _temperature(control.evap_temperature, climate.temperature)
    th = _temperature(control.heat_temperature, climate.temperature)
    fg = np.exp(-.5*((tc-p.t_opt)/p.t_width)**2) if control.temperature_response else 1.
    qr = p.q10_resp**((tr-p.t_ref)/10) if control.temperature_response else 1.
    qh = p.q10_soil**((tr-p.t_ref)/10) if control.temperature_response else 1.
    qe = p.q10_evap**((te-p.t_ref)/10) if control.temperature_response else 1.
    rg = np.sqrt(r[0]*r[2])
    fw = rg**p.photo_power/(rg**p.photo_power+p.photo_half**p.photo_power)
    fh = r[3]/(r[3]+p.soil_half)
    if not control.water_response:
        fw = fh = 1.
    g = p.gpp_max*fg*fw*(-np.expm1(-p.k_leaf*c[0]))
    rp = w[:3].sum()/capacity(c, p).sum()
    mt = p.heat_max*expit(p.heat_slope*(th-p.heat_half))
    mw = p.drought_max*expit(-p.drought_slope*(rp-p.drought_half))
    if control.drought_rate is not None:
        mw = control.drought_rate
    m = (mt if control.mortality in ("heat", "additive") else 0.)
    m += mw if control.mortality in ("water", "additive") else 0.
    mortality = np.full(3, m)
    resp_rate = np.asarray(p.respiration)*qr
    losses = np.asarray(p.turnover)+mortality+resp_rate
    u = np.r_[np.asarray(p.allocation)*g, 0.]
    matrix = np.diag(np.r_[-losses, -p.soil_respiration*qh*fh])
    matrix[3, :3] = np.asarray(p.turnover)+np.asarray(p.eta)*mortality
    q = np.asarray(p.conductance)*positive_flow(r[[3, 2, 1]]-r[[2, 1, 0]], p.flow_smoothing)
    el = p.transpiration*qe*r[0]/(r[0]+p.evap_half)*(-np.expm1(-p.k_area*c[0]))
    eo = p.evaporation*qe*r[3]/(r[3]+p.evap_half)
    release = losses*w[:3]
    incoming = climate.precipitation + np.asarray(p.zeta)@release
    runoff = incoming*(w[3]/p.soil_sat)**p.runoff_power
    return Fluxes(float(g), resp_rate*c[:3], float(-matrix[3, 3]*c[3]),
                  mortality, mortality*c[:3], np.asarray(p.turnover)*c[:3],
                  float(((1-np.asarray(p.eta))*mortality)@c[:3]), q,
                  float(el), float(eo), p.drainage*max(w[3]-p.soil_fc, 0.),
                  float(runoff), release, float((1-np.asarray(p.zeta))@release),
                  float(mt), float(mw), u, matrix)


def rhs(t, state, climate: Climate, p=Parameters(), control=Controls()):
    f = fluxes(state, climate, p, control)
    c = np.maximum(np.asarray(state)[:4], 0.)
    dw = WATER_INCIDENCE@f.transport
    dw[:3] -= f.water_release
    dw[0] -= f.transpiration
    dw[3] += (climate.precipitation + np.asarray(p.zeta)@f.water_release
              - f.evaporation - f.drainage - f.runoff)
    return np.r_[f.carbon_input + f.carbon_matrix@c, dw]


def budget_residuals(state, climate, p=Parameters(), control=Controls()):
    f = fluxes(state, climate, p, control)
    dx = rhs(0., state, climate, p, control)
    return np.array([dx[:4].sum()-(f.nep-f.carbon_export),
                     dx[4:].sum()-(climate.precipitation-f.transpiration
                     -f.evaporation-f.drainage-f.runoff-f.water_export)])


def frozen_carbon_capacity(state, climate, p=Parameters(), control=Controls()):
    """Freeze current GPP AND matrix, including their carbon dependence.

    This is not the solution obtained by re-solving nonlinear carbon feedbacks.
    Singular loss matrices (e.g. zero decomposition at zero water) are rejected.
    """
    f = fluxes(state, climate, p, control)
    return np.linalg.solve(-f.carbon_matrix, f.carbon_input)


def jacobian(state, climate, p=Parameters(), control=Controls(), step=1e-5):
    x = np.asarray(state, dtype=float)
    j = np.empty((8, 8))
    for k in range(8):
        h = step*max(1., abs(x[k]))
        delta = np.zeros(8)
        delta[k] = h
        if x[k] >= h:
            j[:, k] = (rhs(0, x+delta, climate, p, control)-rhs(0, x-delta, climate, p, control))/(2*h)
        else:
            j[:, k] = (rhs(0, x+delta, climate, p, control)-rhs(0, x, climate, p, control))/h
    return j


@dataclass(frozen=True)
class Equilibrium:
    state: np.ndarray
    residual: float
    eigenvalues: np.ndarray
    condition_number: float
    evaluations: int

    @property
    def max_real_part(self):
        return float(self.eigenvalues.real.max())


def equilibrium(climate=Climate(), p=Parameters(), control=Controls(), guess=None):
    """Locally solve one branch; this routine does not prove uniqueness."""
    x0 = np.array([7., 180., 25., 180., 1.6, 13., 4., 290.]) if guess is None else np.asarray(guess)
    fit = least_squares(lambda y: rhs(0, y*STATE_SCALE, climate, p, control)/STATE_SCALE,
                        x0/STATE_SCALE, bounds=(0., np.inf),
                        ftol=1e-13, xtol=1e-13, gtol=None, max_nfev=2000)
    x = fit.x*STATE_SCALE
    residual = float(np.max(np.abs(rhs(0, x, climate, p, control))))
    if not fit.success or residual > 1e-8 or domain_violation(x, p) > 1e-6:
        raise RuntimeError(f"equilibrium not verified: residual={residual}, {fit.message}")
    j = jacobian(x, climate, p, control)
    return Equilibrium(x, residual, np.linalg.eigvals(j), float(np.linalg.cond(j)), fit.nfev)


def domain_violation(state, p=Parameters()):
    x = np.asarray(state)
    return float(max(0., -x.min(), (x[4:7]-capacity(x[:4], p)).max(), x[7]-p.soil_sat))


def integrate(initial, times, climate: Climate | Callable, p=Parameters(), control=Controls(),
              *, method="Radau", max_step=10., rtol=1e-8, atol=1e-10):
    """Forcing must be continuous inside this interval; split monthly jumps.

    Internal steps are <=10 days by default. Observations may be monthly or
    yearly; they do not define the numerical integration step.
    """
    times = np.asarray(times, dtype=float)
    if times.ndim != 1 or len(times) < 2 or np.any(np.diff(times) <= 0):
        raise ValueError("strictly increasing output times required")
    if domain_violation(initial, p) > 1e-8:
        raise ValueError("initial state violates physical domain")
    forcing = climate if callable(climate) else lambda t: climate
    sol = solve_ivp(lambda t, x: rhs(t, x, forcing(t), p, control),
                    (times[0], times[-1]), initial, method=method, t_eval=times,
                    max_step=max_step, rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError(sol.message)
    if max(domain_violation(x, p) for x in sol.y.T) > 1e-6:
        raise RuntimeError("accepted trajectory left physical domain; no clipping applied")
    return sol.y.T


def integrate_months(initial, climates, month_days, p=Parameters(), control=Controls(),
                     **solver):
    """Integrate each piecewise constant month separately, returning boundaries."""
    if len(climates) != len(month_days) or len(climates) == 0:
        raise ValueError("one length per month required")
    states, times = [np.asarray(initial)], [0.]
    for climate, days in zip(climates, month_days):
        if not np.isfinite(days) or days <= 0:
            raise ValueError("positive finite month length required")
        end = times[-1]+days
        states.append(integrate(states[-1], [times[-1], end], climate, p, control, **solver)[-1])
        times.append(end)
    return np.array(times), np.array(states)


def parameter_manifest(p=Parameters()):
    return {f.name: {"value": asdict(p)[f.name], **f.metadata} for f in fields(p)}
