"""Independent, assumed continuous temperature–NSC mass-balance model.

Not a transcription of VISIT or SEIB. See docs/temperature_nsc_equations.md.
Stocks: C_L,C_S,C_R,C_O,N_L,N_S,N_R,N_O [Mg C/ha]; time [day].
"""
from dataclasses import dataclass
import numpy as np

NAMES = ("C_L", "C_S", "C_R", "C_O", "N_L", "N_S", "N_R", "N_O")
PLANT = np.array([0, 1, 2, 4, 5, 6])


@dataclass(frozen=True)
class Parameters:
    gmax: float = 0.05
    leaf_extinction: float = 0.7
    t_opt: float = 20.0
    sigma_g: float = 12.0
    t_ref: float = 20.0
    q10: float = 2.0
    respiration: tuple = (0.001, 0.00004, 0.0003)
    substrate_half: tuple = (0.02, 0.2, 0.05)
    turnover: tuple = (1/365, 1/36500, 1/1095)
    growth: tuple = (0.001, 0.003)
    yield_c: float = 0.75
    leaf_target: float = 3.0
    flush_rate: float = 0.03
    flush_half: float = 0.5
    transfer: tuple = (0.02, 0.01, 0.001, 0.001)
    mu0: float = 1/36500
    mu_max: float = 0.02
    chi_crit: float = 0.03
    hill: float = 4.0
    soil_k: tuple = (1/3650, 1/60)
    soil_q10: float = 2.0
    humification: float = 0.2
    retention: float = 1.0
    # Experimental switches; all defaults are the plan's autonomous candidate.
    leaf_dependence: bool = True
    high_temperature_decline: bool = True
    nsc_mortality: bool = True
    flush_enabled: bool = True
    target_mode: str = "constant"
    support_half: float = 30.0

    def __post_init__(self):
        if not np.isfinite([self.t_opt, self.t_ref]).all():
            raise ValueError("finite reference temperatures required")
        positive = (self.gmax, self.leaf_extinction, self.sigma_g, self.q10,
                    self.leaf_target, self.flush_half, self.chi_crit, self.hill,
                    self.soil_q10, self.support_half, *self.substrate_half,
                    *self.soil_k)
        nonnegative = (self.flush_rate, self.mu0, self.mu_max, *self.respiration,
                       *self.turnover, *self.growth, *self.transfer)
        if not all(np.isfinite(v) and v > 0 for v in positive):
            raise ValueError("positive finite coefficients required")
        if not all(np.isfinite(v) and v >= 0 for v in nonnegative):
            raise ValueError("nonnegative finite rates required")
        if not (0 < self.yield_c <= 1 and 0 <= self.retention <= 1
                and 0 <= self.humification <= 1):
            raise ValueError("invalid efficiency/retention fraction")
        if self.target_mode not in ("constant", "supported"):
            raise ValueError("unknown leaf target closure")
        for name, size in (("respiration", 3), ("substrate_half", 3),
                           ("turnover", 3), ("growth", 2), ("transfer", 4),
                           ("soil_k", 2)):
            if len(getattr(self, name)) != size:
                raise ValueError(f"{name} requires {size} entries")


def concentrations(x):
    """NaN is diagnostic at the genuinely empty plant/organ, never epsilon."""
    x = np.asarray(x, dtype=float)
    c, n = x[:3], x[4:7]
    organ = np.divide(n, c+n, out=np.full(3, np.nan), where=c+n > 0)
    total = c.sum()+n.sum()
    return (float(n.sum()/total) if total > 0 else np.nan), organ


def dry_to_carbon_fraction(q, f_nsc, f_structure):
    """q is NSC dry mass / total dry mass; carbon fractions must be supplied."""
    if not (0 <= q <= 1 and 0 < f_nsc <= 1 and 0 < f_structure <= 1):
        raise ValueError("invalid mass fractions")
    return f_nsc*q/(f_structure*(1-q)+f_nsc*q)


def photosynthesis(leaf, temperature, p):
    thermal = np.exp(-0.5*((temperature-p.t_opt)/p.sigma_g)**2)
    if not p.high_temperature_decline and temperature > p.t_opt:
        thermal = 1.0
    canopy = -np.expm1(-p.leaf_extinction*leaf) if p.leaf_dependence else 1.0
    return p.gmax*thermal*canopy


def maintenance(c, n, temperature, p):
    demand = np.asarray(p.respiration)*c*p.q10**((temperature-p.t_ref)/10)
    paid = demand*n/(np.asarray(p.substrate_half)+n)
    return demand, paid


def mortality(chi, p):
    if np.isnan(chi):
        # Empty plant has zero outgoing flux; no concentration assigned to it.
        return p.mu0
    return p.mu0 + (p.mu_max/(1+(chi/p.chi_crit)**p.hill)
                    if p.nsc_mortality else 0.0)


def flush_draw(c, n, p, phi=1.0):
    """Donor-specific rates, not daily amounts; each vanishes at empty donor.

    Algebraically cancels the singular donor-share expression at N_S+N_R=0.
    An arbitrary explicit Euler step is NOT guaranteed positive.
    """
    target = p.leaf_target
    if p.target_mode == "supported":
        support = c[1]+c[2]
        target *= support/(p.support_half+support)
    deficit = max(target-c[0], 0.0)
    coefficient = p.flush_rate*phi*deficit/(p.flush_half+n[1]+n[2])
    return coefficient*n[1:] if p.flush_enabled else np.zeros(2)


def transfers(n, p):
    """Donor-proportional L→S,L→R,S→L,R→L (assumed, conservative)."""
    return np.asarray(p.transfer)*n[[0, 0, 1, 2]]


def evaluate(x, temperature, p=Parameters(), phi=1.0):
    x = np.asarray(x, dtype=float)
    if x.shape != (8,) or not np.isfinite(x).all() or np.any(x < 0):
        raise ValueError("state must be finite nonnegative length 8")
    if not np.isfinite(temperature) or not np.isfinite(phi) or phi < 0:
        raise ValueError("invalid forcing")
    c, n = x[:3], x[4:7]
    chi, organ = concentrations(x)
    mu = mortality(chi, p)
    g = photosynthesis(c[0], temperature, p)
    demand, paid = maintenance(c, n, temperature, p)
    draw = flush_draw(c, n, p, phi)
    growth = np.asarray(p.growth)*n[1:]
    use = np.r_[draw.sum(), growth]
    j = transfers(n, p)
    loss_rate = np.asarray(p.turnover)+mu
    lc, ln = loss_rate*c, loss_rate*n
    soil = np.asarray(p.soil_k)*p.soil_q10**((temperature-p.t_ref)/10)*x[[3, 7]]
    dx = np.zeros(8)
    dx[:3] = p.yield_c*use-lc
    dx[4:7] = -paid-ln
    dx[4] += g-j[0]-j[1]+j[2]+j[3]
    dx[5] += j[0]-j[2]-draw[0]-growth[0]
    dx[6] += j[1]-j[3]-draw[1]-growth[1]
    dx[3] = p.retention*lc.sum()+p.humification*soil[1]-soil[0]
    dx[7] = p.retention*ln.sum()-soil[1]
    rg = (1-p.yield_c)*use.sum()
    rh = soil[0]+(1-p.humification)*soil[1]
    export = (1-p.retention)*(lc.sum()+ln.sum())
    external = g-paid.sum()-rg-rh-export
    diagnostics = dict(gpp=g, demand=demand, paid=paid, unpaid=demand-paid,
                       flush=draw.sum(), growth=growth, growth_resp=rg,
                       mortality=mu, rh=rh, export=export, chi=chi, organ_chi=organ,
                       external=external, budget=float(dx.sum()-external))
    return dx, diagnostics


def rhs(x, temperature, p=Parameters()):
    return evaluate(x, temperature, p)[0]


def jacobian(x, temperature, p=Parameters()):
    """Interior centered finite difference; bare origin has no unique Jacobian."""
    x = np.asarray(x, dtype=float)
    if x.shape != (8,) or not np.isfinite(x).all() or np.any(x <= 0):
        raise ValueError("interior Jacobian only; analyze bare boundary separately")
    steps = np.minimum(1e-5*np.maximum(1.0, x), 0.1*x)
    return np.column_stack([(rhs(x+np.eye(8)[i]*h, temperature, p)
                             -rhs(x-np.eye(8)[i]*h, temperature, p))/(2*h)
                            for i, h in enumerate(steps)])
