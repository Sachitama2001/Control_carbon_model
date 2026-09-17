"""Model-agnostic C--N benchmark; NOT a VISIT transcription or forest fit.

States c=C/(q*N0), n=N/N0, optionally d=D/(q*N0); time is years.
Andrews (1968), doi:10.1002/bit.260100602, eq. (1) supplies only the
substrate-inhibition shape. Budgets, recycling and parameters are hypotheses.
"""
from dataclasses import dataclass
import numpy as np
from scipy.integrate import solve_ivp


@dataclass(frozen=True)
class Parameters:
    g0: float = 0.4
    m: float = 0.1
    loss: float = 1.0
    recycle: float = 0.5
    inhibition: bool = True
    decomposition: float | None = None

    def __post_init__(self):
        if min(self.g0, self.m, self.loss) <= 0 or not 0 <= self.recycle < 1:
            raise ValueError("Positive rates and 0 <= recycle < 1 required")
        if self.decomposition is not None and self.decomposition <= 0:
            raise ValueError("Decomposition must be positive")


def growth(n, p=Parameters()):
    """NPP per unit live carbon [yr^-1]; g0 is NOT the inhibited maximum."""
    return p.g0 * n / (1 + n + (n*n if p.inhibition else 0))


def growth_derivative(n, p=Parameters()):
    if p.inhibition:
        return p.g0 * (1-n*n) / (1+n+n*n)**2
    return p.g0 / (1+n)**2


def rhs(y, influx, p=Parameters()):
    """Physical normalized states (c,n), or (c,d,n), not log coordinates."""
    c, n = y[0], y[-1]
    g = growth(n, p)
    if p.decomposition is None:
        return np.array([(g-p.m)*c,
                         influx-p.loss*n-(g-p.recycle*p.m)*c])
    d, k = y[1], p.decomposition
    return np.array([(g-p.m)*c, p.m*c-k*d,
                     influx-p.loss*n-g*c+p.recycle*k*d])


def jacobian(y, p=Parameters()):
    c, n = y[0], y[-1]
    g, gp = growth(n, p), growth_derivative(n, p)
    if p.decomposition is None:
        return np.array([[g-p.m, c*gp],
                         [-g+p.recycle*p.m, -p.loss-c*gp]])
    k = p.decomposition
    return np.array([[g-p.m, 0, c*gp], [p.m, -k, 0],
                     [-g, p.recycle*k, -p.loss-c*gp]])


def equilibria(influx, p=Parameters()):
    """Only nonnegative equilibria; forest=low N, saddle=high N root."""
    bare = [0, influx/p.loss] if p.decomposition is None else [0, 0, influx/p.loss]
    result = {"bare": np.array(bare, dtype=float)}
    if p.inhibition:
        b = p.g0/p.m-1
        roots = [] if b < 2 else [("forest", (b-np.sqrt(b*b-4))/2),
                                  ("saddle", (b+np.sqrt(b*b-4))/2)]
    else:
        roots = [] if p.g0 <= p.m else [("forest", p.m/(p.g0-p.m))]
    for name, n in roots:
        c = (influx-p.loss*n)/((1-p.recycle)*p.m)
        if c > 0:
            result[name] = np.array([c, n] if p.decomposition is None
                                    else [c, p.m*c/p.decomposition, n])
    return result


def influx_at(t, duration, initial=3., final=10., smooth=False):
    if duration <= 0:
        raise ValueError("Ramp duration must be positive")
    s = np.clip(np.asarray(t)/duration, 0, 1)
    if smooth:
        s = s*s*(3-2*s)
    return initial+(final-initial)*s


def simulate(duration, p=Parameters(), initial=3., final=10., hold=1000.,
             smooth=False, rtol=1e-9, max_step=np.inf, sample_step=1.):
    """Continuous forcing, split at ramp end; monthly/yearly OUTPUT maps.

    Integrate log(c) so no imposed extinction threshold changes the dynamics.
    This is not a forward-Euler annual update or a seasonally resolved model.
    """
    y0 = equilibria(initial, p)["forest"].copy()
    y0[0] = np.log(y0[0])

    def fun(t, z):
        y = z.copy()
        y[0] = np.exp(z[0])
        f = rhs(y, influx_at(t, duration, initial, final, smooth), p)
        f[0] = growth(y[-1], p)-p.m
        return f

    times, states = [], []
    for a, b in [(0., duration), (duration, duration+hold)]:
        if b <= a:
            continue
        sol = solve_ivp(fun, (a, b), y0, method="Radau", rtol=rtol,
                        atol=rtol*0.01, max_step=max_step, dense_output=True)
        if not sol.success:
            raise RuntimeError(sol.message)
        ts = np.unique(np.r_[np.arange(a, b, sample_step), b])
        if times:
            ts = ts[1:]
        ys = sol.sol(ts).T
        ys[:, 0] = np.exp(ys[:, 0])
        times.extend(ts)
        states.extend(ys)
        y0 = sol.y[:, -1]
    return np.asarray(times), np.asarray(states)


def classify(y, final=10., p=Parameters(), tol=1e-5):
    """Equilibrium convergence, not a carbon-only extinction cutoff."""
    for name in ("forest", "bare"):
        eq = equilibria(final, p).get(name)
        if eq is not None and np.max(np.abs(y-eq)/(1+np.abs(eq))) < tol:
            return name
    return "unresolved"


def critical_duration(p=Parameters(), lo=10., hi=40., iterations=16,
                      rtol=1e-9, hold=3000., smooth=False, max_step=np.inf):
    """Local bracket only: does NOT assert uniqueness across all protocols."""
    def outcome(t):
        return classify(simulate(t, p, hold=hold, rtol=rtol, smooth=smooth,
                                 max_step=max_step)[1][-1], p=p)
    if (outcome(lo), outcome(hi)) != ("bare", "forest"):
        raise ValueError("Endpoints must bracket bare and forest outcomes")
    for _ in range(iterations):
        mid = (lo+hi)/2
        result = outcome(mid)
        if result == "unresolved":
            raise RuntimeError("Longer final hold required near separatrix")
        if result == "bare":
            lo = mid
        else:
            hi = mid
    return lo, hi


def fast_step_certificate():
    """Analytic bootstrap constants for DEFAULTS, instantaneous 3 -> 10."""
    p = Parameters()
    c0, n0 = equilibria(3., p)["forest"]
    cb, nb = 60., 3.
    gmax = p.g0/3
    n_speed_lower = 10-nb-(gmax-p.recycle*p.m)*cb
    time_upper = (nb-n0)/n_speed_lower
    c_upper = c0*np.exp((gmax-p.m)*time_upper)
    inward = 10-nb-(growth(nb)-p.recycle*p.m)*cb
    decay = p.m-growth(nb)
    return dict(c_upper=float(c_upper), c_boundary=cb, n_boundary=nb,
                time_upper=float(time_upper), n_speed_lower=n_speed_lower,
                inward=float(inward), decay=float(decay))
