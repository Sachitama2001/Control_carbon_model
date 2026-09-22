"""Numerical analysis for the independent temperature–NSC track."""
from dataclasses import dataclass
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares, root
from .temperature_nsc import PLANT, evaluate, rhs, jacobian, mortality, photosynthesis


@dataclass
class Trajectory:
    times: np.ndarray
    states: np.ndarray
    carbon_integral: np.ndarray
    min_stock: float
    integrated_budget_error: float
    endpoint_event: bool
    domain_event: bool


def integrate(x0, times, driver, p, max_step=30.0, rtol=1e-8, atol=1e-12):
    """Integrate log stocks, not clipped stocks, plus external carbon integral.

    Strictly positive starts required. Exact bare vegetation is an invariant
    boundary handled analytically, not seeded by an arbitrary epsilon.
    Stop at plant total 1e-12 as a *near-boundary diagnostic*, not an attractor
    certification. Samples include accepted steps in the positivity audit.
    """
    x0, times = np.asarray(x0, float), np.asarray(times, float)
    if x0.shape != (8,) or not np.isfinite(x0).all() or np.any(x0 <= 0):
        raise ValueError("log integration requires strictly positive stocks")
    if len(times) < 2 or not np.isfinite(times).all() or np.any(np.diff(times) <= 0):
        raise ValueError("increasing finite times required")
    def fun(t, z):
        x = np.exp(z[:8])
        dx, d = evaluate(x, driver(t), p)
        return np.r_[dx/x, d['external']]
    def log_jac(t, z):
        # The cumulative-budget column is identically zero. Supplying it avoids
        # adaptive numerical-difference overflow on that unobservable column
        # during very long stationary holds (SciPy's default num_jac).
        h = 1e-5
        j = np.zeros((9, 9))
        for i in range(8):
            e = np.eye(9)[i]*h
            j[:, i] = (fun(t, z+e)-fun(t, z-e))/(2*h)
        return j
    def near_boundary(t, z):
        return np.logaddexp.reduce(z[PLANT])-np.log(1e-12)
    near_boundary.terminal = True
    near_boundary.direction = -1
    def log_domain(t, z):
        return np.min(z[:8])+600
    log_domain.terminal = True
    log_domain.direction = -1
    sol = solve_ivp(fun, (times[0], times[-1]), np.r_[np.log(x0), 0.0],
                    method="Radau", jac=log_jac, dense_output=True, max_step=max_step,
                    rtol=rtol, atol=atol, events=(near_boundary, log_domain))
    if not sol.success:
        raise RuntimeError(sol.message)
    sample = np.unique(np.r_[times[times <= sol.t[-1]], sol.t[-1]])
    z = sol.sol(sample)
    x = np.exp(z[:8].T)
    audit_t = np.unique(np.r_[sol.t, sample,
                              *(sol.t[:-1]+q*np.diff(sol.t) for q in (.25, .5, .75))])
    audit_z = sol.sol(audit_t)
    audit_x = np.exp(audit_z[:8].T)
    if not np.isfinite(audit_x).all() or np.any(audit_x <= 0):
        raise RuntimeError("log integration lost finite positivity")
    err = np.max(np.abs(audit_x.sum(axis=1)-x0.sum()-audit_z[8]))
    return Trajectory(sample, x, z[8], float(audit_x.min()), float(err),
                      bool(len(sol.t_events[0])), bool(len(sol.t_events[1])))


def spinup(seed, temperature, p, years=1000, chunk_years=50, max_step=30,
           rtol=1e-8, atol=1e-12, rhs_tolerance=1e-9, drift_tolerance=1e-7):
    x = np.asarray(seed, float).copy()
    records, converged = [], False
    for start in np.arange(0, years, chunk_years):
        duration = min(chunk_years, years-start)*365
        tr = integrate(x, np.linspace(0, duration, int(duration/365)+1),
                       lambda t: temperature, p, max_step, rtol, atol)
        new = tr.states[-1]
        residual = float(np.max(np.abs(rhs(new, temperature, p))))
        drift = float(np.max(np.abs(new-tr.states[-2])/np.maximum(1, new)))
        records.append(dict(year=float(start+tr.times[-1]/365), state=new.tolist(),
                            sample_days=(start*365+tr.times).tolist(),
                            sample_states=tr.states.tolist(),
                            residual=residual, annual_drift=drift,
                            min_stock=tr.min_stock,
                            integrated_budget_error=tr.integrated_budget_error))
        x = new
        if tr.endpoint_event or tr.domain_event:
            break
        if residual < rhs_tolerance and drift < drift_tolerance:
            converged = True
            break
    return x, dict(converged=converged, near_boundary=tr.endpoint_event,
                   domain_limit=tr.domain_event,
                   records=records, initial=np.asarray(seed).tolist(),
                   temperature=temperature)


def relative_field(logx, temperature, p):
    x = np.exp(logx)
    return rhs(x, temperature, p)/x*365


def equilibrium(seed, temperature, p):
    """Positive direct root. Near-boundary roots are never called positive QSE."""
    fit = least_squares(lambda z: relative_field(z, temperature, p),
                        np.log(np.maximum(seed, 1e-12)), bounds=(-32, 16),
                        xtol=1e-12, ftol=1e-12, gtol=1e-12, max_nfev=500)
    x = np.exp(fit.x)
    if (np.max(np.abs(fit.fun)) > 1e-7 or x[PLANT].sum() < 1e-8):
        return None
    return x


def linear_info(x, temperature, p):
    j = jacobian(x, temperature, p)
    eig = np.linalg.eigvals(j)
    return dict(state=x.tolist(), temperature=float(temperature), jacobian=j.tolist(),
                eigen_real=eig.real.tolist(), eigen_imag=eig.imag.tolist(),
                max_real=float(eig.real.max()),
                residual=float(np.max(np.abs(rhs(x, temperature, p)))))


def arclength(x0, t0, p, direction=1, steps=200, step_size=.06):
    """Positive log-state pseudo-arclength; retains folds/unstable branches.

    T coordinate scaled by 10 K. Failure is reported, not silently bridged.
    Local smooth continuation only; does not prove no disconnected branches.
    """
    z = np.r_[np.log(x0), t0/10]
    def field(v):
        return relative_field(v[:8], v[8]*10, p)
    def tangent(v):
        h = 1e-5
        a = np.column_stack([(field(v+np.eye(9)[i]*h)-field(v-np.eye(9)[i]*h))/(2*h)
                             for i in range(9)])
        return np.linalg.svd(a, full_matrices=True)[2][-1]
    v = tangent(z)
    v *= direction*np.sign(v[-1])
    points = [linear_info(x0, t0, p)]
    reason = "step_limit"
    for _ in range(steps):
        pred = z+step_size*v
        def augmented(y):
            if np.any(y[:8] < -60) or np.any(y[:8] > 20):
                return np.full(9, 1e6)
            return np.r_[field(y), np.dot(y-pred, v)]
        fit = root(augmented, pred, tol=1e-9)
        if not fit.success or np.max(np.abs(augmented(fit.x))) > 1e-7:
            reason = "corrector_failed"
            break
        z = fit.x
        if not (0 < z[-1]*10 < 60) or np.min(z[:8]) < -25:
            reason = "domain_limit"
            break
        nv = tangent(z)
        if np.dot(nv, v) < 0:
            nv = -nv
        v = nv
        row = linear_info(np.exp(z[:8]), z[8]*10, p)
        row['temperature_tangent'] = float(v[-1])
        points.append(row)
    return dict(points=points, stop_reason=reason, step_size=step_size)


def smooth_temperature(t, duration, t0, t1):
    q = min(max(t/duration, 0), 1)
    return t0+(t1-t0)*(6*q**5-15*q**4+10*q**3)


def bare_invasion(temperature, p):
    """Homogeneous small-stock limit, not a spurious Jacobian at chi=0/0.

    Paid maintenance is quadratic near zero for the absolute-pool Michaelis
    constants used here. Common mortality cancels from projective dynamics.
    Supported-target flush is quadratic and has no first-order contribution.
    For constant target + all transfer edges positive, the dominant Perron
    direction yields the asymptotic interior invasion exponent.
    """
    if not p.leaf_dependence:
        return dict(status='bare_not_invariant')
    a = (p.flush_rate*p.leaf_target/p.flush_half
         if p.target_mode == 'constant' and p.flush_enabled else 0.0)
    m = np.zeros((6, 6))
    for i, tau in enumerate(p.turnover):
        m[i, i] -= tau
        m[i+3, i+3] -= tau
    # derivative of G at C_L=0
    m[3, 0] += photosynthesis(1.0, temperature, p)/(-np.expm1(-p.leaf_extinction))*p.leaf_extinction
    m[0, 4:6] += p.yield_c*a
    m[4, 4] -= a; m[5, 5] -= a
    for c, n, k in ((1,4,p.growth[0]), (2,5,p.growth[1])):
        m[c,n] += p.yield_c*k; m[n,n] -= k
    for donor, recipient, k in ((3,4,p.transfer[0]),(3,5,p.transfer[1]),
                                (4,3,p.transfer[2]),(5,3,p.transfer[3])):
        m[donor,donor] -= k; m[recipient,donor] += k
    eig, vec = np.linalg.eig(m)
    i = np.argmax(eig.real)
    v = vec[:,i].real
    if v.sum() < 0:
        v = -v
    v /= v.sum()
    chi = v[3:].sum()
    exponent = eig[i].real-mortality(chi,p)
    return dict(status='homogeneous_limit', matrix=m.tolist(), direction=v.tolist(),
                chi=float(chi), transport_eigenvalue=float(eig[i].real),
                invasion_exponent=float(exponent))
