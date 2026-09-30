"""Deterministic checks of the provisional ten-pool equations, not calibration.

Run: python examples/check_tree_grass_carbon.py
No parameter search, native VISIT transcription, or climate-rate experiment.
All stocks: kg C m^-2; time: yr. Output is JSON to stdout.
"""
import hashlib
import json
from pathlib import Path
import platform

import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

# Illustrative assumptions only; these are not measured Cerrado coefficients.
P = np.array([2.0, 1.8])
a = np.array([2.0, 3.0])
p = np.array([[0.4, 0.35, 0.25], [0.5, 0.2, 0.3]])
d = np.array([[0.5, 0.08, 0.2], [0.8, 0.6, 0.7]])
eta = np.array([[0.8, 0.4, 0.0], [0.9, 0.8, 0.0]])
b = np.array([[0.3, 0.1, 0.0], [0.8, 0.8, 0.0]])
c = np.array([0.4, 0.8])
kD, kH, rho, kappa = 0.3, 0.02, 0.25, 0.15
lam, chi, fuel_half = 0.4, 0.5, 0.2


def npp(x):
    return P * np.array([1.0, np.exp(-kappa*x[0])]) * (-np.expm1(-a*x[[0, 5]]))


def hazard(x):
    fuel = chi*x[8]
    return lam*fuel/(fuel_half+fuel)


def matrix_input(x, fire):
    M, u = np.zeros((10, 10)), np.zeros(10)
    mu = npp(x)
    for i in range(2):
        o = 5*i
        for j in range(3):
            M[o+j, o+j] = -(d[i, j]+eta[i, j]*fire)
            M[o+3, o+j] = d[i, j]+(1-b[i, j])*eta[i, j]*fire
            u[o+j] = p[i, j]*mu[i]
        M[o+3, o+3] = -(kD+c[i]*fire)
        M[o+4, o+3] = rho*kD
        M[o+4, o+4] = -kH
    return M, u


def flux_rhs(x, fire):
    # Construct donor/recipient fluxes independently of matrix assembly.
    z, out, mu = x.reshape(2, 5), np.zeros((2, 5)), npp(x)
    rh = burning = 0.0
    for i in range(2):
        turnover = d[i]*z[i, :3]
        damage = eta[i]*fire*z[i, :3]
        out[i, :3] = p[i]*mu[i]-turnover-damage
        out[i, 3] = sum(turnover+(1-b[i])*damage)-kD*z[i, 3]-c[i]*fire*z[i, 3]
        out[i, 4] = rho*kD*z[i, 3]-kH*z[i, 4]
        rh += (1-rho)*kD*z[i, 3]+kH*z[i, 4]
        burning += sum(b[i]*damage)+c[i]*fire*z[i, 3]
    return out.ravel(), rh, burning


def equilibrium_fixed_fire(fire):
    x = np.zeros(10)
    for i in range(2):
        q = d[i]+eta[i]*fire
        effective = P[i]*np.exp(-kappa*x[0]) if i else P[i]
        assert p[i, 0]*effective*a[i] > q[0]
        def per_leaf(L):
            return p[i, 0]*effective*(-np.expm1(-a[i]*L))/L-q[0]
        L = brentq(per_leaf, 1e-10, p[i, 0]*effective/q[0], xtol=1e-14)
        mu = effective*(-np.expm1(-a[i]*L))
        x[5*i:5*i+3] = p[i]*mu/q
        transfer = d[i]+(1-b[i])*eta[i]*fire
        x[5*i+3] = sum(transfer*x[5*i:5*i+3])/(kD+c[i]*fire)
        x[5*i+4] = rho*kD*x[5*i+3]/kH
    return x


def jacobian(fun, x, step):
    eye = np.eye(10)
    return np.column_stack([(fun(x+step*e)-fun(x-step*e))/(2*step) for e in eye])


def main():
    x = np.linspace(0.1, 1.0, 10)
    f = hazard(x)
    M, u = matrix_input(x, f)
    direct, rh, burning = flux_rhs(x, f)
    matrix_error = np.max(np.abs(M@x+u-direct))
    budget_error = abs(sum(direct)-(sum(npp(x))-rh-burning))
    cap = np.linalg.solve(-M, u)
    cap_error = np.max(np.abs(direct+M@(cap-x)))
    assert max(matrix_error, budget_error, cap_error) < 1e-12
    assert np.all(cap >= 0)
    assert np.max(np.linalg.eigvals(M).real) < 0
    assert np.all(M-np.diag(np.diag(M)) >= 0)
    assert np.max(M.sum(axis=0)) <= 1e-14
    boundary_min = np.inf
    for j in range(10):
        z = x.copy()
        z[j] = 0
        boundary_min = min(boundary_min, flux_rhs(z, hazard(z))[0][j])
    assert boundary_min >= 0
    eq_results = []
    for fixed_fire in [0.0, 0.2]:
        eq = equilibrium_fixed_fire(fixed_fire)
        fun = lambda z: flux_rhs(z, fixed_fire)[0]
        J = jacobian(fun, eq, 1e-5)
        residual = np.max(np.abs(fun(eq)))
        spectral_abscissa = np.max(np.linalg.eigvals(J).real)
        assert residual < 1e-12 and spectral_abscissa < 0
        eq_results.append(dict(fixed_fire=fixed_fire, state=eq.tolist(),
                               residual=float(residual), max_real_eigenvalue=float(spectral_abscissa)))
    # Rank-one identity at an arbitrary positive state (not an equilibrium).
    J0 = jacobian(lambda z: flux_rhs(z, f)[0], x, 1e-5)
    v = (flux_rhs(x, f+1e-5)[0]-flux_rhs(x, f-1e-5)[0])/2e-5
    gamma = lam*chi*fuel_half/(fuel_half+chi*x[8])**2
    expected = J0.copy()
    expected[:, 8] += gamma*v
    fun = lambda z: flux_rhs(z, hazard(z))[0]
    J = jacobian(fun, x, 1e-5)
    Jfine = jacobian(fun, x, 5e-6)
    rank_error = np.max(np.abs(J-expected))
    jac_refinement = np.max(np.abs(J-Jfine))
    assert rank_error < 1e-8 and jac_refinement < 1e-8
    assert J[0, 8] < 0 and J[5, 0] < 0  # Both coupling directions.
    # Integrate net atmosphere loss as an extra diagnostic, not a carbon state.
    def augmented(t, y):
        dx, rh, burning = flux_rhs(y[:10], hazard(y[:10]))
        return np.r_[dx, rh+burning-sum(npp(y[:10]))]
    runs = []
    for max_step in [0.25, 0.125]:
        sol = solve_ivp(augmented, (0, 20), np.r_[x, 0.0],
                        rtol=1e-10, atol=1e-12, max_step=max_step, dense_output=True)
        assert sol.success
        samples = sol.sol(np.linspace(0, 20, 401))
        minimum = min(sol.y[:10].min(), samples[:10].min())
        err = np.max(np.abs(samples[:10].sum(axis=0)+samples[10]-sum(x)))
        assert minimum >= 0 and err < 1e-9
        runs.append((samples, float(err), float(minimum)))
    refinement = np.max(np.abs(runs[0][0]-runs[1][0]))
    assert refinement < 1e-7
    print(json.dumps(dict(
        scope='Uncalibrated algebra/structure checks; no climate or tipping result',
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        versions=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__),
        parameters=dict(P=P.tolist(), a=a.tolist(), allocation=p.tolist(), turnover=d.tolist(),
                        eta=eta.tolist(), b=b.tolist(), c=c.tolist(), kD=kD, kH=kH,
                        rho=rho, kappa=kappa, lam=lam, chi=chi, fuel_half=fuel_half),
        matrix_error=float(matrix_error), budget_error=float(budget_error),
        capacity_identity_error=float(cap_error), boundary_min=float(boundary_min),
        fixed_fire_equilibria=eq_results, rank_one_error=float(rank_error),
        jacobian_refinement=float(jac_refinement), integration_refinement=float(refinement),
        integrated_budget_error=max(r[1] for r in runs),
        minimum_sampled_carbon=min(r[2] for r in runs)), indent=2))


if __name__ == '__main__':
    main()
