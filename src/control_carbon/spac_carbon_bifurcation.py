"""Local B2 analysis of the exact three-state E0--E3 dynamical core.

No quasi-steady water approximation: W/R/S are passive only when wood
feedback is zero. Linear active conductances and fixed B1 parameters.
"""

from itertools import permutations

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm
from scipy.optimize import brentq, least_squares

from .manzoni_hydraulics import conductances
from .spac_carbon import assimilation
from .spac_carbon_two_way import (
    feedback_factors, two_way_rhs, two_way_jacobian, equilibrium_curve,
)


def _full_state(state, feedback):
    state = np.asarray(state, dtype=float)
    if feedback.wood_supply != 0:
        raise ValueError("three-state core requires zero wood feedback")
    if state.shape != (3,) or not np.all(np.isfinite(state)):
        raise ValueError("core state must have three finite entries")
    return np.r_[state, feedback.wood_reference, 1., 1.]


def core_rhs(state, root_potential, vpd, hydraulic, carbon, feedback):
    return two_way_rhs(_full_state(state, feedback), root_potential, vpd, hydraulic, carbon, feedback)[:3]


def core_jacobian(state, root_potential, vpd, hydraulic, carbon, feedback):
    return two_way_jacobian(_full_state(state, feedback), root_potential, vpd, hydraulic, carbon, feedback)[:3, :3]


def cubic_hurwitz(matrix):
    """Coefficients of det(lambda I-A), and a1*a2-a3."""
    matrix = np.asarray(matrix, dtype=float)
    if matrix.shape != (3, 3) or not np.all(np.isfinite(matrix)):
        raise ValueError("finite three-by-three matrix required")
    a1 = -np.trace(matrix)
    a2 = sum(matrix[i, i]*matrix[j, j]-matrix[i, j]*matrix[j, i]
             for i, j in [(0, 1), (0, 2), (1, 2)])
    a3 = -np.linalg.det(matrix)
    return np.array([a1, a2, a3, a1*a2-a3])


def core_tensors(state, root_potential, vpd, hydraulic, carbon, feedback):
    """Exact D^2 f, D^3 f (no factorials) on the open conducting domain."""
    core_jacobian(state, root_potential, vpd, hydraulic, carbon, feedback)
    x, psi, leaf = state
    if psi <= max(2*hydraulic.psi50_x, 2*hydraulic.psi50_s):
        raise ValueError("higher derivatives require both conductances strictly positive")
    gx = hydraulic.gx_max*(1-psi/(2*hydraulic.psi50_x))
    dgx = -hydraulic.gx_max/(2*hydraulic.psi50_x)
    dgs = -hydraulic.gs_max/(2*hydraulic.psi50_s)
    fl, fe, fw, dfl, dfe, dfw = feedback_factors(leaf, feedback.wood_reference, feedback)
    b = np.zeros((3, 3, 3))
    c = np.zeros((3, 3, 3, 3))
    def second(out, i, j, value):
        b[out, i, j] = b[out, j, i] = value
    def third(out, indices, value):
        for perm in set(permutations(indices)):
            c[(out, *perm)] = value
    second(0, 0, 1, -2*dgx*(1+fl)/hydraulic.cx)
    second(0, 0, 2, -2*gx*dfl/hydraulic.cx)
    second(0, 1, 1, 4*dgx*fl/hydraulic.cx)
    second(0, 1, 2, (-2*dgx*dfl*(x-psi)+2*gx*dfl)/hydraulic.cx)
    second(1, 0, 1, 2*dgx*fl/hydraulic.cl)
    second(1, 0, 2, 2*gx*dfl/hydraulic.cl)
    second(1, 1, 1, -4*dgx*fl/hydraulic.cl)
    second(1, 1, 2, (2*dgx*dfl*(x-psi)-2*gx*dfl-dgs*dfe*vpd)/hydraulic.cl)
    third(0, (0, 1, 2), -2*dgx*dfl/hydraulic.cx)
    third(0, (1, 1, 2), 4*dgx*dfl/hydraulic.cx)
    third(1, (0, 1, 2), 2*dgx*dfl/hydraulic.cl)
    third(1, (1, 1, 2), -4*dgx*dfl/hydraulic.cl)
    alpha = carbon.allocation[0]*carbon.assimilation_max
    q = carbon.canopy_coefficient
    h = 1-psi/(2*hydraulic.psi50_s)
    dh = -1/(2*hydraulic.psi50_s)
    decay = np.exp(-q*leaf)
    second(2, 1, 2, alpha*dh*q*decay)
    second(2, 2, 2, -alpha*h*q*q*decay)
    third(2, (1, 2, 2), -alpha*dh*q*q*decay)
    third(2, (2, 2, 2), alpha*h*q**3*decay)
    return b, c


def first_lyapunov(matrix, second, third):
    """Kuznetsov l1 = Re(G21)/(2*omega), <p,q>=1, ||q||2=1.

    Original day time and the original mixed-unit core coordinates.
    Positive l1 means a locally repelling cycle on the stable-equilibrium side.
    """
    matrix = np.asarray(matrix, dtype=float)
    second, third = np.asarray(second), np.asarray(third)
    if matrix.shape != (3, 3) or second.shape != (3, 3, 3) or third.shape != (3, 3, 3, 3):
        raise ValueError("core matrix and derivative tensor shapes do not match")
    if not all(np.all(np.isfinite(a)) for a in [matrix, second, third]):
        raise ValueError("normal-form inputs must be finite")
    eigen, vectors = np.linalg.eig(matrix)
    candidates = np.flatnonzero(eigen.imag > 1e-8)
    if len(candidates) != 1:
        raise ValueError("normal form requires one simple nonzero-frequency conjugate pair")
    index = candidates[0]
    omega = float(eigen[index].imag)
    real_eigen = eigen[abs(eigen.imag) < 1e-8]
    if len(real_eigen) != 1 or real_eigen[0].real >= -1e-8*omega or abs(eigen[index].real) > 1e-7*omega:
        raise ValueError("normal form requires a critical pair with a stable remaining eigenvalue")
    q = vectors[:, index]
    q *= np.exp(-1j*np.angle(q[np.argmax(abs(q))]))
    left_eigen, left_vectors = np.linalg.eig(matrix.T)
    p = left_vectors[:, np.argmin(abs(left_eigen+1j*omega))]
    pairing = np.vdot(p, q)
    if abs(pairing) < 1e-12:
        raise ValueError("critical left/right eigenvectors are ill-defined")
    p = p/np.conj(pairing)
    def b(u, v):
        return np.einsum("ijk,j,k->i", second, u, v)
    c = np.einsum("ijkl,j,k,l->i", third, q, q, q.conj())
    h11 = np.linalg.solve(matrix, b(q, q.conj()))
    h20 = np.linalg.solve(2j*omega*np.eye(3)-matrix, b(q, q))
    g21 = np.vdot(p, c-2*b(q, h11)+b(q.conj(), h20))
    return {
        "omega": omega, "l1": float(g21.real/(2*omega)), "G21": g21,
        "q": q, "p": p, "h11_solve": h11, "h20_solve": h20,
        "eigenvalues": eigen,
    }


def hopf_point(bounds, root_potential, hydraulic, carbon, feedback):
    def determinant(psi):
        y, d, _ = equilibrium_curve(psi, root_potential, hydraulic, carbon, feedback)
        return cubic_hurwitz(core_jacobian(y[:3], root_potential, d, hydraulic, carbon, feedback))[-1]
    psi = brentq(determinant, *bounds, xtol=1e-13)
    y, d, slope = equilibrium_curve(psi, root_potential, hydraulic, carbon, feedback)
    a = core_jacobian(y[:3], root_potential, d, hydraulic, carbon, feedback)
    coefficients = cubic_hurwitz(a)
    if min(coefficients[:3]) <= 0 or abs(slope) < 1e-8:
        raise RuntimeError("Hurwitz crossing is not a nonzero-frequency Hopf candidate")
    normal = first_lyapunov(a, *core_tensors(y[:3], root_potential, d, hydraulic, carbon, feedback))
    derivatives = []
    for step in [1e-5, 5e-6]:
        real_parts = []
        for leaf in [psi-step, psi+step]:
            state, vpd, _ = equilibrium_curve(leaf, root_potential, hydraulic, carbon, feedback)
            eig = np.linalg.eigvals(core_jacobian(state[:3], root_potential, vpd, hydraulic, carbon, feedback))
            real_parts.append(float(eig[np.argmax(eig.imag)].real))
        derivatives.append((real_parts[1]-real_parts[0])/(2*step*slope))
    if abs(derivatives[-1]) < 1e-6 or abs(normal["l1"]) < 1e-8:
        raise RuntimeError("Hopf transversality or cubic nondegeneracy not resolved")
    if abs(derivatives[0]-derivatives[1]) > 1e-4*abs(derivatives[-1]):
        raise RuntimeError("critical eigenvalue slope refinement failed")
    return {"state": y, "D": d, "curve_slope": slope, "hurwitz": coefficients,
            "eigenvalue_slope_D_refinement": derivatives, **normal}


def periodic_orbit(hopf, amplitude, root_potential, hydraulic, carbon, feedback, *,
                   leaf_bracket, rtol=1e-10, steps_per_period=100):
    """Amplitude/phase-fixed shooting, with variational Jacobian and Floquet.

    Unknowns: stable initial coordinate, period, D. Nonzero critical
    coordinate is fixed, excluding the spurious equilibrium solution.
    Passive W/R/S are lifted using their unique periodic forced solution.
    """
    if not np.isfinite(amplitude) or amplitude <= 0 or not np.isfinite(rtol) or not 0 < rtol < 1:
        raise ValueError("positive amplitude and rtol in (0,1) required")
    if steps_per_period < 10:
        raise ValueError("at least ten steps per period required")
    q, omega = hopf["q"], hopf["omega"]
    critical = hopf["state"][:3]
    a = core_jacobian(critical, root_potential, hopf["D"], hydraulic, carbon, feedback)
    eigen, vectors = np.linalg.eig(a)
    stable = vectors[:, np.argmin(eigen.real)].real
    basis = np.column_stack([2*q.real, -2*q.imag, stable])
    inverse = np.linalg.inv(basis)
    period0 = 2*np.pi/omega
    shift = -omega*hopf["l1"]*amplitude**2/hopf["eigenvalue_slope_D_refinement"][-1]
    if not np.isfinite(shift) or shift == 0:
        raise ValueError("Hopf normal form does not provide a resolvable parameter shift")
    curve_d = sorted(equilibrium_curve(psi, root_potential, hydraulic, carbon, feedback)[1]
                     for psi in leaf_bracket)
    parameter_bounds = sorted((d-hopf["D"])/shift for d in curve_d)
    last = {}

    def equilibrium(d):
        psi = brentq(lambda z: equilibrium_curve(z, root_potential, hydraulic, carbon, feedback)[1]-d,
                     *leaf_bracket, xtol=1e-13)
        return equilibrium_curve(psi, root_potential, hydraulic, carbon, feedback)[0][:3]

    def parameter_derivative(y):
        gs = conductances(y[1], hydraulic)[1]
        fe = feedback_factors(y[2], feedback.wood_reference, feedback)[1]
        return np.array([0., -gs*fe/hydraulic.cl, 0.])

    def evaluate(unknown):
        if "unknown" in last and np.array_equal(last["unknown"], unknown):
            return last
        s, period, d = unknown[0]*amplitude**2, unknown[1]*period0, hopf["D"]+unknown[2]*shift
        eq = equilibrium(d)
        y0 = eq+2*amplitude*q.real+s*stable
        kink = max(2*hydraulic.psi50_x, 2*hydraulic.psi50_s)
        if y0[1] <= kink or y0[2] <= feedback.minimum_leaf:
            raise RuntimeError("periodic shooting initial guess left the smooth viable domain")
        def rhs(t, z):
            y = z[:3]
            jac = core_jacobian(y, root_potential, d, hydraulic, carbon, feedback)
            return np.r_[core_rhs(y, root_potential, d, hydraulic, carbon, feedback),
                         (jac@z[3:12].reshape(3, 3)).ravel(),
                         jac@z[12:]+parameter_derivative(y)]
        def boundary(t, z):
            return min(z[1]-kink, z[2]-feedback.minimum_leaf)
        boundary.direction = -1
        boundary.terminal = True
        flow = solve_ivp(
            rhs, (0., period), np.r_[y0, np.eye(3).ravel(), np.zeros(3)],
            method="Radau", rtol=rtol, atol=rtol*.01, max_step=period/steps_per_period,
            dense_output=True, events=boundary,
        )
        if not flow.success:
            raise RuntimeError(flow.message)
        if flow.status == 1:
            raise RuntimeError("trial periodic flow left the smooth viable domain; reduce local amplitude")
        end = flow.y[:3, -1]
        monodromy = flow.y[3:12, -1].reshape(3, 3)
        eq_derivative = -np.linalg.solve(
            core_jacobian(eq, root_potential, d, hydraulic, carbon, feedback), parameter_derivative(eq),
        )
        derivatives = np.column_stack([
            (monodromy-np.eye(3))@stable*amplitude**2,
            core_rhs(end, root_potential, d, hydraulic, carbon, feedback)*period0,
            ((monodromy-np.eye(3))@eq_derivative+flow.y[12:, -1])*shift,
        ])
        last.clear()
        last.update(unknown=unknown.copy(), residual=inverse@(end-y0)/amplitude,
                    jacobian=inverse@derivatives/amplitude, flow=flow, D=d, period=period, initial=y0)
        return last

    solution = least_squares(
        lambda z: evaluate(z)["residual"], [0., 1., 1.],
        jac=lambda z: evaluate(z)["jacobian"],
        bounds=([-np.inf, .7, parameter_bounds[0]], [np.inf, 1.3, parameter_bounds[1]]),
        ftol=1e-11, xtol=1e-11, gtol=1e-11, max_nfev=40,
    )
    data = evaluate(solution.x)
    residual = float(max(abs(data["residual"])))
    if not solution.success or residual > 1e-7:
        raise RuntimeError(f"periodic shooting failed: {solution.message}; normalized residual={residual}")
    flow, period, d = data["flow"], data["period"], data["D"]
    time = np.linspace(0., period, 401)
    active = flow.sol(time)[:3].T
    if np.any(active[:, :2] > 0) or np.any(active[:, 2] <= feedback.minimum_leaf):
        raise RuntimeError("periodic orbit left the declared viable domain")
    if np.ptp(active[:, 1]) < amplitude:
        raise RuntimeError("shooting converged to a spurious equilibrium")
    monodromy = flow.y[3:12, -1].reshape(3, 3)
    multipliers = np.linalg.eigvals(monodromy)
    neutral = int(np.argmin(abs(multipliers-1)))
    if abs(multipliers[neutral]-1) > 1e-7:
        raise RuntimeError("periodic orbit lacks the required neutral Floquet multiplier")

    k = np.array(carbon.losses)
    eta = np.array(carbon.litter_fractions)
    passive_matrix = np.array([[-k[1], 0., 0.], [0., -k[2], 0.],
                               [eta[1]*k[1], eta[2]*k[2], -k[3]]])
    def passive_rhs(t, stocks):
        core = flow.sol(t)[:3]
        stress = conductances(core[1], hydraulic)[1]/hydraulic.gs_max
        gpp = assimilation(core[2], stress, carbon)
        return passive_matrix@stocks+np.array([
            carbon.allocation[1]*gpp, carbon.allocation[2]*gpp, eta[0]*k[0]*core[2],
        ])
    particular = solve_ivp(
        passive_rhs, (0., period), np.zeros(3), method="Radau",
        rtol=rtol, atol=rtol*.01, dense_output=True, max_step=period/steps_per_period,
    )
    if not particular.success:
        raise RuntimeError(particular.message)
    passive0 = np.linalg.solve(np.eye(3)-expm(passive_matrix*period), particular.y[:, -1])
    passive = np.array([expm(passive_matrix*t)@passive0+particular.sol(t) for t in time])
    state = np.column_stack([active, passive])
    closure = float(max(abs(state[-1]-state[0])))
    if np.any(state[:, 2:] <= 0) or not np.all(np.isfinite(state)) or closure > 1e-8:
        raise RuntimeError("full six-state periodic lift failed positivity/closure")
    full_multipliers = np.r_[multipliers, np.exp(-k[1:]*period)]
    return {
        "D": float(d), "period_days": float(period), "amplitude_coordinate": amplitude,
        "predicted_D_shift": float(shift), "shooting_residual": residual,
        "full_state_closure": closure, "solver_message": solution.message, "shooting_nfev": solution.nfev,
        "core_floquet": multipliers, "full_floquet": full_multipliers,
        "radial_floquet": multipliers[np.arange(3) != neutral][np.argmax(abs(multipliers[np.arange(3) != neutral]))],
        "time_days": time, "state": state,
    }
