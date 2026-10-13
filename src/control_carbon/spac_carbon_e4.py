"""MPa-facing nonlinear analysis of the exact E4 four-state active core."""

from itertools import product
import math

import numpy as np
from scipy.integrate import solve_bvp, solve_ivp
from scipy.linalg import expm
from scipy.optimize import brentq

from .spac_carbon_bifurcation import first_lyapunov
from .spac_carbon import assimilation, hydraulic_factor


def _full(state, model):
    y = np.asarray(state, dtype=float)
    if y.shape != (4,) or not np.all(np.isfinite(y)) or y[3] <= 0:
        raise ValueError("E4 core requires four finite entries and positive wood")
    return np.r_[y, 1., 1.]


def core_rhs(state, vpd_mpa, model):
    return model.rhs(_full(state, model), vpd_mpa)[:4]


def core_jacobian(state, vpd_mpa, model):
    return model.jacobian(_full(state, model), vpd_mpa)[:4, :4]


def core_tensors(state, vpd_mpa, model):
    """Exact D2/D3 tensors using the finite Leibniz rule, no factorials."""
    core_jacobian(state, vpd_mpa, model)
    x, psi, leaf, wood = state
    h, c, f = model.hydraulic, model.carbon, model.feedback
    if psi <= max(2*h.psi50_x, 2*h.psi50_s):
        raise ValueError("E4 derivatives require the smooth conducting domain")
    gx, gs = h.gx_max*(1-psi/(2*h.psi50_x)), h.gs_max*(1-psi/(2*h.psi50_s))
    dgx, dgs = -h.gx_max/(2*h.psi50_x), -h.gs_max/(2*h.psi50_s)

    def affine(value, slope):
        return lambda indices: value if not indices else slope[indices[0]] if len(indices) == 1 else 0.

    gx_factor = affine(gx, [0., dgx, 0., 0.])
    gs_factor = affine(gs, [0., dgs, 0., 0.])
    fl = affine(1+f.leaf_supply*((leaf+f.leaf_floor)/(f.leaf_reference+f.leaf_floor)-1),
                [0., 0., f.leaf_supply/(f.leaf_reference+f.leaf_floor), 0.])
    fe = affine(1+f.leaf_demand*(leaf/f.leaf_reference-1),
                [0., 0., f.leaf_demand/f.leaf_reference, 0.])
    def fw(indices):
        if not indices:
            return 1+f.wood_supply*((wood/f.wood_reference)**f.wood_exponent-1)
        if any(i != 3 for i in indices):
            return 0.
        order = len(indices)
        coefficient = math.prod(f.wood_exponent-j for j in range(order))
        return f.wood_supply*coefficient*(wood/f.wood_reference)**(f.wood_exponent-order)/f.wood_reference**order
    q = c.canopy_coefficient
    def canopy(indices):
        if not indices:
            return -np.expm1(-q*leaf)
        return (-1)**(len(indices)+1)*q**len(indices)*np.exp(-q*leaf) if all(i == 2 for i in indices) else 0.
    erx = [gx_factor, fw, affine(model.root_potential_mpa-x, [-1., 0., 0., 0.])]
    exl = [gx_factor, fl, affine(x-psi, [1., -1., 0., 0.])]
    atmosphere = [gs_factor, fe]
    gpp = [affine(gs/h.gs_max, [0., dgs/h.gs_max, 0., 0.]), canopy]
    def derivative(factors, indices):
        value = 0.
        for assignments in product(range(len(factors)), repeat=len(indices)):
            value += math.prod(factor(tuple(i for i, a in zip(indices, assignments) if a == j))
                               for j, factor in enumerate(factors))
        return value
    tensors = []
    d = model.mole_fraction(vpd_mpa)
    for order in [2, 3]:
        tensor = np.empty((4,)*(order+1))
        for indices in product(range(4), repeat=order):
            rx = 2*derivative(erx, indices)
            xl = 2*derivative(exl, indices)
            e = d*derivative(atmosphere, indices)
            p = c.assimilation_max*derivative(gpp, indices)
            tensor[(slice(None), *indices)] = [(rx-xl)/h.cx, (xl-e)/h.cl,
                                               c.allocation[0]*p, c.allocation[1]*p]
        tensors.append(tensor)
    return tuple(tensors)


def hopf_point(leaf_bounds, model):
    def pair_real(leaf):
        state, vpd, _ = model.equilibrium_curve(leaf)
        eig = np.linalg.eigvals(core_jacobian(state[:4], vpd, model))
        return float(eig[np.argmax(eig.imag)].real)
    for z in leaf_bounds:
        state, vpd, _ = model.equilibrium_curve(z)
        eig = np.linalg.eigvals(core_jacobian(state[:4], vpd, model))
        if max(eig.imag) <= 1e-8:
            raise ValueError("E4 Hopf bracket must stay in the complex-pair neighborhood")
    leaf = brentq(pair_real, *leaf_bounds, xtol=1e-13)
    state, vpd, slope = model.equilibrium_curve(leaf)
    normal = first_lyapunov(core_jacobian(state[:4], vpd, model),
                           *core_tensors(state[:4], vpd, model))
    slopes = [(pair_real(leaf+step)-pair_real(leaf-step))/(2*step*slope) for step in [1e-5, 5e-6]]
    if abs(slopes[-1]) < 1e-6 or abs(normal["l1"]) < 1e-8:
        raise RuntimeError("E4 Hopf nondegeneracy unresolved")
    if abs(slopes[0]-slopes[1]) > 1e-4*abs(slopes[-1]):
        raise RuntimeError("E4 Hopf slope refinement failed")
    return {"state": state, "vpd_mpa": vpd, "slope_refinement_per_mpa": slopes,
            "curve_slope": slope, **normal}


def periodic_orbit(hopf, amplitude, model, *, leaf_bounds=(-2.12, -1.9), tol=1e-8):
    """Local periodic collocation; phase/amplitude constraints exclude equilibria.

    Independent of B2 shooting. Four active states retained, including slow W.
    Returns full periodic states and Floquet spectrum; R/S are periodic filters.
    """
    if not np.isfinite(amplitude) or amplitude <= 0 or not 0 < tol < 1:
        raise ValueError("positive finite amplitude and tolerance required")
    q, p = hopf["q"], hopf["p"]
    scale = np.array([1., 1., .1, 3.])
    period0 = 2*np.pi/hopf["omega"]
    shift = -hopf["omega"]*hopf["l1"]*amplitude**2/hopf["slope_refinement_per_mpa"][-1]
    center = hopf["state"][:4]
    mesh = np.linspace(0., 1., 161)
    guess = 2*amplitude*np.real(q[:, None]*np.exp(2j*np.pi*mesh))/scale[:, None]

    def equilibrium(vpd):
        leaf = brentq(lambda z: model.equilibrium_curve(z)[1]-vpd, *leaf_bounds, xtol=1e-13)
        return model.equilibrium_curve(leaf)[0][:4]
    def rhs(t, z, parameters):
        period, vpd = period0*parameters[0], hopf["vpd_mpa"]+shift*parameters[1]
        return np.column_stack([period*core_rhs(center+scale*y, vpd, model)/scale for y in z.T])
    def boundary(left, right, parameters):
        vpd = hopf["vpd_mpa"]+shift*parameters[1]
        projection = np.vdot(p, center+scale*left-equilibrium(vpd))
        return np.r_[(left-right)/amplitude, (projection.real-amplitude)/amplitude, projection.imag/amplitude]
    solution = solve_bvp(rhs, boundary, mesh, guess, p=[1., 1.], tol=tol, bc_tol=tol, max_nodes=5000)
    if not solution.success:
        raise RuntimeError(f"E4 periodic collocation failed: {solution.message}")
    period, vpd = period0*solution.p[0], hopf["vpd_mpa"]+shift*solution.p[1]
    if period <= 0 or vpd <= 0:
        raise RuntimeError("E4 periodic solution has invalid period/forcing")
    phase = np.linspace(0., 1., 401)
    states = center+solution.sol(phase).T*scale
    if np.ptp(states[:, 1]) < amplitude or np.any(states[:, :2] > 0) or np.any(states[:, 2:] <= 0):
        raise RuntimeError("E4 collocation failed nontriviality/physical domain")
    def variation(t, phi):
        y = center+solution.sol(t/period)*scale
        return (core_jacobian(y, vpd, model)@phi.reshape(4, 4)).ravel()
    flow = solve_ivp(variation, (0., period), np.eye(4).ravel(), method="Radau",
                     rtol=tol*.01, atol=tol*.0001, max_step=period/160)
    if not flow.success:
        raise RuntimeError(flow.message)
    multipliers = np.linalg.eigvals(flow.y[:, -1].reshape(4, 4))
    if min(abs(multipliers-1)) > 1e-6:
        raise RuntimeError("E4 orbit lacks neutral Floquet multiplier")
    k = np.array(model.carbon.losses)
    eta = np.array(model.carbon.litter_fractions)
    passive_matrix = np.array([[-k[2], 0.], [eta[2]*k[2], -k[3]]])
    def passive_rhs(t, passive):
        active = center+solution.sol(t/period)*scale
        gpp = assimilation(active[2], hydraulic_factor(active[1], model.hydraulic), model.carbon)
        return passive_matrix@passive+np.array([
            model.carbon.allocation[2]*gpp, eta[0]*k[0]*active[2]+eta[1]*k[1]*active[3],
        ])
    passive_flow = solve_ivp(passive_rhs, (0., period), np.zeros(2), method="Radau",
                             rtol=tol*.01, atol=tol*.0001, dense_output=True, max_step=period/160)
    if not passive_flow.success:
        raise RuntimeError(passive_flow.message)
    passive0 = np.linalg.solve(np.eye(2)-expm(passive_matrix*period), passive_flow.y[:, -1])
    passive = np.array([expm(passive_matrix*t)@passive0+passive_flow.sol(t) for t in phase*period])
    full_states = np.column_stack([states, passive])
    if np.any(full_states[:, 2:] <= 0) or not np.all(np.isfinite(full_states)):
        raise RuntimeError("E4 passive periodic lift failed positivity")
    # Independent IVP checks the collocated core, rather than only its BCs.
    direct = solve_ivp(lambda t, y: core_rhs(y, vpd, model), (0., period), states[0],
                       method="Radau", rtol=tol*.01, atol=tol*.0001,
                       jac=lambda t, y: core_jacobian(y, vpd, model),
                       max_step=period/160, t_eval=phase*period)
    if not direct.success:
        raise RuntimeError(direct.message)
    flow_error = float(max(abs(direct.y.T-states).ravel()))
    if flow_error > 1e-6 or max(abs(full_states[-1]-full_states[0])) > 1e-8:
        raise RuntimeError("E4 periodic solution failed independent flow/closure verification")
    return {"vpd_mpa": float(vpd), "period_days": float(period), "amplitude": amplitude,
            "boundary_residual": float(max(abs(boundary(solution.y[:, 0], solution.y[:, -1], solution.p)))),
            "max_collocation_residual": float(max(solution.rms_residuals)), "nodes": len(solution.x),
            "core_floquet": multipliers,
            "full_floquet": np.r_[multipliers, np.exp(-np.array(model.carbon.losses[2:])*period)],
            "independent_flow_max_state_error": flow_error,
            "full_state_closure": float(max(abs(full_states[-1]-full_states[0]))),
            "phase": phase, "state": full_states}
