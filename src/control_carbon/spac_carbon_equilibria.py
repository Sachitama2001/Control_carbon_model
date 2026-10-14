"""Leaf-carbon elimination of existing linear Manzoni E0/E1/E4 equilibria.

Canonical D=curve(L); public VPD outputs remain MPa through PressureSPAC.
The global shape classification is restricted to checked unsafe baseline
hypotheses. No zero-carbon full vector field or mortality regime is added.
"""

from decimal import Decimal

import numpy as np
from scipy.optimize import brentq
from scipy.special import exprel

from .spac_carbon import positive_carbon_equilibrium
from .spac_carbon_two_way import feedback_factors, effective_conductance
from .spac_carbon_capacity import diagnose_capacity
from .spac_carbon_stability import _Interval


def stock_ratios(model):
    c = model.carbon
    live = np.array(c.allocation)*c.losses[0]/(c.allocation[0]*np.array(c.losses[:3]))
    soil = np.dot(np.array(c.litter_fractions)*np.array(c.losses[:3]), live)/c.losses[3]
    return np.r_[live, soil]


def leaf_opening(leaf, model):
    """h(L), h'(L); the L=0 value is a limit, not a selected full root."""
    if not np.isfinite(leaf) or leaf < 0:
        raise ValueError("leaf carbon must be finite nonnegative")
    c = model.carbon
    r = c.losses[0]/(c.allocation[0]*c.assimilation_max)
    z = c.canopy_coefficient*leaf
    opening = r/c.canopy_coefficient/exprel(-z)
    if z < 1e-3:
        # Bernoulli expansion avoids cancellation near the excluded L=0 limit.
        derivative = r*(.5+z/6-z**3/180+z**5/5040)
    else:
        canopy = -np.expm1(-z)
        derivative = r*(canopy-z*np.exp(-z))/canopy**2
    return float(opening), float(derivative)


def positive_domain(model):
    h = model.hydraulic
    ax, bs = 2*h.psi50_x, 2*h.psi50_s
    h0 = leaf_opening(0., model)[0]
    hx, hr = 1-ax/bs, 1-model.root_potential_mpa/bs
    if model.root_potential_mpa <= ax or hr <= max(h0, hx):
        raise ValueError("no conducting positive-carbon equilibrium domain for these boundaries")
    high = positive_carbon_equilibrium(hr, model.carbon)
    if high is None:
        raise RuntimeError("positive-domain endpoint unexpectedly lost its carbon root")
    low = 0.
    if hx > h0:
        low = brentq(lambda leaf: leaf_opening(leaf, model)[0]-hx, 0., high[0], xtol=1e-14)
    return {"leaf_lower_open": float(low), "leaf_upper_closed_at_zero_vpd": float(high[0]),
            "h_zero_leaf_limit": h0, "h_xylem_boundary": hx, "h_root_boundary": hr,
            "leaf_exit_boundary": model.feedback.minimum_leaf}


def leaf_curve(leaf, model):
    """Reconstruct all six states plus MPa VPD and its derivative per leaf carbon."""
    domain = positive_domain(model)
    if not np.isfinite(leaf) or not domain["leaf_lower_open"] < leaf <= domain["leaf_upper_closed_at_zero_vpd"]:
        raise ValueError("leaf must lie in the conducting positive equilibrium domain")
    h, dh = leaf_opening(leaf, model)
    hydraulic, f = model.hydraulic, model.feedback
    bs, ax = 2*hydraulic.psi50_s, 2*hydraulic.psi50_x
    psi, dpsi = bs*(1-h), -bs*dh
    if leaf == domain["leaf_upper_closed_at_zero_vpd"]:
        psi, h = model.root_potential_mpa, domain["h_root_boundary"]
    stocks = stock_ratios(model)*leaf
    fl, fe, fw, dfl, dfe, dfw = feedback_factors(leaf, stocks[1], f)
    gx, gs = hydraulic.gx_max*(1-psi/ax), hydraulic.gs_max*h
    dgx, dgs = -hydraulic.gx_max*dpsi/ax, hydraulic.gs_max*dh
    gr, gl = 2*gx*fw, 2*gx*fl
    dgr = 2*(dgx*fw+gx*dfw*stock_ratios(model)[1])
    dgl = 2*(dgx*fl+gx*dfl)
    geff = effective_conductance(gr, gl)
    dgeff = (gl*gl*dgr+gr*gr*dgl)/(gr+gl)**2
    drop = model.root_potential_mpa-psi
    base = gx*drop/gs
    dbase = (dgx*drop-gx*dpsi)/gs-base*dgs/gs
    if model.root_potential_mpa == 0. and ax == -3. and bs == -4.:
        factor = 4*hydraulic.gx_max/(3*hydraulic.gs_max)
        base = factor*(5-4*h-1/h)
        dbase = factor*(1/h**2-4)*dh
    strengths = (f.leaf_supply, f.leaf_demand, f.wood_supply)
    if strengths == (0., 0., 0.):
        d, derivative = base, dbase
    elif strengths == (1., 0., 0.):
        denominator = leaf+f.leaf_reference+2*f.leaf_floor
        d = 2*base*(leaf+f.leaf_floor)/denominator
        derivative = 2*(dbase*(leaf+f.leaf_floor)/denominator
                        +base*(f.leaf_reference+f.leaf_floor)/denominator**2)
    elif strengths == (1., 1., 1.) and f.wood_exponent == 1.:
        kappa = stock_ratios(model)[1]*f.leaf_reference/f.wood_reference
        rho = kappa*(f.leaf_reference+f.leaf_floor)/f.leaf_reference
        denominator = (1+rho)*leaf+f.leaf_floor
        supply = (leaf+f.leaf_floor)/denominator
        # Cancel FW/FE algebraically, not as two divergent O(1/L) derivatives.
        d = 2*kappa*base*supply
        derivative = 2*kappa*(dbase*supply-base*rho*f.leaf_floor/denominator**2)
    else:
        d = geff*drop/(gs*fe)
        derivative = ((dgeff*drop-geff*dpsi)/(gs*fe)
                      -d*(dgs/gs+dfe/fe))
    x = (gr*model.root_potential_mpa+gl*psi)/(gr+gl)
    return np.r_[x, psi, stocks], model.pressure(d), derivative*model.atmospheric_pressure_mpa


def baseline_shape_proof(model):
    """Check the hypotheses and outward-rounded margins used by the shape proof."""
    h, c, f = model.hydraulic, model.carbon, model.feedback
    strengths = (f.leaf_supply, f.leaf_demand, f.wood_supply)
    if strengths not in [(0., 0., 0.), (1., 0., 0.), (1., 1., 1.)]:
        raise ValueError("shape proof only covers E0, E1 and full E4")
    if (2*h.psi50_x != -3 or 2*h.psi50_s != -4 or model.root_potential_mpa != 0
            or c.canopy_coefficient != 5 or f.wood_exponent != 1):
        raise ValueError("shape proof requires unsafe thresholds -3/-4, root0, q5 and beta1")
    r = _Interval(c.losses[0])/(_Interval(c.allocation[0])*_Interval(c.assimilation_max))
    h0 = r/_Interval(c.canopy_coefficient)
    if not Decimal(".4444") < h0.lo <= h0.hi < Decimal(".4446"):
        raise ValueError("shape proof h(0+) assumption failed")
    label = {(0., 0., 0.): "E0", (1., 0., 0.): "E1", (1., 1., 1.): "E4"}[strengths]
    evidence = {"model": label, "h_zero_leaf_interval": h0.strings(),
                "positive_shape": "single_maximum" if label != "E4" else "strictly_decreasing",
                "claim": "global_on_positive_conducting_leaf_domain_under_checked_hypotheses"}
    if label == "E4":
        lf, lr, wr = _Interval(f.leaf_floor), _Interval(f.leaf_reference), _Interval(f.wood_reference)
        ratio_w = (_Interval(c.allocation[1])*_Interval(c.losses[0])
                   /(_Interval(c.allocation[0])*_Interval(c.losses[1])))
        kappa = ratio_w*lr/wr
        rho = kappa*(lr+lf)/lr
        a = _Interval(".05")
        decay = (-_Interval(5)*a).exp()
        opening_a = r*a/(1-decay)
        derivative_a = r*(1-decay-_Interval(5)*a*decay)/((1-decay)*(1-decay))
        # h' is increasing; log(D0)_h is decreasing while h<1/2.
        log_base_upper = derivative_a*((1/(h0*h0)-4)/(5-4*h0-1/h0))
        supply_lower = rho*lf/((a+lf)*((1+rho)*a+lf))
        if not (opening_a.lo > Decimal(".5") and derivative_a.hi < Decimal("1.205")
                and log_base_upper.hi < supply_lower.lo):
            raise ValueError("E4 strict monotonicity margin did not certify")
        evidence.update({
            "cut_leaf_carbon": ".05", "opening_at_cut_interval": opening_a.strings(),
            "opening_derivative_at_cut_interval": derivative_a.strings(),
            "kappa_interval": kappa.strings(), "rho_interval": rho.strings(),
            "positive_log_D0_derivative_upper_interval": log_base_upper.strings(),
            "negative_log_supply_derivative_lower_interval": supply_lower.strings(),
            "strict_gap_interval": (supply_lower-log_base_upper).strings(),
            "evidence": "analytic_shape_proof_plus_outward_Decimal_constant_margin_not_stability_certificate",
        })
    return evidence


def baseline_landmarks(model):
    proof = baseline_shape_proof(model)
    domain = positive_domain(model)
    top = domain["leaf_upper_closed_at_zero_vpd"]
    h0 = domain["h_zero_leaf_limit"]
    base0 = 4*model.hydraulic.gx_max/(3*model.hydraulic.gs_max)*(5-4*h0-1/h0)
    label, f = proof["model"], model.feedback
    if label == "E0":
        zero_limit = base0
    elif label == "E1":
        zero_limit = 2*base0*f.leaf_floor/(f.leaf_reference+2*f.leaf_floor)
    else:
        kappa = stock_ratios(model)[1]*f.leaf_reference/f.wood_reference
        zero_limit = 2*kappa*base0
    fold = None
    if label in ["E0", "E1"]:
        hydraulic_fold_leaf = brentq(lambda leaf: leaf_opening(leaf, model)[0]-.5, 0., top, xtol=1e-14)
        fold = (hydraulic_fold_leaf if label == "E0" else
                brentq(lambda leaf: leaf_curve(leaf, model)[2], hydraulic_fold_leaf, top, xtol=1e-14))
    return {"domain": domain, "shape_proof": proof, "zero_leaf_vpd_limit_mpa": model.pressure(zero_limit),
            "fold_leaf": fold, "fold_vpd_mpa": None if fold is None else leaf_curve(fold, model)[1]}


def positive_equilibria(vpd_mpa, model, scales):
    """All positive conducting roots, using proved monotone intervals, not a scan."""
    model.mole_fraction(vpd_mpa)
    landmarks = baseline_landmarks(model)
    top, fold = landmarks["domain"]["leaf_upper_closed_at_zero_vpd"], landmarks["fold_leaf"]
    nodes = [0., top] if fold is None else [0., fold, top]
    def curve(leaf):
        return landmarks["zero_leaf_vpd_limit_mpa"] if leaf == 0 else leaf_curve(leaf, model)[1]
    values = [curve(leaf)-vpd_mpa for leaf in nodes]
    roots = {leaf for leaf, value in zip(nodes, values) if value == 0 and leaf > 0}
    for left, right, fl, fr in zip(nodes[:-1], nodes[1:], values[:-1], values[1:]):
        if fl*fr < 0:
            roots.add(brentq(lambda leaf: curve(leaf)-vpd_mpa, left, right, xtol=np.finfo(float).tiny))
    records = []
    for leaf in sorted(roots, reverse=True):
        y, v, slope = leaf_curve(leaf, model)
        data = diagnose_capacity(y, vpd_mpa, model, scales, equilibrium=y)
        if data["full_residual_norm_scaled_per_day"] > 1e-8 or abs(v-vpd_mpa) > 1e-11:
            raise RuntimeError("scalar equilibrium reconstruction failed full equations")
        records.append({"state": y, "vpd_slope_per_leaf": slope, "diagnostics": data})
    return {"roots": tuple(records), "landmarks": landmarks,
            "zero_leaf_limit_hit": bool(values[0] == 0),
            "status": "positive_roots_reconstructed" if records else "no_positive_conducting_root",
            "excluded_sets": "formal_zero_carbon_nonconducting_continua_and_kinks_not_in_positive_elimination"}
