"""VISIT-informed closures for the independent eight-state NSC model.

Only the explicitly documented algebra is source-transcribed.  Combining it
with the continuous NSC model is a source-grounded reduced compartment model,
not the native sequential VISIT daily map.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
import math
import numpy as np

from .temperature_nsc import Parameters, evaluate

SOURCE_REPOSITORY = "Sachitama2001/VISIT-matrix"
SOURCE_COMMIT = "3285bd8e131a932e338b59892751648fd9edcc7b"


@dataclass(frozen=True)
class VisitTreeParameters:
    """TKY tree column E of Config_TKY.xlsx, with source units."""

    alloc_ass: float = 0.1
    alloc_abg: float = 0.67
    sla: float = 150.0  # cm2 gDM-1
    extinction: float = 0.48
    lue0: float = 0.052
    pmax: float = 13.0  # umol CO2 m-2 s-1
    topt0: float = 18.0
    tmin: float = 2.0
    tmax: float = 35.0
    kmci: float = 40.0
    cmpcd0: float = 50.0
    rgf: float = 0.4
    rgc: float = 0.18
    rgr: float = 0.32
    rmf0: float = 1.3
    rmc_s: float = 0.033
    rmr_s: float = 0.525
    rmc_h: float = 0.004
    rmr_h: float = 0.16
    qtf0: float = 2.0
    qtc0: float = 2.0
    qtr0: float = 2.0
    f_size_stem: float = 40.0
    f_size_root: float = 15.0
    leaf_turnover: float = 8.2e-5
    stem_turnover: float = 7.8e-6
    root_turnover: float = 6.4e-5


@dataclass(frozen=True)
class FixedEnvironment:
    """Fixed non-temperature forcing; values are experiment assumptions."""

    ppfd_top: float = 1000.0  # umol photon m-2 s-1
    day_length: float = 12.0  # hour
    co2: float = 370.0  # ppmv
    ci_fraction: float = 0.7
    water_scalar: float = 1.0
    soil_temperature_offset: float = 0.0


def load_tky_tree_parameters(path):
    """Read the TKY tree column by labels; never silently fill missing cells."""
    import openpyxl
    book = openpyxl.load_workbook(path, data_only=True, read_only=True)
    sheet = book["Parameter-daily"]
    values = {str(row[3].value): row[4].value for row in sheet.iter_rows()
              if row[3].value is not None}
    mapping = {
        "alloc_ass": "alloc_ass", "alloc_abg": "alloc_abg", "sla": "SLA",
        "extinction": "eK0", "lue0": "lueO", "pmax": "pmax",
        "topt0": "topt0", "tmin": "tmin", "tmax": "tmax",
        "kmci": "kmci", "cmpcd0": "cmpcd0", "rgf": "rgf", "rgc": "rgc",
        "rgr": "rgr", "rmf0": "rmf0", "rmc_s": "rmc_s", "rmr_s": "rmr_s",
        "rmc_h": "rmc_h", "rmr_h": "rmr_h", "qtf0": "qTf0",
        "qtc0": "qTc0", "qtr0": "qTr0", "f_size_stem": "f_sz_s",
        "f_size_root": "f_sz_r", "leaf_turnover": "lf0",
        "stem_turnover": "lc0", "root_turnover": "lr0",
    }
    missing = [source for source in mapping.values()
               if source not in values or values[source] is None]
    if missing:
        raise ValueError(f"missing TKY tree parameters: {missing}")
    return VisitTreeParameters(**{field: float(values[source])
                                  for field, source in mapping.items()})


def q10_acclimated(q10_base, temperature):
    """ecophysiology.c::f_q10_ar, lines 303-321."""
    return q10_base * math.exp(-0.009 * (temperature - 15.0))


def leaf_area(leaf_carbon, p):
    """ecophysiology.c::lai_mass, lines 108-126; Mg C ha-1 -> LAI."""
    return max(0.0, p.sla * leaf_carbon * 2.2 / 100.0 / 2.0)


def leaf_carbon_from_lai(lai, p):
    if p.sla <= 0 or lai < 0:
        raise ValueError("positive SLA and nonnegative LAI required")
    return lai * 200.0 / (2.2 * p.sla)


def visit_physiology(temperature, p, env):
    """Fixed-ci version of source temperature/CO2 limitation.

    The native six-pass stomatal-ci iteration is deliberately not claimed here.
    Water is fixed through an explicit scalar, consistent with this experiment.
    """
    ci = env.ci_fraction * env.co2
    topt = p.topt0 + 0.01 * ci
    comp = p.cmpcd0 * max(0.0, 1 + 0.0451*(temperature-20)
                          + 0.000347*(temperature-20)**2)
    numerator = (temperature-p.tmax)*(temperature-p.tmin)
    denominator = numerator-(temperature-topt)**2
    ftemp = numerator/denominator if denominator != 0 else math.nan
    if temperature + env.soil_temperature_offset < 2:
        ftemp = 0.0
    ftemp = min(1.0, max(0.0, ftemp))
    fco2 = min(1.0, max(0.0, 0.30+0.70*(ci-comp)/(p.kmci+ci)))
    psat = p.pmax*ftemp*fco2*env.water_scalar
    lue = p.lue0*(52-temperature)/(3.5+0.75*(52-temperature))*ci/(90+0.6*ci)
    return dict(ci=ci, topt=topt, compensation=comp, ftemp=ftemp,
                fco2=fco2, psat=psat, lue=lue,
                qtf=q10_acclimated(p.qtf0, temperature),
                qtc=q10_acclimated(p.qtc0, temperature),
                qtr=q10_acclimated(p.qtr0, temperature))


def optimum_lai_source(temperature, p, env):
    """Exact valid-branch algebra of ecophysiology.c::f_opt_lai.

    Returns NaN and a reason at a source singularity rather than adding epsilon.
    """
    v = visit_physiology(temperature, p, env)
    if v["psat"] <= 0:
        return 0.0, "psat_nonpositive"
    if p.sla <= 0 or p.extinction <= 0:
        return math.nan, "invalid_sla_or_extinction"
    arm = p.rmf0/1000 * v["qtc"]**((temperature-15)/10)*2.2*10000/p.sla
    arg = p.leaf_turnover*2.2*10000/p.sla*(1+p.rgf)
    cost = arm+arg
    margin = v["psat"]*env.day_length-24*cost
    if cost <= 0 or margin <= 0:
        return math.nan, "source_denominator_nonpositive"
    cc3 = v["psat"]*env.day_length/margin
    c = v["psat"]*(cc3-1)
    if c <= 0:
        return math.nan, "source_cost_nonpositive"
    ratio = p.extinction*v["lue"]*env.ppfd_top/c
    return math.log(max(1.0, ratio))/p.extinction, "valid"


def optimum_leaf_carbon(temperature, p, env):
    lai, status = optimum_lai_source(temperature, p, env)
    return ((leaf_carbon_from_lai(lai, p) if math.isfinite(lai) else math.nan),
            lai, status)


def visit_canopy_gpp(leaf_carbon, temperature, p, env):
    """photosynthesis.c::f_gpp, lines 22-42 [Mg C ha-1 day-1]."""
    v = visit_physiology(temperature, p, env)
    if v["psat"] <= 0:
        return 0.0
    lai = leaf_area(leaf_carbon, p)
    b = p.extinction*v["lue"]*env.ppfd_top/v["psat"]
    top_root = math.sqrt(1+b)
    bottom_root = math.sqrt(1+b*math.exp(-p.extinction*lai))
    # Algebraically identical rationalized root difference. The literal native
    # ratio rounds to one near bare vegetation, falsely setting positive GPP=0.
    root_difference = b*(-math.expm1(-p.extinction*lai))/(top_root+bottom_root)
    conversion = 3600*12/100000000
    return (2*v["psat"]*env.day_length*conversion/p.extinction
            * math.log1p(root_difference/(1+bottom_root)))


def living_tissue(total, size_scale):
    """Source numeric map with an explicit 1 Mg C ha-1 reference quantity."""
    if total < 0 or size_scale <= 0:
        raise ValueError("nonnegative mass and positive size scale required")
    if total == 0:
        return 0.0
    exponent = 1-0.33334*total/(size_scale+total)
    reference = 1.0
    return min(total, reference*(total/reference)**exponent)


def size_dependent_respiration(c, temperature, p):
    sap = living_tissue(c[1], p.f_size_stem)
    fine = living_tissue(c[2], p.f_size_root)
    stem = ((p.rmc_s*sap+p.rmc_h*(c[1]-sap))/(c[1]+1e-5)
            * c[1]/1000*q10_acclimated(p.qtc0, temperature)**((temperature-15)/10))
    root = ((p.rmr_s*fine+p.rmr_h*(c[2]-fine))/(c[2]+1e-5)
            * c[2]/1000*q10_acclimated(p.qtr0, temperature)**((temperature-15)/10))
    leaf = c[0]*p.rmf0/1000*q10_acclimated(p.qtf0, temperature)**((temperature-15)/10)
    return np.array([leaf, stem, root])


def evaluate_variant(x, temperature, base, visit, env, variant="V1t",
                     reference_temperature=20.0):
    """Evaluate V0/V1f/V1t/V2/V3 with NSC mortality disabled in every variant."""
    if base.nsc_mortality or base.mu_max != 0:
        raise ValueError("VISIT allocation experiments require NSC mortality disabled")
    target_temperature = (reference_temperature if variant == "V1f" else temperature)
    if variant == "V0":
        target = base.leaf_target
        target_lai = leaf_area(target, visit)
        status = "baseline_constant"
    elif variant in ("V1f", "V1t", "V2", "V3"):
        target, target_lai, status = optimum_leaf_carbon(target_temperature, visit, env)
        if not math.isfinite(target):
            raise ValueError(f"invalid VISIT optimum LAI branch: {status}")
    else:
        raise ValueError(f"unsupported first-stage variant: {variant}")
    local = replace(base, leaf_target=target, target_mode="constant")
    dx, diagnostics = evaluate(x, temperature, local)
    old_gpp = diagnostics["gpp"]
    if variant in ("V2", "V3"):
        new_gpp = visit_canopy_gpp(x[0], temperature, visit, env)
        dx[4] += new_gpp-old_gpp
        diagnostics["gpp"] = new_gpp
        diagnostics["external"] += new_gpp-old_gpp
        diagnostics["budget"] = float(dx.sum()-diagnostics["external"])
    if variant == "V3":
        old_paid = np.asarray(diagnostics["paid"])
        new_demand = size_dependent_respiration(np.asarray(x[:3]), temperature, visit)
        new_paid = new_demand*np.asarray(x[4:7])/(np.asarray(base.substrate_half)
                                                   + np.asarray(x[4:7]))
        delta_paid = new_paid-old_paid
        dx[4:7] -= delta_paid
        diagnostics["demand"] = new_demand
        diagnostics["paid"] = new_paid
        diagnostics["unpaid"] = new_demand-new_paid
        diagnostics["external"] -= delta_paid.sum()
        diagnostics["budget"] = float(dx.sum()-diagnostics["external"])
    diagnostics.update(target_leaf_carbon=target, target_lai=target_lai,
                       actual_lai=leaf_area(x[0], visit), target_status=status,
                       variant=variant, nsc_mortality_enabled=False)
    return dx, diagnostics
