"""Selectable VISIT-informed closures for the v2 theoretical experiment.

Only the pinned old VISIT-matrix snapshot is used. Native daily scalar algebra
is separated from its normalization/mapping into an eight-state ODE.
"""
from dataclasses import asdict, dataclass
import numpy as np

REPOSITORY = "Sachitama2001/VISIT-matrix"
COMMIT = "3285bd8e131a932e338b59892751648fd9edcc7b"


def visit_temperature(t, tmin, topt, tmax):
    """f_pc_sat temperature scalar on its physiological interval (soil unfrozen).

    Outside [tmin,tmax], return zero explicitly: the literal rational source
    can become positive again far outside this interval. That is not imported.
    """
    if not tmin < topt < tmax:
        raise ValueError("require tmin < topt < tmax")
    if t <= tmin or t >= tmax:
        return 0.
    a = (t-tmax)*(t-tmin)
    return float(a/(a-(t-topt)**2))


def visit_gpp(psat, daylight_hours, extinction, lue, ppfd, lai):
    """Native f_gpp scalar algebra, Mg C/ha/day, from leaf micromol/m2/s.

    daylen is in hours; lTs=3600*12/1e8 converts the ground-area carbon unit.
    No daily updates or stomatal iterations are performed by this function.
    """
    if extinction <= 0 or min(daylight_hours, lue, ppfd, lai) < 0:
        raise ValueError("invalid canopy arguments")
    if psat <= 0:
        return 0.
    b = extinction*lue*ppfd/psat
    numerator = 1+np.sqrt(1+b)
    denominator = 1+np.sqrt(1+b*np.exp(-extinction*lai))
    return float(2*psat*daylight_hours*(3600*12/1e8)/extinction*np.log(numerator/denominator))


def visit_quadratic_supply(supply, demand, curvature=.85):
    """Native quadratic limiter, evaluated by a cancellation-free equivalent.

    Native S and E are one-day amounts [mm]. In the ODE adapter both arguments
    are fluxes [mm/day], with S=W/tau. This tau is a declared physical mapping,
    never the numerical solver timestep.
    """
    if min(supply, demand) < 0 or not 0 < curvature <= 1:
        raise ValueError("invalid supply limiter arguments")
    total = supply+demand
    if total == 0:
        return 0.
    return float(2*supply*demand/(total+np.sqrt((supply-demand)**2+4*(1-curvature)*supply*demand)))


def visit_c3_soil_scalar(water, field_capacity, half):
    """C3 non-stomatal factor only, not the complete stomatal water response."""
    if water < 0 or field_capacity <= 0 or half <= 0:
        raise ValueError("invalid soil scalar inputs")
    return .05+.95*water/(water+field_capacity*half)


@dataclass(frozen=True)
class ProcessOptions:
    photo_temperature: str = "gaussian"
    photo_canopy: str = "beer"
    photo_water: str = "plant_hill"
    evap_supply: str = "relative_hill"
    drainage: str = "excess"
    photo_tmin: float = 0.
    photo_tmax: float = 45.
    canopy_light_load: float = 4.
    soil_half: float = .3
    supply_time_days: float = 1.
    supply_curvature: float = .85
    baseflow_rate: float = .001

    def __post_init__(self):
        choices = {"photo_temperature": ("gaussian", "visit_tem"),
                   "photo_canopy": ("beer", "visit_monsi"),
                   "photo_water": ("plant_hill", "visit_soil_with_plant_gate"),
                   "evap_supply": ("relative_hill", "visit_quadratic"),
                   "drainage": ("excess", "visit_baseflow_plus_excess")}
        for key, allowed in choices.items():
            if getattr(self, key) not in allowed:
                raise ValueError(f"unknown {key}: {getattr(self, key)}; choose {allowed}")
        for key in ("photo_tmin", "photo_tmax", "canopy_light_load", "soil_half",
                    "supply_time_days", "supply_curvature", "baseflow_rate"):
            if not np.isfinite(getattr(self, key)):
                raise ValueError(f"nonfinite option: {key}")
        if not self.photo_tmin < self.photo_tmax:
            raise ValueError("invalid photosynthetic temperature interval")
        if min(self.canopy_light_load, self.soil_half, self.supply_time_days) <= 0:
            raise ValueError("positive closure scales required")
        if not 0 < self.supply_curvature <= 1 or self.baseflow_rate < 0:
            raise ValueError("invalid curvature/baseflow")


def process_provenance(options=ProcessOptions()):
    """Source order, units, mapping, and numerical-parameter status."""
    sources = {
        "visit_tem": ("visit_local/photosynthesis.c", "f_pc_sat", "78-91",
            "dimensionless; T in degC", "tmp_sfc -> air T; topt held fixed; unfrozen soil; restrict to [Tmin,Tmax]"),
        "visit_monsi": ("visit_local/photosynthesis.c", "f_gpp", "22-42",
            "native Mg C ha-1 day-1", "fixed light; psat proportional to T/W factors; normalize saturation to gpp_max; eK*LAI=k_leaf*C_L"),
        "visit_soil_with_plant_gate": ("visit_local/photosynthesis.c", "f_pc_sat", "103-113",
            "dimensionless", "C3 non-stomatal factor mapped to W_O/soil_fc; retain original plant water gate; native stomatal iteration omitted"),
        "visit_quadratic": ("visit_local/hydro_balance.c", "f_hydrology", "127-175",
            "native mm per day-step; adapted mm/day", "W/tau is supply; leaf W supplies transpiration instead of native deep-soil W; simultaneous ODE, not sequential buckets"),
        "visit_baseflow_plus_excess": ("visit_local/hydro_balance.c", "f_hydrology", "182-202,230",
            "0.001 day-1 (native non-QHB unfrozen branch)", "add one D=k*W_O term to current excess drainage; do not duplicate source baseflow debit; no soil-temperature input"),
    }
    selected = {}
    for key, value in asdict(options).items():
        if isinstance(value, str) and value in sources:
            path, function, lines, units, mapping = sources[value]
            selected[key] = {"repository": REPOSITORY, "commit": COMMIT,
                             "path": path, "function": function, "lines": lines,
                             "units": units, "adaptation": mapping,
                             "native_order": "location/ecophysiology prepares scalars; sequential daily hydrology; later daily GPP",
                             "level": "source-grounded reduced compartment model; not native VISIT"}
    return {"options": asdict(options), "sources": selected,
            "parameter_status": "Tmin/Tmax, light load, soil half, supply timescale are uncalibrated assumptions; curvature .85 and baseflow .001 follow cited source defaults",
            "unchanged_transport": "relative-storage O->R->S->L transport; native VISIT has no matching explicit three organ water stores"}
