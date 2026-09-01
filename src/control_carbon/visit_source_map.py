"""Source-grounded mapping from VISIT C code to state-space concepts.

This module does not reimplement VISIT.  It records which native variables,
functions, and source files justify each state/input/output used by the
control-theory adapter.  The C source in Sachitama2001/VISIT-matrix/visit_local
is the authority; the older Python matrix reduction is only a cross-check.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Tuple

from .provenance import SourceRef, provenance_manifest


VISIT_SOURCE_REPOSITORY = "Sachitama2001/VISIT-matrix"
VISIT_SOURCE_COMMIT = "3285bd8e131a932e338b59892751648fd9edcc7b"
VISIT_SOURCE_ROOT = "visit_local"


@dataclass(frozen=True)
class StateVariable:
    """A state variable used by a VISIT state-space representation."""

    name: str
    c_expression: str
    subsystem: str
    meaning: str
    unit: str
    status: str
    sources: Tuple[SourceRef, ...]


@dataclass(frozen=True)
class InputVariable:
    name: str
    c_expression: str
    layer: str
    meaning: str
    sources: Tuple[SourceRef, ...]


@dataclass(frozen=True)
class OutputVariable:
    name: str
    c_expression: str
    meaning: str
    sources: Tuple[SourceRef, ...]


def _visit_source(
    path: str, symbol: str, role: str, *, function: str | None = None
) -> SourceRef:
    return SourceRef(
        path,
        symbol,
        role,
        repository=VISIT_SOURCE_REPOSITORY,
        commit=VISIT_SOURCE_COMMIT,
        function=function,
    )


STRUCTURE = _visit_source(
    "visit_local/structure.h",
    "struct Pmas / struct Smas / struct Mass / struct Flux",
    "native state, parameter, and flux declarations",
)
RESTART = _visit_source(
    "visit_local/experiment.c",
    "f_experiment",
    "restart variables and daily experiment loop",
    function="f_experiment",
)
PLANT_PROCESS = _visit_source(
    "visit_local/plant_proc.c",
    "plant_process",
    "daily plant carbon updates",
    function="plant_process",
)
SOIL_PROCESS = _visit_source(
    "visit_local/soil_proc.c",
    "f_cycle_soil",
    "daily litter/humus mass balance and heterotrophic respiration",
    function="f_cycle_soil",
)
DAILY_SCHEME = _visit_source(
    "visit_local/daily_scheme.c",
    "daily_scheme",
    "sequential ecosystem update and ecosystem diagnostics",
    function="daily_scheme",
)
LOCATION_PROCESS = _visit_source(
    "visit_local/location_proc.c",
    "f_loct_proc",
    "daily meteorological/environmental forcing preparation",
    function="f_loct_proc",
)
PHOTOSYNTHESIS = _visit_source(
    "visit_local/photosynthesis.c",
    "f_gpp / f_pc_sat",
    "daily GPP and photosynthetic environmental response",
    function="f_gpp / f_pc_sat",
)
ECOPHYS = _visit_source(
    "visit_local/ecophysiology.c",
    "f_ecophysiology and helpers",
    "physiological coefficients, stomatal response, LAI relationships",
    function="f_ecophysiology",
)
RESPIRATION = _visit_source(
    "visit_local/respiration.c",
    "f_rfm / f_rcm / f_rrm / f_rfg / f_rcg / f_rrg",
    "autotrophic maintenance and growth respiration",
    function="f_rfm / f_rcm / f_rrm / f_rfg / f_rcg / f_rrg",
)
ALLOCATION = _visit_source(
    "visit_local/allocation.c",
    "f_allocation / reallocation_survival",
    "state- and regime-dependent allocation/reallocation",
    function="f_allocation / reallocation_survival",
)
LITTERFALL = _visit_source(
    "visit_local/litterfall.c",
    "f_lf / f_lc / f_lr",
    "plant-to-litter turnover",
    function="f_lf / f_lc / f_lr",
)
DECOMPOSITION = _visit_source(
    "visit_local/decomposition.c",
    "frl / frh",
    "temperature/moisture environmental scalars for decomposition",
    function="frl / frh",
)


# Structural carbon state used by the first reduced control model.
# These 18 pools are all explicit masses in struct Pmas/Smas and are explicitly
# restored in f_experiment's restart path.
VISIT_CORE_18: Tuple[StateVariable, ...] = (
    StateVariable("tree_fol", "mass->tree.fol", "plant/tree", "tree foliage carbon", "Mg C ha-1", "core", (STRUCTURE, RESTART, PLANT_PROCESS)),
    StateVariable("tree_stm", "mass->tree.stm", "plant/tree", "tree stem and branch carbon", "Mg C ha-1", "core", (STRUCTURE, RESTART, PLANT_PROCESS)),
    StateVariable("tree_rot", "mass->tree.rot", "plant/tree", "tree root carbon", "Mg C ha-1", "core", (STRUCTURE, RESTART, PLANT_PROCESS)),
    StateVariable("c3_fol", "mass->c3.fol", "plant/c3", "C3 foliage carbon", "Mg C ha-1", "core", (STRUCTURE, RESTART, PLANT_PROCESS)),
    StateVariable("c3_stm", "mass->c3.stm", "plant/c3", "C3 stem carbon", "Mg C ha-1", "core", (STRUCTURE, RESTART, PLANT_PROCESS)),
    StateVariable("c3_rot", "mass->c3.rot", "plant/c3", "C3 root carbon", "Mg C ha-1", "core", (STRUCTURE, RESTART, PLANT_PROCESS)),
    StateVariable("c4_fol", "mass->c4.fol", "plant/c4", "C4 foliage carbon", "Mg C ha-1", "core", (STRUCTURE, RESTART, PLANT_PROCESS)),
    StateVariable("c4_stm", "mass->c4.stm", "plant/c4", "C4 stem carbon", "Mg C ha-1", "core", (STRUCTURE, RESTART, PLANT_PROCESS)),
    StateVariable("c4_rot", "mass->c4.rot", "plant/c4", "C4 root carbon", "Mg C ha-1", "core", (STRUCTURE, RESTART, PLANT_PROCESS)),
    StateVariable("ltr_tf", "mass->soil.ltr_tf", "soil/litter", "tree foliage litter carbon", "Mg C ha-1", "core", (STRUCTURE, RESTART, SOIL_PROCESS)),
    StateVariable("ltr_tc", "mass->soil.ltr_tc", "soil/litter", "tree stem litter carbon", "Mg C ha-1", "core", (STRUCTURE, RESTART, SOIL_PROCESS)),
    StateVariable("ltr_tr", "mass->soil.ltr_tr", "soil/litter", "tree root litter carbon", "Mg C ha-1", "core", (STRUCTURE, RESTART, SOIL_PROCESS)),
    StateVariable("ltr_gf", "mass->soil.ltr_gf", "soil/litter", "grass foliage litter carbon", "Mg C ha-1", "core", (STRUCTURE, RESTART, SOIL_PROCESS)),
    StateVariable("ltr_gc", "mass->soil.ltr_gc", "soil/litter", "grass stem litter carbon", "Mg C ha-1", "core", (STRUCTURE, RESTART, SOIL_PROCESS)),
    StateVariable("ltr_gr", "mass->soil.ltr_gr", "soil/litter", "grass root litter carbon", "Mg C ha-1", "core", (STRUCTURE, RESTART, SOIL_PROCESS)),
    StateVariable("msl_a", "mass->soil.msl_a", "soil/humus", "active humus carbon", "Mg C ha-1", "core", (STRUCTURE, RESTART, SOIL_PROCESS)),
    StateVariable("msl_i", "mass->soil.msl_i", "soil/humus", "intermediate humus carbon", "Mg C ha-1", "core", (STRUCTURE, RESTART, SOIL_PROCESS)),
    StateVariable("msl_p", "mass->soil.msl_p", "soil/humus", "passive humus carbon", "Mg C ha-1", "core", (STRUCTURE, RESTART, SOIL_PROCESS)),
)


# Carbon variables that are dynamically modified by native source code and must
# be considered when constructing an exact rather than reduced carbon state.
VISIT_EXTENDED_CARBON: Tuple[StateVariable, ...] = (
    StateVariable("tree_nsch", "mass->tree.nsch_storage", "plant/tree", "non-structural carbohydrate storage", "Mg C ha-1", "extended-required-for-exact-plant-map", (STRUCTURE, PLANT_PROCESS)),
    StateVariable("c3_nsch", "mass->c3.nsch_storage", "plant/c3", "non-structural carbohydrate storage", "Mg C ha-1", "extended-required-for-exact-plant-map", (STRUCTURE, PLANT_PROCESS)),
    StateVariable("c4_nsch", "mass->c4.nsch_storage", "plant/c4", "non-structural carbohydrate storage", "Mg C ha-1", "extended-required-for-exact-plant-map", (STRUCTURE, PLANT_PROCESS)),
    StateVariable("tree_grain", "mass->tree.grain", "plant/tree", "crop grain carbon when applicable", "Mg C ha-1", "configuration-dependent", (STRUCTURE, PLANT_PROCESS)),
    StateVariable("c3_grain", "mass->c3.grain", "plant/c3", "crop grain carbon when applicable", "Mg C ha-1", "configuration-dependent", (STRUCTURE, PLANT_PROCESS)),
    StateVariable("c4_grain", "mass->c4.grain", "plant/c4", "crop grain carbon when applicable", "Mg C ha-1", "configuration-dependent", (STRUCTURE, PLANT_PROCESS)),
)


VISIT_MET_INPUTS: Tuple[InputVariable, ...] = (
    InputVariable("air_temperature", "loct->tmp_2m", "external", "2-m air temperature", (STRUCTURE, LOCATION_PROCESS)),
    InputVariable("surface_temperature", "loct->tmp_sfc", "external", "surface temperature used by photosynthesis and shoot respiration", (STRUCTURE, LOCATION_PROCESS, PHOTOSYNTHESIS, RESPIRATION)),
    InputVariable("soil_temperature_10cm", "loct->tmp10_soil", "external/effective", "upper-soil temperature used by root respiration and litter decomposition", (STRUCTURE, LOCATION_PROCESS, RESPIRATION, DECOMPOSITION)),
    InputVariable("soil_temperature_deep", "loct->tmp200_soil", "external/effective", "deeper-soil temperature used by humus decomposition", (STRUCTURE, LOCATION_PROCESS, DECOMPOSITION)),
    InputVariable("shortwave_radiation", "loct->dswrf_sfc", "external", "downward shortwave radiation", (STRUCTURE, LOCATION_PROCESS, PHOTOSYNTHESIS)),
    InputVariable("precipitation", "loct->prate_sfc", "external", "precipitation forcing; affects carbon mainly through hydrology", (STRUCTURE, LOCATION_PROCESS)),
    InputVariable("atmospheric_co2", "loct->aCO2", "external", "ambient atmospheric CO2 concentration", (STRUCTURE, ECOPHYS, PHOTOSYNTHESIS)),
    InputVariable("soil_water_upper", "loct->soilwtr_l", "effective", "upper-layer soil water controlling litter decomposition", (STRUCTURE, DECOMPOSITION)),
    InputVariable("soil_water_lower", "loct->soilwtr_h", "effective", "lower-layer soil water controlling photosynthesis/humus decomposition", (STRUCTURE, PHOTOSYNTHESIS, DECOMPOSITION)),
)


VISIT_EFFECTIVE_INPUTS: Tuple[InputVariable, ...] = (
    InputVariable("gpp_tree", "flux->tree.gpp", "effective-carbon", "tree gross primary production", (PLANT_PROCESS, PHOTOSYNTHESIS)),
    InputVariable("gpp_c3", "flux->c3.gpp", "effective-carbon", "C3 gross primary production", (PLANT_PROCESS, PHOTOSYNTHESIS)),
    InputVariable("gpp_c4", "flux->c4.gpp", "effective-carbon", "C4 gross primary production", (PLANT_PROCESS, PHOTOSYNTHESIS)),
    InputVariable("epp_tree", "flux->tree.epp", "effective-carbon", "tree GPP minus maintenance respiration", (PLANT_PROCESS, RESPIRATION, ALLOCATION)),
    InputVariable("epp_c3", "flux->c3.epp", "effective-carbon", "C3 GPP minus maintenance respiration", (PLANT_PROCESS, RESPIRATION, ALLOCATION)),
    InputVariable("epp_c4", "flux->c4.epp", "effective-carbon", "C4 GPP minus maintenance respiration", (PLANT_PROCESS, RESPIRATION, ALLOCATION)),
    InputVariable("decomp_scalar_litter", "echar->soil.f_tm_l", "effective-carbon", "litter decomposition environmental scalar", (DECOMPOSITION, SOIL_PROCESS)),
    InputVariable("decomp_scalar_humus", "echar->soil.f_tm_h", "effective-carbon", "humus decomposition environmental scalar", (DECOMPOSITION, SOIL_PROCESS)),
)


VISIT_OUTPUTS: Tuple[OutputVariable, ...] = (
    OutputVariable(
        "NEP",
        "flux->nep",
        "NPP_tree + funder_c3*NPP_c3 + funder_c4*NPP_c4 - soil heterotrophic respiration",
        (DAILY_SCHEME, SOIL_PROCESS),
    ),
    OutputVariable("NPP", "flux->npp", "area-weighted ecosystem net primary production", (DAILY_SCHEME, PLANT_PROCESS)),
    OutputVariable("GPP", "flux->gpp", "area-weighted ecosystem gross primary production", (DAILY_SCHEME, PLANT_PROCESS, PHOTOSYNTHESIS)),
    OutputVariable("NECB", "flux->necb", "NEP with harvest term included in native daily scheme", (DAILY_SCHEME,)),
    OutputVariable("total_carbon", "mass->total_c", "diagnostic ecosystem carbon total", (DAILY_SCHEME, STRUCTURE)),
)


VISIT_SOURCES: Tuple[SourceRef, ...] = (
    STRUCTURE,
    RESTART,
    LOCATION_PROCESS,
    DAILY_SCHEME,
    PLANT_PROCESS,
    PHOTOSYNTHESIS,
    ECOPHYS,
    RESPIRATION,
    ALLOCATION,
    LITTERFALL,
    DECOMPOSITION,
    SOIL_PROCESS,
)


def source_paths() -> tuple[str, ...]:
    """Unique source paths used by the initial VISIT mapping."""
    return tuple(dict.fromkeys(ref.path for ref in VISIT_SOURCES))


def visit_provenance_manifest() -> dict[str, object]:
    """Return provenance for the initial source-grounded VISIT mapping."""
    return provenance_manifest(
        model="VISIT",
        approximation_level="source-map",
        sources=VISIT_SOURCES,
    )
