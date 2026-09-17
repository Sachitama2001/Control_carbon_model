"""Pinned provenance for the current ``visit-manager/VISITc`` source tree.

The older :mod:`control_carbon.visit_source_map` remains the authority for the
existing soil-carbon adapter.  This module deliberately uses a separate name
and commit so equations read from the two VISIT snapshots cannot be mixed
silently.
"""

from __future__ import annotations

from .provenance import SourceRef, provenance_manifest


VISITC_SOURCE_REPOSITORY = "visit-manager/VISITc"
VISITC_SOURCE_COMMIT = "5202debd96df6f88beb7d61f8688fff02ace964a"
VISITC_SOURCE_ROOT = "point"


def _source(path: str, symbol: str, role: str) -> SourceRef:
    return SourceRef(
        path=path,
        symbol=symbol,
        role=role,
        repository=VISITC_SOURCE_REPOSITORY,
        commit=VISITC_SOURCE_COMMIT,
        function=symbol,
    )


VISITC_STRUCTURE_SOURCE = _source(
    "point/structure.h",
    "struct Loct / struct Mass",
    "water-store, flux, and diagnostic declarations",
)
VISITC_HYDROLOGY_SOURCE = _source(
    "point/hydro_balance.c",
    "f_hydrology",
    "source-ordered daily snow, interception, soil-water, transpiration, and runoff update",
)
VISITC_LOCATION_SOURCE = _source(
    "point/location_proc.c",
    "f_loct_proc",
    "hydrology call order and conversion of water stores to carbon-process drivers",
)
VISITC_DAILY_SCHEME_SOURCE = _source(
    "point/daily_scheme.c",
    "daily_scheme",
    "daily carbon update and post-update water diagnostics",
)
VISITC_HYDRO_FLOWS_SOURCE = _source(
    "point/hydro_flows.c",
    "f_airdens/f_vap_pre_sat/f_slope_vps/f_r_aero/pm_incep/pm_evap/pm_transp",
    "atmospheric diagnostics, resistances, and Penman-Monteith potential fluxes",
)
VISITC_RADIATION_SOURCE = _source(
    "point/radiation.c",
    "f_net_rad",
    "LAI-dependent partition of longwave and shortwave net radiation",
)
VISITC_CANOPY_CONDUCTANCE_SOURCE = _source(
    "point/ecophysiology.c",
    "f_canopy_cond",
    "LAI, GPP, CO2, and VPD dependence of canopy conductance",
)

VISITC_SOURCES = (
    VISITC_STRUCTURE_SOURCE,
    VISITC_HYDROLOGY_SOURCE,
    VISITC_LOCATION_SOURCE,
    VISITC_DAILY_SCHEME_SOURCE,
    VISITC_HYDRO_FLOWS_SOURCE,
    VISITC_RADIATION_SOURCE,
    VISITC_CANOPY_CONDUCTANCE_SOURCE,
)


def visitc_provenance_manifest() -> dict[str, object]:
    """Return a JSON-compatible manifest for the inspected VISITc snapshot."""
    return provenance_manifest(
        model="VISITc",
        approximation_level="source-map",
        sources=VISITC_SOURCES,
    )
