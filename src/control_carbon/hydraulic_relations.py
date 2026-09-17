"""Source-labelled storage--water-potential constitutive relations.

The soil relation is a direct transcription of the three texture branches in
``visit-manager/VISITc::point/location_proc.c`` at commit ``5202debd...``.
The plant relation is an aggregation of the symplasmic pressure--volume curve
in Ruffault et al. (2022), SurEau-Ecos v2.0, equations 41--47.  Applying that
relation to leaf, stem, and root stores on a ground-area basis is a reduced
model choice, not native VISITc behaviour.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import IntEnum

import numpy as np

from .visitc_source_map import VISITC_LOCATION_SOURCE


SUREAU_ECOS_DOI = "https://doi.org/10.5194/gmd-15-5593-2022"


class VISITCSoilTexture(IntEnum):
    """Native ``grid->stexture`` selectors used by VISITc."""

    SANDY = 0
    MEDIUM = 1
    FINE = 2


_SOIL_RETENTION_COEFFICIENTS = {
    VISITCSoilTexture.SANDY: (0.121, 4.05),
    VISITCSoilTexture.MEDIUM: (0.478, 5.39),
    VISITCSoilTexture.FINE: (0.405, 11.4),
}


@dataclass(frozen=True)
class PlantPressureVolumeParameters:
    """Parameters of a reduced symplasmic pressure--volume reservoir.

    ``saturated_water`` is in mm per unit ground area, ``osmotic_potential``
    is the negative osmotic potential at full turgor in MPa, and
    ``bulk_modulus`` is the positive bulk modulus of elasticity in MPa.
    ``minimum_relative_water`` is a numerical dry-end regularization; the
    analytical curve tends to minus infinity as RWC tends to zero.
    """

    saturated_water: float
    osmotic_potential: float
    bulk_modulus: float
    minimum_relative_water: float = 1.0e-4

    def __post_init__(self) -> None:
        values = np.asarray(
            (
                self.saturated_water,
                self.osmotic_potential,
                self.bulk_modulus,
                self.minimum_relative_water,
            ),
            dtype=float,
        )
        if not np.all(np.isfinite(values)):
            raise ValueError("pressure-volume parameters must be finite")
        if self.saturated_water <= 0.0:
            raise ValueError("saturated_water must be positive")
        if self.osmotic_potential >= 0.0:
            raise ValueError("osmotic_potential must be negative")
        if self.bulk_modulus <= 0.0:
            raise ValueError("bulk_modulus must be positive")
        if not 0.0 < self.minimum_relative_water < 1.0:
            raise ValueError("minimum_relative_water must lie in (0, 1)")
        if self.relative_water_at_turgor_loss <= 0.0:
            raise ValueError(
                "bulk_modulus must exceed the magnitude of osmotic_potential"
            )

    @property
    def relative_water_at_turgor_loss(self) -> float:
        """Relative water content at the turgor-loss point."""
        return 1.0 + self.osmotic_potential / self.bulk_modulus

    @property
    def turgor_loss_potential(self) -> float:
        """Water potential at turgor loss, in MPa."""
        return 1.0 / (
            1.0 / self.bulk_modulus + 1.0 / self.osmotic_potential
        )


def visitc_soil_matric_potential(
    water: float,
    field_capacity: float,
    texture: VISITCSoilTexture | int,
    *,
    minimum_water: float = 0.2,
) -> float:
    """Return native VISITc matric potential in MPa.

    This reproduces the post-hydrology branches in
    ``point/location_proc.c::f_loct_proc``.  VISITc clips only the value used
    in this diagnostic to 0.2 mm; it does not add that water to the store.
    """
    values = np.asarray((water, field_capacity, minimum_water), dtype=float)
    if not np.all(np.isfinite(values)):
        raise ValueError("soil-retention inputs must be finite")
    if water < 0.0 or field_capacity <= 0.0 or minimum_water <= 0.0:
        raise ValueError(
            "water must be nonnegative and capacities/minimum_water positive"
        )
    try:
        texture_value = VISITCSoilTexture(texture)
    except ValueError as error:
        raise ValueError("texture must be VISITc selector 0, 1, or 2") from error
    coefficient, exponent = _SOIL_RETENTION_COEFFICIENTS[texture_value]
    relative_water = max(water, minimum_water) / field_capacity
    return -coefficient * relative_water ** (-exponent)


def visitc_soil_total_potential(
    water: float,
    field_capacity: float,
    texture: VISITCSoilTexture | int,
    *,
    gravitational_potential: float,
    minimum_water: float = 0.2,
) -> float:
    """Return matric plus source-prescribed gravitational potential, MPa."""
    if not np.isfinite(gravitational_potential):
        raise ValueError("gravitational_potential must be finite")
    return gravitational_potential + visitc_soil_matric_potential(
        water, field_capacity, texture, minimum_water=minimum_water
    )


def plant_pressure_volume_potential(
    water: float, parameters: PlantPressureVolumeParameters
) -> float:
    """Map plant water storage to symplasmic water potential in MPa.

    This is SurEau-Ecos equation 43 with ``RWC = W / W_sat``.  Storage above
    saturation is clipped to RWC=1 because the reduced model does not include
    a separate apoplasmic/overflow reservoir.
    """
    if not np.isfinite(water) or water < 0.0:
        raise ValueError("plant water must be finite and nonnegative")
    relative_water = np.clip(
        water / parameters.saturated_water,
        parameters.minimum_relative_water,
        1.0,
    )
    pi0 = parameters.osmotic_potential
    epsilon = parameters.bulk_modulus
    if relative_water >= parameters.relative_water_at_turgor_loss:
        return float(
            -pi0
            - epsilon * (1.0 - relative_water)
            + pi0 / relative_water
        )
    return float(pi0 / relative_water)


def plant_pressure_volume_capacitance(
    water: float, parameters: PlantPressureVolumeParameters
) -> float:
    """Return ``dW/dpsi`` in mm MPa-1 for the pressure--volume curve.

    The derivative follows SurEau-Ecos equations 42 and 44.  At or above
    saturation the reduced storage-to-potential map is clipped and its local
    derivative is therefore reported as zero.
    """
    if not np.isfinite(water) or water < 0.0:
        raise ValueError("plant water must be finite and nonnegative")
    raw_relative_water = water / parameters.saturated_water
    if (
        raw_relative_water >= 1.0
        or raw_relative_water <= parameters.minimum_relative_water
    ):
        return 0.0
    relative_water = max(raw_relative_water, parameters.minimum_relative_water)
    pi0 = parameters.osmotic_potential
    epsilon = parameters.bulk_modulus
    potential = plant_pressure_volume_potential(water, parameters)
    if relative_water >= parameters.relative_water_at_turgor_loss:
        derivative = relative_water / (
            -pi0 - potential - epsilon + 2.0 * epsilon * relative_water
        )
    else:
        derivative = -pi0 / potential**2
    return float(parameters.saturated_water * derivative)


HYDRAULIC_RELATION_PROVENANCE = {
    "soil": {
        "source": VISITC_LOCATION_SOURCE.to_dict(),
        "equation": "psi_m=-a_texture*(max(W,0.2)/field_capacity)^(-b_texture)",
        "approximation": "exact source diagnostic algebra",
    },
    "plant": {
        "source": SUREAU_ECOS_DOI,
        "equations": "41-47",
        "approximation": (
            "SurEau-Ecos symplasmic pressure-volume curve aggregated to each "
            "leaf/stem/root store on a ground-area basis"
        ),
    },
}


__all__ = [
    "HYDRAULIC_RELATION_PROVENANCE",
    "PlantPressureVolumeParameters",
    "SUREAU_ECOS_DOI",
    "VISITCSoilTexture",
    "plant_pressure_volume_capacitance",
    "plant_pressure_volume_potential",
    "visitc_soil_matric_potential",
    "visitc_soil_total_potential",
]
