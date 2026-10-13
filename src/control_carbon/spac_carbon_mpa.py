"""MPa input/output boundary for the unchanged, canonical Manzoni SPAC model.

VPD is positive vapor-pressure deficit, NOT atmospheric water potential.
Atmospheric pressure is mandatory. Historical paper APIs remain mol/mol.
"""

from dataclasses import dataclass
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING

import numpy as np

from .manzoni_hydraulics import HydraulicParameters, _forcing
from .spac_carbon import CarbonParameters
from . import spac_carbon_two_way as legacy
from . import spac_carbon_bifurcation as bifurcation
from .spac_carbon_stability import certify_stable_path, _Interval, PRECISION
from .visitc_hydrology import visitc_saturated_vapor_pressure


@dataclass(frozen=True)
class PressureSPAC:
    hydraulic: HydraulicParameters
    carbon: CarbonParameters
    feedback: legacy.StructuralFeedback
    root_potential_mpa: float
    atmospheric_pressure_mpa: float

    def __post_init__(self):
        _forcing(self.root_potential_mpa, 0.)
        if not np.isfinite(self.atmospheric_pressure_mpa) or self.atmospheric_pressure_mpa <= 0:
            raise ValueError("atmospheric pressure must be finite positive MPa")

    def mole_fraction(self, vpd_mpa):
        if not np.isfinite(vpd_mpa) or vpd_mpa < 0:
            raise ValueError("VPD must be finite nonnegative MPa")
        result = vpd_mpa/self.atmospheric_pressure_mpa
        if not np.isfinite(result):
            raise ValueError("VPD conversion exceeds representable range")
        return float(result)

    def pressure(self, mole_fraction):
        if not np.isfinite(mole_fraction) or mole_fraction < 0:
            raise ValueError("canonical VPD must be finite nonnegative mol/mol")
        result = mole_fraction*self.atmospheric_pressure_mpa
        if not np.isfinite(result):
            raise ValueError("VPD conversion exceeds representable range")
        return float(result)

    @property
    def gs_max_per_mpa(self):
        """Canopy atmospheric coefficient in m water/day/MPa."""
        return self.hydraulic.gs_max/self.atmospheric_pressure_mpa

    def air_diagnostic(self, vpd_mpa, temperature_c):
        self.mole_fraction(vpd_mpa)
        if not np.isfinite(temperature_c) or temperature_c <= -273.15:
            raise ValueError("diagnostic temperature must be finite and above absolute zero")
        saturation = visitc_saturated_vapor_pressure(temperature_c)*1e-4
        if not np.isfinite(saturation) or saturation <= 0:
            raise ValueError("saturation vapor pressure is not finite positive MPa")
        return {
            "saturation_vapor_pressure_mpa": saturation,
            "implied_relative_humidity": 1-vpd_mpa/saturation,
            "admissible_at_diagnostic_temperature": bool(
                0 <= vpd_mpa <= saturation < self.atmospheric_pressure_mpa
            ),
        }

    @property
    def _parameters(self):
        return self.root_potential_mpa, self.hydraulic, self.carbon, self.feedback

    def rhs(self, state, vpd_mpa):
        r, h, c, f = self._parameters
        return legacy.two_way_rhs(state, r, self.mole_fraction(vpd_mpa), h, c, f)

    def jacobian(self, state, vpd_mpa):
        r, h, c, f = self._parameters
        return legacy.two_way_jacobian(state, r, self.mole_fraction(vpd_mpa), h, c, f)

    def fluxes(self, state, vpd_mpa):
        r, h, c, f = self._parameters
        return legacy.structural_fluxes(state, r, self.mole_fraction(vpd_mpa), h, f)

    def equilibrium_curve(self, leaf_potential_mpa):
        state, d, slope = legacy.equilibrium_curve(leaf_potential_mpa, *self._parameters)
        return state, self.pressure(d), slope*self.atmospheric_pressure_mpa

    def stationary_points(self, bounds, *, samples=300):
        return legacy.stationary_points(bounds, *self._parameters, samples=samples)

    def equilibria(self, vpd_mpa, bounds, stationary):
        return legacy.equilibria_at_vpd(self.mole_fraction(vpd_mpa), bounds, stationary, *self._parameters)

    def continue_curve(self, start_leaf, step, count, bounds, *, vpd_scale_mpa):
        records, status = legacy.continue_curve(
            start_leaf, step, count, bounds, *self._parameters,
            d_scale=self.mole_fraction(vpd_scale_mpa),
        )
        return [(y, self.pressure(d)) for y, d in records], status

    def simulate(self, initial, times, forcing_mpa, *, rtol=1e-9):
        return legacy.simulate_two_way(
            initial, times, lambda t: self.mole_fraction(forcing_mpa(t)), *self._parameters, rtol=rtol,
        )

    def core_rhs(self, state, vpd_mpa):
        r, h, c, f = self._parameters
        return bifurcation.core_rhs(state, r, self.mole_fraction(vpd_mpa), h, c, f)

    def core_jacobian(self, state, vpd_mpa):
        r, h, c, f = self._parameters
        return bifurcation.core_jacobian(state, r, self.mole_fraction(vpd_mpa), h, c, f)

    def core_tensors(self, state, vpd_mpa):
        r, h, c, f = self._parameters
        return bifurcation.core_tensors(state, r, self.mole_fraction(vpd_mpa), h, c, f)

    def hopf_point(self, leaf_bounds):
        data = dict(bifurcation.hopf_point(leaf_bounds, *self._parameters))
        data["vpd_mpa"] = self.pressure(data.pop("D"))
        data["vpd_slope_per_mpa"] = data.pop("curve_slope")*self.atmospheric_pressure_mpa
        data["eigenvalue_slope_per_vpd_mpa_refinement"] = [
            v/self.atmospheric_pressure_mpa for v in data.pop("eigenvalue_slope_D_refinement")
        ]
        return data

    def periodic_orbit(self, hopf_mpa, amplitude, *, leaf_bracket, rtol=1e-10, steps_per_period=100):
        canonical = dict(hopf_mpa)
        canonical["D"] = self.mole_fraction(canonical.pop("vpd_mpa"))
        canonical["curve_slope"] = canonical.pop("vpd_slope_per_mpa")/self.atmospheric_pressure_mpa
        canonical["eigenvalue_slope_D_refinement"] = [
            v*self.atmospheric_pressure_mpa
            for v in canonical.pop("eigenvalue_slope_per_vpd_mpa_refinement")
        ]
        data = dict(bifurcation.periodic_orbit(
            canonical, amplitude, *self._parameters, leaf_bracket=leaf_bracket,
            rtol=rtol, steps_per_period=steps_per_period,
        ))
        data["vpd_mpa"] = self.pressure(data.pop("D"))
        data["predicted_vpd_shift_mpa"] = data.pop("predicted_D_shift")*self.atmospheric_pressure_mpa
        return data

    def certify_stable_path(self, leaf_bounds, vpd_mpa_bounds, **limits):
        low, high = vpd_mpa_bounds
        self.mole_fraction(low)
        self.mole_fraction(high)
        p = Decimal.from_float(float(self.atmospheric_pressure_mpa))
        with localcontext() as ctx:
            ctx.prec = PRECISION
            ctx.rounding = ROUND_FLOOR
            dlow = Decimal.from_float(float(low))/p
            ctx.rounding = ROUND_CEILING
            dhigh = Decimal.from_float(float(high))/p
        report = dict(certify_stable_path(leaf_bounds, (dlow, dhigh), *self._parameters, **limits))
        def scale_interval(values):
            return (_Interval(Decimal(values[0]), Decimal(values[1]))*_Interval(p)).strings()
        report.pop("requested_D_bounds")
        report["requested_vpd_mpa_bounds"] = [str(Decimal.from_float(float(v))) for v in [low, high]]
        report["endpoint_vpd_mpa"] = {k: scale_interval(v) for k, v in report.pop("endpoint_D").items()}
        report["cells"] = [
            {**{k: v for k, v in cell.items() if k not in ["D", "Dprime_L"]},
             "vpd_mpa": scale_interval(cell["D"]),
             "vpd_slope_per_leaf_carbon": scale_interval(cell["Dprime_L"])}
            for cell in report["cells"]
        ]
        upper = report.pop("Dprime_L_upper_bound")
        report["vpd_slope_per_leaf_carbon_upper_bound"] = scale_interval([upper, upper])[1]
        report["claim"] = "unique stable positive equilibrium on the certified branch for every VPD (MPa) in the requested interval"
        report["atmospheric_pressure_mpa"] = self.atmospheric_pressure_mpa
        return report
