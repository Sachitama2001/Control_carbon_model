"""Rational E4 with the Manzoni soil/root pathway, without a soil storage ODE."""

from dataclasses import dataclass

import numpy as np
from numpy.polynomial import Polynomial

from .manzoni_extensions import SoilRootParameters, soil_root_conductance
from .manzoni_hydraulics import conductances
from .spac_carbon_equilibria import stock_ratios
from .spac_carbon_rational import RationalSPAC


@dataclass(frozen=True)
class SoilBoundaryRationalSPAC:
    plant: RationalSPAC
    soil: SoilRootParameters
    psi_soil_mpa: float

    def __post_init__(self):
        if self.plant.label != "E4":
            raise ValueError("Soil-root coupling currently requires E4")
        soil_root_conductance(self.psi_soil_mpa, self.soil)

    def root_potential(self, state):
        y = self.plant._state(state)
        gsr = soil_root_conductance(self.psi_soil_mpa, self.soil)
        gx, _ = conductances(y[1], self.plant.base.hydraulic)
        _, _, fw = self.plant.factors(y)
        return (gsr*self.psi_soil_mpa+2*gx*fw*y[0])/(gsr+2*gx*fw)

    def fluxes(self, state, vpd_mpa):
        y = self.plant._state(state)
        fluxes = self.plant.fluxes(y, vpd_mpa)
        gsr = soil_root_conductance(self.psi_soil_mpa, self.soil)
        gx, _ = conductances(y[1], self.plant.base.hydraulic)
        _, _, fw = self.plant.factors(y)
        root_xylem = 2*gx*fw
        fluxes[0] = gsr*root_xylem/(gsr+root_xylem)*(self.psi_soil_mpa-y[0])
        return fluxes

    def rhs(self, state, vpd_mpa):
        rhs = self.plant.rhs(state, vpd_mpa)
        erx, exl, _ = self.fluxes(state, vpd_mpa)
        rhs[0] = (erx-exl)/self.plant.base.hydraulic.cx
        return rhs

    def capacity(self, state):
        return self.plant.capacity(state)

    def jacobian(self, state, vpd_mpa):
        y = self.plant._state(state)
        jac = self.plant.jacobian(y, vpd_mpa)
        h, f = self.plant.base.hydraulic, self.plant.base.feedback
        gx, _ = conductances(y[1], h)
        dgx = -h.gx_max/(2*h.psi50_x) if gx > 0 else 0.
        fl, _, fw = self.plant.factors(y)
        gsr = soil_root_conductance(self.psi_soil_mpa, self.soil)
        root_xylem = 2*gx*fw
        effective = gsr*root_xylem/(gsr+root_xylem)
        sensitivity = (gsr/(gsr+root_xylem))**2
        drop = self.psi_soil_mpa-y[0]
        jac[0] = 0.
        jac[0, 0] = (-effective-2*gx*fl)/h.cx
        jac[0, 1] = (sensitivity*2*dgx*fw*drop
                      -2*dgx*fl*(y[0]-y[1])+2*gx*fl)/h.cx
        jac[0, 2] = -2*gx*(y[0]-y[1])/(h.cx*(f.leaf_reference+f.leaf_floor))
        jac[0, 3] = sensitivity*2*gx*drop/(h.cx*f.wood_reference)
        return jac

    @property
    def wet_leaf(self):
        c = self.plant.base.carbon
        opening = 1-self.psi_soil_mpa/(2*self.plant.base.hydraulic.psi50_s)
        return c.allocation[0]*c.assimilation_max*opening/c.losses[0]-1/c.canopy_coefficient

    def curve(self, leaf):
        if not np.isfinite(leaf) or not 0 < leaf <= self.wet_leaf:
            raise ValueError("Soil-boundary curve requires a positive physical leaf carbon")
        y, _, _ = self.plant.curve(leaf)
        if leaf == self.wet_leaf:
            y[1] = self.psi_soil_mpa
        gx, gs = conductances(y[1], self.plant.base.hydraulic)
        fl, fe, fw = self.plant.factors(y)
        gsr = soil_root_conductance(self.psi_soil_mpa, self.soil)
        supply = gsr*(2*gx*fw)/(gsr+2*gx*fw)
        leaf_supply = 2*gx*fl
        y[0] = (supply*self.psi_soil_mpa+leaf_supply*y[1])/(supply+leaf_supply)
        flux = supply*leaf_supply/(supply+leaf_supply)*(self.psi_soil_mpa-y[1])
        vpd = flux/(gs*fe)*self.plant.base.atmospheric_pressure_mpa
        return y, vpd

    def equilibrium_polynomial(self, vpd_mpa):
        self.plant.base.mole_fraction(vpd_mpa)
        c, h, f = self.plant.base.carbon, self.plant.base.hydraulic, self.plant.base.feedback
        leaf = Polynomial([0., 1.])
        opening = Polynomial([self.plant.opening(0.),
                               c.losses[0]/(c.allocation[0]*c.assimilation_max)])
        psi = -2*h.psi50_s*(opening-1)
        gx, gs = h.gx_max*(1-psi/(2*h.psi50_x)), h.gs_max*opening
        fl = (leaf+f.leaf_floor)/(f.leaf_reference+f.leaf_floor)
        ratio = stock_ratios(self.plant.base)[1]
        fw = ratio*leaf/f.wood_reference
        gsr = soil_root_conductance(self.psi_soil_mpa, self.soil)
        numerator = (2*gx*gsr*(ratio*f.leaf_reference/f.wood_reference)
                     *fl*(self.psi_soil_mpa-psi))
        denominator = gs*(gsr*(fw+fl)+2*gx*fw*fl)
        return self.plant.base.mole_fraction(vpd_mpa)*denominator-numerator

    def equilibria(self, vpd_mpa):
        polynomial = self.equilibrium_polynomial(vpd_mpa)
        if vpd_mpa == 0 and self.wet_leaf > 0:
            return [self.curve(self.wet_leaf)[0]]
        roots = []
        for root in polynomial.roots():
            if abs(root.imag) <= 1e-9*max(1., abs(root.real)) and 0 < root.real <= self.wet_leaf:
                state, _ = self.curve(float(root.real))
                if np.max(np.abs(self.rhs(state, vpd_mpa))) > 1e-9:
                    raise RuntimeError("Soil-boundary root failed the full-system residual")
                if not any(abs(state[2]-other[2]) < 1e-10 for other in roots):
                    roots.append(state)
        return sorted(roots, key=lambda state: -state[2])
