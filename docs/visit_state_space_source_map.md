# VISIT state-space source map

This note records the first source-grounded state-space decomposition of the VISIT site model stored in `Sachitama2001/VISIT-matrix/visit_local`.

## Scope and source revision

The source headers identify the model as **VISIT (Vegetation Integrative Simulator for Trace gases)**, formerly Sim-CYCLE, developed by A. Ito. The source snapshot states a January 24, 2013 version, with later inline revisions in several routines.

Primary source files inspected:

- `visit_local/structure.h`: state, parameter, and flux structures.
- `visit_local/experiment.c`: experiment loop and restart variables.
- `visit_local/location_proc.c`: meteorological forcing and local environmental preprocessing.
- `visit_local/daily_scheme.c`: order of the daily ecosystem update and ecosystem-level outputs.
- `visit_local/plant_proc.c`: plant carbon state updates.
- `visit_local/photosynthesis.c`: daily GPP and photosynthetic temperature/water/CO2 response.
- `visit_local/ecophysiology.c`: photosynthetic/stomatal/respiration parameter updates and diagnostic LAI relationships.
- `visit_local/respiration.c`: plant maintenance and growth respiration.
- `visit_local/allocation.c`: state-dependent allocation of effective primary production.
- `visit_local/litterfall.c`: plant-to-litter turnover.
- `visit_local/decomposition.c`: litter/humus environmental response factors.
- `visit_local/soil_proc.c`: litter and humus decomposition, humification, heterotrophic respiration, and soil carbon mass balance.

An existing reduced matrix implementation in `VISIT-matrix/visit_matrix/` is treated as a useful prior implementation, but the C source above is the authority for this project.

## 1. Native model is a discrete-time nonlinear state transition

The experiment driver executes, for each day,

1. `f_loct_proc(...)` to construct the daily environment,
2. disturbance handling,
3. `daily_scheme(...)` for biological and soil carbon updates.

Thus the most faithful state-space representation is initially discrete:

```text
x[k+1] = F(x[k], u[k], theta, q[k])
y[k]   = H(x[k], u[k], theta, q[k])
```

where `q[k]` denotes discrete regime variables such as phenological stage, crop stage, and switches. A continuous-time ODE representation is a later approximation/localization, not the native numerical form.

## 2. Carbon state hierarchy

### 2.1 Core 18-pool structural carbon state

The first model-independent carbon state used for control/IRF work is

```text
x_C = [
  tree.fol, tree.stm, tree.rot,
  c3.fol,   c3.stm,   c3.rot,
  c4.fol,   c4.stm,   c4.rot,
  soil.ltr_tf, soil.ltr_tc, soil.ltr_tr,
  soil.ltr_gf, soil.ltr_gc, soil.ltr_gr,
  soil.msl_a, soil.msl_i, soil.msl_p
]^T
```

All 18 variables are explicit carbon masses in `struct Pmas` or `struct Smas`. They are also the carbon pools explicitly written/read in the restart path in `experiment.c`, making this a defensible minimal dynamic carbon state.

The existing `VISIT-matrix/visit_matrix/matrices.py` independently uses the same 18-pool ordering. This agreement is useful validation, but the state definition here is grounded in `structure.h`, `experiment.c`, `plant_proc.c`, and `soil_proc.c`.

### 2.2 Extended carbon states

`struct Pmas` additionally contains:

- `nsch_storage`: non-structural carbohydrate storage,
- `grain`: crop grain carbon,
- `stm_sp`, `stm_ht`: sapwood/heartwood sub-compartments,
- `rot_fn`, `rot_tp`: fine/coarse-root sub-compartments.

`plant_proc.c` explicitly changes `nsch_storage` during emergence and allocation, and changes `grain` through `flux->tpg` and harvest. Therefore, for an **exact full carbon state-space representation**, these variables cannot automatically be treated as parameters or diagnostics. The 18-pool system is a reduced structural-carbon subsystem.

`struct Smas` also declares `mcrb`, `doc`, and `sdw`. Their inclusion will be decided only after tracing all update routines that alter them; declaration alone is not sufficient evidence that they belong to the active state of the selected model configuration.

### 2.3 Diagnostics, not independent states

At the carbon-only level:

- `mass->tree.plant`, `mass->c3.plant`, `mass->c4.plant` are sums of foliage/stem/root after the daily plant update.
- `mass->soil.soil` is the sum of the six litter and three humus pools in `f_cycle_soil()`.
- `mass->total_c` is assembled in `daily_scheme()` from plant and soil totals.
- `LAI` is normally recalculated by `lai_mass()` from foliage carbon and SLA, although fixed-phenology modes may overwrite it.

Therefore these are outputs/algebraic variables for the default reduced carbon state, not additional independent carbon states.

## 3. Inputs and forcing hierarchy

The project will distinguish three input layers.

### Layer A: meteorological/external forcing

From `Loct` and `f_loct_proc()`:

- air/surface/soil temperature (`tmp_2m`, `tmp_sfc`, `tmp10_soil`, `tmp200_soil`),
- shortwave radiation (`dswrf_sfc`),
- precipitation (`prate_sfc`),
- humidity / VPD,
- wind,
- atmospheric CO2 (`aCO2` / scenario variables).

These are the scientifically desirable external inputs for transfer functions.

### Layer B: effective physiological/environmental drivers

Examples:

- PFT GPP / EPP,
- `psat`, `lue`, `ci`, stomatal response,
- litter decomposition scalar `f_tm_l = frl(...)`,
- humus decomposition scalar(s) from `frh(...)`,
- phenological stage and allocation fractions.

These are useful intermediate inputs because the carbon subsystem is much closer to linear/affine in them.

### Layer C: carbon flux inputs into pools

For a frozen regime, plant allocation fluxes `tpf`, `tpc`, `tpr` and litter inputs `li_*` can be treated as direct inputs to carbon boxes. This is the most classical compartmental representation, but it hides physiological feedbacks and is therefore a reduced input-output model rather than the full VISIT system.

## 4. Process-to-equation mapping

### Photosynthesis

Source: `photosynthesis.c::f_gpp()` and `f_pc_sat()`, called in the plant/ecophysiology sequence.

`f_gpp()` uses the Monsi-Saeki/Kuroiwa canopy-light equation and depends on `mass->lai`, `pchar->psat`, `pchar->lue`, PPFD, day length, and light attenuation. `f_pc_sat()` makes `psat` nonlinear in surface temperature, intercellular CO2, and soil water; it also imposes a frozen-soil threshold when `tmp10_soil < 2 C`.

Consequence: meteorology-to-GPP is nonlinear and piecewise smooth. A local Jacobian or empirical IRF is appropriate; a single global LTI transfer function is not.

### Autotrophic respiration

Source: `respiration.c::{f_rfm,f_rcm,f_rrm,f_rfg,f_rcg,f_rrg}`.

Maintenance respiration is proportional to the relevant plant carbon stock with an exponential Q10-type temperature factor. Growth respiration is proportional to allocated carbon flux. These terms provide analytically differentiable contributions to the local state Jacobian and direct-feedthrough terms.

### Allocation

Source: `allocation.c::f_allocation()`.

Allocation depends on EPP sign, LAI relative to optimum LAI, phenological season, and crop stage. Therefore allocation is state dependent and piecewise defined. This is a major reason the full plant system is nonlinear/hybrid even before considering photosynthesis.

### Plant turnover and litter input

Sources: `litterfall.c::{f_lf,f_lc,f_lr}`, `plant_proc.c`, and `daily_scheme.c`.

Outside special crop/phenology branches, litterfall/mortality is first order in foliage/stem/root carbon. `daily_scheme()` then maps tree litter fluxes directly to tree litter pools and combines C3/C4 litter fluxes using `funder_c3` and `funder_c4` before the soil update.

### Soil decomposition and humification

Sources: `decomposition.c::{frl,frh}` and `soil_proc.c::f_cycle_soil()`.

For each litter/humus pool the degraded amount is stock times a base decomposition rate times an environmental scalar. Litter decomposition is partitioned between microbial respiration and transfers to active/intermediate/passive humus. Hence, conditional on environmental scalars, the nine-pool soil subsystem is linear compartmental.

The environmental scalars are nonlinear in soil temperature and moisture, including Lloyd-Taylor-type temperature expressions and saturating moisture/aperture limitations.

## 5. Output equation: NEP

Source: `daily_scheme.c`.

VISIT computes ecosystem NEP explicitly as

```text
NEP = NPP_tree + funder_c3*NPP_c3 + funder_c4*NPP_c4 - soil.hr
```

where `soil.hr` is the sum of heterotrophic respiration terms assembled in `f_cycle_soil()`.

Thus NEP is an output, not a state. Around an operating point `(x*, u*)`, the local output equation is

```text
delta_NEP = C_NEP delta_x + D_NEP delta_u
```

with both terms generally nonzero because:

- respiration depends directly on state stocks,
- GPP/NPP and decomposition depend directly on meteorological/effective forcing,
- GPP/NPP also depend indirectly on foliage and other state-dependent physiological variables.

For the reduced frozen-coefficient compartment model, NEP can also be expressed through the total-carbon derivative, subject to careful treatment of harvest, BVOC, DOC export, disturbances, and other NECB terms.

## 6. First state-space levels to implement

We will maintain explicitly different approximation levels rather than call all of them "VISIT".

### Level 0: native nonlinear daily map

```text
x[k+1] = F_VISIT(x[k], z[k])
NEP[k] = H_VISIT(x[k], z[k])
```

Reference implementation: the C source itself.

### Level 1: 18-pool frozen-driver compartment model

```text
dx/dt = B(t) mu(t) + M(t) x
```

This is close to the existing `visit_matrix` implementation and is exact for the selected soil equations given fixed environmental scalars, but only an approximation to the complete plant update because allocation, respiration parameters, phenology, NSC storage, and GPP depend on state/regime.

### Level 2: local linearization of the native carbon map

For a daily discrete map:

```text
delta_x[k+1] = A_d delta_x[k] + B_d delta_z[k]
delta_y[k]   = C_d delta_x[k] + D_d delta_z[k]
```

where the Jacobians are evaluated at a specified equilibrium, periodic orbit, or climatological operating point.

For continuous-time IRF notation, a consistent continuous realization may then be obtained when justified.

### Level 3: periodic linear time-varying system

Because phenology and climate are strongly seasonal, an annual periodic linearization is likely more faithful than one annual-mean LTI model. Its monodromy matrix will provide annual modes and an annual IRF/transfer representation.

## 7. Important source-code features for control analysis

1. **Sequential update:** plant pools are updated before their litter fluxes are passed into `f_cycle_soil()` on the same day. The exact daily model is therefore a composed discrete map, not a simultaneous Euler update of one ODE.
2. **Hybrid/piecewise behavior:** phenological stages, crop stages, frozen-soil photosynthesis, allocation branches, nonnegative clipping, and disturbance switches introduce nonsmooth boundaries.
3. **State-dependent coefficients:** allocation, LAI-related photosynthesis, size-dependent respiration parameters, NSC storage constraints, and hydrology prevent a globally constant matrix representation.
4. **Soil subsystem is especially amenable to classical control/IRF analysis:** once temperature/moisture scalars are supplied, its carbon equations are first-order linear transfers.
5. **NEP is not NECB:** `daily_scheme.c` adds harvest separately in `necb`, and optional BVOC/DOC losses can alter total carbon. Source/sink statements must specify whether the target output is NEP, NEE, or NECB.

## 8. Immediate implementation tasks

- Build a provenance-aware `VISITSourceMap` in `Control_carbon_model`.
- Implement generic LTI/LTV state-space and IRF utilities independent of VISIT.
- Reproduce the 18-pool reduced VISIT model without silently copying assumptions from the older matrix code.
- Add tests that compare one-day reduced-model fluxes against algebraic calculations extracted from the source equations.
- Then construct Jacobians of the native daily map and compare its local IRF with the 18-pool reduced IRF.

The key design rule is: every matrix element or nonlinear response included in the adapter must point back to a named VISIT source file and function.