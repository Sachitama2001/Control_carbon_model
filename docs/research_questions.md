# Research Questions and Analysis Matrix

This document keeps the implementation tied to the scientific questions.

## RQ1. Can terrestrial ecosystem models be compared through common input-output dynamics despite different internal structures?

### Hypothesis

Different models may have different state dimensions and pool definitions but still admit comparable quantities:

- dominant poles / time scales,
- total-carbon and NEP transfer functions,
- low-frequency QSE sensitivities,
- forcing-rate response,
- modal residues and source/sink reversal behavior.

### Minimum demonstration

Run the same analysis pipeline for VISIT and at least one second terrestrial ecosystem model.

---

## RQ2. Can NEP source/sink sign be inferred analytically from QSE disequilibrium and forcing rate?

### Core objects

```text
x_star(z)
p = x_star - x
NEP = h(x,z)
```

Local linear theory should quantify how `p`, `dz/dt`, and output sensitivities combine.

### Test

Compare analytical sign boundaries against explicit time-domain simulation under ramps in temperature and/or other forcing.

### Deliverable

A formula or computable criterion giving the local source/sink boundary in forcing-rate space, plus parameter sensitivities of that boundary.

---

## RQ3. What additional information does IRF/modal analysis reveal beyond carbon storage capacity and equilibrium diagnostics?

### Candidate findings

- multiple time scales around the same QSE;
- transient source-to-sink or sink-to-source reversal;
- hidden internal modes with weak NEP projection;
- forcing-frequency-dependent sensitivity;
- cancellation between photosynthetic and decomposition pathways;
- reduced-order representations for fast mapping.

### Test

Construct examples where total QSE displacement alone cannot predict transient NEP trajectory, but pole/residue information can.

---

## RQ4. Is nonlinearity necessary for interesting non-monotone source/sink trajectories?

### Working expectation

No. A stable linear multi-pool system can show sign reversal when multiple modes contribute opposite-signed residues to NEP.

For a two-mode response:

```text
NEP(t) = a1 exp(lambda1 t) + a2 exp(lambda2 t)
```

opposite signs in `a1`, `a2` may produce a crossing.

### Test

Identify source-grounded VISIT operating points where local linearized modes have opposing NEP residues. Then verify the predicted zero crossing with small perturbations.

---

## RQ5. When does rate-induced behavior require genuinely nonlinear/global analysis rather than local IRF theory?

### Strategy

Separate:

1. local lag under moving QSE;
2. local source/sink reversal due to modal competition;
3. loss of tracking / large departure;
4. true rate-induced tipping between attractors or basins.

Do not use "tipping" for ordinary transient sign reversal.

### Test

Increase forcing amplitude/rate until local linear predictions fail systematically, then inspect whether the native nonlinear model contains multiple attracting regimes or switching mechanisms relevant to R-tipping.

---

## RQ6. Which forcing should be treated as the primary external input?

### Initial priority

Temperature, because it directly affects:

- photosynthetic capacity,
- autotrophic respiration,
- litter decomposition,
- soil organic matter decomposition,
- phenology indirectly.

### But temperature alone may be insufficient

Add in stages:

```text
T -> (T, soil water) -> (T, soil water, radiation) -> (+ CO2/VPD)
```

Use effective-driver transfer functions as intermediate diagnostics when direct meteorology-to-carbon dynamics are too nonlinear.

---

## RQ7. How much state is required for a faithful local input-output representation?

### Candidate hierarchy

```text
9 soil pools
18 structural carbon pools
18 + NSC
18 + NSC + hydrology/phenology memory
full relevant native state
```

### Experiment

At matched structural carbon state, vary omitted state variables and measure whether the next-state/output differs.

### Interpretation

If it differs, the reduced state is not Markov. This can motivate state augmentation or reduced memory kernels.

---

## RQ8. Can IRFs reduce computational cost for global source/sink mapping?

### Candidate workflow

At representative operating points/grid cells:

1. solve QSE/reference trajectory;
2. compute local linear model or identify IRF;
3. reduce modes;
4. use convolution rather than rerunning full process model for many small perturbation scenarios;
5. quantify approximation error.

### Metrics

```text
runtime speedup
NEP RMSE
Ctot RMSE
sign classification accuracy
critical-rate error
```

Do not claim speedup until end-to-end benchmarks include model setup and linearization cost.

---

## RQ9. How do parameters control source/sink boundaries?

For parameter `theta`, examine derivatives of:

```text
poles
residues
QSE sensitivity
DC gain
critical forcing rate
```

Potential parameter classes:

- decomposition temperature sensitivity;
- base turnover rates;
- photosynthetic optimum temperature;
- allocation parameters;
- maintenance respiration Q10;
- soil moisture limitation parameters.

Map direction and magnitude, not merely rank correlation.

---

## RQ10. Is an annual-mean LTI system adequate?

### Concern

VISIT contains strong seasonal phenology, climate, and LAI variation.

### Comparison

Evaluate:

```text
annual-mean LTI
season-specific LTI
periodic LTV/Floquet representation
native perturbation response
```

### Scientific criterion

Annual-mean LTI is inadequate if it mispredicts important NEP sign changes, dominant time scales, or forcing-phase dependence.

---

# Experiment matrix

| Experiment | Model level | Input | Output | Main purpose |
|---|---|---|---|---|
| E01 | 9-pool VISIT soil | litter impulse | pool C, Rh | validate compartment IRF |
| E02 | 9-pool VISIT soil | decomposition-scalar step | Rh, total soil C | test transfer response |
| E03 | 9-pool VISIT soil | temperature step via `frl/frh` | Rh | test local nonlinear-to-linear mapping |
| E04 | reduced18 | EPP impulse per PFT | NEP, Ctot | identify modal pathways |
| E05 | reduced18 | small temperature step | NEP | find opposite-signed residues |
| E06 | reduced18 | temperature ramp | NEP | test forcing-rate/QSE relation |
| E07 | native VISIT | same small perturbations | NEP, pools | validate local model |
| E08 | native vs reduced | varied perturbation amplitude | errors | linearity radius |
| E09 | periodic VISIT | perturbation by day-of-year | NEP | phase-dependent IRF |
| E10 | cross-model | matched forcing | normalized transfer metrics | generality |

# Publication-quality outputs to target eventually

1. Diagram translating source-code process topology into state-space representation.
2. Pole/time-scale spectrum by model and biome/site.
3. NEP modal residue decomposition.
4. Example analytical source/sink reversal with simulation validation.
5. QSE sensitivity and critical forcing-rate maps.
6. Error/speed comparison of full simulation vs IRF emulator.
7. Cross-model comparison showing robust vs model-specific modes.

# Interpretation discipline

- A pole is not automatically a named ecological process; use participation factors/provenance before assigning interpretation.
- A source/sink sign reversal is not automatically tipping.
- A local QSE is not necessarily globally attracting.
- Fast computational classification is useful only if its validity domain is mapped.
- Differences between models are scientific information, not errors to be calibrated away by default.
