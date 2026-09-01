# Implementation Roadmap

## Guiding principle

Build from **source fidelity -> mathematical abstraction -> local linear analysis -> IRF/control diagnostics -> cross-model generalization**.

The roadmap is intentionally staged so that each layer can be validated independently.

---

## Phase A — Stabilize the core mathematical library

### A1. Continuous LTI API

Current file: `src/control_carbon/state_space.py`.

Tasks:

- add explicit stability checks;
- add DC gain with singularity handling;
- add step response;
- add forced response by convolution;
- add modal decomposition utilities;
- expose residues for SISO/MIMO channels where numerically sensible;
- add controllability/observability matrices and ranks as diagnostics;
- add balanced-truncation hooks later, after validation.

Tests:

- first-order scalar system with analytical impulse/step/transfer function;
- two-pole system with known modal sum;
- system with nonzero D;
- unstable system should not silently expose finite DC gain as if meaningful.

### A2. Discrete-time LTI API

Native VISIT is daily/discrete. Implement:

```python
@dataclass(frozen=True)
class DiscreteLTI:
    A: np.ndarray
    B: np.ndarray
    C: np.ndarray
    D: np.ndarray
    dt: float
```

Required methods:

- poles/eigenvalues;
- stability (`abs(lambda)<1`);
- transfer function in z-domain;
- impulse response `C A^k B`;
- step response;
- forced response;
- conversion to/from continuous LTI using documented assumptions (`scipy.signal.cont2discrete`, matrix logarithm only with care).

Tests should emphasize that discrete and continuous representations are not interchangeable by notation alone.

### A3. LTV/periodic API

Do not implement full Floquet theory prematurely. Start with:

```text
x[k+1] = A[k] x[k] + B[k] u[k]
y[k]   = C[k] x[k] + D[k] u[k]
```

Required utilities:

- state transition product Phi(k,j),
- finite-horizon impulse response,
- monodromy matrix for periodic sequences,
- Floquet multipliers = eigenvalues of annual monodromy matrix.

This will later support seasonal VISIT linearization.

---

## Phase B — Provenance infrastructure

### B1. General source reference schema

Refactor or wrap `visit_source_map.py` so provenance is not VISIT-specific.

Recommended dataclasses:

```python
SourceRef
EquationRef
StateRef
FluxRef
ApproximationRef
```

Minimum metadata:

- repository;
- commit/ref if known;
- path;
- function;
- symbol;
- units;
- source role;
- assumptions;
- approximation level.

### B2. Machine-readable provenance export

Support export to JSON/YAML-like dictionaries so every model adapter can produce a reproducibility manifest.

Desired artifact:

```text
provenance/visit_reduced18.json
```

Generated, not hand-maintained.

### B3. Source commit pinning

Determine and record the exact commit SHA of the source snapshot in `Sachitama2001/VISIT-matrix` used for derivation. Every published analysis should be reproducible against that SHA.

---

## Phase C — VISIT soil subsystem: exact algebraic reconstruction

This is the best first scientific implementation target because it is close to a linear compartment system.

### C1. State vector

```text
x_soil = [
  ltr_tf, ltr_tc, ltr_tr,
  ltr_gf, ltr_gc, ltr_gr,
  msl_a, msl_i, msl_p
]^T
```

### C2. Inputs

Direct litter input:

```text
u_litter = [li_tf, li_tc, li_tr, li_gf, li_gc, li_gr]^T
```

Environmental effective forcing:

```text
f_tm_l,
f_tm_h
```

or active/intermediate/passive-specific humus scalars when that code path is enabled.

### C3. Exact source equations

Implement pure functions reproducing `f_cycle_soil()`:

```text
degrade_litter_i = C_i * sr_i/1000 * f_tm_l
degrade_humus_j  = C_j * sr_j/1000 * f_tm_h_j
```

Litter degradation is split into:

```text
CO2 fraction = f_co2_type
humified fraction = 1 - f_co2_type
```

Humified carbon enters active/intermediate/passive pools according to `f_hm_a/i/p`.

### C4. Matrix form

Construct a matrix satisfying

```text
x_soil[k+1] = A_soil[k] x_soil[k] + B_soil u_litter[k]
```

for the native one-day update, and optionally a continuous approximation

```text
dx/dt = M_soil x + B u.
```

Do not claim continuous exactness where the native code is discrete.

### C5. Tests

- each decomposition term compared against a hand-calculated scalar value;
- each humification path tested separately;
- total soil carbon loss equals heterotrophic respiration when no litter input and no optional exports;
- matrix update equals direct algebraic function for random positive states/parameters;
- environmental scalar changes only intended turnover pathways.

### C6. Temperature response

Implement `frl()` and `frh()` separately from the box update.

Tests:

- known temperatures;
- threshold below -20 C;
- monotonic behavior in the locally relevant regime;
- soil moisture saturation behavior;
- branch behavior for `EX_DECTMP` modes.

Flag suspicious source branches rather than silently fixing them. In particular, any repeated conditional branch that looks like a source bug should be documented and tested as-source before a corrected alternative is proposed.

---

## Phase D — VISIT plant structural-carbon subsystem

### D1. Structural states

Per PFT:

```text
fol, stm, rot
```

Start with tree only, then parameterize C3/C4 through the same class.

### D2. Litterfall

Implement `f_lf`, `f_lc`, `f_lr` exactly, including regime switches.

### D3. Respiration

Implement:

- maintenance respiration as state x temperature-dependent rate;
- growth respiration as allocation-flux-dependent loss.

Keep Q10 acclimation (`f_q10_ar`) and size-dependent specific respiration (`f_spcfc_resp`) as separate coefficient-update functions rather than burying them inside one opaque routine.

### D4. GPP

Implement only after the carbon bookkeeping layer is stable.

Recommended progression:

1. treat GPP/EPP as prescribed effective input;
2. implement source-derived GPP equation;
3. connect temperature/water/CO2/radiation responses;
4. validate against native VISIT outputs if executable data are available.

### D5. Allocation

Allocation is state-dependent and piecewise. Implement as a nonlinear pure function returning allocation fractions and fluxes, with explicit regime labels.

Tests should cover:

- EPP > 0 and EPP <= 0;
- LAI above/below optimum;
- dormant vs growing season;
- crop-stage branch;
- conservation of allocated carbon before growth-respiration subtraction.

### D6. NSC decision

Determine whether `nsch_storage` must enter the minimal Markov state for the analysis target.

Required experiment:

- construct two native states with identical 18 structural pools but different NSC;
- show whether their next-day structural-carbon updates differ under identical forcing.

If yes, document formally that the 18-pool state is non-Markov for native plant dynamics.

---

## Phase E — Assemble reduced VISIT carbon adapter

Recommended interface:

```python
class VISITReduced18Adapter(TerrestrialCarbonAdapter):
    state_schema()
    input_schema()
    output_schema()
    provenance()
    rhs_or_step(...)
    linearize(...)
    qse(...)
```

Key metadata:

```text
model_name = "VISIT"
approximation = "reduced18"
time_domain = "discrete_daily" or explicit continuous approximation
```

The adapter must expose source-grounded matrices/functions rather than import assumptions from the old `visit_matrix` package without verification.

---

## Phase F — QSE and local linearization

### F1. QSE solver

For a continuous autonomous reduced model:

```text
f(x_star, z) = 0
```

For the discrete daily map:

```text
F(x_star, z) - x_star = 0.
```

Implement both and label them separately.

### F2. Jacobians

Support:

- analytical Jacobian where easy;
- finite differences as reference implementation;
- optional automatic differentiation only if architecture permits without rewriting source-grounded equations unnaturally.

Compute:

```text
A = df/dx
B = df/dz
C = dh/dx
D = dh/dz
```

or discrete equivalents.

Tests:

- analytical and finite-difference Jacobians agree for soil subsystem;
- perturbation prediction error scales quadratically with perturbation size where smooth.

### F3. QSE sensitivity

From implicit differentiation:

```text
dx_star/dz = -A^{-1} B
```

for continuous systems, with discrete analogue based on `I - A_d`.

Add direct numerical comparison to finite perturbations of solved QSE.

---

## Phase G — NEP/IRF/transfer-function theory

### G1. Output operators

Implement explicit output channels:

```text
GPP
NPP
Rh
NEP
Ctot
```

with provenance.

### G2. Core identities

For each locally linear stable model, test identities connecting:

- forcing -> total carbon transfer function;
- forcing -> NEP transfer function;
- forcing-rate -> NEP low-frequency gain;
- QSE sensitivity.

Do not hard-code identities unless the output definition and carbon closure assumptions actually justify them.

### G3. Modal analysis

For each input-output channel compute:

- poles;
- modal time scales;
- residues / modal amplitudes;
- sign of each residue;
- pole-zero cancellation indicators.

Scientific target: detect source/sink reversals generated by competition between modes of opposite-signed residues.

### G4. Frequency analysis

Produce Bode-like data without requiring plotting in the core library:

```text
omega
magnitude
phase
```

Plotting should be a separate example/experiment layer.

### G5. Reduced-order models

Only after full local systems are validated:

- residue-based mode pruning;
- balanced truncation;
- error metrics in NEP and total carbon channels.

---

## Phase H — Native VISIT perturbation bridge

Goal: compare local IRFs derived from reduced/Jacobian models against perturbations of the original C model.

Possible approaches, in preferred order:

1. compile native VISIT and drive a controlled site experiment;
2. add a minimal wrapper/output mode that exposes state and flux trajectories;
3. perturb one forcing component or effective driver by a small impulse/step;
4. subtract baseline trajectory;
5. compare with linear IRF prediction.

Do not modify the scientific source more than necessary. Wrapper/instrumentation code should be separate and documented.

Validation metrics:

- normalized RMSE by state/output;
- error vs perturbation amplitude;
- horizon-dependent error;
- mode/time-scale agreement.

---

## Phase I — Seasonal/periodic dynamics

A constant-QSE picture may be inadequate for deciduous/seasonal systems.

Construct a periodic reference trajectory under repeated climatology and linearize daily:

```text
A_1,...,A_365
B_1,...,B_365
```

Then evaluate:

- annual monodromy matrix;
- Floquet multipliers;
- phase-dependent IRFs;
- whether annual mean LTI loses important source/sink sign changes.

This is likely a strong scientific extension rather than merely numerical refinement.

---

## Phase J — Forcing-rate and source/sink criteria

Once QSE sensitivity and IRF are stable:

1. define moving QSE `x_star(z(t))`;
2. define disequilibrium `p = x_star - x`;
3. derive local forcing-rate response;
4. identify sufficient/necessary local conditions for NEP sign;
5. compare analytical boundaries against simulation;
6. map sensitivity of the boundary to parameters.

Potential outputs:

```text
critical dT/dt
critical dP/dt
critical multidimensional forcing direction
parameter derivative of critical rate
```

This phase is the bridge back to rate-induced tipping / QSE-tracking questions.

---

## Phase K — Cross-model abstraction

Only after VISIT workflow works end-to-end, define the common adapter protocol around what was actually needed.

Likely concepts:

```python
StateSchema
ForcingSchema
OutputSchema
ModelAdapter
OperatingPoint
Linearization
IRFResult
QSEReport
ProvenanceManifest
```

Target future adapters may include CABLE, CLM/CTSM, ORCHIDEE, LPJ-family, CARDAMOM-like box systems, etc., but no adapter should be designed solely from assumed model structures.

Cross-model comparison should allow different state dimensions while comparing common quantities:

- total-carbon transfer function;
- NEP transfer function;
- dominant poles/time scales;
- DC/QSE sensitivities;
- mode participation by broad pool class;
- source/sink reversal conditions.

---

## Milestones and acceptance criteria

### M0 — repository engineering baseline

- tests run with one command;
- lint/type-check strategy documented;
- package imports cleanly;
- continuous and discrete LTI basics tested.

### M1 — VISIT reduced carbon model

- source-grounded 18-pool model;
- soil exact algebra tests;
- reduced plant bookkeeping;
- NEP output;
- QSE and local IRF.

### M2 — native-vs-reduced validation

- native VISIT perturbation experiments;
- quantitative IRF comparison;
- documented failure modes of 18-pool reduction.

### M3 — moving QSE/source-sink analysis

- QSE sensitivities;
- forcing-rate transfer function;
- analytical NEP sign/reversal conditions;
- simulation validation.

### M4 — seasonal/LTV analysis

- periodic linearization;
- Floquet modes;
- phase-dependent IRFs.

### M5 — second ecosystem-model adapter

- same analysis pipeline works for another terrestrial ecosystem model;
- cross-model comparison demonstrates that the framework is not VISIT-specific.

---

## Priority order for the next coding agent

If only one continuous work session is available, do this sequence:

1. implement discrete LTI core + tests;
2. refactor provenance into generic schema;
3. implement VISIT 9-pool soil subsystem exactly;
4. add matrix-vs-direct-algebra tests;
5. implement QSE for soil subsystem;
6. compute soil IRF and temperature-scalar sensitivity;
7. only then proceed to plant structural carbon.

This sequence yields a scientifically defensible end-to-end vertical slice early.