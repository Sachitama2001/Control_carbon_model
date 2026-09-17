# Carbon-water matrix research plan

## One-sentence objective

陸域炭素循環の行列アプローチに水の貯留・輸送を動的状態として結合し、水理過程が炭素の貯留容量、非平衡性、応答時間を変える経路を、少数の状態変数と共通診断量で説明する。

## Main claim and non-claims

The scientific contribution should be the **diagnosis of coupled carbon-water
dynamics**, not merely the statement that water can be written as differential
equations. SPAC and plant-hydraulic models already represent water storage and
transport explicitly, and matrix land models already include soil moisture as
an environmental modifier.

The defensible target claim is:

> A water-storage and transport subsystem can be coupled to the land-carbon
> matrix framework so that productivity effects, residence-time effects,
> hydraulic memory, and transient disequilibrium can be separated and tested.

Do not claim, without a systematic review, that this is the first matrix
representation of water cycling. Do not make R-tipping a required outcome.
Tipping analysis is a later diagnostic if multiple attractors or loss of
tracking can first be established.

## Research questions

### RQ1 - Representation

Can an existing land-model water balance be rewritten with explicit stores,
internal flux topology, and a closed mass budget without changing its native
processes or update order?

### RQ2 - Reduction

Which water stores must remain dynamic to reproduce the carbon response of the
reference model? In particular, does finite plant water storage add behavior
that a soil-moisture modifier or quasi-steady hydraulic closure cannot retain?

### RQ3 - Non-equilibrium diagnosis

How does water change total carbon storage capacity
\(X_c=\mu\tau_E\), current stock \(X\), and signed storage potential
\(X_p=X_c-X\)? Does drought primarily change carbon input \(\mu\), ecosystem
transit time \(\tau_E\), or their interaction?

### RQ4 - Response times and memory

Which coupled modes govern recovery after precipitation or atmospheric-demand
perturbations? Over what forcing periods does finite water storage alter phase,
amplitude, or recovery relative to a quasi-steady water approximation?

### RQ5 - Strong nonlinear response (conditional)

After a normal tracking target is defined, can a sufficiently rapid forcing
change move the coupled trajectory into a different basin or persistent
attractor? If not, report the nonlinear transient range without calling it
R-tipping.

## Two model tracks that must remain separate

| Track | Purpose | Time semantics | Current implementation |
|---|---|---|---|
| Faithful VISITc audit | Reproduce and understand native water/carbon process order | Sequential daily path in inspected source; executable path still to be pinned | `visitc_source_map.py`, `visitc_hydrology.py`, source ledger |
| Four-component synthesis | Derive and test low-dimensional carbon-water theory | Continuous ODE in days | `coupled_carbon_water.py`, equations document, idealized experiment |

The faithful track validates representation. The synthesis track tests theory.
A mismatch between them cannot be attributed to “matrix form” until differences
in process content, aggregation, and time discretization are controlled.

## Minimal state and scope

Begin with four shared components:

\[
\boldsymbol C=(C_{leaf},C_{stem},C_{root},C_{soil})^\top,
\qquad
\boldsymbol W=(W_{leaf},W_{stem},W_{root},W_{soil})^\top.
\]

Initial restrictions:

- one vegetation type and one idealized site;
- no snow, fire, harvest, land-use change, or crop state in the reduced model;
- no dynamic N or P state in the first comparison;
- one aggregated soil-carbon pool and one soil-water store;
- idealized forcing before global meteorology;
- every omitted process listed rather than absorbed into an unnamed parameter.

Add a state only when one of these criteria is met:

1. a required mass balance cannot be expressed;
2. the reference model contains memory that materially changes the target
   response;
3. aggregation fails a predeclared validation metric;
4. the state distinguishes a mechanism central to the paper's question.

Likely first additions are a second soil-water layer, interception/snow storage,
or a litter/soil carbon split - not more plant organs by default.

## Hypotheses

These are alternatives to test, not conclusions:

- H1: Water limitation lowers \(\mu\) and therefore reduces the increase in
  carbon storage capacity under otherwise favorable forcing.
- H2: The same drying also suppresses decomposition and lengthens
  \(\tau_E\), opposing the input effect.
- H3: Finite plant/soil water stores create phase lag and hysteresis-like
  transient loops even when the frozen system has a unique equilibrium.
- H4: A quasi-steady elimination of water is accurate only when water modes are
  sufficiently faster than carbon modes and forcing.
- H5: Large or fast drying can invalidate the local linear and quasi-steady
  approximations; this does not by itself imply tipping.

## Work packages and acceptance criteria

### WP0 - Provenance and reproducibility

- Pin both repositories and never call them interchangeable VISIT versions.
- Record source file, function, variable names, units, branch assumptions, and
  approximation level for every derived equation.
- Keep literature PDFs in `docs/` when redistribution permits; otherwise store
  DOI/publisher/PDF links in `docs/literature_index.md`.
- Keep generated experiment data out of git unless deliberately selected.

Done when another agent can identify the authority for every implemented
carbon-water relation without consulting chat history.

### WP1 - Faithful VISITc water ledger

1. Inventory native stores, diagnostics, fluxes, and call order.
2. Transcribe store-update algebra with potential evaporation/transpiration as
   effective inputs.
3. Build a minimal C bridge for exact one-day comparison.
4. Audit water balance, clipping, and the apparent baseflow double subtraction.
5. Trace water-to-carbon and carbon-to-water dependencies.

Done when Python and native C agree for all exposed one-day outputs across
normal and edge cases, and every residual has a named source branch.

### WP2 - Four-component coupled equations

1. Prove carbon and water budget identities.
2. Establish nonnegative-domain behavior and solver checks.
3. Solve a frozen-forcing coupled equilibrium.
4. Compute instantaneous carbon capacity separately.
5. Partition the Jacobian and test the water Schur complement.

The first executable version of these items is now in the repository. Soil
retention and plant pressure-volume relations are source-labelled; their
parameters and the remaining conductance/external-flux closures still require calibration.

### WP3 - Constitutive grounding and calibration

- Calibrate and sensitivity-test the implemented soil water-retention and
  plant pressure-volume relations.
- Select conductance/vulnerability and stomatal response formulations from a
  declared model family (e.g. SurEau-Ecos, FETCH, CLM5 hydraulics).
- Reconcile water units and area bases explicitly.
- Estimate or bound parameters from literature or native VISITc settings.

Done when each parameter has a unit, source/range, and sensitivity plan, and
the system passes dimensional checks.

### WP4 - Controlled reduction experiments

Compare in this order:

| Comparison | Changes held fixed | Scientific question |
|---|---|---|
| Native source vs faithful transcription | All process equations and order | Is representation correct? |
| Full water stores vs aggregated stores | Carbon equations and forcing | What memory is lost by aggregation? |
| Dynamic water vs quasi-steady water | Constitutive laws and parameters | When does finite storage matter? |
| Coupled vs one-way water->carbon | External forcing and initial state | Does carbon feedback materially alter water/carbon response? |
| Nonlinear vs local tangent model | Operating point and forcing protocol | Where does linear diagnosis fail? |

Do not compare a reduced continuous model directly with native VISITc and label
the entire difference “the effect of water.”

### WP5 - Idealized forcing experiments

Use forcing that maps to a specific question:

- small precipitation or potential-transpiration steps for local modes;
- periodic forcing over a frequency sweep for attenuation and phase;
- equal-total-rainfall pulse trains for storage/memory effects;
- dry-down and re-wetting for signed capacity and recovery;
- ramps with identical endpoints but different rates only after cumulative
  exposure is also controlled.

Report state trajectories, every external flux, budget residuals,
\(X, X_c, X_p\), the three \(\Delta X_c\) terms, eigenvalues/response times,
and numerical refinement error.

### WP6 - Periodic reference and conditional tipping tests

For seasonal forcing, compute or approximate a stable periodic reference orbit
rather than using an annual-mean fixed point. Test phase-dependent perturbation
responses before invoking rate-induced tipping from a periodic attractor.

Only use “R-tipping” when:

1. a family of quasistatic attractors/reference trajectories is defined;
2. tracking and tipping outcomes have operational criteria;
3. the outcome changes with forcing rate while endpoints are fixed;
4. basin or connecting-orbit evidence distinguishes tipping from slow recovery;
5. time-step and protocol refinements preserve the threshold.

## Experiment ladder

| Stage | Input | Model pair | Primary diagnostic | Pass condition |
|---|---|---|---|---|
| A | Constant forcing | faithful source/Python | one-day stores and fluxes | numerical agreement and explained budget |
| B | Small step | dynamic vs tangent | Jacobian/IRF | error scales quadratically with perturbation |
| C | Frequency sweep | dynamic vs quasi-steady water | gain and phase | crossover tied to water eigenmodes |
| D | Dry-down/re-wet | coupled vs one-way | \(X,X_c,X_p\) decomposition | budget and refinement checks pass |
| E | Equal-total pulse trains | alternate storage resolutions | memory and recovery | differences survive matched totals |
| F | Rate ramps | nonlinear attractor tests | tracking/basin outcome | tipping criteria above satisfied |

## Current checkpoint

Implemented:

- pinned source map for current `visit-manager/VISITc`;
- native-validated source-order hydrology, atmospheric, radiation,
  canopy-conductance, and Penman-Monteith transcription;
- executable eight-state coupled ODE and flux diagnostics;
- instantaneous carbon capacity and Wei-style finite-change decomposition;
- coupled equilibrium, Jacobian blocks, and Schur complement;
- dry-down/re-wetting example with numerical refinement;
- tests of budgets, provenance, coupling, equilibrium, and integration.

Next highest-value actions:

1. implement dynamic-vs-quasi-steady water experiments under matched processes;
2. connect the PM chain to controlled meteorological forcing;
3. calibrate/bound plant saturated water, pressure-volume, and conductance parameters;
4. test alternate compartment resolution and carbon-dependent allometry;
5. only then expand toward seasonal or tipping analyses.
