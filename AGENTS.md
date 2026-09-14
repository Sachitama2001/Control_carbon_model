# AGENTS.md

## Mission

Continue this repository as a rigorous research-code project for adding
explicit water storage and SPAC transport to matrix land-carbon models, then
using QSE, disequilibrium, IRF, modal, and forcing-rate analysis on the coupled
system.

VISIT is the first validation target, not the final scope.

The active scientific question is how dynamic water changes carbon storage
capacity, signed disequilibrium, and response time. R-tipping is conditional
and downstream; do not tune the model merely to produce tipping.

## Read first

Before carbon-water coding, read in this order:

1. `docs/carbon_water_research_plan.md`
2. `docs/coupled_carbon_water_equations.md`
3. `docs/visitc_carbon_water_source_ledger.md`
4. `docs/literature_index.md`
5. `HANDOFF.md`
6. `docs/current_status.md`
7. `docs/architecture_and_conventions.md`

## Two VISIT source boundaries

The existing nine-pool soil-carbon adapter uses:

- repository: `Sachitama2001/VISIT-matrix`
- commit: `3285bd8e131a932e338b59892751648fd9edcc7b`
- source directory: `visit_local/`

The Python code under `VISIT-matrix/visit_matrix/` is useful prior work but is not authoritative when it disagrees with, simplifies, or omits behavior in the C source.

The active water audit uses:

- repository: `visit-manager/VISITc`
- commit: `5202debd96df6f88beb7d61f8688fff02ace964a`
- source directory: `point/`

Never merge assumptions or coefficients across these snapshots without an
explicit comparison. Use `visit_source_map.py` for the old snapshot and
`visitc_source_map.py` for the current one.

## Non-negotiable provenance rule

Every VISIT-derived equation, coefficient, state mapping, or output mapping
must identify repository, commit, C-source path, function, units, update order,
and whether the implementation is exact or approximate.

Do not add undocumented matrix coefficients.

## Approximation levels

Always label model objects/results as one of:

- native nonlinear discrete VISIT map;
- source-grounded reduced compartment model;
- local linearization;
- periodic/LTV linearization;
- model-agnostic derived analysis.

Never label a reduced 18-pool model as simply "the VISIT model" without qualification.

## Active implementation order

The old nine-pool soil target is implemented. Continue in this order:

1. validate `visitc_hydrology.py` against a minimal native C bridge;
2. transcribe `pm_incep`, `pm_evap`, and `pm_transp` dependencies;
3. finish the bidirectional carbon-water dependency ledger;
4. replace placeholder normalized storage in the eight-state model with
   sourced soil retention and plant pressure-volume relations;
5. compare dynamic and quasi-steady water with matched forcing/processes;
6. build a periodic reference trajectory before any seasonal R-tipping work.

Do not expand pool count until a mass-balance, memory, validation, or research
criterion requires it.

## Native and reduced time semantics

The inspected VISITc path includes a sequential daily update, but the selected
executable and compile flags must be traced before declaring one universal
native timestep. `coupled_carbon_water.py` is a separate continuous model in
days. Do not describe it as a native VISITc ODE.

## State caution

The 18 structural carbon pools are a reduced state. `nsch_storage`, phenology/hydrology memory, and other variables may be required for a Markov representation of native VISIT. Test this; do not assume it away.

## Output caution

Keep NEP distinct from NECB and from observational NEE sign conventions.

## Testing requirement

For every source-derived process:

1. pure function reproducing source algebra;
2. hand-computable unit test;
3. mass-balance/structural test where applicable;
4. provenance metadata;
5. matrix or linearized integration only afterward.

## Carbon-water testing requirement

For every source-derived water process:

1. reproduce source order in a pure function;
2. add a hand-computable test;
3. test water mass balance and any clipping correction;
4. compare with native C before calling the transcription exact;
5. record provenance and approximation level;
6. only then insert it into a reduced or linearized model.

For every reduced coupled equation, test the carbon and water budgets,
nonnegative domain, equilibrium residual, both Jacobian coupling directions,
and numerical refinement. Keep instantaneous carbon capacity distinct from the
coupled equilibrium.

## Definition of scientific success

The current milestone succeeds when a faithful VISITc water audit and a
source-grounded low-dimensional carbon-water model can be compared through
storage capacity, signed disequilibrium, modes, IRFs, and hydraulic memory,
while every result retains traceability and mass balance.
