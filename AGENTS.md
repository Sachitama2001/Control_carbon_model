# AGENTS.md

## Mission

Continue this repository as a rigorous research-code project for applying state-space, IRF, transfer-function, QSE, and forcing-rate analysis to terrestrial ecosystem carbon-cycle models.

VISIT is the first validation target, not the final scope.

## Read first

Before coding, read in this order:

1. `HANDOFF.md`
2. `docs/implementation_roadmap.md`
3. `docs/architecture_and_conventions.md`
4. `docs/visit_state_space_source_map.md`
5. `docs/research_questions.md`

## Source of truth for VISIT

Use:

- repository: `Sachitama2001/VISIT-matrix`
- source directory: `visit_local/`

The Python code under `VISIT-matrix/visit_matrix/` is useful prior work but is not authoritative when it disagrees with, simplifies, or omits behavior in the C source.

## Non-negotiable provenance rule

Every VISIT-derived equation, coefficient, state mapping, or output mapping must identify its C-source path and function and state whether the implementation is exact or approximate.

Do not add undocumented matrix coefficients.

## Approximation levels

Always label model objects/results as one of:

- native nonlinear discrete VISIT map;
- source-grounded reduced compartment model;
- local linearization;
- periodic/LTV linearization;
- model-agnostic derived analysis.

Never label a reduced 18-pool model as simply "the VISIT model" without qualification.

## First implementation target

Implement the 9-pool VISIT soil carbon subsystem directly from:

- `visit_local/soil_proc.c::f_cycle_soil`
- `visit_local/decomposition.c::frl`
- `visit_local/decomposition.c::frh`

Add direct-algebra vs matrix-equivalence tests before expanding the plant subsystem.

## Native time semantics

VISIT's native update is daily and sequential. Build discrete-time tools before treating continuous-time LTI notation as the primary representation.

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

## Suggested first coding session

1. add `DiscreteLTI` + tests;
2. generalize provenance dataclasses;
3. implement exact daily VISIT soil subsystem;
4. test matrix equivalence;
5. add soil fixed point/QSE;
6. compute soil IRFs;
7. then begin plant structural-carbon reconstruction.

## Definition of scientific success

The project succeeds when the same pipeline can compare multiple terrestrial ecosystem models through input-output dynamical quantities such as poles, residues, IRFs, QSE sensitivities, NEP transfer functions, and forcing-rate source/sink boundaries while retaining traceability to each native model's equations.
