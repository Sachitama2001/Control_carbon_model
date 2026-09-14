# Control Carbon Model

A research project for analyzing terrestrial ecosystem carbon-cycle models using control theory, impulse response functions (IRFs), transfer functions, quasi-static equilibria (QSEs), and source-sink diagnostics.

The current focus is deterministic nonlinear state-space dynamics: which
state feedbacks create multiple equilibrium branches and B/R-tipping?
Only selected environmental drivers vary; stochastic mechanisms are deferred.
See the [active plan](docs/deterministic_tipping_plan.md) for equations,
assumptions, implemented benchmarks, and the next source-grounding steps.

The long-term goal is **not to reproduce one specific land model**, but to build a model-agnostic framework that can compare many terrestrial ecosystem models through common state-space and input-output dynamical quantities.

## Current first testbed

VISIT is the first source-grounded implementation target. The authoritative source snapshot currently used for derivation is in `Sachitama2001/VISIT-matrix/visit_local`.

The existing Python matrix implementation in `VISIT-matrix/visit_matrix` is useful prior work, but the VISIT C source is treated as the authority for state definitions, update order, process equations, and provenance.

## Start here

For coding agents and new contributors:

1. `AGENTS.md` — concise coding-agent rules and first tasks.
2. `HANDOFF.md` — full scientific and implementation handoff.
3. `docs/implementation_roadmap.md` — phased engineering/research roadmap and acceptance criteria.
4. `docs/current_status.md` — completed implementation, reproducible checks, and current blockers.
5. `docs/progress_and_next_plan.md` — Japanese progress summary and historical P0-P5 work.
6. `docs/architecture_and_conventions.md` — package boundaries, time/sign/unit conventions, testing rules.
7. `docs/visit_state_space_source_map.md` — source-grounded VISIT state/input/output decomposition.
8. `docs/research_questions.md` — scientific questions and planned experiment matrix.

## Initial goals

- Develop model-agnostic continuous, discrete, and eventually periodic/LTV state-space tools.
- Connect QSE disequilibrium to NEP/NEE/NECB source-sink behavior with explicit sign conventions.
- Use IRFs, transfer functions, modal analysis, and frequency response to characterize carbon-cycle dynamics.
- Derive forcing-rate and QSE-sensitivity relations that can support analytical source/sink criteria.
- Build provenance-preserving adapters for process-based terrestrial ecosystem models, beginning with VISIT and later extending to additional models.
- Compare full simulation with IRF/reduced-order approaches for speed, interpretability, and validity range.

## Current code

The new [nonlinear core](src/control_carbon/nonlinear.py) provides single-driver
continuous integration, frozen-equilibrium diagnostics, and signed compartment
algebra. [Deterministic benchmarks](src/control_carbon/tipping.py) add analytic
saturating-feedback branches, a fold experiment, and finite-ramp basin tests.
These hypothetical models are not native VISIT equations or TKY predictions.

Reproduce the temperature-only VISIT soil control and synthetic B/R figures:

```bash
python -m pip install -e '.[dev,analysis]'
PYTHONPATH=src python examples/plot_deterministic_tipping.py
```

The script writes a PNG, numerical NPZ arrays, and a JSON experiment manifest
under `artifacts/`. No CDS access or native compiler is needed for this example.

`src/control_carbon/state_space.py` contains generic continuous- and discrete-time linear state-space utilities, including impulse responses, forced simulation, and constant-input fixed points.

`src/control_carbon/provenance.py` and `src/control_carbon/visit_source_map.py` provide machine-readable source provenance pinned to the authoritative VISIT source revision.

`src/control_carbon/visit_decomposition.py` reconstructs the VISIT litter and humus temperature/moisture multipliers. The undefined `frh()` branches in the source snapshot are rejected by default; the apparent `stype=0/1/2` correction is available only as an explicitly labelled inferred behavior.

`src/control_carbon/visit_soil.py` implements the source-derived daily 9-pool soil carbon balance, constant-environment fixed points, heterotrophic-respiration impulse responses, and common-condition trajectory data for comparing IRF predictions with direct daily simulation.

`src/control_carbon/visit_native.py` builds a minimal harness against the pinned VISIT `soil_proc.c` and `decomposition.c` sources and compares native C baseline/perturbation trajectories with the Python 9-pool implementation and its IRF. See `docs/native_visit_soil_bridge.md` for scope and exclusions.

Grouped discrete soil modes, basis-independent participation, Rh residues, and implicit fixed-point environmental sensitivities are documented in `docs/visit_soil_modal_qse.md`.

`src/control_carbon/visit_plant.py` begins the P5 plant reconstruction with source-grounded foliage/stem/root turnover, Q10 and size-dependent respiration, and piecewise allocation. Each slice has hand calculations and direct native C checks, but they are not yet integrated into a complete plant daily update.

`src/control_carbon/era5.py` provides point-scale ERA5 acquisition, strict UTC daily aggregation, an intended VISIT `WMODE=1` NCEP/NCAR climate writer, and monthly request construction for direct Copernicus CDS retrieval. See `docs/era5_visit_weather.md` for source-format differences, known compatibility defects in the pinned VISIT snapshot, transformation assumptions, and the tested command-line workflow.

The native soil bridge supports day-varying environments, local environmental
Jacobians, modal analysis, QSE sensitivity, and a native/nonlinear/tangent
comparison figure. Plant turnover, respiration, and allocation slices have
also been individually matched to native C. Their full daily integration
remains pending; immediate research priorities follow the
[deterministic tipping plan](docs/deterministic_tipping_plan.md).
