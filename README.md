# Control Carbon Model

A research project for analyzing terrestrial ecosystem carbon-cycle models using control theory, impulse response functions (IRFs), transfer functions, quasi-static equilibria (QSEs), and source-sink diagnostics.

The long-term goal is **not to reproduce one specific land model**, but to build a model-agnostic framework that can compare many terrestrial ecosystem models through common state-space and input-output dynamical quantities.

## Current first testbed

VISIT is the first source-grounded implementation target. The authoritative source snapshot currently used for derivation is in `Sachitama2001/VISIT-matrix/visit_local`.

The existing Python matrix implementation in `VISIT-matrix/visit_matrix` is useful prior work, but the VISIT C source is treated as the authority for state definitions, update order, process equations, and provenance.

## Start here

For coding agents and new contributors:

1. `AGENTS.md` — concise coding-agent rules and first tasks.
2. `HANDOFF.md` — full scientific and implementation handoff.
3. `docs/implementation_roadmap.md` — phased engineering/research roadmap and acceptance criteria.
4. `docs/architecture_and_conventions.md` — package boundaries, time/sign/unit conventions, testing rules.
5. `docs/visit_state_space_source_map.md` — source-grounded VISIT state/input/output decomposition.
6. `docs/research_questions.md` — scientific questions and planned experiment matrix.

## Initial goals

- Develop model-agnostic continuous, discrete, and eventually periodic/LTV state-space tools.
- Connect QSE disequilibrium to NEP/NEE/NECB source-sink behavior with explicit sign conventions.
- Use IRFs, transfer functions, modal analysis, and frequency response to characterize carbon-cycle dynamics.
- Derive forcing-rate and QSE-sensitivity relations that can support analytical source/sink criteria.
- Build provenance-preserving adapters for process-based terrestrial ecosystem models, beginning with VISIT and later extending to additional models.
- Compare full simulation with IRF/reduced-order approaches for speed, interpretability, and validity range.

## Current code

`src/control_carbon/state_space.py` contains the first generic continuous-time LTI/IRF utilities.

`src/control_carbon/visit_source_map.py` contains the first machine-readable VISIT provenance/state map.

The next high-priority implementation is a discrete-time systems core followed by an exact source-derived reconstruction of the 9-pool VISIT soil carbon subsystem.
