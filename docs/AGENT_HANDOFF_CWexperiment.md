# Coding-agent handoff — SPAC–carbon multistability master's research

**Context / goal.** Extend the checked, synthetic, **one-way** plant hydraulics + 4-pool carbon PoC in this directory into (i) a scientifically consistent, dimensionally explicit Leuning–photosynthesis–water balance closure; (ii) two-way carbon structural feedback; (iii) nonlinear equilibrium continuation and stability analyses; and (iv) rate-induced tipping tests. Distinguish mathematical results, computational numerics, and empirical claims. Do not overclaim prior novelty.

## First actions

1. Read `README.md`, `model.py`, `run_experiments.py`, `test_model.py` and the accompanying Japanese research plan DOCX.
2. Run `pytest -q` and `python run_experiments.py`; confirm existing outputs. Commit this baseline without changing it.
3. Keep **model variants** separated rather than overwriting the Manzoni-style baseline. Suggested files: `model_leuning.py`, `model_two_way.py`, `continuation.py`, `experiments_leuning.py`, `test_leuning.py`, `test_continuation.py`.
4. Add parameter files (`synthetic_parameters.yaml` and later site-calibrated cases) with explicit units and source provenance.

## Primary scientific questions

- Under what mathematical conditions does the frozen-time carbon capacity become set-valued when the water subsystem is multistable?
- When is zero plant-carbon equilibrium locally stable, when does a positive plant-carbon equilibrium exist, and does the complete system (including soil carbon and hydraulic state) really have two **isolated stable** equilibria?
- Does two-way leaf/root/wood allocation alter the **location** of a hydraulic fold (where the one-way model cannot)?
- For smooth, physically regularized models, can temporal forcing induce basin escape at rates where all frozen-time equilibria of interest persist (genuine R-tipping), as distinct from B-tipping with delayed departure?

## Literature: read the originals

1. Luo et al. 2017, `10.5194/bg-14-145-2017` (moving carbon storage capacity and potential).
2. Luo et al. 2022, `10.1029/2022MS003008` (matrix approach; state of art; environmental modifiers).
3. Sierra & Müller 2015, `10.1890/15-0361.1`; Müller & Sierra 2017, `10.1007/s12080-017-0342-3` (reservoir stability and nonlinear ISS).
4. Manzoni et al. 2014, `10.1002/2013WR015236` (equations 1–5; piecewise linear hydraulic equilibrium; fold).
5. Leuning 1995, `10.1111/j.1365-3040.1995.tb00370.x`; Tuzet, Perrier & Leuning 2003, `10.1046/j.1365-3040.2003.01035.x` (already couples stomata, leaf psi, photosynthesis and SPAC; central prior art).
6. Parolari et al. 2014, `10.1002/2013JG002592`; Sperry et al. 2017, `10.1111/pce.12852`; Johnson et al. 2022 `10.1111/pce.14327`; Ashwin et al. 2012 `10.1098/rsta.2011.0306`.
7. FETCH3: `10.5194/gmd-15-2619-2022`. Consider root hydraulic network incidence-matrix precedents.

## Model A: existing baseline; never alter without regression tests

State `(psi_x,psi_l, C_leaf,C_wood,C_root,C_soil)`, all rates and conductances fully specified in `model.py`.

- Piecewise Manzoni-conductance: `gx= gxmax*max(0,1-psi_l/(2*p50x))`, analogous `gs`.
- Carbon `dC/dt=b*Pmax*h(psi_l)*(1-exp(-kappa*LAI)) + A*C`, `LAI=SLA*C_leaf/f_C`.
- `A` is a compartmental matrix: off-diagonal nonnegative, column sums <= 0, negative eigenvalue real parts; carefully distinguish GPP/NPP/respiration to avoid double counting.
- Frozen-state carbon threshold: `R0=b_leaf*Pmax*kappa*SLA*h/(k_leaf*f_C)`.
- If `R0<1`, all-carbon zero is locally stable **for frozen water**, if `R0>1`, zero is unstable and a unique positive leaf equilibrium exists.
- The water-to-carbon Jacobian is lower block triangular, so the hydraulic fold does not move.

## Model B: Leuning + C3 photosynthesis + hydraulic feedback

1. Decide explicit convention: Leuning1995 **CO2** stomatal conductance `g_sc=g0+a1*A_n/[(C_s-Gamma)*(1+D_s/D0)]` (molar units). Stomatal H2O conductance is approximately `1.6*g_sc`. Check water vapour, CO2 and VPD units and pressure normalization in `E=g_sw*(D_s/P_atm)`; do not pretend the old synthetic `gs*D` shares molar units.
2. Leuning A_n and g_sc must be solved **simultaneously**: specify FvCB `A_n(C_i,I,T)` or a clearly labeled surrogate and use `A_n=g_sc*(C_s-C_i)` in consistent molar units. Solve for C_i in physically allowed interval with Brent root solver; reject negative or nonphysical solutions.
3. To enable hydraulic closure under drought, use a separate water stress function linked to psi_l, e.g. Tuzet-style response; justify whether it multiplies a1, photosynthesis capacity, or g_sc; avoid double-applying the same decline as both Leuning and Manzoni g_s. Option 1: replace Manzoni gs with a consistent Leuning–Tuzet closure. Option 2: retain Manzoni gs as a baseline only and use Leuning to compute carbon assimilation, clearly labelling any partial coupling.
4. `g0>0` denotes residual **stomatal** conductance. Cuticular leakage must be separately handled and may destroy exact `A_n=0` and `E=0` equilibria; avoid assuming `gs=0` on a dead branch if g0 positive.
5. `LAI=(SLA/f_C)*C_leaf` only if C_leaf is carbon mass per **ground area** and SLA is m² leaf/kg **dry mass**. Add a sun/shade canopy treatment if upgrading beyond simple Beer–Lambert radiation saturation; the integrated canopy flux is not generally `LAI * A_leaf` when light is vertically varying.
6. If C_leaf -> 0, the notion of plant leaf psi and fixed leaf capacitance becomes ill-posed; regularize via minimal viable biomass, a vegetation occupancy variable, or explicit model regime switch. Be explicit about biological meaning.
7. Make xylem hydraulic conductance depend on wood allocation only using an independently sourced scaling/allometry, and soil-root conductance depend on root allocation; avoid assigning unrealistic dimensions.

## Model C: equilbria, continuation, eigenanalysis

1. For fixed `(D,T,psi_r,I,Ca)`, solve `F(y;mu)=0` from multiple initial guesses. Detect all **isolated** physically admissible equilibria, classification via eigvals of Jacobian (symbolic/autodiff or checked finite differences).
2. Track branches with pseudo-arclength continuation (scipy custom or AUTO-07p / PyDSTool / BifurcationKit in Julia), especially folds and local extinction boundary bifurcations. Do not infer number of equilibria from solver single convergence.
3. Verify Jacobian with central finite differences and optionally JAX. Check conservation and positivity.
4. Define linear frozen carbon capacity `-A(psi,T)^(-1)u(psi,D,T)` only for affine carbon subsystem with Hurwitz A. For nonlinear assimilation/allometry, report **set of admissible stable roots** as a generalized instantaneous capacity, not a single vector automatically.
5. At `gx=0` original dry equilibrium is typically not isolated and xylem water potential can freeze. A zero eigenvalue exists from degeneracy, *not necessarily from a fold*. Rework the dry branch with a physiologically justified smooth residual transport/recovery if rigorous full-state stability desired.
6. Decompose eigenvectors into water/carbon coordinates and quantify any slow eigenmodes. For a block-triangular one-way case, eigenvalue union is exact. In two-way case, analyze Schur complement / singular perturbations.

## Model D: climate forcing and tipping classification

- Start with equilibrium continuation vs constant D; add T independently and with physically consistent atmospheric humidity relationship. Keep radiation, psi_r, CO2 separate unless explicitly coupled.
- For B-tipping, parametric variation crosses a fold of the frozen system; identify its location and any finite-rate delay.
- For **genuine R-tipping**, design trajectories with endpoints and entire parameter path retaining a feasible attracting frozen branch. Demonstrate existence of a rate threshold for basin escape; run slow vs fast forcing with the SAME parameter endpoints. Do not mislabel fold-crossing runs as R-tipping.
- For `x_c(t)` moving equilibrium, compare `x(t)` vs `x_c(t)` and quantify lag; note `x_c(t)` is not generally a solution of the nonautonomous ODE.
- Add temperature only after base equations are validated; choose measured or source-backed functions for Vcmax, Jmax, maintenance respiration, tissue turnover, hydraulic vulnerability, VPD.

## Recommended deliverables and acceptance tests

1. `tests`: Baseline 7 tests plus >15 new tests covering C_i gas-exchange closure, stomata/CO2 H2O conversion, flow conservation, Jacobian finite differences, physical admissibility, stable/unstable classification, zero state, pseudo-arclength step convergence, reproducibility, R-tipping classification.
2. `figures`: at least hydraulic + carbon bifurcation plots, phase portraits/basins, conditional capacity vs actual stocks, VPD/temperature sensitivity, slow-vs-fast *non-bifurcation-crossing* forcing.
3. `results/`: CSV with parameter units, root residuals, Jacobian eigenvalues, stability, model configuration hash.
4. `docs`: equations with parameter units and literature provenance, a corrected derivation of zero equilibrium threshold, an explicit statement of the dry-branch degeneracy, and a 'not yet validated' section.
5. `run` CLI: `python run_experiments.py`, `pytest -q`, plus separate `python experiments_leuning.py`.

## Non-negotiable scientific guardrails

- Land absorbs a **fraction**, not a majority, of anthropogenic CO2; Global Carbon Budget 2025 puts the decadal terrestrial fraction at ~21% (methodology revised).
- `x=0` is **not** an equilibrium of forced affine C models with nonzero u.
- Failure to prove ISS is NOT proof that trajectories do not track; ISS applies to nonlinear systems too.
- A moving capacity is a **frozen-time** root unless shown otherwise.
- A saddle-water branch may have a conditional carbon equilibrium but is NOT a stable attractor of the full water-carbon model.
- The Leuning model and Manzoni stomatal relation are alternatives/closures requiring reconciliation; do not multiply without a stated physical interpretation.
- Near zero leaf biomass the hydraulic leaf state and capacitance assumptions may fail.
- Avoid empirical claims before parameter calibration and independent validation.
