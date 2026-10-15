# Control Carbon Model — IDE Coding Agent Handoff

> 2026-10-16 rational E4 diagnostic figures: [Japanese figure guide](docs/e4_rational_diagnostics_20261016.md),
> [corrected plots and trajectories](artifacts/spac_carbon/rational_e4_diagnostics_20261016/figures_v4/).
> Follow-up: [equilibrium correction and conditional decay proof](docs/e4_rational_equilibria_and_decay_20261016.md),
> [synchronized trajectory movies](artifacts/spac_carbon/rational_e4_diagnostics_20261016/movies_v4/index.html).
> IMPORTANT: v3's branch plot spuriously joined two disjoint stable intervals.
> E4 has at most one positive equilibrium at fixed VPD, not the two suggested by
> that chord. v4 masks the gaps and colors unstable equilibria red.
> Both forcing plots use 0-26 years; day-scale details resolve the unchanged C1
> 10-day smoothstep. Movies share physical frame times (nonuniform playback).
> Includes frozen-slice leaf/water nullclines; E4 equilibrium/stability branch;
> 10-day versus 10-year bounded 2.44->3.00->2.44 kPa VPD pulses; psi_x/psi_l
> tracks, projected phase paths, and original K=-M^-1u against all four C pools.
> K is state-evaluated frozen-input capacity, not a direct VPD function or
> moving full equilibrium. Fast pulse crosses psi_l=-3 at day 11.094697707.
> Thereafter the invariant dry region gives C_L' <= -0.35/year*C_L and all
> carbon pools/capacities tend to zero conditionally on continuing this ODE.
> Fast integration stops at C_L=1e-6 after ~13.94 years; not mortality.
> Numerical gs=K=0 near psi_l=-4 is not proof of exact finite-time closure.
> Zero-C equilibria form a continuum, not an isolated stable point.
> Slow pulse completes 400-year post-pulse hold and returns within ~0.03% of
> initial stocks for this one path; neither result establishes global basins
> or rate-induced tipping. Capacity identity residual <1.5e-19; source/output
> hashes in results.json/provenance.json. Run
> `python examples/analyze_spac_carbon_rational_e4_diagnostics.py --output NEW_PATH`.

> 2026-10-15 rational-canopy simplification: [Japanese TeX](docs/manzoni_carbon_rational_ja.tex)
> and [25-page PDF](artifacts/papers/manzoni_carbon_rational_20261015/manzoni_carbon_rational_ja.pdf).
> Independent `spac_carbon_rational.py` uses p*h*q*C_L/(1+q*C_L);
> original exponential ODEs, coefficients and structural reference stocks unchanged.
> Positive E0/E2 equilibria are quadratic; E1/E3/E4 are cubic.
> Baseline E0/E1 have one fold; E2/E3/E4 are strictly decreasing
> (outward E4 shape margin >1.21032, NOT a stability certificate).
> E1 fold moves to 1.100569383 kPa; E4 imaginary-pair crossings are
> 2.845318750 and 4.791972651 kPa. No new l1/cycle/global-basin claim.
> Rational RHS/Jacobian accept the formal zero-carbon boundary; simulation
> still stops at C_L_min and adds no mortality continuation.
> Linear canopy is an analytic control: fixed-water isolated positive carbon
> equilibria are lost; this is NOT a blanket claim about full coupled E1-E4.
> 69 focused/regression tests pass; 15 panels, 14 positive roots, 58 trajectories,
> 15 tolerance checks in [verified evidence](artifacts/spac_carbon/rational_20261015/verified/).
> Run `python examples/analyze_spac_carbon_rational.py --output NEW_PATH`.

> 2026-10-15 manuscript revision: [revised PDF and evidence](artifacts/papers/manzoni_carbon_capacity_20261015_revised/)
> repositions linear time-varying capacity versus nonlinear time-varying systems.
> Compartment carbon notation C_L/C_W/C_R/C_S; full six-state ODE, forcing,
> 6x6 Jacobian, and all E0-E4 positive equilibrium vectors are displayed.
> E2/E3 positive curves are globally decreasing under the baseline assumptions:
> one positive conducting equilibrium for every finite VPD, not globally stable.
> New projected full-six-state phase portraits cover 15 fixed-VPD regimes;
> no mortality, global basin or new rate-tipping claim. Original ODEs unchanged.
> Explicit numerical equilibrium vectors/Jacobians are in the manuscript appendix.

> 2026-10-15 Japanese manuscript: [TeX source](docs/manzoni_carbon_capacity_paper_ja.tex)
> consolidates the Manzoni reproduction, carbon coupling, scalar proofs,
> stability/cycles, capacity diagnostics and bounded trajectories.
> Analytic deductions, numerical checks and unresolved outcomes are separated.
> [PDF and validation artifacts](artifacts/papers/manzoni_carbon_capacity_20261015/)
> are built with LuaLaTeX; existing models, parameters and PowerPoint files unchanged.
> H1 hydraulic rerun matches archived numerical results exactly; historical
> source-hash differences for H1/B2 are disclosed, not silently overwritten.

> 2026-10-14 scalar-equilibrium milestone: latest user direction is
> [equilibrium/capacity plan](docs/research_plan_hydraulic_failure_capacity_20261014.md);
> it supersedes earlier "Next" ramp/cycle/global-basin priorities below.
> Read [scalar derivation and proofs](docs/spac_carbon_scalar_equilibrium_results.md).
> Baseline629314706aaac5f2f78ff37f6e19f01a01d58b42; A/B/C ODEs/configs unchanged.
> Independent spac_carbon_equilibria.py eliminates W/R/S/h/psi/x into D_curve(L).
> E0/E1 each have one analytic fold; E4 beta1 is globally decreasing,
> with an outward Decimal constant gap >.22718 (shape, NOT stability proof).
> L=0 excluded; v=0 wet endpoint included; L_min separate from root domain.
> E0/E1 high-L roots stable, low-L unstable; E1 uses transformed Metzler core.
> E4 spectral scan finds .002664752904905 and .004791823126468 MPa crossings.
> Second has L=.00003059249>L_min but is outside routine/25 C forcing.
> Do not hide far-branch restabilization or claim a global stability certificate.
> Second crossing has no new l1/cycle analysis. Near-zero tolerance is not
> proof of a genuinely zero eigenvalue. All positive roots satisfy K=C*,
> including unstable ones; formal zero-carbon continua counted separately.
> gs=0 gives immediate K=0; both closed give frozen water and exp(Mt)C ->0.
> gx-only closed set invariant, baseline R0<=.5625 gives C/K ->0.
> E4 in that set retains positive limiting h if initially h>0.
> Eight fixed-.00244 MPa examples distinguish recovery-neighborhood approach
> from L_min exit, NOT mortality/global-basin classification.
> Exact conductance kink: report flux/K/C, explicitly mark Jacobian undefined.
> 24 new tests;236 related regression passed. 262 roots;15 source hashes.
> Adjacent-float zero-leaf endpoint cases preserve positive roots without
> artificial branch/endpoint discontinuities or absolute-tolerance clipping.
> [Final artifacts](artifacts/spac_carbon/scalar_equilibria_final_v5_20261014/).
> Run `python examples/analyze_spac_carbon_equilibria.py --output NEW_PATH`.
> STOP here per current plan; no ramps/global basins/new mortality without
> further direction. Remaining: E4 continuous stability and empirical assumptions.

> 2026-10-14 WP-C2 / bounded C3 milestone: read
> `docs/spac_carbon_stage_c23_results.md`. Baseline757f359; A/B/C1 unchanged.
> New `spac_carbon_tracking.py`: F_v, y*'=-J^-1F_v, exact/linearized error,
> fixed-scale nonnormality and frozen Lyapunov Q/Q_v. Q_v*v' matters;
> sampled metric drift fails the sufficient contraction test in all six
> big ramps, NOT proof of necessary instability. No uniform certificate.
> E4 smoothstep .00101325 -> .0024 or .00266 MPa, 10day/10year/100year,
> 300year hold. .00266/10day exits L_min around15.78year; other five remain
> positive but fail strict recovery-neighborhood criteria: finite-time unresolved.
> Large-ramp linearization inaccurate; tiny ramps verify quadratic error.
> New `spac_carbon_section.py`: phase-fixed E4 Poincare DP, left unstable
> covector and stable tangent plane, independent nonlinear one-return check.
> At r=.002, VPD~.00266474448, radial e-folding~474year; DIFFERENT slice
> from rate endpoint .00266. Codim1 tangent candidate, NOT computed nonlinear
> stable manifold or global recovery/exit boundary. No true E4 R-tipping.
> Positive scalar benchmark20/120year reproduces low/high branch transition;
> exact linear periodic-response benchmark is nonautonomous, not autonomous cycle.
> New13 tests;212 related regression passed. Final artifacts/16source hashes:
> `artifacts/spac_carbon/stage_c23_tracking_final_v2_20261014/`.
> Run `python examples/analyze_spac_carbon_tracking.py --output NEW_PATH`.
> Next: same-VPD cycle continuation, nonlinear local stable-manifold geometry,
> bounded outcomes with exit/local-leaving/unresolved kept separate. C4 new
> laws and C5 manuscript deferred; no mortality/new allometry mechanism.

> 2026-10-13 Stage C capacity-first milestone: newest user direction is
> `docs/SPAC_nonlinear_carbon_storage_capacity_research_plan.md`.
> Read `docs/spac_carbon_capacity_definition_ledger.md` and
> `docs/spac_carbon_stage_c_results.md` before further tipping exploration.
> Baseline 2ef9438c080d8183edc1efeff6e3a3b9fad62eb1; original A/B unchanged.
> WP-C0/C1 completed with independent `spac_carbon_capacity.py`, 27 new cases.
> K=-M^-1u, P=K-C=-M^-1 Cdot: maps/conditional roots/full roots/attracting
> sets/nonautonomous trajectories must stay distinct. P=0 need not mean H=0;
> sum(P)=0 need not mean P=0. Formal zero carbon is NOT a defined full root.
> 232 bounded E0-E4 root records, 15 counterexamples, E4 frozen transient and
> repelling-cycle diagnostics; potential norm is nonmonotone on the cycle.
> Numerical Floquet modes imply local stable-manifold dimension5/codimension1,
> NOT a computed global recovery/exit boundary. Radial e-folding ~1657 years
> at r=.001 makes short-time outcome labels particularly unreliable.
> Conditional tracking bound derived under UNVERIFIED uniform transition-matrix,
> nonlinear-remainder and positive-local-domain assumptions; no E4 guarantee.
> Results: `artifacts/spac_carbon/stage_c_capacity_final_v2_20261013/`.
> `python examples/analyze_spac_carbon_capacity.py --output NEW_PATH`.
> Directly related regression:199 passed. 18 current source hashes saved.
> Old A/B1/B3 artifact hashes match; B2/pressure archive differs only in
> bifurcation helper generalized during B3, not changed during Stage C.
> Next WP-C2/C3: tracking and bounded Floquet-directed recovery/exit sections
> within 0--.003 MPa, not a new mortality/allometry/tipping mechanism.

> 2026-10-13 B3 E4 nonlinear milestone: read
> `docs/spac_carbon_stage_b3_e4_results.md`. Latest MPa/scale review takes
> precedence over the older E3 high-VPD global-continuation priority.
> Exact four-state E4 tensors + generalized existing l1 helper verify
> NUMERICAL subcritical Hopf at VPD=.002664752905 MPa, l1=+39.69534297,
> slope ~419.869/day/MPa, with two stable noncritical active modes.
> Independent periodic collocation at r=.0005/.001/.002, full R/S lift,
> independent IVP, Floquet and tolerance refinement verify repelling cycles.
> Seven predefined sensitivities: capacities x.5/x2 retain type; wood
> exponent .75 or loss x.5 moves crossing outside 25 C RH range.
> Exponent1.25 has a fold .002128793 MPa and a complex-pair crossing on an
> ALREADY UNSTABLE branch; no attracting-side Hopf/l1/cycle claim there.
> Results + ten input hashes: `artifacts/spac_carbon/stage_b3_e4_final_v2_20261013/`.
> Run `python examples/analyze_spac_carbon_e4.py --output NEW_PATH`.
> Targeted/neighboring tests:163 passed. Original physiology/config unchanged.
> Next: bounded E4 cycle/stable-manifold vs recovery/exit boundary analysis
> in routine VPD, coefficient provenance; consult before post-exit physiology.
> No global basin, interval-certified Hopf, mortality, or R-tipping claim.

> 2026-10-13 unit/scale review: start with `docs/spac_carbon_units_and_scales.md`.
> User confirmed VPD remains positive vapor-pressure deficit, NOT psi_air.
> Current input/output boundary `spac_carbon_mpa.PressureSPAC` uses MPa and
> requires explicit atmospheric pressure. Original equations/coefficients,
> Stage A/B1/B2 APIs/config/outputs remain canonical mol/mol for provenance.
> New `configs/spac_carbon_pressure_review.json`: p_atm=.101325 MPa;
> 25 C is RH diagnosis ONLY; routine VPD 0--.003 MPa, formal 0--.02 MPa.
> Preserve gs scaling gs/p_atm, slope/parameter-shift chain rules, and
> outward interval inclusion. 27 new unit-equivalence tests + 127 related
> baseline tests passed. Current runner `examples/review_spac_carbon_units.py`.
> Final `artifacts/spac_carbon/pressure_review_final_v2_20261013/` includes
> MPa curves/thresholds, E2/E3 11/26-cell certificates and E4 local check.
> E2/E3 Hopf at .015337756/.004863191 MPa exceed saturation at 25 C
> (.003167489 MPa). Do not treat them as main climate-range phenomena.
> E4 candidate .002664753 MPa remains within routine range: exact
> (psi_x,psi_leaf,L,W) four-state core, passive R/S. Quartic Delta3 and
> critical pair cross zero locally; slope ~419.869 day^-1/MPa refines.
> NOT nondegenerate-Hopf/l1/cycle verification; coefficients uncalibrated.
> Updated 3D MPa figures in `artifacts/spac_carbon/phase3d_mpa_verified_20261013/`.
> Next: synthetic-coefficient sensitivity and E4 local nonlinear analysis,
> not high-VPD cycles as ordinary climate exploration. Consult user before
> post-exit physiology. No commitments/empirical mortality claims added.
> Historical sections below retain canonical D; do not pass MPa directly.

> 2026-10-13 explanatory synthesis:
> `docs/spac_carbon_research_logic_report.md` explains question -> reduction ->
> method -> evidence -> scope, with analytic/interval/numerical/hypothesis labels.
> Read it before presenting results to the user. Standard local bifurcation
> methods are not in themselves a mathematical frontier; global basin geometry
> and post-exit physiology remain unresolved. The D=.03 rate comparison and
> D~.048 Hopf portrait are different parameter slices; do not assert that the
> plotted repelling cycle has been proved to bound the rate-exit basin.
> No new experiment, code change or stronger mortality/R-tipping claim.

> 2026-10-12 Stage B2: short plan `docs/SPAC_carbon_stage_B2_plan.md`,
> results `docs/spac_carbon_stage_b2_results.md`; baseline commit 05a7609.
> Stage A/B1 equations, synthetic coefficients and prior outputs are unchanged.
> E2/E3 have an exact three-state dynamic core (water x2 + leaf), not a
> quasi-steady approximation; W/R/S contribute stable passive eigenvalues.
> Outward-rounded 40-digit Decimal intervals certify Hurwitz inequalities,
> monotone D(L), and endpoint coverage over the whole D=.01-.03 branch
> (E2:11 cells, E3:26). Claim is stable frozen branch, NOT tracking/no-R.
> Analytic second/third tensors + Lyapunov formula + transversality resolve
> E2/E3 candidates as NUMERICALLY verified subcritical Hopf:
> D_H=.15137188728/.04799596087, l1=+11.53120/+22.67162.
> Nonzero-amplitude shooting and full periodic W/R/S lift verify repelling
> cycles at r=.0005/.001/.002; Floquet radial multiplier >1, neutral ~1.
> Unstable cycles are NOT alternate attractors or proof of R-tipping.
> Formal L=0/E2/E3 demand vanishes: gx=0 suffices for a zero-carbon water
> continuum even if gs>0. No isolated dead equilibrium or physical mortality
> regime was added. Conserved H_i(psi,organ C), tissue water/carbon transfers
> and regeneration inputs must be specified before a post-exit variant.
> Tests: 129 targeted/neighboring passed. Final artifacts + 9 source hashes:
> `artifacts/spac_carbon/stage_b2_verified_20261012/`.
> Run `python examples/analyze_spac_carbon_bifurcation.py --output NEW_PATH`.
> Next: bounded global cycle/stable-manifold vs recovery/exit boundary analysis;
> consult user before choosing post-exit physiology. E4 nonlinear analysis,
> empirical forest mortality, NSC, dynamic soil and Leuning remain unimplemented.

> 2026-10-09 Stage B1: newest user specification is
> `docs/SPAC_carbon_stage_B_plan.md`; results in `docs/spac_carbon_stage_b_results.md`.
> It supersedes the older "Leuning first" next-step suggestion below.
> Independent `spac_carbon_two_way.py` implements E0–E4 structural flux feedback
> with fixed capacities/references and unchanged Stage A GPP/carbon budgets.
> Correct series conductance is Gr*Gl/(Gr+Gl), NOT gx*FL.
> Positive equilibrium reduction + reduced pseudo-arclength + full 6x6 Jacobian
> + Schur/nondegeneracy checks + full multistart roots are verified.
> E0 fold D=.02439716312; E1 leaf-supply-only fold D=.01526557528.
> E2/E3/E4 have no detected fold in the stated viable-domain scan but have
> oscillatory stability candidates (.1513719/.0479960/.0262991); no Hopf proof.
> Same D=.01 -> .03 path: 10-day E2/E3 exits L_min, 10-year recovers positive
> equilibrium after 500-year hold. This is rate-dependent DOMAIN EXIT, NOT
> verified R-tipping/mortality: no physical alternate attractor beyond L_min.
> E2 invariant-box non-exit bound D<=.01859886673 is analytic and rate-independent,
> but not a global convergence/no-R theorem. Keep the dry continuum unresolved.
> Targeted + neighboring tests: 115 passed; original Stage A 88 revalidated.
> Final results/source hashes: `artifacts/spac_carbon/stage_b1_final_20261009/`.
> Reproduce: `python examples/analyze_spac_carbon_two_way.py --output NEW_PATH`.
> Next: rigorous continuous-path stability, nonlinear oscillatory bifurcation
> analysis, physically justified domain extension BEFORE alternate-attractor
> R-tipping claims. E4 transient, Leuning, soil/precipitation, variable capacities,
> NSC/damage/recovery remain deferred. No native VISIT or Stage A code changed.

> 2026-10-09 user-directed SPAC–carbon plan: read
> `docs/SPAC_carbon_research_proposal_2026.docx`,
> `docs/AGENT_HANDOFF_CWexperiment.md`, and `docs/spac_carbon_stage_a.md`.
> User waived separate PoC Python files. Independent six-state one-way
> `spac_carbon.py` reuses checked Manzoni, adds four carbon pools and concave
> leaf-dependent gross assimilation. GPP + explicit respiration/litter partition
> replaces proposal's effective-NPP convention; all carbon/T coefficients synthetic.
> Proven: R0 zero/positive threshold, unique positive root/global conditional
> convergence, input-off zero capacity/decay, block-triangular stability lifting,
> and no finite-temperature-pulse R-extinction if final R0>1 and water stays safe.
> Sigmoid three-water-root example lifts to positive stable / saddle / zero stable
> selected carbon roots, but zero-leaf water states are only mathematical extension.
> Linear dry continuum is NOT regularized. No two-way/Leuning/real mortality claim.
> VPD-only fold crossing and same-endpoint temperature recovery controls saved;
> final outputs: `artifacts/spac_carbon/stage_a_final_20261009/`.
> Run `python examples/analyze_spac_carbon.py --output NEW_PATH`.
> Dedicated + neighboring tests: 88 passed. Next: plan Stage B, consistent
> Leuning/assimilation closure and physically meaningful zero-biomass regime,
> before structural feedback or genuine R-tipping. No native VISIT changes.

> 2026-10-07 Manzoni remaining-figure follow-up: read
> `docs/manzoni2014_remaining_results.md`. Fig. 3 median/trait thresholds,
> Appendix A sigmoids/all finite roots/full Jacobian, and soil-root continuity
> are implemented in `manzoni_extensions.py`. Root counts use log-domain
> exponential-sum isolation with dominant-tail bounds, not a finite dry cutoff;
> tangent/equal-exponent/extreme-root checks included. Use log flux residuals
> and positive-gx-scaled equilibrium Jacobians when absolute fluxes underflow.
> Appendix maps exhibit 0/1/2/3 roots; 3364 common coarse/fine counts agree.
> User approved assumed-soil Fig. 4 sensitivity (NOT original numeric replication):
> psi_sat=-0.1 MPa, gsr_sat=0.03/0.1/0.3 m/day/MPa, fixed damage fraction=0.6.
> Fig. 3 observational percentiles/points and original soil/capacitance values
> remain unavailable. A2 low/high D and unsafe traits are explicitly chosen.
> Reproduce: `python examples/reproduce_manzoni2014_remaining.py --shape-step .25 --output NEW_PATH`.
> Final artifacts: `artifacts/manzoni2014/remaining_final_20261007/`.
> Dedicated + neighboring regressions: 66 passed; no full suite run.
> No carbon/cover coupling, mortality, recovery dynamics or R-tipping claim.

> 2026-10-07 user-directed Manzoni (2014) reproduction: read
> `docs/manzoni2014_reproduction.md`. Local paper main text and Appendix A
> were read. Independent `manzoni_hydraulics.py` implements equations (1)--(4),
> admissible active equilibria, the full Jacobian and the Fig. 2 fold.
> Wet-root zero-limit / -1 MPa critical D: 0.02439716312 / 0.01307440031 mol/mol.
> D is NOT kPa; capacitances are numerical fixtures, not paper-calibrated.
> The nonconducting dry state is an equilibrium continuum in the full 2D model,
> not an isolated hyperbolic sink. No carbon coupling or R-tipping claim.
> Reproduce with `python examples/reproduce_manzoni2014.py --output NEW_PATH`.
> First plots/results/hashes: `artifacts/manzoni2014/figure2_20261007/`.
> Dedicated + neighboring regressions: 38 passed; no full-suite run.
> Fig. 3/Appendix/soil continuity follow-up is recorded above.
> Original soil coefficients and observational trait data remain unresolved.
> Existing VISIT and cover/carbon tracks unchanged.

> 2026-10-02 current direction: read docs/forest_grass_cover_climate_research_plan.md
> first. Derive an 11-state baseline (grass cover + ten carbon stores), including
> area conversion budgets and the block-triangular stability proof. Then source
> temperature/relative-humidity response laws and test the scalar no-R barrier.
> Use fuel-moisture memory or carbon feedback only if justified. Cover bistability
> does not establish climate-driven R-tipping or different total carbon stocks.
> The plan itself is not a climate experiment or a tipping result.

> 2026-10-02 Stage A baseline implemented in
> `src/control_carbon/cover_carbon_baseline.py`; assumptions and the equilibrium
> lifting proof are in `docs/cover_carbon_baseline.md`. Carbon coefficients in
> tests are uncalibrated fixtures. Distinct cover equilibria imply distinct
> grass living-carbon pools, but do not determine total-carbon ordering. No
> temperature/RH response or climate-driven R-tipping is implemented. Next:
> source-constrain one climate-to-process response before rate experiments.
> The phase/function/capacity visualization is in
> `docs/figures/cover_carbon_baseline/cover_carbon_baseline_atlas.pdf`; regenerate
> with `python examples/plot_cover_carbon_baseline.py`.

> 2026-10-02 climate-response literature ledger: see
> `docs/cover_climate_response_ledger.md`. The ledger separates VPD, plant
> productivity, seedling recruitment, fire/fuel moisture, and soil decomposition.
> Only the physical T/R-to-air-VPD transform is presently code-ready. Do not
> use RH as soil moisture or map fire intensity directly to the cover rate.

> 2026-09-30 first tree–grass carbon exploration: read
> docs/tree_grass_process_ledger.md, tree_grass_carbon_equations.md and
> tree_grass_bistability_analysis.md. Cerrado shade–grass fuel–fire damage is
> selected as a source-supported interaction, with uncalibrated closures.
> Ten-pool budgets, positivity and boundedness are derived; shade-only positive
> equilibrium uniqueness, full-system stability and conditional no-R are proved.
> Endogenous-fire equilibria satisfy f=lambda*h(chi*D_G(f)); their full Jacobian
> is a rank-one update of the fixed-fire Jacobian. Do not infer stability from
> the scalar root slope alone, or eliminate fuel litter from transient dynamics.
> Climate-to-fuel-moisture closure, surface/dead-root litter separation and
> resprouting after complete topkill remain unresolved. No fire bistability or
> climate R-tipping claim. Reproduce checks with
> python examples/check_tree_grass_carbon.py; saved JSON is in docs/.
> Next: constrain fuel response and topkill/regrowth before rate experiments.

> 2026-09-30 next research plan: for the tree–grass carbon track, read
> docs/tree_grass_carbon_climate_research_plan.md first. Start with ten carbon
> pools and a source-backed reciprocal interaction. Derive carbon budgets,
> equilibrium reductions and no-tipping/tipping conditions before numerical
> campaigns. Constrain temperature/moisture response laws before selecting a
> forcing path. Prove or delimit the connection from reduced dynamics to the
> full carbon system. No new experiment or climate threshold is claimed by this plan.

> 2026-09-30 forest–grass theory: read docs/forest_grass_nbr_tipping.md
> and the Japanese TeX/PDF first for the new track. Two supplied papers were
> read; only Kumar K & Dutta (2026) is a forest–grass model, with fire feedback.
> Proven: two interior stable cover states, B-folds, reflected-noise N escape,
> conditional no-R theorem for alpha-only forcing, and existence/unique rate
> threshold on an explicitly constructed alpha/phi1 path. The last path is
> mathematical, not empirical; no native VISIT tipping or temperature threshold
> is claimed. VISIT canopy shading is one-way in the audited light-only slice;
> N/water/other couplings are not excluded. Native code unchanged.
> New module forest_grass_tipping.py; 7 dedicated tests pass.
> Reproduce with examples/analyze_forest_grass_tipping.py and a fresh output
> directory; saved artifacts are git-ignored. Use forward/backward matching,
> not direct forward endpoint shooting near the unstable critical trajectory.
> Next: biological feedback/forcing constraints, then carbon pool accounting.

> 2026-09-22 NSC first campaign: independent 4 structural + 4 NSC/labile
> continuous model, no water coupling. Read `docs/temperature_nsc_first_results.md`
> and `temperature_nsc_equations.md` (Japanese TeX also supplied). Standard
> constant leaf target: 42/42 multiseed spinups converge; supported leaf target:
> stable positive + saddle + bare, numerical fold at 37.50696 C. No R-tipping
> found in 8 core ramps or 24 follow-up controls. Near-fold 500-year deviations
> disappear by 2000 years. Bistability remains without NSC mortality: do NOT
> attribute it solely to the NSC threshold. All numeric coefficients assumed,
> NOT literature-calibrated. Next constrain canopy rebuilding and absolute-pool
> versus concentration-based maintenance before adding states or seasonality.
> Final artifacts under `artifacts/temperature_nsc/`: `autonomous_v1_final`,
> `supported_followup_v1`, `controls_v1`, `audit_v2`, `report_v1`. Earlier failed
> runs retained separately; fixed log-Jacobian avoids num_jac overflow. No clipping.
> Tests: 204 passed / 7 skipped, including 14 dedicated NSC tests.

> 2026-09-22 T/P v2: explicit constant-climate spinup (375–400 years), JSON
> process selection, accepted-step/dense-sample nonnegativity guards, and
> native C checks of old-VISIT formula components are implemented separately.
> See `docs/tp_v2_results.md`, `docs/tp_v2_spinup_and_functions.md`, and
> `configs/tp_v2_spinup.json`. The original 4 variants reproduce v1 trajectories;
> negative I_C was never a negative stock. Ten configurations/40 factorial runs
> passed; equation choice can change interaction sign. Artifacts live only in
> `artifacts/temperature_precipitation/v2/spinup_comparison/`. Do not overwrite
> v1 outputs. Previous plan is archived verbatim; revision 1.1 adds only §4.4.
> Tests: 190 passed / 7 skipped. No R-tipping claim or parameter calibration.

> 2026-09-22 initial T/P implementation: WP1 passes and first WP2 comparisons
> are in `docs/temperature_precipitation_first_results.md` and companion TeX.
> New module `temperature_precipitation.py` is an independent theoretical
> eight-state model with explicitly assumed parameters. Mortality hierarchy,
> two rainfall endpoints, one-driver continuation, and 200-year hold are saved
> under `artifacts/temperature_precipitation/`. Water-mortality I_C changes
> from -9.901 at year 20 to +37.380 at year 200; no R-tipping claim.
> Fixed absolute-water clamps violate shrinking storage capacities and are
> rejected after violation (audit.json). Next: valid clamp definition,
> boundary invasion stability and two-parameter branch/basin analysis.
> Do not infer uniqueness from stable sampled branches or add states yet.
> Tests: 176 passed, 7 skipped; dedicated tests: 13.

> 2026-09-22 temperature–precipitation plan: see
> `docs/temperature_precipitation_tipping_research_plan.md`. The new theoretical
> track starts from matching leaf/stem/root/soil carbon and water stores, monthly
> temperature/precipitation forcing, and a four-run factorial design. It keeps
> quantitative interaction, frozen-system bifurcation, and R-tipping claims
> separate. It complements rather than replaces the native 37-state workbook,
> VISITc audit, and eight-state explicit-SPAC model. No implementation or
> calibration is claimed at this checkpoint.

> 2026-09-21 native-state follow-up: pinned 2013 VISIT was compiled and run
> with supplied TKY forcing/deposition/GHG and reconstructed text parameters.
> `VISIT_TKY_matrix_equations_native_20000715.xlsx` now has all 37 start states,
> 15 histories and post-location forcing for 2000-07-15. Missing inputs fell
> 70 -> 2 (only two source-undefined diagnostics). All 37 one-day endpoints
> match native C, maximum scaled error 5.44e-16. Runtime Saxton capacities
> replace raw Config capacities (104.8103 / 704.0597 mm). Native execution also
> exposed a missing `==0` in a `strcmp` vapor-pressure branch; actual TKY vp is
> retained and flagged. See `docs/visit_native_snapshot_20000715.md`.

> 2026-09-17 actual-data follow-up: summer snapshot workbook
> `artifacts/visit_matrix_workbook/VISIT_TKY_matrix_equations_20000715.xlsx`.
> `visit_workbook_data.py` imports supplied visitb data read-only, retaining the
> pinned equations/Config. 24 of 70 missing inputs resolved (14 forcing/aggregate
> inputs, 10 histories reconstructed with initially unknown memory); 46 remain.
> All 37 summer states remain unknown. The 21 spinup-end restart stocks are a
> separate reference sheet, not substituted into summer X. New sheets 06–08
> record conversion, provenance, restart and SHA256 comparisons. Older template
> remains unchanged. Source-version differences in mean-temperature period,
> annual-rain accumulation and VPD branch are documented; pinned rules win.

> 2026-09-17 TKY matrix workbook: see `docs/visit_matrix_workbook.md`.
> User-requested Excel audit covers the selected old `visit_local` daily
> C/N/water path with 37 states, 132 flux channels and a 37×78 B matrix.
> Config_TKY defaults are imported; actual states/weather/history remain blank.
> `visit_workbook.py`, `_export.py`, `_validation.py` generate real Excel formulas,
> cached values and source/unknown inventories. Seven artificial cases match
> 1,379 native-C values to rounding error; 13 dedicated tests pass. It is a
> source-order daily effective-rate representation, not a continuous VISIT ODE
> or the 30-minute pathway. Native baseflow/N-unit/undefined-gas issues are
> exposed, not silently fixed. Excel-application recalculation remains untested.

> 2026-09-17 capacity/QSE exploration: `docs/capacity_qse_attractor_exploration.*`
> derives the state-dependent capacity map, its self-consistency equation, and
> the full-equilibrium Jacobian `J=M(I-D capacity)`. A model-agnostic positive
> scalar counterexample has two nonzero stable capacity branches and gives
> different final branches for 20-year and 120-year ramps; it is not a VISIT
> result or calibration. `docs/visit_r_tipping_process_screening.*` audits the
> 2020 explanatory slides and the pinned old VISIT source. The source-grounded
> next target is the complete daily plant map, especially LAI--GPP--allocation--NSC.
> New module/tests/example: `capacity_feedback`; nine dedicated tests pass.
> Full suite at this checkpoint: 146 passed, 7 source-dependent skips.

> 2026-09-15 carbon–nitrogen follow-up: `docs/nitrogen_tipping_minimal.md` and TeX/PDF.
> Independent 2-state fixed-stoichiometry benchmark: monotone N limitation has no
> permanent R-tipping when the final positive equilibrium exists. An explicitly
> hypothetical excess-N inhibition law admits rate-dependent basin changes;
> the fast-step basin entry has an analytic sufficient proof. Equal-dose centered
> ramps and 3-state delayed recycling are tested. No VISIT/water coefficients changed.
> New module/tests/example: `minimal_nitrogen_tipping`; full suite: 144 passed.

> 2026-09-15 theory follow-up: see `docs/water_tipping_minimal.md` and its TeX/PDF.
> User-requested minimal R-tipping analysis is a separate uncalibrated benchmark
> (`minimal_water_tipping.py`), not a completion of the native calibration path.
> It proves a scalar rainfall-only no-tipping result and a finite-rate threshold
> for a paired environmental path, plus a conditional dynamic-water example.
> VISITa was inspected at `5c513196f21c1b1efb9ead540e3d40865b15e07e` but not run.
> Existing hydraulic unit/boundary caveats remain unresolved. Full suite: 129 passed.

> Active priority (2026-09-14): explicit water storage and SPAC transport in
> the land-carbon matrix framework. Read `docs/carbon_water_research_plan.md`,
> `docs/coupled_carbon_water_equations.md`, and
> `docs/visitc_carbon_water_source_ledger.md` before using the historical
> skeleton below. Deterministic tipping code remains a tested, deferred track.

## 0. Purpose of this document

This repository is intended to become a **model-agnostic control-theoretic framework for terrestrial ecosystem carbon-cycle models**. VISIT is the first source-grounded implementation and validation target, not the final scope.

The immediate coding-agent task is **not** to produce a polished VISIT emulator as quickly as possible. The task is to build a rigorous chain from process-model source code to state-space/IRF representations while preserving provenance and clearly separating exact equations from reductions and local approximations.

The active chain is now:

```text
current VISITc water source
  -> source-order water ledger and native comparison
  -> four-component carbon-water balances
  -> coupled equilibrium/capacity/Jacobian diagnostics
  -> dynamic-vs-quasi-steady water experiments
  -> periodic tracking and conditional tipping analysis
```

Implemented in the current working tree:

- `visitc_source_map.py`: current `visit-manager/VISITc` provenance at
  `5202debd96df6f88beb7d61f8688fff02ace964a`;
- `visitc_hydrology.py`: `f_hydrology`, atmospheric/resistance, LAI/radiation,
  canopy-conductance, and PM transcriptions with an explicit baseflow audit;
- `visitc_native.py`: direct native-C checks for all exposed one-day hydrology,
  radiation, PM, and canopy-conductance outputs;
- `hydraulic_relations.py`: native VISITc soil retention and reduced
  SurEau-Ecos plant pressure-volume/capacitance relations;
- `coupled_carbon_water.py`: eight-state reduced ODE, budgets, carbon capacity,
  capacity-change decomposition, coupled equilibrium, Jacobian blocks, Schur
  reduction, and simulation;
- `examples/run_coupled_carbon_water.py`: idealized dry-down/re-wetting
  experiment with refinement and provenance outputs;
- dedicated tests and carbon-water documentation.

Do not call the eight-state continuous model “VISITc.” It is a synthesis. The
implemented hydrology/PM slices are native-validated, but they are not a
complete source-faithful VISITc model: initialization, all site branches, and
the complete carbon system remain outside the bridge.

The research objective is to make it possible to compare many terrestrial ecosystem models in a common language:

- state variables and flux topology,
- quasi-static equilibrium (QSE),
- disequilibrium from QSE,
- impulse response functions (IRFs),
- poles / modal time scales / residues,
- transfer functions and frequency responses,
- NEP/NEE/NECB input-output behavior,
- local linearization and, later, periodic/LTV response,
- forcing-rate response and analytical source/sink criteria.

## 1. Authoritative repositories

### Primary implementation repository

- `Sachitama2001/Control_carbon_model`

### Source/reference repository for the first model adapter

- `Sachitama2001/VISIT-matrix`
- pinned commit `3285bd8e131a932e338b59892751648fd9edcc7b`
- authoritative VISIT site-model source is under `visit_local/`

Important: `VISIT-matrix/visit_matrix/` already contains a reduced 18-pool matrix implementation. Use it as prior work and a cross-check, **not as the authority**. The C source in `visit_local/` is the authority when deciding whether an equation, state, input, or parameter is actually present in VISIT.

### Source/reference repository for the active water audit

- `visit-manager/VISITc`
- pinned commit `5202debd96df6f88beb7d61f8688fff02ace964a`
- inspected source is under `point/`

Keep this provenance separate from the older VISIT-matrix adapter. A formula
found in one snapshot is not automatically authoritative for the other.

## 2. Current repository status

> Status note: this section describes the original skeleton. The implemented
> checkpoint has advanced substantially; read `docs/current_status.md` before
> using the older task ordering below.

Current package skeleton:

```text
Control_carbon_model/
├── README.md
├── HANDOFF.md
├── pyproject.toml
├── docs/
│   └── visit_state_space_source_map.md
├── src/control_carbon/
│   ├── __init__.py
│   ├── state_space.py
│   └── visit_source_map.py
└── tests/
    └── test_state_space.py
```

Existing functionality:

1. `state_space.py`
   - generic continuous-time LTI representation,
   - poles,
   - transfer-function evaluation,
   - frequency response,
   - output IRF `C exp(A t) B`,
   - state IRF `exp(A t) B`.

2. `visit_source_map.py`
   - source-provenance registry for VISIT variables/processes.

3. `docs/visit_state_space_source_map.md`
   - first source-grounded decomposition of VISIT into state/input/output/process layers.

4. `tests/test_state_space.py`
   - basic analytical consistency checks for generic LTI utilities.

## 3. Scientific design rule

Every VISIT-derived model object must belong to one explicit approximation level.

### Level 0 — native nonlinear discrete map

```text
x[k+1] = F_VISIT(x[k], z[k], q[k]; theta)
y[k]   = H_VISIT(x[k], z[k], q[k]; theta)
```

This is the conceptual reference. VISIT is run daily and contains sequential updates and hybrid/piecewise rules.

### Level 1 — source-grounded reduced compartment model

```text
dx/dt = M(t) x + B(t) mu(t)
```

or, where the native daily discretization matters,

```text
x[k+1] = A_d[k] x[k] + B_d[k] mu[k].
```

The soil carbon subsystem is especially close to this form because, conditional on environmental scalars, decomposition is first-order in carbon stocks.

### Level 2 — local linearization of the native model

```text
delta x[k+1] = A_d delta x[k] + B_d delta z[k]
delta y[k]   = C_d delta x[k] + D_d delta z[k]
```

with Jacobians evaluated at an explicitly defined operating point.

### Level 3 — periodic linear time-varying representation

Because phenology and meteorology are strongly seasonal, an annual periodic trajectory is likely a more physically meaningful reference than a single annual-mean equilibrium.

### Level 4 — model-agnostic comparison layer

Different terrestrial ecosystem models should expose the same abstract interfaces even when their native state dimensions and process structures differ.

Never collapse these levels into a single class called simply `VISITModel` without metadata indicating which approximation is being used.

## 4. Central mathematical targets

For a nonlinear model

```text
x_dot = f(x, z, theta)
y     = h(x, z, theta)
```

or its discrete equivalent, implement tools to obtain a local state-space system

```text
delta x_dot = A delta x + B delta z
delta y     = C delta x + D delta z
```

and then

```text
G(s) = C (sI - A)^(-1) B + D.
```

For carbon stock total

```text
Ctot = w^T x
```

and NEP-like outputs, the long-term target is to connect:

```text
QSE geometry -> disequilibrium -> forcing rate -> NEP sign
```

with IRF / transfer-function quantities.

A recurring object will be the moving QSE

```text
x_star(z) : f(x_star, z) = 0
```

and its sensitivity

```text
dx_star/dz.
```

In an LTI/local-linear setting, the low-frequency response to forcing rate should recover the QSE sensitivity. This identity should eventually become a tested theorem/example in the repository rather than remain only in notes.

## 5. First model: VISIT

### 5.1 Minimal structural-carbon state

Use the following 18-pool ordering as the first reduced carbon state:

```text
0  tree_fol
1  tree_stm
2  tree_rot
3  c3_fol
4  c3_stm
5  c3_rot
6  c4_fol
7  c4_stm
8  c4_rot
9  ltr_tf
10 ltr_tc
11 ltr_tr
12 ltr_gf
13 ltr_gc
14 ltr_gr
15 msl_a
16 msl_i
17 msl_p
```

Do not describe this as the complete native VISIT state. `nsch_storage`, crop `grain`, sub-compartments, phenology memory, water stores, and potentially other variables can carry dynamical memory.

### 5.2 External forcing candidates

At minimum track:

```text
Tair, Tsurface, Tsoil10, Tsoil200,
SW radiation, precipitation, humidity/VPD, wind, atmospheric CO2
```

For early IRF work, effective inputs are allowed:

```text
GPP or EPP by PFT,
f_tm_l,
f_tm_h or pool-specific humus scalars,
allocation fractions,
phenology/regime variables.
```

But label them explicitly as **effective/internal drivers**, not meteorological forcing.

### 5.3 Outputs

At minimum:

```text
GPP, NPP, heterotrophic respiration, NEP, total carbon
```

Keep `NEP` distinct from `NECB`. Harvest, BVOC, DOC export, disturbance, and other lateral/non-respiratory losses can cause `dCtot/dt`, NEP, and NECB to differ.

## 6. Source-provenance requirement

Every VISIT adapter equation or matrix element must have provenance metadata with:

- repository,
- source path,
- function name,
- source variable names,
- equation/rule summary,
- approximation notes,
- regime assumptions,
- unit conversion.

A matrix coefficient with no provenance should not be merged.

Preferred pattern:

```python
SourceRef(
    repository="Sachitama2001/VISIT-matrix",
    path="visit_local/soil_proc.c",
    function="f_cycle_soil",
    symbol="degrade_tf",
    note="ltr_tf * sr_lf / 1000 * f_tm_l",
)
```

## 7. Implementation order

Follow the roadmap in `docs/implementation_roadmap.md`. The short version is:

1. Harden generic state-space API.
2. Complete provenance representation.
3. Implement VISIT soil subsystem exactly from C equations.
4. Verify one-step mass balance algebraically.
5. Implement reduced plant structural-carbon subsystem.
6. Make every approximation explicit.
7. Build discrete-time local-linear tools.
8. Build QSE solver and QSE sensitivity.
9. Implement IRF/step/ramp/frequency diagnostics.
10. Derive and test NEP transfer-function identities.
11. Compare reduced model IRFs against native VISIT perturbation experiments.
12. Add periodic/LTV treatment.
13. Generalize adapters for additional ecosystem models.

## 8. Definition of done for the first major milestone

Milestone M1 is complete when all of the following are true:

- the 18-pool reduced VISIT carbon model can be built from documented parameters and drivers;
- each nonzero matrix/flux relation is traceable to VISIT C source;
- one-day soil and plant-litter mass-balance unit tests pass;
- NEP from the reduced model is explicitly computed and tested;
- constant-forcing QSE can be solved;
- the local LTI/discrete state-space system around QSE is constructible;
- impulse and step responses are available for at least EPP, temperature/decomposition-scalar perturbations;
- direct time integration agrees with convolution/IRF in the linearized regime;
- README contains a reproducible example;
- no result is labelled "VISIT" without its approximation level.

## 9. Things the coding agent must not do

- Do not silently copy numerical defaults from `VISIT-matrix/visit_matrix` and call them VISIT parameters unless verified against C parameter initialization/input files.
- Do not assume `LAI` is always an independent state; in normal code it is recalculated from foliage mass, but fixed-phenology modes can overwrite it.
- Do not assume the native code is a simultaneous continuous ODE; it performs sequential daily updates.
- Do not treat all forcing as temperature-only. Temperature is the first scientific focus, but soil water, radiation, CO2, and phenology can dominate particular paths.
- Do not infer exact native Jacobians from the reduced 18-pool model.
- Do not optimize performance before establishing source fidelity and tests.
- Do not merge large undocumented refactors.

## 10. Recommended working style

For each new VISIT process:

1. locate the exact C routine;
2. write the equation in a source note;
3. identify state/input/parameter/output roles;
4. implement the smallest pure Python function reproducing the algebra;
5. add a unit test with hand-computable values;
6. only then insert it into a matrix or nonlinear adapter;
7. add provenance metadata;
8. document any approximation or branch excluded.

## 11. Expected eventual architecture

```text
src/control_carbon/
├── core/
│   ├── continuous.py
│   ├── discrete.py
│   ├── ltv.py
│   ├── irf.py
│   ├── qse.py
│   ├── modal.py
│   └── reduction.py
├── provenance/
│   └── source.py
├── adapters/
│   ├── base.py
│   └── visit/
│       ├── schema.py
│       ├── source_map.py
│       ├── soil.py
│       ├── plant.py
│       ├── outputs.py
│       ├── reduced18.py
│       ├── linearize.py
│       └── native_bridge.py
└── experiments/
    ├── irf.py
    ├── qse_tracking.py
    └── forcing_rate.py
```

Do not move files merely to match this tree immediately; migrate incrementally with tests.

## 12. Long-term research direction

The final scientific contribution should not be "a matrix version of VISIT". The intended contribution is a framework able to answer, across many terrestrial ecosystem models:

- what are the dominant carbon-cycle modes and time scales?
- which forcing frequencies project onto which modes?
- which internal modes are visible in NEP and which are hidden by cancellation?
- when does NEP change sign because of modal competition?
- how does QSE motion map to forcing-rate thresholds?
- how do model parameters shift poles, residues, QSE sensitivities, and source/sink boundaries?
- which features are robust across model structures?

VISIT is the first testbed for proving that this program can be carried out rigorously from real ecosystem-model source code.
