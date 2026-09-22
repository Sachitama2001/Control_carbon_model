# Literature index: matrix land-carbon models and explicit water dynamics

Last verified: 2026-09-22. Prefer DOI and publisher links because local PDFs
may differ by accepted-manuscript or typeset version. A local copy is listed
only when it is present in `docs/` and redistribution appears permitted.

## Local PDFs

| Local file | Identity and use | Pages | SHA-256 |
|---|---|---:|---|
| [`Luo2022MatrixApproach.pdf`](Luo2022MatrixApproach.pdf) | Luo et al. (2022), core review and notation; public full-text copy | 15 | `b0f879cb2025eccb63e2d930a3aff26697656437baac5457358c1ea3cb1b6465` |
| [`Land Carbon Cycle Modeling_26_03_09_16_10_15.pdf`](Land%20Carbon%20Cycle%20Modeling_26_03_09_16_10_15.pdf) | *Land Carbon Cycle Modeling: Matrix Approach, Data Assimilation, & Ecological Forecasting*; chapter-level background | 399 | `cc77b8af9a214a48b68209741ed138a308d8f19aa4060f872c8bafe424dfc295` |
| [`Ashwin2018FormulateTipping.pdf`](Ashwin2018FormulateTipping.pdf) | Ashwin et al., corrected manuscript of the B/N/R-tipping classification paper; use only after a tracking problem is defined | 22 | `b96de9acf2745212321b044fdf947d1b22b2827aa1bfd150bcdc86e0f2579bd5` |
| [`Feudel2023RateInducedTipping.pdf`](Feudel2023RateInducedTipping.pdf) | Feudel (2023), unstable states and basin-boundary mechanism in R-tipping | 22 | `8e8c1eaadb3ec779b2578eea3f1c4735c59524112dd22f71c7291477dfc16874` |

The Luo review PDF was retrieved from the
[NSF Public Access Repository](https://par.nsf.gov/servlets/purl/10337983);
the canonical article is at
[doi:10.1029/2022MS003008](https://doi.org/10.1029/2022MS003008).

## Core matrix and disequilibrium papers

| Reference | Primary link | What this project should take from it | What it does not establish for us |
|---|---|---|---|
| Luo et al. (2022), “Matrix Approach to Land Carbon Cycle Modeling” | [Article/PDF](https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2022MS003008) | Common form \(\dot X=B\mu+A\xi KX\); carbon input, transfer, environmental modifier and turnover separation; storage capacity and potential | It includes moisture effects in \(\xi\), but does not by that fact make water storage a joint state vector |
| Wei et al. (2022), “Nutrient Limitations Lead to a Reduced Magnitude of Disequilibrium…” | [Article/PDF](https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2021JG006764) | Signed \(X_p=X_c-X\); exact decomposition \(\Delta X_c=NPP_0\Delta\tau_E+\tau_{E,0}\Delta NPP+\Delta NPP\Delta\tau_E\); controlled process comparisons | Water need not have the same sign or dominant term as nutrient limitation |
| Lu et al. (2020), “Full Implementation of Matrix Approach to Biogeochemistry Module of CLM5” | [Article/PDF](https://agupubs.onlinelibrary.wiley.com/doi/pdf/10.1029/2020MS002105) | Methodological standard: reorganize existing equations without changing processes; run matrix and original implementations as replicas; verify compatibility | Our four-pool ODE is a reduction, so it cannot serve as this fidelity test |
| Sierra & Müller (2015), “A general mathematical framework for representing soil organic matter dynamics” | [doi:10.1890/15-0361.1](https://doi.org/10.1890/15-0361.1) | Compartmental-system conventions and distinction between linear/nonlinear kinetics | Does not select a water hydraulic closure |
| Metzler et al. (2018), “Transit-time and age distributions for nonlinear time-dependent compartmental systems” | [doi:10.1073/pnas.1705296115](https://doi.org/10.1073/pnas.1705296115) | Nonautonomous nonlinear compartment systems and trajectory-dependent transit-time/age diagnostics | A representation along one solution does not preserve state-feedback behavior away from that trajectory |

## Explicit SPAC storage and transport

| Reference | Primary link | State/flux lesson for the reduced model | Planned use |
|---|---|---|---|
| Mirfenderesgi et al. (2016), FETCH2 | [PDF](https://agupubs.onlinelibrary.wiley.com/doi/pdfdirect/10.1002/2016JG003467) | Above-ground water storage and tree hydraulic strategy must be dynamic when sub-daily buffering matters | Candidate reference for stem/leaf storage and validation experiments |
| Silva et al. (2022), FETCH3 | [Open article/PDF](https://gmd.copernicus.org/articles/15/2619/2022/) | Coupled soil-root-stem transport, variable capacitance, continuity of potential, mass-conserving PDE formulation | Reference for deriving a low-order compartment closure from a spatial model |
| Ruffault et al. (2022), SurEau-Ecos v2.0 | [Open article](https://gmd.copernicus.org/articles/15/5593/2022/) | Compartment balances use conductance times water-potential gradients; soil retention differs from plant capacitance/pressure-volume relations; capacitance buffers short-term potential variation | Primary reference for replacing normalized storage with \(W\leftrightarrow\psi\) relations and for “remove capacitance” sensitivity tests |
| Kennedy et al. (2019), plant hydraulics in CLM5 | [PDF](https://agupubs.onlinelibrary.wiley.com/doi/pdf/10.1029/2018MS001500) | Example of hydraulics embedded in an Earth-system land model; connects soil layers, xylem states, stomata, and carbon exchange | Bridge between idealized SPAC equations and production land-model architecture |
| Sperry et al. (2017), hydraulic-cost stomatal optimization | [PDF](https://onlinelibrary.wiley.com/doi/pdf/10.1111/pce.12852) | Stomatal response couples photosynthetic gain to hydraulic risk rather than using only an empirical moisture scalar | Later alternative to the first bounded transpiration/photosynthesis functions |
| Eller et al. (2020), SOX in a land-surface model | [PDF](https://nph.onlinelibrary.wiley.com/doi/pdfdirect/10.1111/nph.16419) | Evaluates hydraulically informed stomatal optimization within a land model | Candidate model hierarchy for carbon-water feedback experiments |
| Anderegg & Venturas (2020), “Plant hydraulics play a critical role in Earth system fluxes” | [PDF](https://nph.onlinelibrary.wiley.com/doi/pdf/10.1111/nph.16548) | Concise synthesis of why stomatal/hydraulic representation matters for carbon and water fluxes | Framing, not an equation source by itself |

## Tipping and nonautonomous tracking

The temperature-only NSC track is specified in
[`temperature_nsc_r_tipping_research_plan.md`](temperature_nsc_r_tipping_research_plan.md).
It adapts organ-level NSC bookkeeping to an area-based theoretical model and keeps
frozen-system multistability, B-tipping, and R-tipping as separate claims.

The new temperature–precipitation factorial track is specified in
[`temperature_precipitation_tipping_research_plan.md`](temperature_precipitation_tipping_research_plan.md).
It uses a four-carbon/four-water theoretical model and must remain distinct from
the native VISIT audits and the eight-state explicit-SPAC synthesis.

The current capacity-focused synthesis is in
[`capacity_qse_attractor_exploration.tex`](capacity_qse_attractor_exploration.tex),
with the VISIT equation screen in
[`visit_r_tipping_process_screening.tex`](visit_r_tipping_process_screening.tex).
The central distinction is now explicit: a second stable equilibrium is another
self-consistent QSE, a seasonally forced orbit differs from instantaneous QSE,
and a genuine non-QSE attractor under constant forcing would require a limit
cycle or more complex invariant set.

For the 2026-09-15 user-requested minimal-model analysis, see
[`water_tipping_minimal.tex`](water_tipping_minimal.tex) and
[`water_tipping_minimal.md`](water_tipping_minimal.md).
The new benchmark is independent of the native VISIT adapters.

| Reference | Primary link | Use here |
|---|---|---|
| Ashwin et al. (2012), “Tipping points in open systems…” | [Royal Society article](https://royalsocietypublishing.org/rsta/article/370/1962/1166/114607/Tipping-points-in-open-systems-bifurcation-noise), [doi:10.1098/rsta.2011.0306](https://doi.org/10.1098/rsta.2011.0306) | Definitions separating B-, N-, and R-tipping; need a nonautonomous tracking formulation |
| Alkhayuon & Ashwin (2018), periodic attractors | [doi:10.1063/1.5000418](https://doi.org/10.1063/1.5000418) | Phase-dependent partial/total tipping when the reference attractor is periodic |
| Feudel (2023), basin boundaries and transients | [Open article](https://npg.copernicus.org/articles/30/481/2023/) | Basin-boundary and unstable-state evidence needed to distinguish R-tipping from a long transient |
| Ashwin, Perryman & Wieczorek (2017), parameter shifts | [Author version](https://arxiv.org/abs/1506.07734), [DOI](https://doi.org/10.1088/1361-6544/aa675b) | Low-dimensional tracking and basin-stability background; the new finite-ramp integral is derived separately in the TeX |
| Siero et al. (2019), Grazing Away the Resilience of Patterned Ecosystems | [Primary article](https://doi.org/10.1086/701669), Appendix A1/A6 | Source for the Klausmeier-family water uptake structure, with grazing and spatial terms removed; not forest parameter calibration |
| Kiers (2020), Rate-Induced Tipping in Discrete-Time Dynamical Systems | [Author version](https://arxiv.org/abs/1907.11601), [DOI](https://doi.org/10.1137/19M1276297) | Why a discrete update must be specified rather than assumed equivalent to the ODE |

## NSC, leaf flush, and carbon-starvation mortality

| Reference | Primary link | Use here | Constraint |
|---|---|---|---|
| Ninomiya et al. (2023), SEIB-DGVM-NSC v1.0 | [Open article](https://gmd.copernicus.org/articles/16/4155/2023/), [DOI](https://doi.org/10.5194/gmd-16-4155-2023) | Organ-level leaf/trunk/root NSC pools; NSC use for metabolism and bud flush; model-structure reference | Individual-based DGVM; the reduced area-based ODE is an adaptation, not a transcription |
| Dietze et al. (2014), “Nonstructural Carbon in Woody Plants” | [DOI](https://doi.org/10.1146/annurev-arplant-050213-040054), [author PDF](https://people.bu.edu/dietze/manuscripts/annurev-arplant-050213-040054.pdf) | NSC as a buffer for maintenance and growth; storage/allocation model taxonomy | NSC regulation and mortality linkage remain uncertain |
| Klein et al. (2016), growth, phenology, and storage | [DOI](https://doi.org/10.1093/treephys/tpw030) | Empirical basis for stored-carbon consumption during spring leaf flush | Three temperate deciduous species; not a universal tropical parameterization |
| Barker Plotkin et al. (2021), defoliation and stored-carbon threshold | [DOI](https://doi.org/10.1111/1365-2435.13891), [author PDF](https://harvardforest1.fas.harvard.edu/publications/pdfs/BarkerPlotkin_FunctionEcol_2021.pdf) | Evidence linking low TNC concentration after defoliation with high mortality risk | Oak/defoliation-specific threshold; necessary-like rather than universally sufficient condition |

## Temperature, drought, and tropical-forest mortality

| Reference | Primary link | Use here | Constraint |
|---|---|---|---|
| Huntingford et al. (2017), improved plant respiration | [Open article](https://doi.org/10.1038/s41467-017-01774-z) | Instantaneous respiration response and acclimation alternative to fixed Q10 | Leaf-to-whole-plant aggregation remains model dependent |
| O'Sullivan et al. (2017), thermal limits across biomes | [DOI](https://doi.org/10.1111/gcb.13477) | Tcrit/Tmax variation and thermal-safety-margin context | Short heat assays are not monthly mortality functions |
| Faber et al. (2024), thermal death time in plants | [DOI](https://doi.org/10.1093/jxb/erae096) | Temperature-duration damage accumulation for a later daily extension | Do not reconstruct acute dose from monthly means |
| Aleixo et al. (2019), Amazon tree mortality | [Article](https://doi.org/10.1038/s41558-019-0458-0) | Observational evidence that drought and heat are associated with mortality and delayed effects | Does not uniquely identify an ODE mortality law |
| Nepstad et al. (2002), partial throughfall exclusion | [DOI](https://doi.org/10.1029/2001JD000360) | Eastern-Amazon experimental context for rainfall reduction | Early response paper, not a complete mortality calibration |
| Rowland et al. (2015), tropical drought death | [Article](https://doi.org/10.1038/nature15539) | Hydraulic constraint on drought mortality and warning against a carbon-starvation-only mechanism | A one-bucket water state remains phenomenological |
| Sanchez-Martinez et al. (2025), long-term experimental drought | [Article](https://doi.org/10.1038/s41559-025-02702-x) | Long-term adjustment under Amazon rainfall exclusion | Site response is not a universal tropical parameter set |

The additional sources support `src/control_carbon/minimal_water_tipping.py`,
`tests/test_minimal_water_tipping.py`, and `examples/analyze_minimal_water_tipping.py`.
The TeX records units (years and ground-area stocks), all new assumptions,
monthly/yearly interval maps, and the distinction between analytical examples
and source-faithful equations. VISITa was inspected separately at commit
`5c513196f21c1b1efb9ead540e3d40865b15e07e`; no coefficients were imported.

## Evidence-to-equation routing

| Model choice | Minimum sources to read before changing code | Current status |
|---|---|---|
| Carbon compartment signs and capacity | Luo et al.; Sierra & Müller; Wei et al. | Implemented with tests |
| Native VISITc water bookkeeping | `point/hydro_balance.c`, `structure.h`, `location_proc.c`, `daily_scheme.c` | Partial source-order transcription implemented |
| Soil \(W\leftrightarrow\psi\) | Current VISITc `location_proc.c` texture branches; SurEau-Ecos as structural context | Native VISITc diagnostic relation implemented |
| Plant \(W\leftrightarrow\psi\) and capacitance | SurEau-Ecos; FETCH2/FETCH3 | Reduced SurEau-Ecos symplasmic P-V relation implemented; parameters uncalibrated |
| Stomata/photosynthesis-hydraulic coupling | VISITc PM/ecophysiology source; Kennedy; Sperry or SOX | Native `f_canopy_cond` and PM dependencies implemented; reduced-ODE closure remains idealized |
| Dynamic vs quasi-steady water comparison | SurEau-Ecos capacitance-removal design; Jacobian Schur complement | Analysis tools implemented; experiment pending |
| Periodic reference and R-tipping | Ashwin; Alkhayuon & Ashwin; Feudel | Deferred until ordinary dynamics are validated |

## Reading protocol for additions

For each new paper, record:

1. full citation and stable DOI/publisher URL;
2. which state, flux, constitutive relation, parameter range, or validation
   method is being used;
3. time and area units;
4. whether the equation is copied, adapted, inferred, or used only as context;
5. the corresponding code and test file;
6. a local PDF only when access and redistribution terms permit it.
