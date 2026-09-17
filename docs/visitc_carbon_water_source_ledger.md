# VISITc carbon-water source ledger

## Provenance boundary

This ledger refers to [`visit-manager/VISITc`](https://github.com/visit-manager/VISITc)
at commit
[`5202debd96df6f88beb7d61f8688fff02ace964a`](https://github.com/visit-manager/VISITc/tree/5202debd96df6f88beb7d61f8688fff02ace964a).
It is distinct from the older `Sachitama2001/VISIT-matrix/visit_local`
snapshot used by the existing nine-pool soil-carbon adapter. Do not combine
coefficients from the two snapshots without a documented comparison.

`src/control_carbon/visitc_hydrology.py` transcribes the source-order store
update below and the Penman-Monteith dependency chain. The effective-input
boundary remains available so store algebra and PM algebra can be tested independently.
`src/control_carbon/coupled_carbon_water.py` is a separate four-component
continuous synthesis and is not expected to reproduce a VISITc day.

## Native stores and diagnostics

| Quantity | C expression | Units/comment in source | Declared/read at | Role in this repository |
|---|---|---|---|---|
| Snow water | `mass->snwa` | mm water equivalent | `point/structure.h::struct Mass`; `hydro_balance.c` | Native daily state in hydrology transcription; excluded from first reduced ODE |
| Upper soil water | `mass->sw30` | mm, above 30 cm | same | Native daily state; source of soil evaporation and C3/C4 transpiration |
| Whole soil water | `mass->sww` | mm, source comment says whole soil | same | Native daily state; source of tree transpiration and baseflow |
| Upper water diagnostic | `loct->soilwtr_l` | mm water | assigned after hydrology in `location_proc.c` | Carbon-process/effective driver |
| Whole water diagnostic | `loct->soilwtr_h` | mm water | assigned after hydrology in `location_proc.c` | Carbon-process/effective driver |
| Upper/whole aperture | `soilappr_l`, `soilappr_w` | fraction | `location_proc.c` | Moisture stress/decomposition drivers; not independent stores |
| Water-filled pore space | `loct->wfps` | fraction | `location_proc.c` | GHG/decomposition-related diagnostic |
| Matric/total potential | `pot_matric_l/h`, `pot_total_l/h` | MPa | `location_proc.c` | Constitutive relation candidate for later SPAC reconstruction |
| Previous moisture | `m_casa_pre`, `vmc_pre` | source-specific | end of `daily_scheme.c` | Memory variables required for an exact Markov state where used |

The names `sw30` and `sww` should be retained in exact audits. Although `sww`
is commented as whole-soil storage, the update algebra transfers `ro1` into it
and applies an exchange `retran` between it and `sw30`. Treating the two as
strictly disjoint layers is therefore an interpretation, not a source fact.

## `f_hydrology` process ledger

| Source order | Flux/update | Main source algebra | Donor -> receiver/outside | Current status |
|---:|---|---|---|---|
| 1 | Snow fraction | `1/(1+exp(0.75*(tmp_2m-2)))` | precipitation partition | Native-C validated |
| 2 | Thaw | logistic temperature factor times snow store; separate small-snow branch | snow -> liquid input | Native-C validated |
| 3 | Canopy interception | LAI capacity plus quadratic supply-potential limitation, separately for tree/C3/C4 | rain -> atmosphere | Store update exact; PM potentials are effective inputs |
| 4 | Liquid soil input | `(rain - incep) + thaw` | atmosphere/snow -> upper bucket | Exact |
| 5 | `ro1` | cubic bucket overflow using `fieldcap30 - sw30` | upper -> `sww` path | Exact, including exponent `0.33333` |
| 6 | Soil evaporation | quadratic limitation of `pm_evpr` by updated `sw30` | upper -> atmosphere | Exact given effective potential |
| 7 | C3/C4 transpiration | each potential limited by the same current `sw30`; sum then removed | upper -> atmosphere | Exact given effective potentials |
| 8 | Tree transpiration | tree potential limited by `sww` | `sww` -> atmosphere | Exact given effective potential |
| 9 | Baseflow | `0.001*sww`, `0.003*sww` at `QHB`, zero for frozen deep soil | `sww` -> outside | Exact and separately exposed |
| 10 | `ro2` | cubic bucket overflow of `ro1`, then baseflow added | `sww` -> outside | Exact |
| 11 | `retran` | capacity-ratio imbalance, asymmetric factor, hydraulic-conductivity cap | exchange between stores | Exact |
| 12 | Final slow balance | `sww += ro1 - ro2 - trnsp_tree` | mixed | Exact; audit warning below |

The source function is
[`point/hydro_balance.c::f_hydrology`](https://github.com/visit-manager/VISITc/blob/5202debd96df6f88beb7d61f8688fff02ace964a/point/hydro_balance.c#L21-L238).

## Confirmed accounting issue to preserve and test

At lines 184-202, positive `baseflow` is subtracted once from `mass->sww` and
then added to `loct->ro2`. At line 230, all of `ro2` is subtracted again. In a
no-clipping case, direct algebra gives

\[
\Delta W_{\rm stores}
-\left(P-AET-ro2\right)=-baseflow.
\]

Therefore baseflow contributes an additional source-level loss. The Python
transcription deliberately reproduces this behavior and reports
`water_budget_residual` and `expected_source_residual`. It does **not** assert
that the source is wrong in intent, and it does not silently repair the update.
Any corrected experiment must be a named alternative with native-vs-corrected
tests.

Clipping at lines 222-226 also introduces water when an intermediate store is
negative. Those clip corrections are reported separately.

## Call order and carbon coupling

In `point/location_proc.c`, the observed order is:

1. radiation and net-radiation calculations;
2. `f_hydrology(...)`;
3. copy `snwa`, `sw30`, and `sww` into location diagnostics;
4. calculate aperture, WFPS, and soil water potentials;
5. call `f_ecophysiology(...)` for tree, C3, and C4;
6. later photosynthesis/decomposition processes use the prepared diagnostics.

Thus water affects carbon through state-dependent physiological and
environmental coefficients, not through carbon-water mass transfers. Carbon
already affects water through LAI and canopy conductance entering interception
and transpiration. A coupled Jacobian should therefore have both
\(J_{CW}\ne0\) and \(J_{WC}\ne0\).

The repository README calls VISITc a 30-minute model, while the inspected
execution path contains `daily_scheme()` and day-unit flux comments. Do not
infer one universal timestep from the README or function name alone. Trace the
chosen executable/compile flags and record its actual call frequency before
converting coefficients to continuous rates.

## What remains to read before a faithful carbon-water adapter

| Priority | Source work | Acceptance criterion |
|---|---|---|
| complete | Transcribe `pm_incep`, `pm_evap`, and `pm_transp` plus resistance/radiation/LAI/canopy-conductance dependencies | Native C and Python agree for hand-picked daily inputs |
| complete | Compile minimal `f_hydrology` and PM bridges | One-day states and all exposed fluxes agree, including normal, frozen, QHB, and clipping branches |
| P1 | Trace `gc`, LAI, photosynthesis moisture response, and decomposition scalars | Each carbon-water arrow has file, function, variable, unit, and update order |
| P1 | Determine runtime meaning of `sww` from initialization, outputs, and parameter generation | Layer interpretation is evidence-backed |
| P1 | Add snow-free and snow-enabled mass-balance tests | Residuals are explained by named source branches only |
| P2 | Compare original and explicitly corrected baseflow variants | Scientific results are not conditioned on an unnoticed bookkeeping choice |
| P2 | Map native daily memory variables | Exact Markov state list is documented |
