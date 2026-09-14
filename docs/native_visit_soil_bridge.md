# Native VISIT soil perturbation bridge

## Purpose

This bridge validates the Python 9-pool soil subsystem against execution of the
authoritative VISIT C functions, rather than against another Python equation.
It compiles these files directly from the pinned checkout:

- `visit_local/soil_proc.c::f_cycle_soil`;
- `visit_local/decomposition.c::frl`;
- `visit_local/decomposition.c::frh`.

The checkout HEAD must equal
`3285bd8e131a932e338b59892751648fd9edcc7b`, and the compiled C/header files
must have no tracked differences from that commit; otherwise the build stops.

## Scope

The bridge controls and records:

- initial values of the nine litter/humus carbon pools;
- six daily litter carbon inputs;
- the 12 carbon turnover/partition parameters in `Schar`;
- `kml`, `kmh`, `kmsl`, and `kmsh`;
- upper/deep soil temperature, water, aperture, and field capacity;
- daily updated pools, Rh, `f_tm_l`, and `f_tm_h`.

The bridge accepts all eight decomposition environment values independently
for every day. It covers baseline `EX_DECTMP=0`, where `frl` and `frh` both
receive `stype=0`. It does not execute the undefined `frh(stype=1/2)` source
path.

The pinned source has `N_CYCLE=1`, so `f_cycle_soil` references nitrogen
routines after completing its carbon update. The bridge links those five
nitrogen routines as explicit no-op stubs. Stable-isotope code is disabled by
the pinned `SCI_SCHEME=0`. Consequently this is a native-source carbon-subsystem
validation, not a full soil C-N-isotope validation.

## Parameter fixture

Regression tests use the Deciduous Broadleaf Forest (MOD12 class 4) Soil
column for TKY from the supplementary workbook
`visit_local/INPUT/parameter_VISITc_16.xlsx`, SHA-256
`f8748b9dca3a7e7e38a2aa7ce93fd54c22653948f2467fd2ed7c1f6181edbad3`.
The C routine `parameter.c::set_parameter` remains authoritative for symbol
interpretation; the test passes the extracted values directly to `Schar`.

## Validation contracts

`tests/test_visit_native.py` compiles the bridge during the test session and
checks:

1. one native C day against the Python source-order update;
2. a 40-day litter-pulse trajectory against repeated Python updates;
3. a 100-day native baseline-minus-perturbation trajectory against the exact
   fixed-coefficient soil IRF;
4. a 120-day synthetic seasonal temperature/water trajectory;
5. native and nonlinear Python temperature-pulse responses against a daily
   tangent linearization;
6. explicit rejection of linearization at the `-20 C` source threshold.

The litter pulse must have zero Rh effect on its input day and a positive Rh
effect on the following day, preserving VISIT's pre-input decomposition order.
The 100-day state-difference tolerance is `2e-10`; the observed error is about
`8.03e-11 Mg C ha-1`, dominated by subtraction of separately evolved native trajectories
whose fixed-point stocks are around `1e3`. The Rh-difference tolerance is
`1e-13`; the observed maximum is about `2.94e-15 Mg C ha-1 day-1`.

For the 120-day varying environment, the maximum native-minus-Python state
difference is `1.42e-13 Mg C ha-1` and the maximum Rh difference is
`5.55e-17 Mg C ha-1 day-1`.

For an upper-soil temperature pulse, the maximum state errors of the tangent
prediction are:

```text
2.0 K pulse: 4.09e-4 Mg C ha-1
1.0 K pulse: 1.02e-4 Mg C ha-1
0.5 K pulse: 2.54e-5 Mg C ha-1
```

Halving the perturbation reduces error by approximately four, as expected for
the remainder after a first-order approximation in a smooth regime.

## Run

```bash
pytest -q tests/test_visit_native.py
```

The comparison figure is generated with:

```bash
PYTHONPATH=src python examples/plot_visit_soil_temperature_pulse.py
```

It writes `artifacts/visit_soil_temperature_pulse.png` with forcing, native C /
nonlinear Python / tangent Rh responses, the total-litter versus total-humus
projected perturbation trajectory, and approximation errors.

The test requires a C compiler and the sibling checkout
`/mnt/d/CT/VISIT-matrix/visit_local`. The reusable Python entry points are:

- `build_native_visit_soil_bridge`;
- `run_native_visit_soil`;
- `compare_native_visit_soil_to_python`;
- `compare_native_visit_soil_perturbation_to_irf`.

## Next extension

The next extension is modal/QSE analysis and then a controlled native
meteorological/hydrological trajectory. Full executable validation remains
blocked on the native input bundle and known initialization/format defects in
`current_status.md`.