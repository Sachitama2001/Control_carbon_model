# VISIT soil modal and QSE analysis

## Scope

This analysis uses the validated fixed-environment 9-pool daily soil system for
TKY, with the MOD12 class 4 Deciduous Broadleaf Forest soil defaults. It is a
source-grounded carbon-subsystem analysis, not a full VISIT ecosystem result.

## Grouped discrete modes

The daily matrix has repeated tree/grass litter poles. Individual eigenvectors
inside a repeated eigenspace are not unique, so the implementation groups equal
poles and computes a spectral projector for each group. The diagonal of that
projector is used as a basis-independent state-participation measure. Rh
residues are computed as `C P_i B` for each grouped projector `P_i`.

For a stable discrete pole `lambda` with daily timestep, the reported decay
time is

```text
tau = -1 / log(abs(lambda)) days.
```

At the current test environment (upper/deep soil temperatures 10/8 C), the TKY
soil modes are approximately:

| Pole | Multiplicity | Decay time (years) |
| ---: | ---: | ---: |
| 0.99925948 | 2 | 3.70 |
| 0.99949510 | 2 | 5.42 |
| 0.99979804 | 2 | 13.56 |
| 0.99990698 | 1 | 29.45 |
| 0.99995866 | 1 | 66.27 |
| 0.99998708 | 1 | 212.07 |

The grouped residue sum reconstructs the native-order Rh impulse response, and
`D + sum(R_i / (1-lambda_i))` reconstructs the DC gain. The measured matrix and
projector reconstruction error is about `6.2e-14`.

These are mathematical modes of the frozen 9-pool system. The repeated litter
groups can be associated with leaf, stem, and root turnover coefficients, but
the humus-coupled eigenvectors should not be assigned a unique ecological label
without inspecting their grouped participation.

## Fixed-point environmental sensitivity

For constant litter input, the fixed point satisfies

```text
(I - A(e)) x*(e) = B u.
```

Within one smooth decomposition branch, implicit differentiation gives

```text
(I - A) dx*/de = (dA/de) x*.
```

The implementation obtains `dA/de` and `dC/de` by central differences of the
source-grounded daily system and solves the implicit equation. It independently
re-solves the upper/lower perturbed fixed points as a validation reference.

For the numerical fixture with one unit/day in each litter channel:

| Driver | d(total soil C*)/d(driver) |
| --- | ---: |
| Upper-soil temperature | -1612.6 Mg C ha-1 K-1 |
| Deep-soil temperature | -6721.1 Mg C ha-1 K-1 |
| Upper-soil aperture | -429.8 Mg C ha-1 per fraction |
| Whole-soil aperture | -214.1 Mg C ha-1 per fraction |

The absolute carbon values are not a calibrated TKY estimate: the litter input
of six Mg C ha-1 day-1 is deliberately simple and unrealistically large for a
hand-checkable numerical fixture. The signs and numerical identities are the
validated result, not the ecological magnitude.

Equilibrium Rh sensitivity is numerically zero. At steady state and with no
other carbon export, all constant litter input must eventually leave as Rh;
temperature changes equilibrium stocks and turnover times, but not that
steady-state carbon-throughput identity.

The fixture's moisture response is aperture-limited in both layers. Therefore
small local changes in water amount or field capacity do not change the active
`min(water_factor, aperture_factor)` branch, and their local sensitivities are
zero. This is branch-specific and must not be generalized to water-limited
conditions.

## Numerical safeguards

- repeated poles are analyzed as groups;
- non-diagonalizable Jordan systems are rejected;
- modal residues must reconstruct both IRF and DC gain;
- implicit QSE sensitivity is checked against re-solved fixed points;
- finite differences crossing the `-20 C` threshold or a moisture/aperture
  branch switch are rejected as smooth Jacobians.