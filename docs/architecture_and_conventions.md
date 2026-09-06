# Architecture and Conventions

## 1. Package boundary

The package should separate **generic systems mathematics** from **ecosystem-model adapters**.

Generic code must not import VISIT-specific modules.

Recommended dependency direction:

```text
core/provenance  <- adapters/visit <- experiments
```

Never allow:

```text
core -> adapters/visit
```

## 2. Terminology

Use these terms consistently.

### Native model

Original ecosystem-model implementation and update order.

### Reduced model

A lower-dimensional or coefficient-frozen representation derived from the native model.

### Operating point

The state/forcing/regime at which a local linearization is evaluated.

### QSE

An equilibrium corresponding to frozen forcing. For discrete models use fixed-point language where appropriate.

### IRF

Impulse response of an explicitly identified linearized/reduced input-output system. Do not call arbitrary perturbation trajectories an IRF unless linearity/locality is specified.

### Effective forcing

An internal driver such as EPP or a decomposition scalar. Keep distinct from meteorological forcing such as temperature or precipitation.

## 3. Time-domain conventions

Every system object must carry time-domain metadata:

```text
continuous
or
discrete with dt and unit
```

VISIT's native carbon update is daily and sequential. Continuous-time forms are mathematical approximations/reformulations and must be labelled as such.

## 4. Sign conventions

Recommended default:

- carbon stocks positive;
- GPP/NPP into ecosystem positive;
- heterotrophic/autotrophic respiration to atmosphere positive as flux magnitudes;
- NEP positive = ecosystem sink / net carbon gain;
- matrix donor columns and recipient rows when using compartment transfer matrices.

Document any adapter that differs.

Do not assume NEE sign convention from observational datasets; make NEE convention explicit whenever added.

## 5. Units

Core mathematics should be unit-agnostic numerically but metadata must record units.

VISIT source commonly uses:

```text
carbon stock: Mg C ha^-1
carbon flux: Mg C ha^-1 day^-1
```

Do not silently convert to kg C m^-2. If conversion is desired, do it in an explicit adapter/conversion layer and test the factor.

## 6. Data structures

Prefer immutable dataclasses for mathematical systems and schemas.

Suggested concepts:

```python
@dataclass(frozen=True)
class StateVariable:
    name: str
    unit: str
    category: str
    source_refs: tuple[SourceRef, ...]

@dataclass(frozen=True)
class OperatingPoint:
    x: np.ndarray
    u: np.ndarray
    metadata: Mapping[str, Any]
```

Avoid passing unlabelled arrays across model boundaries when state ordering matters. An array may be used internally for performance, but schemas must define its ordering.

## 7. Provenance

Each source-grounded model component must provide provenance.

Provenance should answer:

1. Where did this equation come from?
2. What source variable/function does it correspond to?
3. Was it copied exactly or approximated?
4. What branch/configuration was assumed?
5. What units were used?

Where possible, pin source repository commit SHA.

## 8. Exact vs approximate code

Use names that advertise approximation.

Good:

```text
visit_soil_daily_exact_from_source
VISITReduced18Adapter
linearize_native_daily_map
annual_mean_lti
```

Bad:

```text
visit_model
exact_visit_matrix
```

unless the implementation is genuinely exact under documented configuration.

## 9. Testing hierarchy

### Unit tests

Pure process equations against hand-calculated values.

### Structural tests

Mass conservation, dimensions, positivity where expected, topology.

### Equivalence tests

Direct algebra vs matrix form.

### Linearization tests

Finite perturbation vs Jacobian prediction.

### Native validation tests

Reduced/local model vs original VISIT executable perturbations.

### Scientific regression tests

Known poles, QSE, NEP signs, IRF features for frozen fixtures.

## 10. Numerical differentiation

Implement a central finite-difference reference with configurable relative/absolute step.

For state component `x_i`, choose perturbation from scale metadata where possible. Avoid a single absolute epsilon across carbon pools with very different magnitudes.

Tests should examine derivative convergence across step sizes.

## 11. QSE solving

QSE routines must report diagnostics rather than only return a vector:

```text
converged
residual norm
iterations/method
stability of linearization
conditioning
active regime/branch
```

For hybrid/piecewise models, a root near a switching surface requires a warning.

## 12. Modal analysis

Do not equate a raw matrix diagonal element with ecosystem residence time.

Keep distinct:

- pool turnover rate;
- eigenvalue/pole;
- modal time scale;
- mean transit/residence quantities;
- Luo-style traceability residence-time matrices.

These can be related but are not interchangeable.

For continuous stable pole `lambda`, modal decay time may be reported as

```text
tau = -1 / Re(lambda)
```

when meaningful. For discrete poles, convert using the known timestep only with an explicitly documented mapping.

## 13. IRF conventions

For continuous LTI:

```text
h(t) = C exp(A t) B, t > 0
```

with direct feedthrough `D delta(t)` treated separately.

For discrete LTI:

```text
h[0] = D
h[k] = C A^(k-1) B, k >= 1
```

if the convention is input at k affecting x[k+1]. Lock this convention in tests and docstrings.

## 14. Plotting and notebooks

Core library should return arrays/data objects, not create figures implicitly.

Use `examples/` or `notebooks/` for exploratory figures after numerical core is tested.

Every published figure should be reproducible from a script with fixed inputs/configuration.

## 15. Cross-model adapter contract

Do not force all models into 18 pools. The common interface should tolerate arbitrary dimensions.

An adapter should eventually expose:

```text
state schema
forcing schema
output schema
step/rhs
QSE/fixed-point operation where meaningful
linearization
provenance
```

Common comparison occurs at derived input-output quantities, not identical state dimension.

## 16. Git conventions

Prefer focused commits such as:

```text
Add discrete LTI impulse response
Reconstruct VISIT soil decomposition matrix
Add source provenance for frl/frh
Validate VISIT soil mass balance
```

Avoid mixing scientific-equation changes and unrelated refactors in one commit.

Any change to a source-derived equation should update its provenance note/test in the same PR.

## 17. Scientific failure is an acceptable result

If the 18-pool reduction cannot reproduce native VISIT IRFs because omitted NSC/phenology/hydrology states carry memory, document and quantify the failure instead of tuning coefficients to hide it. Identifying the minimal state needed for input-output equivalence is itself a valuable result.
