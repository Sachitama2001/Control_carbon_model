# Four-component carbon-water equations

## Status and intended use

This document specifies the first executable carbon-water model in
`src/control_carbon/coupled_carbon_water.py`. It is a **source-grounded reduced
continuous model**, not an exact rewrite of VISITc. The exact/partial VISITc
water transcription is kept separately in `visitc_hydrology.py`.

The purpose of this model is to connect four things in one inspectable system:

1. water mass storage and internal SPAC transport;
2. water limitation of carbon input and decomposition;
3. carbon control of hydraulic transport and transpiring area;
4. matrix diagnostics of carbon capacity, disequilibrium, and response time.

Default values are numerical examples only. They are not calibrated VISITc,
species, biome, or site parameters.

## State, units, and forcing

The eight-state vector is

\[
\boldsymbol{x}
=
\begin{pmatrix}\boldsymbol C\\\boldsymbol W\end{pmatrix},\qquad
\boldsymbol C=
\begin{pmatrix}C_l\\C_t\\C_r\\C_s\end{pmatrix},\qquad
\boldsymbol W=
\begin{pmatrix}W_l\\W_t\\W_r\\W_s\end{pmatrix}.
\]

Subscripts denote leaf, stem, root, and soil. Carbon is in Mg C ha\(^{-1}\),
water is in mm, and time is in days. The forcing object contains precipitation,
a photosynthesis multiplier, a decomposition multiplier, potential
transpiration, and potential soil evaporation.

The first experiment excludes snow, interception storage, separate soil
layers, nitrogen and phosphorus states, fire, and harvest. This scope is a
testable starting point rather than a claim that those processes are minor.

## Carbon balance

At a fixed water state and forcing, carbon is written as

\[
\dot{\boldsymbol C}
=\boldsymbol u_C(\boldsymbol W,t)
+M_C(\boldsymbol W,t)\boldsymbol C.
\]

With allocation fractions \(a_l+a_t+a_r=1\), plant turnover rates
\(k_l,k_t,k_r\), and effective soil loss rate \(k_s f_d\),

\[
\boldsymbol u_C=
\begin{pmatrix}a_l\mu\\a_t\mu\\a_r\mu\\0\end{pmatrix},\qquad
M_C=
\begin{pmatrix}
-k_l&0&0&0\\
0&-k_t&0&0\\
0&0&-k_r&0\\
k_l&+k_t&+k_r&-k_s f_d
\end{pmatrix}.
\]

Plant turnover is an internal transfer to the aggregated soil pool; the only
carbon loss in this first model is heterotrophic respiration
\(R_h=k_s f_d C_s\). Therefore

\[
\frac{d}{dt}\boldsymbol 1^\top\boldsymbol C=\mu-R_h.
\]

The executable model tests this equality at every right-hand-side evaluation.

Water limitation is represented by bounded half-saturation functions. Let
\(\theta_i=W_i/W_i^{\rm ref}\) and \(h(z;a)=z/(z+a)\) for \(z\ge0\). Then

\[
\mu=\mu_{\max} f_{\rm photo}(t)
\sqrt{h(\theta_l;a_p)h(\theta_s;a_p)},
\]

\[
f_d=f_{\rm decomp}(t)h(\theta_s;a_d).
\]

These smooth choices are placeholders whose replacements must be accompanied
by provenance and tests.

## Water balance and incidence structure

The internal path is soil \(\leftrightarrow\) root \(\leftrightarrow\) stem
\(\leftrightarrow\) leaf. With

\[
\boldsymbol q=
\begin{pmatrix}q_{sr}\\q_{rt}\\q_{tl}\end{pmatrix},
\quad
S_W=
\begin{pmatrix}
0&0&1\\
0&1&-1\\
1&-1&0\\
-1&0&0
\end{pmatrix},
\]

the water equation is

\[
\dot{\boldsymbol W}
=S_W\boldsymbol q+\boldsymbol p-\boldsymbol e.
\]

Because \(\boldsymbol 1^\top S_W=0\), internal transport cancels exactly.
The external vectors are

\[
\boldsymbol p=(0,0,0,P)^\top,
\]

\[
\boldsymbol e=(E_l,0,0,E_s+D+R)^\top,
\]

so the model must satisfy

\[
\frac{d}{dt}\boldsymbol 1^\top\boldsymbol W
=P-E_l-E_s-D-R.
\]

The code reports the residual of this identity and tests it near machine
precision.

Internal flows use differences in normalized storage as the first, deliberately
simple proxy for water-potential differences:

\[
\begin{aligned}
q_{sr}&=G_{sr}(\boldsymbol C)(\theta_s-\theta_r),\\
q_{rt}&=G_{rt}(\boldsymbol C)(\theta_r-\theta_t),\\
q_{tl}&=G_{tl}(\boldsymbol C)(\theta_t-\theta_l).
\end{aligned}
\]

The conductances are multiplied by bounded functions of root, stem, and leaf
carbon. Thus carbon changes water transport without pretending that carbon is
physically transferred into water. Negative \(q\) is allowed and represents
reverse equilibration; it is not clipped.

For research-grade SPAC interpretation, normalized storage must eventually be
replaced by explicit constitutive relations

\[
\psi_s=\psi_s(W_s),\qquad
\psi_i=\psi_i(W_i,C_i),\qquad
q_{ij}=K_{ij}(\boldsymbol C,\boldsymbol W)(\psi_i-\psi_j),
\]

using a soil retention curve and plant pressure-volume/capacitance relations.
SurEau-Ecos and FETCH models are implementation references; see
`literature_index.md`.

## Carbon capacity is not the coupled equilibrium

The frozen-coefficient diagnostic is

\[
\boldsymbol C_{\rm cap}(t)=-M_C(t)^{-1}\boldsymbol u_C(t).
\]

Define total capacity, current stock, and signed storage potential as

\[
X_c=\boldsymbol1^\top\boldsymbol C_{\rm cap},\qquad
X=\boldsymbol1^\top\boldsymbol C,\qquad
X_p=X_c-X.
\]

`instantaneous_carbon_capacity` computes these frozen-water quantities. They
are diagnostic targets, not necessarily fixed points of the full system. A
coupled frozen-forcing equilibrium instead solves

\[
F_C(\boldsymbol C^*,\boldsymbol W^*;\boldsymbol u)=0,
\qquad
F_W(\boldsymbol C^*,\boldsymbol W^*;\boldsymbol u)=0
\]

simultaneously. With periodic forcing, neither object is the full tracking
target; a periodic attracting trajectory is the relevant reference.

If \(X_c=\mu\tau_E\), its finite change is decomposed exactly as

\[
\Delta X_c
=\tau_{E,0}\Delta\mu
+\mu_0\Delta\tau_E
+\Delta\mu\,\Delta\tau_E.
\]

This Wei et al.-style decomposition is implemented and tested without assuming
in advance which term dominates under drought.

## Local coupled dynamics

At a frozen-forcing equilibrium, partition the full Jacobian as

\[
J=
\begin{pmatrix}
J_{CC}&J_{CW}\\
J_{WC}&J_{WW}
\end{pmatrix}.
\]

`coupled_jacobian_blocks` evaluates these blocks by centered finite differences.
Both off-diagonal blocks should be nonzero:

- \(J_{CW}\): water changes photosynthesis and decomposition;
- \(J_{WC}\): plant carbon changes conductance and transpiration.

If water anomalies are stable and fast enough to approximate
\(\delta\dot{\boldsymbol W}=0\), then

\[
\delta\boldsymbol W\simeq-J_{WW}^{-1}J_{WC}\delta\boldsymbol C,
\]

and the reduced carbon Jacobian is the Schur complement

\[
J_{\rm eff}=J_{CC}-J_{CW}J_{WW}^{-1}J_{WC}.
\]

This relation creates a direct test of whether explicit water storage changes
carbon response times beyond a quasi-steady water approximation.

## Numerical and scientific guardrails

- The equilibrium solver enforces nonnegative states and reports residuals,
  eigenvalues, stability, and conditioning.
- The integrator uses DOP853; `max_step` must resolve forcing discontinuities.
- `examples/run_coupled_carbon_water.py` compares two maximum step sizes and
  writes the difference to its manifest.
- A transient pulse, mortality event, or long recovery time is not by itself
  evidence of rate-induced tipping. A tracking target, basin distinction, and
  forcing-rate threshold are required.
- Before interpreting fitted parameter values, replace placeholder hydraulic
  relations, add dimension checks, and validate against a source model or data.

