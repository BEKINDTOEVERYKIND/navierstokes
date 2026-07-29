# Cheverry audit of the single-carrier material phase

Date: 2026-07-29

## Verdict

There are two different issues in Cheverry's phase cascade, and they have
different answers for the proposed single-carrier construction.

1. A \(K\)-dependent phase transported by the **complete \(K\)-dependent
   zero-charge velocity** exactly resums Cheverry's *geometrical* phase
   cascade, provided that the coupled phase/profile system is actually
   constructed.  Expanding that one material phase reproduces Cheverry's
   recursive eikonal equations coefficient by coefficient.
2. Cheverry's *adjusting* phases are a representation gauge.  They can be
   moved exactly from the phase into a \(K\)-dependent translation of the
   periodic profile.  They therefore do not require genuinely new fast
   directions.  They do, however, still have to be fixed, explicitly or
   implicitly, to triangularize a coefficient recursion.
3. This does **not** close the normalized-amplitude, order-one-time
   construction.  After the material phase removes low--high \(O(K)\)
   transport and incompressibility removes \(K\) from high--high
   transport, the leading slow modulation contains the hydrostatic Euler
   system.  That system loses one slow derivative and has unbounded
   unstable spectra for some analytic shears.  Thus the present
   Gevrey-2 interaction-order majorant is not a Cauchy estimate for the
   profile PDE.

The phase repair is real, but the current proposed theorem must add an
**analytic or spectrally stable hydrostatic-profile hypothesis**, or find
an invariant subclass which avoids hydrostatic Euler.  The material
pressure parametrix is not the only missing estimate.

Throughout, write

\[
\delta=K^{-1}
\]

for wavelength and \(\nu_K\) for viscosity, so that the regime in the
working note is

\[
\nu_K K^2\longrightarrow0.
\tag{C.1}
\]

This avoids confusing Cheverry's small parameter with viscosity.

---

## 1. The exact one-phase system

Use the lifted ansatz

\[
\begin{aligned}
u^K(t,x)
&=
U^K(t,x)
+W^K(t,x,K\phi^K(t,x)),\\
p^K(t,x)
&=
P^K(t,x)
+\Pi^K(t,x,K\phi^K(t,x)),
\end{aligned}
\tag{C.2}
\]

with

\[
\langle W^K\rangle_\vartheta
=\langle\Pi^K\rangle_\vartheta=0.
\tag{C.3}
\]

Put

\[
\xi^K=\nabla\phi^K,\qquad
{\cal D}_K=\nabla_x+K\xi^K\partial_\vartheta
\tag{C.4}
\]

and transport the phase by the complete zero-charge field:

\[
\partial_t\phi^K+U^K\cdot\nabla\phi^K=0.
\tag{C.5}
\]

Exact physical incompressibility is equivalent to

\[
\operatorname{div}U^K=0,\qquad
\operatorname{div}_xW^K
+K\xi^K\cdot\partial_\vartheta W^K=0.
\tag{C.6}
\]

For a mean-zero periodic function, let
\(\partial_\vartheta^{-1}\) denote the mean-zero inverse.  Equation
(C.6) gives the exact identity

\[
K\,\xi^K\cdot W^K
=-\partial_\vartheta^{-1}\operatorname{div}_xW^K.
\tag{C.7}
\]

Consequently define the \(K\)-free quadratic transport

\[
{\cal Q}(W)
:=
(W\cdot\nabla_x)W
-
\big(
\partial_\vartheta^{-1}\operatorname{div}_xW
\big)\partial_\vartheta W.
\tag{C.8}
\]

Substitution into forced Navier--Stokes and use of (C.5)--(C.7) give the
exact mean equation

\[
\begin{aligned}
\partial_tU^K
+(U^K\cdot\nabla)U^K
+\nabla P^K
+\langle{\cal Q}(W^K)\rangle
=\nu_K\Delta U^K+\langle F^K\rangle,
\end{aligned}
\tag{C.9}
\]

and the exact oscillatory equation

\[
\begin{aligned}
(\partial_t+U^K\cdot\nabla)W^K
&+(W^K\cdot\nabla)U^K\\
&+{\cal Q}(W^K)-\langle{\cal Q}(W^K)\rangle\\
&+\nabla_x\Pi^K+K\xi^K\partial_\vartheta\Pi^K\\
&=\nu_K{\cal L}_{K,\phi^K}W^K+F^{K,*},
\end{aligned}
\tag{C.10}
\]

where

\[
\begin{aligned}
{\cal L}_{K,\phi}W
=\Delta_xW
&+2K\nabla\phi\cdot\nabla_x\partial_\vartheta W\\
&+K(\Delta\phi)\partial_\vartheta W
+K^2|\nabla\phi|^2\partial_\vartheta^2W.
\end{aligned}
\tag{C.11}
\]

Equations (C.5), (C.6), (C.9), and (C.10) are the actual coupled closure
problem.  They verify the valuable algebraic claim in the working note:
there is no naked \(K\) in the nonlinear transport after the constraint
is used.  They do not by themselves give a \(K\)-uniform Cauchy theory.

---

## 2. Cheverry's geometrical phases are resummed

Cheverry uses, for fixed finite \(l\), a geometrical phase

\[
\phi_g^\delta
=
\phi_0+\sum_{n=1}^{l-1}\delta^{n/l}\phi_n.
\tag{C.12}
\]

His equations (3.2) and (4.9) are the coefficient equations for an
approximate material eikonal.

This follows directly from (C.5).  Suppose formally that

\[
U^\delta=\sum_{m\ge0}\delta^{m/l}\bar U_m,\qquad
\phi^\delta=\sum_{n\ge0}\delta^{n/l}\phi_n.
\tag{C.13}
\]

The coefficient of \(\delta^{n/l}\) in
\(\partial_t\phi^\delta+U^\delta\cdot\nabla\phi^\delta=0\) is

\[
(\partial_t+\bar U_0\cdot\nabla)\phi_n
+
\sum_{m=1}^{n}
\bar U_m\cdot\nabla\phi_{n-m}
=0.
\tag{C.14}
\]

This is exactly Cheverry's geometrical hierarchy.  In particular, at the
second nontrivial level it contains

\[
\partial_t\phi_2
+(\bar U_0\cdot\nabla)\phi_2
+(\bar U_1\cdot\nabla)\phi_1
+(\bar U_2\cdot\nabla)\phi_0
=0,
\tag{C.15}
\]

the displayed mechanism in his equation (3.10).

Thus one should not freeze \(\phi\) at the phase transported only by the
initial or affine low field.  But if (C.5) is solved with the complete
generated mean \(U^K\), one \(K\)-dependent phase contains every
geometrical correction in (C.12).  Cheverry's infinitely many
coefficient phases are then the asymptotic expansion of one function;
they are not independent fast angles.

There is an important logical qualifier: (C.14) proves equivalence of
the formal hierarchies.  It does not prove that (C.13) converges, is
Gevrey-summable with the required constants, or remains nondegenerate.

---

## 3. The adjusting phases are gauge, not extra geometry

For any smooth scalar \(a^K(t,x)\), set

\[
\widetilde\phi^K
=\phi^K+\frac{a^K}{K},
\qquad
\widetilde W^K(t,x,\vartheta)
=W^K(t,x,\vartheta-a^K(t,x)).
\tag{C.16}
\]

Then, exactly,

\[
\widetilde W^K(t,x,K\widetilde\phi^K(t,x))
=W^K(t,x,K\phi^K(t,x)).
\tag{C.17}
\]

The physical field has not changed.  Cheverry's complete phase is

\[
\phi_\natural^\delta
=\phi_g^\delta+\delta\,\phi_a^\delta.
\tag{C.18}
\]

Since \(K=\delta^{-1}\), (C.18) changes the fast argument by precisely
\(\phi_a^\delta\), which is the gauge transformation (C.16).  Cheverry
explicitly rewrites the complete-phase profiles in terms of the
geometrical phase by a Taylor-expanded angular translation in his
formulas (2.6)--(2.7).

Therefore:

* adjusting phases need not be added as new material covectors;
* an exact \(K\)-dependent profile can absorb them;
* but a formal recursive proof still needs a gauge choice.  Saying only
  “transport the material phase” does not select the adjusting gauge
  which Cheverry uses to make his hierarchy triangular.

There is a small but essential qualification.  The change (C.16)
preserves the material eikonal (C.5) only when

\[
(\partial_t+U^K\cdot\nabla)a^K=0.
\tag{C.19a}
\]

For a general adjusting gauge, moving \(a^K\) into the profile creates
the compensating term
\[
-
\big((\partial_t+U^K\cdot\nabla)a^K\big)
\partial_\vartheta\widetilde W^K
\tag{C.19b}
\]
in its equation.  Thus “gauge” means equivalence of physical
representations, not that the adjusting equations may be deleted.  A
viable proof may either carry the \(\phi_a^\delta\)'s explicitly or
impose an equivalent normalization on \(W^K\) and retain (C.19b) at
every order.

---

## 4. The surviving equation is hydrostatic Euler

The exact cancellation (C.7) removes the large coefficient but leaves a
derivative-loss operator.  This can be seen without any estimate.

Take the flat phase

\[
\phi=x_3,\qquad \xi=e_3,
\tag{C.19}
\]

and profiles independent of the slow variable \(x_3\).  Write the
order-one horizontal profile as

\[
u=u(t,x_h,\vartheta)\in\mathbb R^2
\tag{C.20}
\]

and the physical oscillatory vertical component as \(K^{-1}v\).
Incompressibility gives

\[
\operatorname{div}_h u+\partial_\vartheta v=0.
\tag{C.21}
\]

The fast pressure equation gives

\[
\partial_\vartheta p=0
\tag{C.22}
\]

at leading order.  The horizontal part of (C.9)--(C.10), in the regime
(C.1), is

\[
\boxed{
\begin{aligned}
\partial_tu
+(u\cdot\nabla_h)u
+v\partial_\vartheta u
+\nabla_h p
&=f_h,\\
\operatorname{div}_h u+\partial_\vartheta v&=0,\\
\partial_\vartheta p&=0.
\end{aligned}}
\tag{C.23}
\]

The zero-charge horizontal velocity is already included in the
\(\vartheta\)-mean of \(u\).  Transporting \(\phi\) by that mean does not
remove \(v\partial_\vartheta u\), because \(v\) depends on
\(\vartheta\).

Restricting (C.23) to one horizontal coordinate and one horizontal
velocity component gives the standard two-dimensional hydrostatic Euler
system

\[
\partial_tu+u\partial_xu+v\partial_\vartheta u+\partial_xp=0,
\qquad
\partial_xu+\partial_\vartheta v=0,
\qquad
\partial_\vartheta p=0.
\tag{C.24}
\]

Equivalently, if \(\omega=\partial_\vartheta u\), then

\[
\partial_t\omega+u\partial_x\omega+v\partial_\vartheta\omega=0,
\qquad
v=-\partial_x\partial_\vartheta^{-2}\omega
\tag{C.25}
\]

after fixing the usual means.  The second formula loses one \(x\)
derivative.

This is the concrete large-amplitude closure problem behind Cheverry's
warning.  His paper states that the naive order-one one-phase modulation
equations are nonhyperbolic and Hadamard ill posed, and that an
order-one time description formally requires an infinite cascade of
phases.  The material phase resums the latter bookkeeping, but (C.23)
shows that it does not remove the nonhyperbolic profile dynamics.

Han-Kwan and Nguyen later proved nonlinear Hadamard ill-posedness for
hydrostatic Euler near some analytic stationary shears, based on
unbounded unstable spectra.  Their displayed theorem uses a bounded
vertical channel.  The same derivative-loss issue is present for the
periodic three-dimensional inviscid primitive equations, which are
exactly (C.23) on a torus; the primitive-equation literature records
ill-posedness in Sobolev spaces and in deterministic Gevrey classes of
order strictly greater than one.  The unstable growth is proportional
to slow tangential frequency.  In an unstable sector it has the form

\[
\|e^{tL_\eta}\|\gtrsim e^{\gamma|\eta|t}.
\tag{C.26}
\]

An analytic weight \(e^{\rho|\eta|}\) can pay for (C.26) for a short
time by shrinking \(\rho\) linearly.  A Gevrey-2 weight
\(e^{\rho|\eta|^{1/2}}\) cannot pay for (C.26) on any fixed positive
time.

The vanishing viscosity in (C.1) does not repair this retained band.  If
\(|\eta|\le M\ll K\), then

\[
\nu_K|\eta|^2\le \nu_KM^2
\ll \nu_KK^2=o(1),
\tag{C.26a}
\]

while the fast viscous coefficient is itself
\(\nu_KK^2=o(1)\).  Thus neither horizontal nor fast dissipation gives a
uniform hydrostatic stability estimate in the proposed scaling.

The elementary convolution estimate

\[
r^2\sum_{p=1}^{r-1}\binom rp^{-2}\lesssim1
\tag{C.27}
\]

in the working note controls **interaction order**.  It does not control
the slow-frequency semigroup in (C.26), so it is not a proof of
Gevrey-2 well-posedness for (C.10).

The statement here is adversarial but limited: hydrostatic Euler is not
ill posed near every profile, and its analytic Cauchy problem has local
solutions.  The cited ill-posedness theorem does not prove that the
particular three-pulse trajectory is unstable.  It proves that a generic
Gevrey-2 theorem with no stability hypothesis cannot follow from the
\(K\)-null identity alone.

---

## 5. What the forced approximate setting does and does not buy

The desired object is a forced approximate trajectory.  It need not be
stable to arbitrary perturbations and need not be shadowed by an exact
unforced Euler solution.  Therefore Hadamard ill-posedness is not, by
itself, a no-go theorem for the construction.

It does rule out the following proof shortcut:

> solve a generic Gevrey-2 one-phase profile Cauchy problem for
> order-one time, use only (C.27), and conclude a \(K\)-uniform
> exponentially accurate parametrix.

For a chosen trajectory, there are three plausible routes.

### Route A: analytic active transition

Work in a slow analytic norm

\[
\sum_{\eta,h}
e^{\rho(t)|\eta|}
e^{\sigma|h|}
|\widehat W_{\eta,h}(t)|,
\qquad
\rho(t)=\rho_0-\Lambda t>0.
\tag{C.28}
\]

Construct (C.5), (C.9), and (C.10) on a stage shorter than the analytic
lifespan.  Use the conic separation

\[
|\eta|\le M\ll K|h|,
\qquad
|\nabla\phi|\ge c>0
\tag{C.29}
\]

for the high-charge pressure inverse.  The inverse-\(K\) coefficient
hierarchy may still have Gevrey-2 growth in its *order*, in which case
the existing choice \(K_j=j^A,\ A>4\), remains compatible with optimal
truncation.

This route moves, rather than solves, the localization problem: a
nonzero compactly supported analytic profile does not exist.  Analytic
tails, inter-cell leakage, and the annular wake seam would need a new
quantitative design.

#### A viable near-identity version

There is a concrete way in which this route can genuinely pass the
hydrostatic regularity gate.

Choose a fixed cascade ratio

\[
r=1+\delta,\qquad 0<\delta\ll1,
\tag{C.31}
\]

and ask each normalized stage to produce only logarithmic gain
\(\log r=O(\delta)\).  If the selected analytic profile path has analytic
radius \(\rho_0>0\), an analytic Cauchy estimate of the form

\[
\rho(t)\ge \rho_0-\Lambda t
\tag{C.32}
\]

survives a normalized stage of duration

\[
T_{\rm st}=c_0\log r=O(\delta)
\tag{C.33}
\]

once \(\delta<\rho_0/(2c_0\Lambda)\).  This is exactly how analytic
regularity pays for the unstable factor
\(e^{\gamma|\eta|T_{\rm st}}\).  No Gevrey-2 Cauchy theory is being
asserted.

Let the uncut rescaled packet have Gaussian or periodic heat-kernel
tails.  Cut it only at

\[
R_j=L j\,\ell_j,\qquad \ell_j=r^{-j},
\tag{C.34}
\]

where \(L\) is a sufficiently large fixed constant.  For every fixed
physical derivative order \(m\), the seam then has size

\[
C_m\ell_j^{-m}e^{-cL^2j^2}
\le
C_m\exp\big(-c'L^2j^2+mj\log r\big)
=O_m(e^{-c''j^2}).
\tag{C.35}
\]

The same conclusion holds for fixed time derivatives.  With
\(a_j=\ell_j^{-\gamma}K_j^\gamma\), the physical stage time is the
normalized duration times the inertial unit
\[
\tau_j=\frac{\ell_j}{a_j}
=\ell_j^{1+\gamma}K_j^{-\gamma},
\]
not \(O(\delta\ell_j^2)\).  Every fixed power of \(\tau_j^{-1}\) still
costs only \(e^{O(j+\log j)}\), so (C.35) is unchanged.  If derivatives
up to the truncation order \(M_j\asymp j^2/\log j\) are used internally,
their factorial/Gevrey constants contribute \(e^{Cj^2}\); increasing
\(L\) can still leave a negative quadratic exponent.  The exact shifted
carrier ledger is recorded in
`2026-07-29-near-identity-cascade-ledger.md`.

The physical cutoff radii do eventually shrink:

\[
\frac{R_{j+1}}{R_j}
=r^{-1}\left(1+\frac1j\right)<1
\quad\text{for }j>\delta^{-1}+O(1).
\tag{C.36}
\]

The finitely many initial stages can be arranged separately.  Moreover,
the small gain does not destroy the usual geometric budgets.  For
example,

\[
\sum_j \delta\ell_j^2
=\frac{\delta}{1-r^{-2}}=O(1),
\qquad
\sum_j \delta\ell_j
=\frac{\delta}{1-r^{-1}}=O(1)
\tag{C.37}
\]

as \(\delta\downarrow0\).

Therefore the proposal

> short analytic transition \(+\) near-identity scale gain \(+\)
> cutoff at \(Lj\) Gaussian radii

**does evade the hydrostatic Hadamard/Gevrey-2 obstruction
conditionally**.  Cheverry's large-amplitude warning does not rule it
out in the forced approximate setting.

It is not yet a transition theorem, for three precise reasons.

1. One needs a \(K\)-uniform analytic construction for the full
   constrained system (C.5)--(C.10), including the material pressure
   inverse, not only analytic existence for the limiting hydrostatic
   system.
2. Hydrostatic pressure is horizontally nonlocal.  A Gaussian velocity
   does not automatically have Gaussian pressure tails: on
   \(\mathbb R^2\), the inverse horizontal Laplacian generally creates
   algebraic multipole tails.  Cutting those tails at \(Lj\) would not
   give (C.35).  The annular wake must cancel pressure multipoles through
   the required order, or an exact pressure-localized subclass must be
   found.
3. Every stage must land, with a uniformly nondegenerate rank-five
   endpoint map, in the same bounded analytic packet class needed to
   start the next stage.  Resetting the analytic radius between stages is
   legitimate only after that recurrence statement is proved.

These are separate endpoint/pressure gates.  They do not reinstate the
generic hydrostatic no-go.

### Route B: prove stability of the actual pulse path

Derive the linearized hydrostatic operator along the three sequential
single-parent pulses and prove a slow-frequency estimate weaker than
(C.26), strong enough for the desired Gevrey class.  It is not enough to
check only the initial sinusoidal parent: partner sidebands and the
generated charge-two and higher modes change the vertical shear.

A targeted spectral computation can falsify this route, but cannot prove
the required nonlinear estimate.

### Route C: find an invariant transparent subclass

Impose algebraic conditions for which

\[
\big(\partial_\vartheta^{-1}\operatorname{div}_xW\big)
\partial_\vartheta W
\tag{C.30}
\]

either vanishes or closes in a finite-dimensional stable family, while
retaining the rank-five zero-charge strain map.  The pure
slow-independent one-phase shear has this property but produces no
zero-charge child.  The missing object is a nontrivial modulated
rank-five family with the same closure.

This is the sharpest long-shot route because it would remove both
Cheverry's ill-posed modulation and the localization regularity mismatch
at once.

---

## 6. Revised theorem gate

The phase part of the proposed single-parent theorem may be replaced by
the following exact statement.

> **Material-phase resummation lemma.**  If a smooth lifted solution of
> (C.5), (C.6), (C.9), and (C.10) exists with
> \(|\nabla\phi^K|\ge c\), then all low--high \(O(K)\) transport is
> cancelled, all high--high \(O(K)\) transport reduces exactly to
> (C.8), Cheverry's geometrical phase hierarchy is the coefficient
> expansion of \(\phi^K\), and all adjusting phases can be absorbed by
> (C.16).

What is still unproved is:

> **Hydrostatic-profile gate.**  Construct the selected sequential
> rank-five trajectory for order-one normalized time with either
>
> 1. a slow analytic radius that stays positive, or
> 2. a verified spectral/structural stability condition yielding the
>    needed Gevrey estimate,
>
> and then prove a no-repeated-loss material pressure estimate and an
> exponentially accurate localization/wake coupling.

Until this gate is closed, the material phase is a correct algebraic
repair, not a complete transition theorem.

---

## Primary references

* Christophe Cheverry,
  [*Cascade of phases in turbulent flows*](https://arxiv.org/abs/math/0402408),
  especially formulas (2.3), (2.6)--(2.7), (3.2), (3.7)--(3.10), and
  (4.9), and the discussion of large-amplitude modulation.
* Daniel Han-Kwan and Toan T. Nguyen,
  [*Ill-posedness of the hydrostatic Euler and singular Vlasov
  equations*](https://arxiv.org/abs/1507.01813), especially the
  hydrostatic Euler system (1.1)--(1.3) and Theorem 1.1.
* Ruimeng Hu, Quyuan Lin, and Rongchang Liu,
  [*Regularization by noise for the inviscid primitive
  equations*](https://arxiv.org/abs/2407.21336), especially
  (1.1)--(1.2) for the periodic three-dimensional system and the
  deterministic regularity discussion in its introduction.
