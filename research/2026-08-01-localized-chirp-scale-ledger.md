# Localized-chirp scale ledger after the affine-capture lemma

## Status

This note records an exact envelope obstruction, an asymptotic tube repair, and
the remaining phase-localization route for the Palasek--Gavrilov transition.
It is not a Navier--Stokes singularity proof.  Its purpose is to prevent the
finite-dimensional capture lemma from being mistaken for a localized PDE
module.

## 1. Stage and carrier scales

At parent frequency `N`, use the intermittent Gavrilov scaling

\[
 r_N=N^{-1},\qquad R_N=N^{4-2\alpha},\qquad
 |\Omega_N|\asymp N^{-2(\alpha-1)},
\]

and the stage amplitudes

\[
 X_N=\|u_N\|_2=N^{\beta-\alpha},\qquad
 Y_N=\|u_N\|_\infty=N^{\beta-1},\qquad
 g_N\asymp N Y_N=N^\beta .
\]

Let the child frequency be `Q=N^b`.  Write

\[
 \epsilon=b-1,\qquad \Delta=\alpha-\beta,
 \qquad \mu=\beta-2b.
\]

The basic region is

\[
 \epsilon,\Delta,\mu>0,\qquad
 2\epsilon+\mu+\Delta<\frac12,                    \tag{1}
\]

equivalent to `1<b`, `2b<beta<alpha<5/2`.  The child/parent
`L^2` ratio and carrier frequency are

\[
 \rho_N=\frac{X_Q}{X_N}=N^{-\epsilon\Delta},
 \qquad \Lambda\asymp N^{\beta/2},
 \qquad \frac Q\Lambda=N^{-\mu/2}.                \tag{2}
\]

For normalized spatial modes, the parent--carrier coupling is
`d~N^alpha`.  If the carrier occupies volume `V_a` and overlaps a child
of volume `V_c=Q^{-2(alpha-1)}`, the child--carrier coupling is

\[
 c\asymp Q^\alpha f,\qquad f=V_c/V_a.              \tag{3}
\]

The affine capture lemma needs `c/d~rho_N^{-1}`.  Hence

\[
 f=N^{-\beta(b-1)},\qquad
 V_a=N^{-A},\qquad
 A=2b(\alpha-1)-\beta(b-1).                        \tag{4}
\]

This is the origin of the small-overlap localization problem; it is not an
optional numerical tuning.

## 2. Exact scalar-envelope obstruction

For a one-dimensional chirp put

\[
 \eta=t\cdot x+\frac{\psi(s)}{2\Lambda},\qquad
 \xi=t\cdot x-\frac{\psi(s)}{2\Lambda},
 \qquad s=e_r\cdot x,
\]

and `a=w-delta t`, `b=w+delta t`, with
`w=e_r+H e_h` and `delta=psi'/(2 Lambda)`.  Scalar envelopes make

\[
 u_+=A a e^{i\Lambda\eta},\qquad
 u_-=B b e^{i\Lambda\xi}.
\]

Exact divergence forces

\[
 A=F(\eta,\zeta),\qquad B=G(\xi,\zeta),
 \qquad \zeta=(e_h-H e_r)\cdot x.                  \tag{5}
\]

The exact high-envelope darkness equation then reduces to

\[
 G F_\eta=F G_\xi .                               \tag{6}
\]

Where `delta` is nonzero, `(s,t dot x)` and `(eta,xi)` are local
coordinates.  On every connected nonzero overlap, (6) implies

\[
 F=C_A(\zeta)e^{c(\zeta)\eta},\qquad
 G=C_B(\zeta)e^{c(\zeta)\xi}.                     \tag{7}
\]

Thus no nonzero smooth scalar envelopes are at once compact, exactly
divergence-free, and exactly high-dark.  Real periodic envelopes force
`c=0`, leaving unlocalized slabs.  This is an exact obstruction, not a
power-counting concern.

## 3. Characteristic tubes and their geometric limit

The own-characteristic choice (5) keeps each carrier exactly divergence-free
and kills its self-envelope derivative.  The cross-envelope defect is

\[
 2\delta[-F G_\xi b+G F_\eta a]e^{i2\Lambda t\cdot x},       \tag{8}
\]

so a taper of width `L` costs `1/(Lambda L)` relative to the desired
low beat.  Flattening `psi` while `F=G` and tapering only after
`delta'=0` also confines the zero-frequency curvature term to the wake.

There is, however, a geometric price for obtaining the overlap (4).  The two
characteristic directions meet at angle

\[
 \sin\theta\asymp Q/\Lambda .                     \tag{9}
\]

For tubes of radius `R_a` and length `L_a`,

\[
 V_a\asymp R_a^2L_a,\qquad
 f\asymp\min\left(1,\frac{R_a}{L_a\sin\theta}\right).
\]

Prescribing the small value in (4) therefore gives

\[
 R_a\asymp f L_a Q/\Lambda,
 \qquad \frac1{\Lambda R_a}\asymp\frac1{fQ L_a}.  \tag{10}
\]

Even with the maximal parent length `L_a~R_N`, a small tube wake requires

\[
 b+4-2\alpha-\beta(b-1)>0.                         \tag{11}
\]

If a complete next Gavrilov torus must fit inside the parent minor radius,
then

\[
 b(2\alpha-4)>1.                                   \tag{12}
\]

Conditions (11)--(12) are incompatible with `beta>2b`.  Indeed, putting
`p=2alpha-4`, (11) has left margin strictly below

\[
 b-p-2b(b-1)<3b-2b^2-1/b<0\qquad(b>1),            \tag{13}
\]

where (12) gives `p>1/b`.  Hence two nearly parallel characteristic tubes
cannot simultaneously realize the required overlap, have a perturbative
taper, and geometrically contain the next self-similar torus.

This does not rule out non-tubular or phase-localized carriers.

## 4. The surviving phase-localization route

Keep both carriers on the larger volume `V_a`, but make their relative phase
constant outside the child core.  The low interaction is proportional to a
derivative of the relative phase and is therefore supported only where that
phase varies.  This separates carrier energy support from child-writing
support without using the forbidden compact exact-dark scalar envelope.

On a one-dimensional core, the clean choice is the small quadratic chirp

\[
 \psi(s)=a s^2,\qquad \chi=|a|R_c^2\ll1.           \tag{14}
\]

It writes

\[
 8Hca\,s\cos(as^2)=\gamma s[1+O(\chi^2)]          \tag{15}
\]

and has three advantages.

* `psi'''=0`, so the omitted third-derivative self residual vanishes on the
  exact core.
* The real-field magnitude-imbalance term is a spatially constant force,
  hence a removable Galilean acceleration with zero strain.
* Choosing `chi=G_N^{-1}` costs only a polylogarithmic carrier factor and gives
  accumulated affine error `G_N chi^2=G_N^{-1}`.

The chirped high and self harmonics are smaller than the written child by

\[
 (Q/\Lambda)^2=N^{-\mu}.                           \tag{16}
\]

An unmatched nonconstant zero-mode term would instead be `O(Q/Lambda)`.
The exact affine-capture estimate bounds its integrated coefficient by less
than `0.48 Q/Lambda` for damping ratio `0.9` and `r>=10`; (14) removes its
strain on the core.

For a generalized phase localized in all slow coordinates, exact
transversality can still be imposed, but divergence introduces WKB
polarization corrections.  The missing lemma is now precise:

> Construct two solenoidal chirped carriers on `V_a` whose relative phase is
> constant off a set of volume `V_c`, whose low beat agrees with an affine (or
> quadratic) child jet on that set, and whose full Euler--viscous residual is
> `o(rho_N X_N)` through a gain action `G_N`.

Unlike the finite tube proposal, this lemma is not contradicted by (13): the
two carriers may have the same support, and localization is carried by their
relative phase rather than by a small geometric intersection.

## 5. Power window not yet refuted

For the bare chirp, curvature, and viscosity errors, a conservative sufficient
ledger is

\[
 5-2\alpha>\epsilon\Delta,
 \qquad \mu>2\epsilon\Delta .                      \tag{17}
\]

The first controls thin-torus curvature relative to the child ratio; the
second even controls an uncancelled `Q/Lambda` term.  It is compatible with
(1).  For example

\[
 b=1.05,\qquad \beta=2.20,\qquad \alpha=2.35       \tag{18}
\]

gives

\[
 \rho_N=N^{-0.0075},\quad
 (Q/\Lambda)^2=N^{-0.10},\quad
 Q/\Lambda=N^{-0.05},\quad
 r_N/R_N=N^{-0.30}.
\]

This example is only a residual window; it does not satisfy the nested-torus
condition (12).  A nested geometry, phase-localization construction, and the
transported BAS estimate must be selected together rather than by this bare
ledger.

## 6. Remaining theorem gap

The finite-dimensional component is now quantitative: the affine core has a
captured orbit with a fixed terminal spectral gap.  The decisive missing PDE
result is no longer “find a favorable triad.”  It is a localized, transported
geometric-optics lemma that simultaneously proves:

1. phase localization with the support ratio (4);
2. invariance up to summable error under the actual curved Gavrilov parent;
3. robustness of viscous capture under those errors;
4. conversion of the written jet into the full next Gavrilov block; and
5. summability through the infinite cascade.

No GPU search resolves these analytic requirements.
