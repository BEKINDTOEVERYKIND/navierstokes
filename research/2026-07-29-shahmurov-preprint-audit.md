# Audit of Shahmurov's interior-quadrupole Euler preprint

Date: 2026-07-29

## Scope and verdict

This note audits the load-bearing frozen hyperbolic calculation in R. Shahmurov,
*Euler Singularities II: Interior Quadrupole Blow-Up for Smooth Axisymmetric
Euler with Swirl in \(\mathbb R^3\)*,
[arXiv:2605.04526v1](https://arxiv.org/pdf/2605.04526).

The claimed blow-up proof has a leading-order sign inconsistency. With the
paper's definitions and its exact hyperbolic model,

\[
  \lambda'=-\sigma\lambda,\qquad b'=\sigma b,\qquad
  C=\lambda^2b,
\]

one obtains

\[
  C'=-\sigma C,
\]

not \(C'\gtrsim Q C>0\). This directly contradicts equation (37) in Lemma
5.2 and consequently breaks the Riccati comparison (66). Independently, the
isotropically shrinking cutoff does not follow the anisotropic hyperbolic
flow, so the moving-projection commutator asserted in equation (35) does not
vanish in the exact model.

This establishes that the stated preprint does **not** prove its claimed
smooth whole-space Euler blow-up theorem. It does not prove that the theorem
itself is false.

## 1. Definitions used by the paper

In local coordinates \(x=r-r_*(t)\), \(y=z\), equations (8), (15)--(17), and
(42) use

\[
 (U,V)=(\sigma x,-\sigma y),\qquad
 \lambda'=-\sigma\lambda,
\]

\[
 \Gamma=\Gamma_*+\frac12 bxy^2+R_\Gamma,\qquad
 C=\lambda^2b.
\]

The projected definitions (31)--(33) are

\[
 b_\lambda
 =\frac{2\langle\Gamma-\Gamma_*,xy^2\rangle_\lambda}
        {\langle xy^2,xy^2\rangle_\lambda},
 \qquad C_\lambda=\lambda^2b_\lambda.
\]

For the exact profile \(R_\Gamma=0\), this projection gives
\(b_\lambda=b\) for every value of \(\lambda\). Thus the issue below cannot
be removed by distinguishing the Taylor coefficient from the projected
coefficient.

## 2. Exact transport calculation

The circulation equation is equation (3):

\[
 D_t\Gamma=0.
\]

Insert the exact hyperbolic field and the active jet:

\[
\begin{aligned}
0
&=\partial_t\!\left(\frac12bxy^2\right)
  +\sigma x\,\partial_x\!\left(\frac12bxy^2\right)
  -\sigma y\,\partial_y\!\left(\frac12bxy^2\right)\\
&=\frac12(b'-\sigma b)xy^2.
\end{aligned}
\]

Therefore

\[
 b'=\sigma b.
\]

This agrees with the paper's Lemma 6.3, equation (46), and with its detailed
jet calculation in Section 15. But differentiating the paper's normalized
amplitude gives

\[
\begin{aligned}
C'
  &=2\lambda\lambda'b+\lambda^2b'\\
  &=2\lambda(-\sigma\lambda)b+\lambda^2(\sigma b)\\
  &=-\sigma\lambda^2b=-\sigma C.
\end{aligned}
\]

More generally, the paper itself repeatedly uses
\(b'=\sigma b+O(E\sigma b)\). Hence in its advertised small-error regime,

\[
 C'=(-1+O(E))\sigma C,
\]

which is strictly negative when \(\sigma,C>0\) and \(E\) is sufficiently
small.

By contrast, Lemma 5.2 asserts in equation (37)

\[
 D^+C_\lambda
 \ge cQ_\lambda C_\lambda-CE\,Q_\lambda C_\lambda,
\]

and the strain estimate used there identifies
\(\sigma\gtrsim Q_\lambda>0\). In the exact regime \(E=0\), the right-hand
side is positive while the exact left-hand side is
\(-\sigma C_\lambda<0\). This is an opposite-sign, leading-order
contradiction, not a missing constant or lower-order error.

The proof of equation (37) tracks only the growth of \(b_\lambda\); it omits
the contribution

\[
 2\lambda\lambda'b_\lambda=-2\sigma C_\lambda
\]

from differentiating \(C_\lambda=\lambda^2b_\lambda\).

Consequently Corollary 7.2, which propagates
\(C\ge\kappa Q^2\) by using \(C'\ge cQC\), and the comparison system (66)
do not follow.

## 3. A second manifestation: the moving cutoff is not material

The same problem appears in Lemma 5.1. Let

\[
 w_\lambda(x,y)=w(x/\lambda,y/\lambda)
\]

be the paper's isotropically scaled cutoff. At fixed \((x,y)\), the choice
\(\lambda'=-\sigma\lambda\) gives

\[
 \partial_t w_\lambda
 =\sigma(x\partial_x+y\partial_y)w_\lambda.
\]

For the frozen hyperbolic velocity \(v=(\sigma x,-\sigma y)\),

\[
 v\cdot\nabla w_\lambda
 =\sigma(x\partial_x-y\partial_y)w_\lambda.
\]

Therefore

\[
 (\partial_t+v\cdot\nabla)w_\lambda
 =2\sigma x\partial_xw_\lambda,
\]

which is generally nonzero on the cutoff collar.

Equation (35) bounds the commutator by
\[
 C\sigma(\varepsilon_{\rm strain}+\rho)
 \int_{P_{2\lambda}}|F||\Phi|.
\]
In the exact flat model
\(\varepsilon_{\rm strain}=\rho=0\), this bound sets the commutator to zero,
but the displayed cutoff flux remains. The statement in the proof that
\(\lambda'=-\sigma\lambda\) cancels the leading dilation of the cutoff is
false for a scalar, isotropically scaled packet: the flow expands in \(x\)
while compressing in \(y\).

There is a parallel omission in equation (36). For \(G=axy\), the projected
coefficient is exactly \(a_\lambda=a\) and
\[
 Q_\lambda=c_Qa\lambda^2.
\]
The monomial \(xy\) has zero transport growth in
\((\sigma x,-\sigma y)\), so if its source coefficient is \(s\), then
\[
 Q_\lambda'=c_Qs\lambda^2-2\sigma Q_\lambda.
\]
Equation (36) retains the source term but not the leading
\(-2\sigma Q_\lambda\) term.

## 4. Geometric meaning

The circulation jet coefficient \(b\) grows because the \(y\)-direction is
compressed more strongly than the \(x\)-direction is expanded for the
monomial \(xy^2\). But incompressibility gives simultaneous expansion in
\(x\). An isotropically shrinking square discards that expanding material.
The discarded flux is precisely what appears as the missing cutoff term and
as the negative derivative of \(\lambda^2\).

A material hyperbolic packet would instead use anisotropic lengths
\[
 \lambda_x'=\sigma\lambda_x,\qquad
 \lambda_y'=-\sigma\lambda_y.
\]
Its area is constant. Such a formulation preserves the useful jet-growth
identity, but it does not yield the paper's isotropic concentration or its
Riccati system. Any repaired mechanism would have to recapture the expanding
direction through a genuinely new folding, reset, or multi-stage return
mechanism.

## 5. What remains potentially reusable

The following local algebra survives this audit:

1. The exact source compatibility
   \[
   r^{-4}\partial_y(\Gamma^2)
   =2r_*^{-4}\Gamma_*bxy+\text{lower-order terms}
   \]
   for \(\Gamma=\Gamma_*+\frac12bxy^2+\cdots\).
2. Hyperbolic transport amplifies the **differential coefficient**
   \(b\) by \(b'=\sigma b\).
3. A four-quadrant \(xy\) vorticity pattern is a plausible sign design for
   an interior strain kernel, although the paper's global kernel and
   bootstrap estimates were not independently verified here.

What does not survive is the closed, isotropically concentrating feedback
loop. For our cascade program this reinforces the need for an anisotropic
Kelvin amplifier plus a separate nonlinear return/wake mechanism; a single
shrinking square cannot perform both jobs.

## 6. Relation to the author's Navier--Stokes claims

The Euler preprint explicitly says its theorem is logically independent of
the author's Navier--Stokes manuscripts, so the calculation above is not by
itself a refutation of those manuscripts.

However, the claimed full three-dimensional proof
[arXiv:2605.01873v2](https://arxiv.org/pdf/2605.01873) is not independent
confirmation: its Theorem 3.3 explicitly takes the axisymmetric-with-swirl
global theorem from
[arXiv:2605.01875v2](https://arxiv.org/pdf/2605.01875) as an endpoint input.
The later 99-page manuscript
[arXiv:2606.07869v1](https://arxiv.org/pdf/2606.07869) again claims that
axisymmetric-with-swirl theorem by a different expanded framework.

A limited check of the later manuscript finds an unresolved uniformity step
in the terminal argument. Lemma 34.6 uses absolute continuity of an
individual \(L^{14/5}\) norm to choose a smaller cylinder, whereas Theorem
34.7 requires a **universal** scale \(\theta\) for the normalized endpoint
class. Boundedness in \(L^{14/5}\) alone does not give uniform absolute
continuity over a concentrating family, and the adjacent compactness Lemma
34.4 states strong compactness only in \(L^2_{\rm loc}\), not in
\(L^{14/5}\). A separate uniform-integrability or stronger compactness
argument would be needed. This is an audit flag rather than a complete
line-by-line refutation of that much longer manuscript.

Accordingly, none of these preprints should presently be treated as a
verified Navier--Stokes breakthrough. The clean conclusion of this note is
narrower and fully explicit: the Euler proof's load-bearing comparison
system contradicts its own exact model.

## Primary sources

- R. Shahmurov, *Euler Singularities II: Interior Quadrupole Blow-Up for
  Smooth Axisymmetric Euler with Swirl in \(\mathbb R^3\)*,
  [arXiv:2605.04526v1](https://arxiv.org/pdf/2605.04526).
- R. Shahmurov, *Large-Data Global Regularity for Three-Dimensional
  Navier--Stokes I*,
  [arXiv:2605.01875v2](https://arxiv.org/pdf/2605.01875).
- R. Shahmurov, *Large-Data Global Regularity for Three-Dimensional
  Navier--Stokes II*,
  [arXiv:2605.01873v2](https://arxiv.org/pdf/2605.01873).
- R. Shahmurov, *Global Regularity for Axisymmetric Navier--Stokes Flows
  with Swirl*,
  [arXiv:2606.07869v1](https://arxiv.org/pdf/2606.07869).
