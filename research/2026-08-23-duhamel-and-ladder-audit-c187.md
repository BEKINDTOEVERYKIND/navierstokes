# C187: an explicit finite-stage Duhamel constant and withdrawal of unlanded ladder numbers

**Date:** 2026-08-23

**Status:** exact \(H^3\) Fourier-algebra constant and explicit
linearized-Navier--Stokes Duhamel estimate; C136's abstract (5.1) folded into
the full renormalization trapping premise; two session-only spectral numbers
withdrawn

**Checker:**
[checks/explicit_h3_duhamel_c187.py](../checks/explicit_h3_duhamel_c187.py)

## 0. Audit verdict

The C136 conditional theorem displayed

\[
 \left\|\int {\cal U}_j\mathbb P\nabla\!\cdot(f\otimes g)\right\|_{{\cal Y}_j}
 \le K\|f\|_{{\cal Y}_j}\|g\|_{{\cal Y}_j}.             \tag{0.1}
\]

Neither \({\cal Y}_j\) nor its active/wake component spaces and moving
projections were defined.  Therefore no explicit uniform \(K\) can be
proved from the landed statement.  Equation (0.1) is now a clause of the
single active-plus-wake renormalization trapping hypothesis, not a proved or
separate named gate.

There is an ordinary positive estimate with completely explicit constants.
On the normalized three-torus, the forced linearized Navier--Stokes
Duhamel map about a background bounded by \(V\) in \(H^3\) satisfies

\[
 \boxed{
 \|z\|_{L^\infty(0,T;H^3)}
 \le44\sqrt{2T\over\nu}
       \exp\!\left({7744V^2T\over\nu}\right)
 \|f\|_{L^\infty H^3}\|g\|_{L^\infty H^3}.}            \tag{0.2}
\]

This proves finite-stage continuity for an actual viscous PDE operator.
Its constant diverges as \(\nu\downarrow0\), and therefore it cannot be
substituted for the scale-uniform backward-weighted estimate required by
the old C136 architecture.

The session-only infinite-ladder enclosures \(0.66855\ldots\) and
\(2.63707\ldots\) are also withdrawn.  No operator, normalization, rational
interval certificate, tail resolvent bound, or checker was landed for
them.  They carry zero evidentiary weight.

## 1. Normalization and an explicit \(H^3\) algebra constant

Use Fourier coefficients normalized so that Parseval has no volume factor,
and define

\[
 \|u\|_{H^3}^2=\sum_{k\in\mathbb Z^3}\langle k\rangle^6
                    |\widehat u(k)|^2,qquad
 \langle k\rangle=(1+|k|^2)^{1/2}.                       \tag{1.1}
\]

Set

\[
                         A_3^2=\sum_{k\in\mathbb Z^3}
                                      \langle k\rangle^{-6}. \tag{1.2}
\]

The max-norm shell \(\|k\|_\infty=r\) has exactly

\[
             (2r+1)^3-(2r-1)^3=24r^2+2\le26r^2         \tag{1.3}
\]

lattice points.  Since

\[
 \sum_{r=1}^{\infty}{1\over r^4}
 \le\sum_{r=1}^{5}{1\over r^4}+\int_5^\infty{x^{-4}}\,dx
 <{13\over12},                                          \tag{1.4}
\]

one obtains

\[
 A_3^2<1+26{13\over12}={175\over6}<{121\over4},
 \qquad A_3<{11\over2}.                                 \tag{1.5}
\]

For \(k=\ell+(k-\ell)\),

\[
 \langle k\rangle^3
 \le4\{\langle\ell\rangle^3+
              \langle k-\ell\rangle^3\}.               \tag{1.6}
\]

Fourier convolution, Young's inequality
\(\|a*b\|_{\ell^2}\le\|a\|_{\ell^2}\|b\|_{\ell^1}\),
and
\(\|\widehat u\|_{\ell^1}\le A_3\|u\|_{H^3}\) give

\[
 \boxed{\|f\otimes g\|_{H^3}
       \le8A_3\|f\|_{H^3}\|g\|_{H^3}
       <44\|f\|_{H^3}\|g\|_{H^3}.}                    \tag{1.7}
\]

The same constant applies to scalar, vector, and Frobenius tensor norms by
taking Euclidean magnitude inside the convolution estimate.

## 2. Explicit Duhamel estimate

Fix \(\nu>0\) and \(T>0\).  Let \(v\) be smooth and divergence free on
\([0,T]\times\mathbb T^3\),
with

\[
                         \sup_{0\le t\le T}\|v(t)\|_{H^3}\le V. \tag{2.1}
\]

For smooth \(f,g\), let \(z\) be the zero-initial-data solution of

\[
 \partial_tz-\nu\Delta z
 +\mathbb P\nabla\!\cdot(v\otimes z+z\otimes v)
 =-\mathbb P\nabla\!\cdot(f\otimes g),qquad z(0)=0.    \tag{2.2}
\]

Put

\[
 Z=\|z\|_{H^3},\quad D=\|\nabla z\|_{H^3},\quad
 F=\|f\|_{H^3},\quad G=\|g\|_{H^3}.                     \tag{2.3}
\]

The Leray projection is an orthogonal Fourier multiplier.  Pairing (2.2)
in \(H^3\), integrating the divergence by parts, and applying (1.7) to the
two background terms and the source gives

\[
 {1\over2}(Z^2)'+\nu D^2
       \le88VZD+44FGD.                                  \tag{2.4}
\]

Applying Young's inequality separately to the two terms, each with
\(\nu D^2/4\), and multiplying by two yields

\[
 (Z^2)'+\nu D^2
 \le {15488\over\nu}V^2Z^2+{3872\over\nu}F^2G^2.       \tag{2.5}
\]

Gronwall's inequality and \(Z(0)=0\) now give

\[
 Z(t)^2\le {3872T\over\nu}
  \exp\!\left({15488V^2T\over\nu}\right)
  \|f\|_{L^\infty H^3}^2\|g\|_{L^\infty H^3}^2.       \tag{2.6}
\]

Since \(3872=2\cdot44^2\), taking square roots proves (0.2).

All constants in (1.3)--(2.6) are explicit.  They are deliberately crude;
their role is to settle the audit, not to claim scale uniformity.

## 3. Why this does not prove the old (5.1)

The C127 normalized viscosity is \(\mu_j=\nu(j!)^{-2}\).  Replacing
\(\nu\) by \(\mu_j\) in (0.2) produces the explicit factor

\[
 44\sqrt{2T\over\mu_j}
 \exp\!\left({7744V^2T\over\mu_j}\right),               \tag{3.1}
\]

which is not stage-uniform.  Moreover, ordinary \(H^3\) does not contain
the active \(n^2\) backward-response weight, the retained-wake projection,
or the moving exit chart.  Thus (0.2) proves neither BAFL nor continuity of
the proposed renormalization map in its needed norm.

The architecture-correct obligation is now:

> specify the complete structured-state Banach space, prove a coercive
> stable-complement or validated finite-rank splitting for
> \(D{\cal R}(X_*)\), and give one explicit nonlinear trapping constant.

Until that space is fixed, a symbol \(K\) in the old (5.1) is a conditional
placeholder, not a theorem.

## 4. Infinite-ladder audit

An exhaustive current-tree and fetched-history search found no mathematical
occurrence of either withdrawn number.  Four raw CSV cells happen to begin
with \(0.66855\), but they are unnamed cascade/optimization columns and are
not spectral intervals.  They remain untouched.

The actual landed spectral record is:

1. C120 certifies only the finite weighted \(6\times6\) hexagon polynomial
   \(\lambda^6-9\lambda^4+18\lambda^2-9\) and
   \(633/250<\sigma_*<2533/1000\).
2. C148 defines an infinite variable-fiber nearest-neighbor operator and
   checks its algebra, but proves no spectral enclosure.
3. C168 defines the complete integer charge ladder and proves norm/tail
   identities, but no spectral enclosure.

Any future infinite-ladder enclosure must land, in one auditable chain, the
Hilbert space and domain, exact normalized operator, interval matrix
enclosure, residual, tail resolvent or Schur-complement bound, isolated
Riesz gap, and checker.  Retrofitting the withdrawn numbers to C148 or C168
would be fabrication.

## 5. Verification boundary

The checker verifies with exact rational arithmetic the shell count, the
zeta-tail bound in (1.4), the constants \(44,88,15488,3872,7744\), and the
Young/Gronwall coefficient ledger.  The Fourier-algebra and energy
arguments are written explicitly above and are not replaced by finite
tests.

The checker does not define the old \({\cal Y}_j\), prove a stage-uniform
trapping estimate, validate an infinite-ladder spectrum, or prove UVSR or a
Navier--Stokes singularity.
