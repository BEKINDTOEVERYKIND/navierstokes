# Audit of arXiv:2604.09949: the displayed lifted system is not 3D Navier--Stokes

Date: 2026-07-29

## Verdict

The preprint

> R. Shahmurov, *Stable Finite-Time Singularity Formation for 3D
> Navier--Stokes via 5D-Lifted Axisymmetric Reductions*,
> [arXiv:2604.09949](https://arxiv.org/abs/2604.09949)

in its current v1 version dated 10 April 2026 does **not** establish
blow-up for the three-dimensional Navier--Stokes
equations.  Its central scalar profile equation is built from a recovery
law that is not the axisymmetric Navier--Stokes recovery law.  Two exact
kinematic discrepancies are already fatal, before any numerical constants
are considered:

1. the meridional velocity is made divergence-free with respect to the
   five-dimensional weight \(r^3\,dr\,dz\), whereas a physical
   three-dimensional axisymmetric velocity is divergence-free with respect
   to \(r\,dr\,dz\);
2. the independently evolving azimuthal vorticity
   \(\omega^\theta/r\), which determines the meridional velocity, is
   omitted and replaced by an elliptic equation whose right-hand side is
   only the swirl forcing \(\partial_z(u^\theta/r)^2\).

Consequently, even an exact zero of the paper's operator \(G\) would solve
a different reduced system.  The Newton--Kantorovich and torus-transfer
claims cannot repair an incorrect governing equation.  This is not a
matter of sign convention: substituting the paper's displayed velocity
formula into the physical divergence gives a nonzero scalar multiple of
its displayed radial velocity.

There are also independent reproducibility and functional-analytic gaps:
the advertised coefficient/certificate/code package is absent; the
displayed tail coercivity starts at mode \(1200\) while the finite inverse
ends at mode \(450\); the stated bilinear estimate truncates both inputs at
mode \(450\); and a small periodic-image error is treated as though it
preserved an exact time-dependent self-similar solution.

No conclusion about the actual Millennium problem follows from this
preprint.

## 1. The exact axisymmetric variables

For a physical axisymmetric velocity

\[
 u=u^r e_r+u^\theta e_\theta+u^z e_z ,
\]

three-dimensional incompressibility is

\[
 \boxed{\partial_r u^r+\frac1r u^r+\partial_z u^z=0.}
\tag{1.1}
\]

Introduce

\[
 U=\frac{u^\theta}{r},
 \qquad
 W=\frac{\omega^\theta}{r},
 \qquad
 B=\partial_{rr}+\frac3r\partial_r+\partial_{zz}.
\tag{1.2}
\]

The standard axisymmetric Navier--Stokes system contains the two distinct
evolution equations

\[
\begin{aligned}
(\partial_t+u^r\partial_r+u^z\partial_z-\nu B)U
+2\frac{u^r}{r}U&=0,\\
(\partial_t+u^r\partial_r+u^z\partial_z-\nu B)W
-\partial_z(U^2)&=0,
\end{aligned}
\tag{1.3}
\]

up to the harmless sign convention for \(\omega^\theta\).  The meridional
velocity \((u^r,u^z)\) is recovered by the axisymmetric Biot--Savart law
from \(W\), not directly from \(\partial_z(U^2)\).

For example, with the standard reduced streamfunction \(\psi_1\), one may
write

\[
 u^r=-r\partial_z\psi_1,\qquad
 u^z=2\psi_1+r\partial_r\psi_1,\qquad
 W=-B\psi_1,
\tag{1.4}
\]

where the sign in the last identity changes if the opposite convention
for \(\omega^\theta\) is used.  The coefficient \(2\) in \(u^z\) is forced
by three-dimensional incompressibility.

These equations are displayed, for example, in Chen--Fang--Zhang,
[*Regularity of 3D axisymmetric Navier--Stokes
equations*](https://arxiv.org/abs/1505.00905), equations (1.2), (1.4),
and (1.5).  Their notation
\(\Theta=ru^\theta\), \(\Gamma=\omega^\theta/r\) is equivalent to
(1.1)--(1.3).

## 2. First fatal mismatch: five-dimensional versus physical divergence

The audited preprint defines \(b=(u^r,u^z)\) as the meridional velocity but
imposes, in its equation (5),

\[
 \partial_r u^r+\frac3r u^r+\partial_z u^z=0.
\tag{2.1}
\]

Its streamfunction formula (9) is explicitly chosen to satisfy (2.1).
Comparing (2.1) with the physical identity (1.1) gives

\[
 \operatorname{div}_{\mathbb R^3}u
 =-\frac{2u^r}{r}.
\tag{2.2}
\]

Thus the reconstructed three-dimensional velocity is not divergence-free
unless \(u^r\equiv0\).  The coefficient \(3/r\) is correct for the lifted
diffusion operator acting on the scalar \(U=u^\theta/r\); it is not the
physical incompressibility coefficient for the meridional velocity.
Treating a useful five-dimensional scalar lift as though the velocity
itself lived in five dimensions changes the equation.

The contradiction can also be read directly from the paper's
streamfunction formula.  Set

\[
 \phi=\frac{\bar\psi}{\rho^4}.
\]

Its equation (9) then gives

\[
 \bar u^\rho=-\rho\phi_\zeta,\qquad
 \bar u^\zeta=4\phi+\rho\phi_\rho.
\tag{2.3}
\]

The physical streamfunction formula (1.4), with the same radial
component, would instead have

\[
 \bar u^\zeta=2\phi+\rho\phi_\rho.
\]

Indeed, (2.3) gives exactly

\[
 \partial_\rho\bar u^\rho+\frac1\rho\bar u^\rho
 +\partial_\zeta\bar u^\zeta
 =2\phi_\zeta=-2\frac{\bar u^\rho}{\rho}.
\tag{2.4}
\]

Nor can \(u^\rho\equiv0\) rescue the claimed nontrivial profile.  The
paper imposes that \(\bar\Omega\) is odd in \(\zeta\).  If
\(\partial_\zeta(\bar\Omega^2)\equiv0\), oddness forces
\(\bar\Omega\equiv0\).  Conversely, if \(u^\rho\equiv0\), then
\(\phi_\zeta\equiv0\), so the left side of the paper's elliptic recovery
equation is independent of \(\zeta\); it cannot equal the nonzero odd
source \(\partial_\zeta(\bar\Omega^2)\).  Thus the intended nonzero
profile necessarily has a nonzero physical divergence.

There is a useful way to locate the exact conceptual slip.  For a
*physical* divergence-free meridional velocity,

\[
 \operatorname{div}_{5}b
 :=\partial_r u^r+\frac3r u^r+\partial_z u^z
 =2\frac{u^r}{r}.
\tag{2.5}
\]

Consequently the correct lifted swirl equation can be written in
five-dimensional conservative form

\[
 \partial_t U+\operatorname{div}_{5}(bU)=\nu BU.
\tag{2.6}
\]

The preprint instead sets \(\operatorname{div}_5 b=0\) while retaining the
separate \(2(u^r/r)U\) stretching term.  It has promoted the scalar lift
to a false kinematic condition on the physical velocity.

This is a purely symbolic contradiction between displayed equations.  It
does not depend on interval arithmetic, a norm choice, or a numerical
tolerance.

## 3. Second fatal mismatch: the meridional-vorticity equation is missing

The preprint's equation (8) asserts an elliptic recovery of its
streamfunction from

\[
 B\!\left(\frac{\bar\psi}{\rho^4}\right)
 =\partial_\zeta(\bar U^2).
\tag{3.1}
\]

But \(\partial_z(U^2)\) is the **forcing term** in the evolution equation
for \(W=\omega^\theta/r\); it is not \(W\) itself and is not an
instantaneous Biot--Savart source.

Under backward self-similar scaling, the omitted equation does not
disappear.  It becomes schematically

\[
\frac32\bar W+\frac12 y\cdot\nabla\bar W
 +\bar b\cdot\nabla\bar W-\nu B\bar W
 =\partial_\zeta(\bar U^2),
\tag{3.2}
\]

while the separate elliptic identity
\(\bar W=-B\bar\psi_1\) recovers \(\bar b\) from \(\bar W\) through
(1.4).  The coefficient \(3/2\) follows from the scaling
\(W(r,z,t)=(T-t)^{-3/2}\bar W(r/\sqrt{T-t},z/\sqrt{T-t})\).
Equation (3.1) drops the scaling drift, transport, and viscous terms in
(3.2), then identifies the remaining forcing with the velocity recovery
source.  The resulting scalar fixed-point operator \(G(U,\nu)\) is
therefore not equivalent to the axisymmetric Navier--Stokes system.

This also explains why the paper appears to close a three-component
velocity problem using only one scalar \(U\): it has removed one of the two
active scalar degrees of freedom.

## 4. Conflict with backward self-similar rigidity

The preprint claims an exactly backward self-similar velocity

\[
 u(x,t)=(T-t)^{-1/2}\bar u(x/\sqrt{T-t}).
\tag{4.1}
\]

The manuscript later calls this core Gaussian-localized.  If its claimed
physical profile had that decay, then \(\bar u\in L^3(\mathbb R^3)\)
(indeed it would lie in much stronger decay classes).  Nontrivial Leray
backward self-similar Navier--Stokes profiles in \(L^3\) were ruled out by
Nečas--Růžička--Šverák,
[*On Leray's self-similar solutions of the Navier--Stokes
equations*](https://doi.org/10.1007/BF02551584).  Tsai,
[*On Leray's self-similar solutions of the Navier--Stokes equations
satisfying local energy estimates*](https://doi.org/10.1007/s002050050099),
proved extensions under local-energy and \(L^p\) hypotheses.

This classical theorem is not merely an external reason for skepticism.
It is consistent with the direct equation audit above: the paper evades
the rigidity result only because its recovered profile is not a
three-dimensional incompressible Navier--Stokes profile.

This rigidity check is deliberately secondary.  The manuscript does not
actually construct the advertised map from its abstract compactified
\(\mathcal X\) norm to a physical Gaussian-decaying velocity, so the
decay hypothesis itself is not certified.  The divergence contradiction
in Section 2 requires no such functional-analytic assumption.

## 5. The numerical certificate is not a reproducible proof

Even for the altered scalar operator, the manuscript does not supply the
objects needed to verify its claimed computer-assisted proof.

* Appendix A prints only five selected coefficients out of the claimed
  \(450\), and says the full data are “intended to be provided.”
* Appendix E says the exact source files are not reproduced because
  previously uploaded files expired.
* Appendix F says a complete release **should** provide the coefficient
  export, Jacobian generator, and quadrature/projection routines.
* The displayed “audit log” is text containing asserted constants, not an
  independently checkable interval certificate.

There are structural gaps in the written estimates as well.

1. The finite Jacobian inverse is for modes \(1,\dots,450\), while the tail
   coercivity constant is defined using an infimum over \(j\ge1200\).
   Modes \(451,\dots,1199\) are not covered by either displayed argument.
   In particular, an infimum over \(j\ge1200\) cannot imply the paper's
   stated estimate for \(Q_{>450}\).
2. Proposition 10.2 defines the quadratic interaction with
   \(k,l\le450\), but then states a bilinear bound for general inputs in
   the infinite-dimensional spaces.  High--low and high--high input
   interactions are absent from the displayed left-hand side.
3. A diagonal tail lower bound and a finite inverse do not by themselves
   invert the full block operator; the finite--tail off-diagonal Schur
   couplings must be bounded.
4. The linearized nonlinear term couples spectral modes.  The manuscript
   passes from a scalar quantity \(\mathrm{Inter}_j\) to a mode-by-mode
   lower bound proportional to \(\|h_j\|\), without a displayed
   row-sum/operator-norm estimate controlling the contributions of all
   other \(h_k\).  No coefficient data are supplied from which that
   missing estimate could be checked.

Any one of these would prevent the written Newton--Kantorovich argument
from being referee-grade.  They are secondary, however, because Sections
2--3 already show that the operator being certified is the wrong one.

## 6. Why the torus transfer does not restore a solution

The paper periodizes a localized profile, applies the Leray projector, and
estimates the interaction with periodic images by a tiny nonzero number.
This gives an approximate profile, not an exact singular solution.

There is an even more elementary inconsistency in its description of the
projection step.  Divergence is linear.  If the whole-space core were
divergence-free, then its Poincaré sum

\[
 \sum_{n\in\mathbb Z^3}\bar u(x+2\pi n)
\]

would already be exactly divergence-free term by term; image
*interactions* cannot generate a divergence defect.  For the actual core,
Section 2 gives a nonzero central divergence defect.  Summing distant
Gaussian images does not turn that central defect into an exponentially
small overlap effect.  Thus the asserted \(10^{-20}\) Leray correction is
not supported by the displayed construction.

A Newton correction could in principle turn a small profile defect into an
exact zero only for a correctly defined, fixed operator on a correctly
defined function space.  Here:

* the underlying operator already fails the physical equations;
* a physical periodized self-similar path would have the form
  \[
  u_{\rm per}(x,t)=(T-t)^{-1/2}
  \sum_{n\in\mathbb Z^3}
  \bar u\!\left(\frac{x+2\pi n}{\sqrt{T-t}}\right),
  \]
  so in similarity variables the image spacing is
  \(2\pi/\sqrt{T-t}\); it is not a fixed stationary torus profile;
* Leray projection enforces divergence-free data at one time but does not
  make the periodized self-similar path solve Navier--Stokes at all later
  times; and
* “exponentially small” is not “zero,” especially for an unforced exact
  theorem.

The torus step therefore cannot support Theorem 2.1.

## Bottom line

The smallest indisputable fatal error is already the first pair of
identities:

\[
\boxed{\partial_r u^r+\frac1r u^r+\partial_z u^z=0}
\qquad\hbox{versus the paper's}\qquad
\boxed{\partial_r u^r+\frac3r u^r+\partial_z u^z=0}.
\]

For the paper's own streamfunction formula their difference is
\(-2u^r/r\), and the paper's odd nontrivial source forces
\(u^r\not\equiv0\).  Independently, the exact physical system also
requires

\[
\boxed{
(\partial_t+b\cdot\nabla-\nu B)
\frac{\omega^\theta}{r}
=\partial_z\left(\frac{u^\theta}{r}\right)^2 }.
\]

The preprint replaces the physical divergence by the lifted one and omits
this second evolution equation.  Its fixed point is consequently not a
three-dimensional Navier--Stokes solution.  This claimed shortcut is
closed; it supplies no usable core for the present cascade program.
