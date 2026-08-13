# C162: an exact obstruction for the tuned circular pure-normal gate block

**Date:** 2026-08-05

**Status:** exact one-edge forward bound and exact counterexample to the
specified tuned circular one-polarization block; broader gate bundles remain
open

**Checker:**
[checks/circular_normal_gate_star_no_go_c162.py](../checks/circular_normal_gate_star_no_go_c162.py)

## 0. Claim boundary

C161 requires a normalized, energy-preserving source-to-bright-daughter
rotation.  This note tests the most direct physical candidate on the same
C149 limiting elliptic ring: one pure-normal gate frequency and one circular
polarization.

There are two exact conclusions.

1. The forward Leray edge never vanishes.  For every source angle and every
   positive normal height its norm is at least half the source radius; at the
   tuned height used below its norm is exactly the source radius divided by
   \(\sqrt2\).  The same formula extends algebraically to zero height, but a
   complex circular coefficient at the zero Fourier mode would violate the
   reality condition and is not treated as a physical gate.
2. At the phase-zero source and that tuned height, the normalized forward
   and physical reverse coefficients have product

   \[
   -{1\over4}+{i\over4}.
   \]

   This nonreal product is invariant under source/daughter phase choices and
   reciprocal positive coordinate rescalings.  It therefore cannot be a
   two-mode skew rotation, whose off-diagonal product is negative real.

The second statement is deliberately narrow.  It rules out the specified
autonomous tuned circular block as an exact realization of the C161 star.
The real gate also sends the source toward \(p-g\), but an exact star in the
natural Fourier daughter coordinates must satisfy the skew reverse relation
on every individual edge; failure on the \(p+g\) edge is therefore already a
necessary-condition obstruction.  A non-diagonal mixing of daughter or
polarization coordinates is a larger block and is not classified here.
It does **not** prove a no-go for other gate heights, other single
polarizations, paired polarizations, time-modulated phases, gate-depleting
trajectories, or full multi-frequency bundles.  Nor does the one-edge lower
bound prove a uniform multi-source star realization.

## 1. Symbol and selected source line

Use the C155 orthonormal coordinates, put \(r=1/\sqrt3\), and write

\[
 p_\phi=A(\cos\phi,\sin\phi,r),\qquad
 g=G e_3,\qquad A>0,\quad G>0.
 \tag{1.1}
\]

The unit C149 selected line is

\[
 a_\phi={1\over2\sqrt2}
 \begin{pmatrix}
 \cos ^2\phi+2\sin ^2\phi-\sin\phi\cos\phi\\
 \sin ^2\phi+2\cos ^2\phi-\sin\phi\cos\phi\\
 -\sqrt3(\sin\phi+\cos\phi)
 \end{pmatrix},
 \qquad p_\phi\cdot a_\phi=0.
 \tag{1.2}
\]

Choose the unit circular gate polarization

\[
 E={1\over\sqrt2}(1,i,0),\qquad E\cdot g=0.
 \tag{1.3}
\]

For distinct modes use the same symmetric projected Euler symbol as C155,

\[
 \mathcal S(p,u;q,v)
 =P_{p+q}\{(u\cdot q)v+(v\cdot p)u\}.
 \tag{1.4}
\]

The forward source-to-daughter vector is

\[
 F_\phi=\mathcal S(p_\phi,a_\phi;g,E).
 \tag{1.5}
\]

This is an algebraic calculation in
\(\mathbb Q(\sqrt2,\sqrt3,i)\) after rational-circle
parameterization of \(\phi\).

## 2. Exact forward norm and a uniform one-edge bound

Set

\[
 x={G\over A},\qquad y=\sqrt3x.
 \tag{2.1}
\]

Direct projection in (1.5) gives

\[
 \boxed{
 |F_\phi|^2={A^2\over2}+A^2 f(x)|(a_\phi)_3|^2,}
 \tag{2.2}
\]

where

\[
 \begin{aligned}
 f(x)
 &=x^2-{x\over\sqrt3}
   -{2x^2\over1+(x+1/\sqrt3)^2}\\
 &={x(3x+\sqrt3)(3x^2-4)
       \over3(3x^2+2\sqrt3x+4)}\\
 &={y(y-2)(y+1)(y+2)
       \over3(y^2+2y+4)}.
 \end{aligned}
 \tag{2.3}
\]

For \(y>0\) (and algebraically also at \(y=0\)),

\[
 f+{1\over3}
 ={P(y)\over3(y^2+2y+4)},
 \qquad
 P(y)=y^4+y^3-3y^2-2y+4.
 \tag{2.4}
\]

The numerator is strictly positive by the following two exact
decompositions:

\[
 P(y)=
 \begin{cases}
 (y-1)^2(y^2+3y+2)+(2-y),&0\leq y\leq2,\\
 y^2(y^2-3)+y(y^2-2)+4,&y\geq2.
 \end{cases}
 \tag{2.5}
\]

Thus \(f(x)\geq-1/3\).  The vertical component in (1.2) obeys

\[
 |(a_\phi)_3|^2
 ={3\over8}(\sin\phi+\cos\phi)^2
 \leq {3\over4},
 \tag{2.6}
\]

because

\[
 2-(\sin\phi+\cos\phi)^2
 =(\sin\phi-\cos\phi)^2\geq0.
 \tag{2.7}
\]

Combining (2.2), (2.4), and (2.6) proves the exact lower bound

\[
 \boxed{|F_\phi|^2\geq {A^2\over4},\qquad
        |F_\phi|\geq {A\over2}.}
 \tag{2.8}
\]

At the tuned height

\[
 x={2\over\sqrt3},\qquad G={2A\over\sqrt3},
 \tag{2.9}
\]

one has \(f=0\), so the angular dependence drops out:

\[
 \boxed{|F_\phi|^2={A^2\over2}\quad\hbox{for every }\phi.}
 \tag{2.10}
\]

Equations (2.8)--(2.10) show that a missing forward edge is not the defect
of this candidate.  They say nothing about simultaneous daughter
orthogonality, collision control, or the reverse edge.

## 3. The tuned circular block is not a skew rotation

Normalize \(A=1\) and take the representative source \(\phi=0\).  Then

\[
 \begin{gathered}
 p=(1,0,1/\sqrt3),\qquad
 g=(0,0,2/\sqrt3),\\
 a={1\over2\sqrt2}(1,2,-\sqrt3),qquad
 E={1\over\sqrt2}(1,i,0).
 \end{gathered}
 \tag{3.1}
\]

The exact forward vector and its unit daughter line are

\[
 F=\mathcal S(p,a;g,E)
   =\left(0,{1-i\over2},0\right),
 \qquad
 v={F\over|F|}
   =\left(0,{1-i\over\sqrt2},0\right).
 \tag{3.2}
\]

Reality supplies the reverse gate \((-g,\overline E)\).  Returning the
daughter at \(p+g\) to the original source gives

\[
 R=\mathcal S(p+g,v;-g,\overline E)
   =\left(0,{1-i\over2},0\right)=F.
 \tag{3.3}
\]

With the Fourier Euler convention \(\dot u_k=-i\sum\mathcal S\), the two
normalized off-diagonal ODE coefficients are therefore

\[
 \kappa_{\rm f}=-i\,v^*F=-{i\over\sqrt2},
 \qquad
 \kappa_{\rm r}=-i\,a^*R=-{1+i\over2\sqrt2}.
 \tag{3.4}
\]

Their exact product is

\[
 \boxed{
 \kappa_{\rm f}\kappa_{\rm r}
 =-{1\over4}+{i\over4}.}
 \tag{3.5}
\]

If source and daughter coordinates are independently rephased and
positively rescaled, the two off-diagonal entries acquire reciprocal
factors, so their product is unchanged.  If the physical gate coefficient
is multiplied by \(\gamma\), reality multiplies the reverse gate by
\(\overline\gamma\), and (3.5) is multiplied only by
\(|\gamma|^2\).  Hence its nonreal argument cannot be removed by any of
these freedoms.

In contrast, a two-mode skew-Hermitian block satisfies

\[
 \kappa_{\rm r}=-\overline{\kappa_{\rm f}},
 \qquad
 \kappa_{\rm f}\kappa_{\rm r}
 =-|\kappa_{\rm f}|^2\in\mathbb R_{\leq0}.
 \tag{3.6}
\]

The same negative-real condition holds after any positive diagonal energy
weight.  Thus the tuned block (3.1) cannot be conjugated by phase and
positive coordinate rescaling into the C161 skew star.  This is only a
necessary two-coordinate test: the reverse vector also has components away
from the selected source line, and allowing those extra polarization
coordinates leads to a larger block not classified here.

## 4. What remains open

The obstruction (3.5) identifies reverse-edge geometry, not forward rank,
as the defect of this particular block.  The natural same-stage repairs
are still open:

* use two polarizations whose reverse defects cancel while their forward
  bright directions add;
* modulate the gate phase or polarization in time so that the integrated
  two-edge propagator is unitary even though the instantaneous block is
  not;
* retain gate depletion and solve a full finite-dimensional triad/star
  trajectory instead of freezing an autonomous gate;
* realize the construction uniformly on the balanced \(q^2\) source
  packet while controlling charge collisions, collars, pressure, and
  finite-\(\epsilon\) errors.

No claim here rules out any of those repairs, and no BAFL or Navier--Stokes
stage map follows from C162 alone.
