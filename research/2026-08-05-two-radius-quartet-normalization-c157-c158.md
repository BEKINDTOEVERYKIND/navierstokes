# C157--C158: a two-radius additive-quartet rank test and the exact packet normalization

**Date:** 2026-08-05
**Status:** exact limiting-cone algebra and exact conditional normalization;
finite \(\epsilon\), a full radial slab, localization, and the nonlinear
stage map remain open
**Checker:**
[checks/two_radius_quartet_normalization_c157_c158.py](../checks/two_radius_quartet_normalization_c157_c158.py)

## 0. Claim boundary

This note stays on the limiting C149 elliptic ray of the existing C121
\(A_2\) pump.  It takes the first step beyond C156's single resonant ring by
putting two distinct radii on the same resonant cone.

There are two exact conclusions.

1. One explicit cross-radius additive quartet has two shared-output wake
   vectors which are linearly independent for every pair of distinct
   positive radii.  Even complex scalar Fourier phases therefore cannot
   make those two channels cancel each other.  The determinant degenerates
   linearly as the radii coalesce, so this is not a uniform radial-slab
   coercivity theorem.
2. If a diagonal cubic normal form is bounded by an order-one kernel, the
   exact upper-control parameter in the reduced positive-frequency
   coefficient ledger is \(Q^2\sum|z_i|^2\).  For \(M=q^3\) equal
   independent coefficients this is \(Q^2A_{\ell^1}^2/M\), where
   \(A_{\ell^1}=\sum|z_i|\).  The scale \(\sqrt M/Q=q^{1/2}\) is only
   where this upper parameter is order one.  This is not a claim that the
   physical thick-packet kernel has a sign, a lower bound, or a dynamical
   onset threshold.

The correction to C156 is essential here: its former aggregate secular
coercivity is not used.  At the limiting rigid-rotation problem, the
physical shared-wake directed matrix is antisymmetric under its precise
bounded-wake hypotheses and can have positive balanced states.

## 1. C157: an exact two-radius additive quartet

Write

\[
 \eta_\phi=(\cos\phi,\sin\phi,1/\sqrt3),
 \qquad k_{R,\phi}=R\eta_\phi,
 \tag{1.1}
\]

and use the limiting selected unit line

\[
 a(\phi)={1\over2\sqrt2}
 \begin{pmatrix}
 \cos ^2\phi+2\sin ^2\phi-\sin\phi\cos\phi\\
 \sin ^2\phi+2\cos ^2\phi-\sin\phi\cos\phi\\
 -\sqrt3(\sin\phi+\cos\phi)
 \end{pmatrix}.
 \tag{1.2}
\]

For distinct \(A,B>0\), the first decomposition is

\[
 p=k_{A,0},\qquad q=k_{B,\pi/2}.
 \tag{1.3}
\]

Put \(D=A^2+B^2\) and define a reflected decomposition by

\[
 \begin{aligned}
 \cos\alpha&={A^2-B^2\over D},&
 \sin\alpha&={2AB\over D},\\
 \cos\beta&={2AB\over D},&
 \sin\beta&={B^2-A^2\over D},
 \end{aligned}
 \tag{1.4}
\]

\[
 p_*=k_{A,\alpha},\qquad q_*=k_{B,\beta}.
 \tag{1.5}
\]

Reflection across the horizontal output direction \((A,B)\) gives the
exact collision

\[
 p+q=p_*+q_*
   =\left(A,B,{A+B\over\sqrt3}\right).
 \tag{1.6}
\]

The decompositions are distinct when \(A\ne B\).  With the symmetric
projected Euler symbol

\[
 \mathcal S(p,u;q,v)
 =P_{p+q}\{(u\cdot q)v+(v\cdot p)u\},
 \tag{1.7}
\]

set

\[
 F=\mathcal S(p,a(0);q,a(\pi/2)),\qquad
 F_*=\mathcal S(p_*,a(\alpha);q_*,a(\beta)),
 \quad K=p+q.
 \tag{1.8}
\]

Exact projection and simplification give

\[
 \boxed{
 K\cdot(F\times F_*)=
 -{\sqrt3\,AB(A-B)(A+B)^2\over12(A^2+B^2)^3}
 \mathcal P(A,B),}
 \tag{1.9}
\]

where

\[
 \mathcal P(A,B)=
 3A^4-4A^3B+14A^2B^2-4AB^3+3B^4.
 \tag{1.10}
\]

This quartic is strictly positive.  Indeed, with
\(y=A/B+B/A\geq2\),

\[
 \mathcal P(A,B)=A^2B^2(3y^2-4y+8)>0.
 \tag{1.11}
\]

For example,
\(3y^2-4y+8=3(y-2/3)^2+20/3\), so no sign inference is being
made from numerical samples.

Thus \(F,F_*\in K^\perp\) are linearly independent whenever \(A\ne B\).
They are real vectors, so their real independence is also independence
over \(\mathbb C\): apply the real system separately to the real and
imaginary parts.  For complex scalar products \(c=z_pz_q\) and
\(c_*=z_{p_*}z_{q_*}\),

\[
 cF+c_*F_*=0\quad\Longrightarrow\quad c=c_*=0.
 \tag{1.12}
\]

The common Fourier multiplier \(i\), and the common factor caused by
summing the two orderings of each unordered pair, do not affect (1.12).
Within the four displayed positive-cone modes and their four reality
partners, the only unordered pairs producing \(+K\) for \(A\ne B\) are
\(\{p,q\}\) and \(\{p_*,q_*\}\); the conjugate pairs produce \(-K\).
Thus reality does not add a hidden cancellation at this output.  This is
nevertheless a deliberately narrow additive-quartet obstruction.  A
larger radial slab can supply more decompositions of the same output, as
well as different polarization channels, and (1.12) does not classify
those higher-multiplicity cancellations.

The loss of uniformity is exact.  At \(A=B=a\),

\[
 {K\cdot(F\times F_*)\over A-B}
 \longrightarrow -{\sqrt3\over2}a^2.
 \tag{1.13}
\]

Consequently, with the common radius restricted to a compact subset of
\((0,\infty)\), this isolated two-column inverse loses exactly one power
of the radial separation.  Adjacent normalized charge/radius layers
separated by \(\Theta(q^{-1})\) therefore have a \(\Theta(q)\) condition
loss in this two-channel solve, before finite-\(\epsilon\) or localization
errors are included.  This statement is not an upper or lower bound for
an inverse using all channels in a radial slab.

## 2. C158: exact \(M\)-coordinate diagonal normalization

Choose one representative from each reality pair and let \(z_i\) be its
complex coefficient on a unit real polarization line; the negative-mode
coefficient is \(\overline{z_i}\).  Suppose only for this reduced
positive-frequency ledger that a diagonal cubic coordinate has the form

\[
 D_i=Q^2 z_i\sum_{j\ne i}\kappa_{ij}|z_j|^2,
 \qquad |\kappa_{ij}|\leq\kappa_*.
 \tag{2.1}
\]

Writing the reduced half-lattice energy

\[
 E_+=\sum_i|z_i|^2,
 \tag{2.2}
\]

gives the uniform upper bound

\[
 \|D\|_{\ell^2}\leq Q^2\kappa_*E_+\|z\|_{\ell^2}.
 \tag{2.3}
\]

For \(M\) equal-modulus coefficients \(|z_i|=c\), define their exact
coefficient \(\ell^1\) scale by \(A_{\ell^1}=\sum_i|z_i|\).  Then

\[
 A_{\ell^1}=Mc,
 \qquad E_+=Mc^2={A_{\ell^1}^2\over M}.
 \tag{2.4}
\]

Hence the diagonal upper-control parameter is

\[
 \boxed{\Lambda_{\rm up}={Q^2A_{\ell^1}^2\over M}.}
 \tag{2.5}
\]

This convention matters.  With normalized Fourier measure and unit
polarizations, the corresponding real field has Fourier energy
\(2E_+\), not \(E_+\), and its point amplitude is at most
\(2A_{\ell^1}\), not identically \(A_{\ell^1}\).  A physical point scale
is comparable to \(A_{\ell^1}\) only when a uniform coherence projection
is supplied.  The selected cone line (1.2) does supply such a kinematic
projection: for \(e=(1,1,0)/\sqrt2\),

\[
 e\mathbin\cdot a(\phi)={3-2\sin\phi\cos\phi\over4}
 \in[1/2,1].
 \tag{2.6}
\]

Thus positive real coefficients and their reality partners, aligned at a
coherence point, give a physical \(e\)-component between
\(A_{\ell^1}\) and \(2A_{\ell^1}\).  This is comparability, not equality.

For \(M=q^3\) independent conjugate pairs and \(Q=q\), the formal
coefficient scale where \(\Lambda_{\rm up}=1\) is

\[
 A_{\ell^1}={\sqrt M\over Q}=q^{1/2}.
 \tag{2.7}
\]

On the C147 schedule \(q=n^8\), identifying its point-amplitude exponent
with the coefficient \(\ell^1\) exponent through an order-one coherence
factor, the exact exponent ledger is

\[
 \begin{array}{c|c|c}
 &A_{\ell^1}&\Lambda_{\rm up}\\ \hline
 \text{seed}&n^{-16}=q^{-2}&n^{-40}=q^{-5}\\
 \text{formal unit-upper-parameter scale}&n^4=q^{1/2}&1\\
 \text{target}&n^{10}=q^{5/4}&n^{12}=q^{3/2}.
 \end{array}
 \tag{2.8}
\]

The seed-to-unit-upper-parameter gain is \(n^{20}\), whereas C147 asks
for \(n^{26}\); its target coefficient scale is larger by \(n^6\).
Because (2.3) is only an upper bound, crossing this scale proves neither
growth nor obstruction.

This does not reinstate a cubic no-go.  Equation (2.3) is an upper bound,
not a lower bound.  The corrected fixed-ring normal form has exact
antisymmetric/balanced cancellations, and radial additive quartets create
off-diagonal phase-dependent terms not represented in (2.1).  Conversely,
no claim is made here that those cancellations persist across \(q\) radial
layers.

## 3. Exact remaining obstruction

C157 rules out the simplest hope that this cross-radius shared output can
be cancelled by its reflected partner alone with scalar phases.  It does
not rule out cancellation after other decompositions are added.  C158
shows the coefficient normalization any surviving diagonal upper-bound
argument must use.  Neither result classifies the full \(q^3\)-packet
additive energy.

The next exact target on this geometry is therefore:

> classify the multi-decomposition quadratic wake map on a finite radial
> slab together with the finite-\(\epsilon\) physical shared-wake kernel,
> and either construct a positive balanced invariant distribution that
> survives the slab coupling or prove that one named output retains a
> uniform component after all decompositions are summed.

Finite \(\epsilon\), the \(O(\log q)\) gain interval, the neutral-wake
invariant graph, lattice discretization, and physical localization remain
separate obligations.  No one-cell Navier--Stokes stage or Millennium
conclusion is claimed.
