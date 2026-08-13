# C161: a normalized terminal charge-star ledger from \(q^2\) to \(q^3\) modes

**Date:** 2026-08-05
**Status:** exact Hilbert-space rotation, support/cardinality, gate-darkness,
and exponent ledger; the Navier--Stokes realization and BAFL remain open
**Checker:**
[checks/normalized_terminal_charge_star_c161.py](../checks/normalized_terminal_charge_star_c161.py)

## 0. Claim boundary

This note stays on the C121--C160 \(A_2\) stage.  It tests a just-in-time
repair suggested by the corrected DACR algebra:

1. carry only a balanced \(q^2\)-mode resonant-cone packet through the long
   amplifier, obtaining coefficient-\(\ell^1\)/\(\ell^2\) gain \(q\);
2. after the gain, apply one energy-preserving \(q\)-way charge splitter to
   create \(q^3\) endpoint modes and the missing \(\sqrt q\) coefficient
   gain.

The exact arithmetic works, and the splitter has a useful structured
support realization: \(q\) pure-normal charge shifts serve all \(q^2\)
source modes.  The correct norm is an **abstract collective star** with
per-daughter coefficient \(q^{-1/2}\).  Treating the daughters as \(q\)
independent gates gives the wrong energy and time ledger.  The abstract
star coefficient is not yet a physical gate amplitude.
Likewise, the coefficient \(\ell^1\) scales below are only upper scales for
a physical point value.  Turning them into a lower bound for one component
of the real velocity requires a uniform polarization/phase coherence chart,
which C161 does not supply.

C161 does **not** prove that the Leray-projected source--gate coefficients
form this star uniformly over the source sheet.  It also assumes, without
proving, that a balanced \(q^2\)-mode packet survives the long gain.  The
only exact positive balance currently exhibited by the corrected C156
calculation is finite-dimensional; extending it to \(q^2\) modes remains
an obligation.

The cardinality also matches C154's geometry.  Its shear-free
periodic-covector set is the two-dimensional plane \(u^\perp\), so a
\(q\)-scale lattice window in that plane has the dimensional capacity
\(O(q^2)\).  The terminal charge star supplies the third Fourier coordinate
only after the long gain.  The fixed-charge product sheet displayed in
Section 2 is a support model and generally is **not** the tilted plane
\(u^\perp\).  A uniform lattice/microlocal approximation of that periodic
plane, compatible with the charge translations below, is part of the open
PDE realization.

## 1. The exact normalized star

Let \(e_0\) denote one source coordinate and
\(e_1,\ldots,e_q\) its desired charge-shift daughters.  Define the unit
bright daughter

\[
 B_q={1\over\sqrt q}\sum_{r=1}^q e_r.
 \tag{1.1}
\]

The skew star generator \(G_q\) is

\[
 G_qe_0=B_q,qquad G_qB_q=-e_0,
 \qquad G_q|_{\{e_0,B_q\}^\perp}=0.
 \tag{1.2}
\]

Thus

\[
 G_q^*=-G_q,qquad
 G_q^2=-I\quad\hbox{on }\operatorname{span}\{e_0,B_q\},
 \tag{1.3}
\]

and the exact quarter rotation satisfies

\[
 e^{(\pi/2)G_q}e_0=B_q.
 \tag{1.4}
\]

It preserves energy and sends a source coefficient \(c\) to \(q\)
coherent daughter coefficients \(c/\sqrt q\).  Every daughter-dark
combination, characterized by coefficient sum zero, is annihilated by
\(G_q\).

The normalization is load-bearing.  If the collective pulse strength is
\(\theta\), then

\[
 q\left({\theta\over\sqrt q}\right)^2=\theta^2.
 \tag{1.5}
\]

Hence the operator strength is \(\theta\), independent of \(q\).  Giving
each daughter the unnormalized strength \(\theta\) instead would have
squared row norm \(q\theta^2\).

## 2. A charge-translation support chart

Use the C144 product coordinates, where \(N,r_1,r_2\) are linearly
independent,

\[
 s_{bc}=m_0N+b r_1+c r_2,qquad 0\leq b,c<q,
 \tag{2.1}
\]

for a fixed-normal-charge source sheet \(S\), so \(|S|=q^2\).  Choose
\(q\) distinct nonzero pure-normal shifts

\[
 g_a=aN,qquad a\in\mathcal A,qquad |\mathcal A|=q.
 \tag{2.2}
\]

Then

\[
 T=\{s_{bc}+g_a:(b,c,a)\in[0,q)^2\times\mathcal A\}
 \tag{2.3}
\]

has exactly

\[
 |T|=q^3.
 \tag{2.4}
\]

Choose \(m_0\) and the shift band so that

\[
 T\cap(S\cup(-S)\cup(-T))=\varnothing.
 \tag{2.5}
\]

For example, \(|a|\le q/2\) and \(m_0=4q\) have this separation in the
product chart.  This matches the source/daughter decomposition in (1.1)
after reality completion.

Each single shift \(g_a\) serves all \(q^2\) sources.  Therefore the
first-step splitter requires only \(q\) structured gate frequencies, not
one frequency for every source--daughter pair.  This does not contradict
C146: its \(M/6\) lower bound concerns collapsing \(M\) sources to six
fixed targets, whereas (2.3) preserves the two source labels and creates
translated target sheets.

For a real gate bundle the shift set must be closed under \(a\mapsto-a\).
This is compatible with exactly \(q\) distinct nonzero shifts when \(q\) is
even, for example

\[
 \mathcal A=\{-q/2,\ldots,-1,1,\ldots,q/2\}.
 \tag{2.6}
\]

We therefore take the harmless even subsequence of the stage schedule below.
Without this parity restriction, starting from \(q\) one-sided shifts and
then imposing reality produces \(2q\) shifts and a \(2q\)-daughter star; its
power exponents are unchanged, but (2.4) is no longer the literal
cardinality.  Boundary layers and repeated charge sums still require
control.

There is also an exact gate-gate darkness.  Fix \(E\perp N\) and let

\[
 W(x)=E F(N\cdot x)
 \tag{2.7}
\]

be any real finite Fourier bundle supported on the shifts (2.2).  Then

\[
 \nabla\cdot W=0,
 \qquad
 (W\cdot\nabla)W=(E\cdot N)E F F'=0.
 \tag{2.8}
\]

Thus the pure-normal gate bundle is an exact nonlinear shear.  Equation
(2.8) eliminates gate--gate self-interaction, but not source--gate symbol
variation, daughter--gate repeated shifts, Leray pressure, or localization
collars.

## 3. Exact C147 cardinality and exponent ledger

Retain the schedule along even integers \(n\) (hence even \(q\))

\[
 q=n^8,qquad b=n^{-2},qquad H=n^{26},qquad
\varepsilon=n^{-28}.
\tag{3.1}
\]

Restricting to this cofinal subsequence changes no exponent.  The long
packet has \(M_0\) equal independent half-lattice coordinates, where

\[
 M_0=q^2=n^{16}
 \tag{3.2}
\]

Here and below \(M_0,M_1\) and the \(\ell^2\) norms use independent
half-lattice Fourier coordinates.  Reality adds the \(-S\) and \(-T\)
companions.  In the full signed lattice this doubles the mode counts and
changes the displayed real-field norms by fixed convention-dependent
factors only; none of the \(n,q,b\) exponents changes.  The seed coefficient,
seed coefficient-\(\ell^1\) scale, and post-gain coefficient are

\[
 c_{\rm seed}={\varepsilon\over q}=n^{-36},
 \qquad
 A_{\ell^1,{\rm seed}}=q^2c_{\rm seed}=n^{-20},
 \qquad
 c_0=Hc_{\rm seed}={b\over q}=n^{-10}.
 \tag{3.3}
\]

After the long gain,

\[
 \|z\|_{\ell^2}=q c_0=b=n^{-2},
 \qquad
 A_{\ell^1,0}=q^2c_0=bq=n^6.
 \tag{3.4}
\]

The normalized star produces

\[
 M_1=q^3=n^{24},
 \qquad
 c_1={c_0\over\sqrt q}=n^{-14}.
 \tag{3.5}
\]

Therefore

\[
 \sqrt{M_1}c_1=n^{12}n^{-14}=n^{-2}=b,
 \tag{3.6}
\]

while the coefficient \(\ell^1\) scale is

\[
 M_1c_1=n^{24}n^{-14}=n^{10}.
 \tag{3.7}
\]

The terminal star preserves the scheduled \(L^2\) size and supplies
exactly the missing \(\sqrt q=n^4\) coefficient-\(\ell^1\) gain.  With the
half-lattice convention above, a physical real-field point value is at most
twice this scale by the triangle inequality.  Equality up to that fixed
reality factor, or even a uniform lower comparison, is conditional on the
open coherence chart stated above; varying Leray polarizations can otherwise
cancel at a point.

## 4. Compatibility with the C146 unfolded powers

Put

\[
 J=b^{-1}=n^2,
 \qquad \theta=b=n^{-2}.
 \tag{4.1}
\]

Use normalized splitter time \(\tau\in[0,1]\) and generator
\((\pi/2)G_q\).  Split the exact propagator into \(J\) identical factors,

\[
 \left(e^{(\pi/2J)G_q}\right)^J=e^{(\pi/2)G_q}.
 \tag{4.2}
\]

Thus \(J\theta=1\), with \(\theta=J^{-1}=b\), is a normalized-time
ledger; the actual star angle of one factor is \(\pi/(2J)\), not \(b\).
Its leading per-daughter increment is

\[
 {\pi\over2}{\theta c_0\over\sqrt q}
 +O\!\left({\theta^2c_0\over\sqrt q}\right)
 \asymp n^{-16}.
 \tag{4.3}
\]

The exact product (4.2), rather than literal addition of the leading
increments, produces the full daughter coefficient \(c_1\).  Across all
\(q^3\) targets, the leading one-factor target scale is

\[
 \sqrt{q^3}\,{\theta c_0\over\sqrt q}
 \asymp b\theta=b^2,
 \tag{4.4}
\]

with the fixed \(\pi/2\) factor suppressed in the power ledger.  Hence

\[
 Jb\theta=b
 \tag{4.5}
\]

at the exponent level.  Under C146's displayed **collective per-factor**
hypotheses, the wake and active-return power sums are

\[
 Jb\theta^2=b^2=n^{-4},
 \qquad
 Jb\theta^3=b^3=n^{-6}.
 \tag{4.6}
\]

Thus the normalized star does not exceed the retained wake or active
return budgets at the scalar-power level.

The checker verifies the exact monomials in \(n,q,b,J\) after stripping
the fixed transfer-angle constants.  It does not assert that \(J\) literal
Euler increments of angle \(b\) sum to a quarter turn, and it does not
identify the abstract star angle with a physical Fourier gate amplitude.

This calculation cannot be replaced by \(q\) independent daughter gates.
First, \(q/J=n^6\), so they do not fit the \(J=b^{-1}\) sequential pulse
budget.  Second, assigning unnormalized strength \(\theta=b\) to every
daughter gives

\[
 q\theta^2=qb^2=n^4>1,
 \tag{4.7}
\]

which overfills the source energy.  The \(q^{-1/2}\) bright normalization
and a leading source-depletion rotation are essential.

## 5. The exact remaining theorem

C161 shows that the endpoint cardinality and all scheduled exponents are
compatible with a collective charge star.  It does not establish the PDE
block needed to use it.  The load-bearing statement is:

> **Uniform normalized charge-star realization.** On the localized
> \(A_2\) parent, construct source and pure-normal gate polarizations such
> that the Leray-projected source--daughter block is, uniformly over the
> \(q^2\) source sheet, the same skew star (1.2) up to a factorial-safe
> error; incorporate order-one source depletion into the leading
> trajectory; and bound all repeated charge collisions, boundary layers,
> pressure/collar terms, and the C125 backward-weighted residual by the
> C146 powers independently of \(q\).  The chart must also give a fixed
> physical component whose endpoint point value is uniformly comparable
> from below to the coefficient \(\ell^1\) scale.

The following remain explicitly open:

1. an exactly balanced \(q^2\)-mode long-gain packet;
2. a uniform lattice/microlocal approximation of C154's tilted
   shear-free plane by the source sheet;
3. the uniform Leray skew-star coefficient and its polarization chart;
4. a uniform physical coherence projection converting the coefficient
   \(\ell^1\) scale into an endpoint point-value lower bound;
5. the leading source-depletion trajectory through all \(J\) factors;
6. daughter--gate cross-charge collisions and radial-band boundaries;
7. localization collars, global pressure, C125 relative control, and BAFL.

No one-cell Navier--Stokes stage or Millennium conclusion is claimed.
