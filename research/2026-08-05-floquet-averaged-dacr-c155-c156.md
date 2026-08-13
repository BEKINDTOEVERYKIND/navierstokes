# C155--C156: Floquet-averaged DACR and the fixed-ring cancellation boundary

**Date:** 2026-08-05
**Status:** exact directed limiting-cocycle Picard coefficient and exact
fixed-ring quadratic-null algebra.  The formerly claimed pair-symmetric
cubic coercivity has been **retracted by exact countercalculation**; the
finite-\(\epsilon\) localized thick-packet stage remains open.
**Checker:**
[checks/floquet_averaged_dacr_c155_c156.py](../checks/floquet_averaged_dacr_c155_c156.py)

## 0. Claim boundary

This note stays on the C149 elliptic ray of the C121 \(A_2\) pump.  It
audits a gap in C151: a nonzero *instantaneous* cubic symbol need not have a
nonzero projection after a complete Kelvin period.  Here the complete
order-two wake is evolved from zero, returned once against the parent
modes, and paired with the selected left Floquet solution.

The directed answer is nonzero.  For the quarter-period pair the
first-period coefficient returning to the phase-zero parent is an explicit
number less than \(-21/8\).  For frozen parents over many limiting periods,
its directed secular kernel can also be computed for every angular
separation.  Exact cancellation of the whole quadratic complement is
possible only on one Fourier ray.

An adversarial reverse-return calculation changes the cubic conclusion.
At quarter separation the secular coefficient returning to the first
parent is \(-9\pi/8\), while the coefficient returning to the second parent
is \(+9\pi/8\).  More generally, the order-four energy identity forces the
two physical directed coefficients of every fixed interior pair to be
exact negatives.  The earlier instruction to obtain the reverse
coefficient by merely replacing \(t\) with \(-t\) held the target left line
fixed and was not the physical reverse return.  Consequently the former
pair-symmetric formula and the cubic phase-cancellation no-go have been
removed.

These are not a no-go theorem for the proposed \(q^3\)-mode child packet.
C145 already forces that packet to have radial/normal thickness.  Different
rings admit additive quartets that the calculation below does not control;
finite \(\epsilon\), evolving amplitudes, localization, and the
neutral-wake invariant graph remain open.

## 1. Limiting resonant modes

Use the orthonormal coordinates of C149 and put

\[
 J=\begin{pmatrix}0&-1&0\\1&0&0\\0&0&0\end{pmatrix},
 \qquad r={1\over\sqrt3},
 \qquad \eta_\phi=(\cos\phi,\sin\phi,r).
 \tag{1.1}
\]

At \(\epsilon=0\), the rotating-frame Kelvin operator at a fixed covector
\(\eta\) is

\[
 C_\eta=-2J+2{\eta\eta^T\over|\eta|^2}J.
 \tag{1.2}
\]

The expanding line selected by the \(\epsilon>0\) splitting of C149 tends
to the
unit periodic solution

\[
 a(\phi)={1\over2\sqrt2}
 \begin{pmatrix}
 \cos ^2\phi+2\sin ^2\phi-\sin\phi\cos\phi\\
 \sin ^2\phi+2\cos ^2\phi-\sin\phi\cos\phi\\
 -\sqrt3(\sin\phi+\cos\phi)
 \end{pmatrix}.
 \tag{1.3}
\]

Thus the two phase-shifted rotating-frame modes are

\[
 b_1(s)=R(-s)a(s),\qquad
 b_2(s)=R(-s)a(s+\phi),
 \tag{1.4}
\]

at the fixed covectors η_0 and η_φ.  For the first mode the selected
left solution may be represented by

\[
 \psi(s)=
 \begin{pmatrix}
 2\sin(s+\pi/4)-3\sqrt2/4\\
 \cos(s+\pi/4)\\
 -\sqrt6/4
 \end{pmatrix},
 \qquad \psi(s)\cdot b_1(s)=1.
 \tag{1.5}
\]

It differs from \(b_1\) only by a multiple of \(\eta_0\), so both give the same
pairing against divergence-free forcing.

For later use, write the symmetric projected Euler symbol as

\[
 \mathcal S(p,A;q,B)
 =P_{p+q}\{(A\cdot q)B+(B\cdot p)A\}.
 \tag{1.6}
\]

In the common rotating frame the Fourier coefficients obey

\[
 u_p'=C_pu_p-i\sum_{q+r=p}
 P_p\{(u_q\cdot r)u_r\}.
\]

For distinct \(q,r\), the two ordered convolution terms are exactly the
single \(\mathcal S(q,u_q;r,u_r)\) in (1.6); no additional factor \(2\)
is missing. We use ordinary complex Taylor coefficients rather than raw
Fréchet derivatives, so no factorial is inserted in the two successive
Duhamel integrations. Reality gives
\(u_{-\eta_j}=\overline{u_{\eta_j}}\). The negative-frequency paths
produce the conjugate return at \(-\eta_0\); they do not double the
coefficient at \(+\eta_0\).

## 2. C155: exact full-period return for the quarter pair

Take φ=π/2.  The complete order-two real interaction has a sum wake at
η_0+η_{π/2} and a difference wake at
η_0-η_{π/2}.  In the respective transverse bases

\[
 E_\Sigma=((1,-1,0),(1,1,-\sqrt3)),\qquad
 E_\Delta=((1,1,0),(0,0,1)),
 \tag{2.1}
\]

their exact equations, with zero initial wake, are

\[
 w_\Sigma'=A_\Sigma w_\Sigma+g_\Sigma,
 \quad
 A_\Sigma=\begin{pmatrix}0&2\\-4/5&0\end{pmatrix},
 \quad g_\Sigma=\binom0{3/10},
 \tag{2.2}
\]

and

\[
 w_\Delta'=g_\Delta(s),\qquad
 g_\Delta(s)=\binom{\frac12\sin2s}
 {-\frac{\sqrt3}{2}\sin2s}.
 \tag{2.3}
\]

Put \(\Omega=\sqrt{8/5}\). Direct integration gives

\[
 w_\Sigma(s)=
 \binom{\frac38(1-\cos\Omega s)}
 {\frac{3\Omega}{16}\sin\Omega s},
 \qquad
 w_\Delta(s)=
 \binom{\frac14(1-\cos2s)}
 {-\frac{\sqrt3}{4}(1-\cos2s)}.
 \tag{2.4}
\]

The selected growing-line rows for the two returns are

\[
 h_\Sigma(s)=
 \left({\cos2s\over4}-1,-{5\sin2s\over4}-{3\over4}\right),
 \tag{2.5}
\]

\[
 h_\Delta(s)=
 \left(-{\sin2s\over2}-{3\cos2s\over4}-{1\over2},
 \sqrt3\left({1\over12}+{\sin2s\over4}-{\cos2s\over4}\right)
 \right).
 \tag{2.6}
\]

These rows already include the two Euler factors \((-i)^2=-1\).  Hence the
two contributions over one covector period are

\[
 I_\Sigma={3\over64}
 \left[-16\pi-3+3\cos\Theta+4\sqrt{10}\sin\Theta\right],
 \qquad
 I_\Delta=-{3\pi\over8},
 \tag{2.7}
\]

where

\[
 \Theta=2\pi\Omega={4\sqrt{10}\pi\over5}.
 \tag{2.8}
\]

Therefore the complete causal first-period coefficient is

\[
 \boxed{
 \mathfrak C_0={3\over64}
 \left[-24\pi-3+3\cos\Theta+4\sqrt{10}\sin\Theta\right]
 =-3.097745763580765\ldots .}
 \tag{2.9}
\]

This is rigorously nonzero without numerical transcendental evaluation:
using \(\pi>3\), \(\sqrt{10}<4\), and the
\(\Theta\)-independent bounds \(\cos\Theta\le1\),
\(\sin\Theta\le1\), gives

\[
 \boxed{\mathfrak C_0<-{21\over8}<0.}
 \tag{2.10}
\]

If the positive Fourier coefficients of the two rays are \(z_1,z_2\), the
sum/difference wakes carry \(z_1z_2\) and
\(z_1\overline{z_2}\). Returning
them with the appropriate reality partner gives exactly

\[
 \mathfrak C_0 z_1|z_2|^2.
 \tag{2.11}
\]

Thus spatial scalar phases cannot remove this two-ray coefficient.

The tuned monodromy has the simple first-order splitting computed in C149.
After choosing continuous unit right and dual-left branches for
\(\epsilon>0\) with the displayed \(\epsilon\downarrow0\) limits, analytic
ODE dependence gives
\(\mathfrak C_\epsilon=\mathfrak C_0+O(\epsilon)\).
Consequently (2.9) persists with the same sign for every sufficiently
small positive \(\epsilon\). This is an existence statement for a small
interval in the frozen principal cocycle, not a validated numerical
endpoint or a finite-\(\epsilon\) PDE packet.

After wave-number scaling by \(Q\), and with zero incoming complement, the
parent projection of the one-period map has the local form

\[
 z_1^{\rm out}=\rho_\epsilon z_1
   +Q^2\mathfrak C_\epsilon z_1|z_2|^2+O(Q^3|z|^4),
 \tag{2.12}
\]

for \(Q|z|\ll1\), with \(|z|=|z_1|+|z_2|\) and fixed
\(\epsilon\). The endpoint complement is already \(O(Q|z|^2)\) and is
generally nonzero. In fact (2.4) gives
\(w_\Sigma(2\pi)\ne0\), because
\(\Omega=\sqrt{8/5}\notin\mathbb Z\). Comparing the cubic term with the
whole retained amplitude gives the single-pair scale
\(Q|z|\asymp1\). When \(|z_1|\) and \(|z_2|\) are comparable, comparison
with the small one-period linear increment
\(\rho_\epsilon-1=O(\epsilon)\) gives the sharper small-\(\epsilon\)
balance

\[
 |z|\asymp {\sqrt\epsilon\over Q}.
 \tag{2.13}
\]

For fixed small \(\epsilon\), this is still \(Q^{-1}\) up to a constant.
This is a threshold for one Fourier coefficient in the two-ray slice.
The C147 quantity \(q^{-2}\) is instead a total coherent point seed; for
\(q^3\) comparable modes its per-mode normalization and the sum of
cross-channel coefficients must be computed separately. No C147 packet
gain bound follows directly from (2.13).

Equation (2.12) is not by itself an iterable scalar map. The sum wake in
(2.4) does not reset at \(2\pi\), so the next period depends on the
incoming complement. A repeated stage must evolve or slave that wake as
part of the state; neither a reset nor an invariant graph is proved here.

## 3. C156: exact quadratic null equations on one resonant ring

Fix one positive normal charge and one limiting resonant radius,

\[
 \Gamma=\{\eta_\phi:\phi\in\mathbb R/2\pi\mathbb Z\}.
 \tag{3.1}
\]

For a finite packet

\[
 v_1=\sum_j z_j a(\phi_j)e^{i\eta_{\phi_j}\cdot x}+\text{c.c.},
 \tag{3.2}
\]

exact cancellation of every quadratic complement coefficient means

\[
 \sum_{i<j:\,\eta_i+\eta_j=k}
 z_i z_j\,\mathcal S(\eta_i,a_i;\eta_j,a_j)=0,
 \tag{3.3}
\]

\[
 \sum_{i,j:\,\eta_i-\eta_j=k}
 z_i\overline{z_j}\,
 \mathcal S(\eta_i,a_i;-\eta_j,a_j)=0
 \tag{3.4}
\]

for every output \(k\). Self interactions vanish.

After a common phase shift take \(\phi_i=0\), \(\phi_j=\phi\), and put
\(t=\tan(\phi/2)\). Exact projection gives

\[
 \mathcal S(\eta_0,a_0;\eta_\phi,a_\phi)
 ={2t^2(t^2+2)\over(1+t^2)^2(t^2+4)}
 (1,t,-\sqrt3),
 \tag{3.5}
\]

\[
 \mathcal S(\eta_0,a_0;-\eta_\phi,a_\phi)
 ={t(t^2-1)\over(1+t^2)^2}
 (1,t,-\sqrt3).
 \tag{3.6}
\]

For two unit horizontal vectors, a nonzero horizontal sum determines the
unordered pair uniquely: if \(s=u+v\ne0\) and \(d=u-v\), then
\(d\perp s\) and \(|d|^2=4-|s|^2\), so in two dimensions \(d\) is fixed
up to sign, which only swaps \(u,v\). Therefore every non-antipodal pair
in (3.3) has a unique output and (3.5) is nonzero. An antipodal pair has
zero sum
symbol, but the limit of (3.6) at \(\phi=\pi\) is a nonzero horizontal
vector, explicitly \((0,1,0)\) in the normalized coordinates, and its
difference output is unique once only those two phases remain. It follows
that

\[
 \boxed{\mathbb P(v_1\cdot\nabla v_1)=0\quad\Longrightarrow\quad
 \#\{\phi_j:z_j\ne0\}\le1}
 \tag{3.7}
\]

on the fixed ring (3.1).  A one-ray reality pair is the familiar nonlinear
plane wave.  It supplies only one Fourier direction, not \(q^3\) coherent
degrees. This is an exact fixed-ring no-go for cancelling *all* of the
\(QA^2\) complement; it is not a thick-packet no-go.

## 4. The directed secular cubic kernel

There is a second fixed-ring result even when the wake is retained. Let
\(0<\phi<\pi\), \(t=\tan(\phi/2)>0\), and keep the scalar coefficients of
the limiting linear parent solutions fixed. The sum wake is a constant
equilibrium plus a bounded oscillator. In physical rotating coordinates
its equilibrium is

\[
 w_\Sigma^{\rm eq}=
 \left({t^3(t^2+2)\over(1+t^2)^3},
       -{t^2(t^2+2)\over(1+t^2)^3},0\right).
 \tag{4.1}
\]

The zero-initial difference wake is periodic and equals

\[
 w_\Delta(s)=F_t(s)(t,t^2,-\sqrt3t),
 \tag{4.2}
\]

\[
 F_t(s)={(t^2-1)\sin2s-2t\cos2s+2t
          \over2(1+t^2)^2}.
 \tag{4.3}
\]

There is no per-period reset in this calculation.  For all \(s\geq0\), the
causal sum wake is exactly

\[
 w_\Sigma(s)=w_\Sigma^{\rm eq}
              -e^{sC_{\eta_0+\eta_\phi}}w_\Sigma^{\rm eq},
 \tag{4.4}
\]

while (4.2) is periodic and satisfies \(w_\Delta(2\pi N)=0\).  Thus the
persistent full complement state, rather than a freshly zeroed wake, is
used in every later period.

The sum oscillator has frequency

\[
 \Omega_t={2\over\sqrt{1+3\cos^2(\phi/2)}}\in(1,2).
 \tag{4.5}
\]

It is therefore nonresonant with the zero and second harmonics for every
fixed interior separation. The return row contains only those harmonics,
so integrating its product with the sum oscillator gives a bounded
transient. Over \(N\) full periods, direct trigonometric integration gives

\[
 \Pi_1v_3(2\pi N)
 =\{N\mathfrak K(t)+O_t(1)\}z_1|z_2|^2,
 \tag{4.6}
\]

where, for each fixed \(t\in(0,\infty)\), the bounded term is the
nonresonant sum-wake transient. The \(O_t(1)\) constant is not uniform as
\(t\to\infty\), where \(\Omega_t\to2\); no endpoint-uniform estimate as
\(t\to0\) is claimed either. The directed secular coefficient is

\[
 \boxed{
 \mathfrak K(t)=
 -{\pi t^3(t^2+2)(t^4+2t^3+3t^2+2t+4)
       \over(1+t^2)^5}<0.}
 \tag{4.7}
\]

At the quarter separation this is \(-9\pi/8\) per period; the first-period
number (2.9) differs because the nonresonant sum wake starts from zero.
Equation (4.7) is a **directed** return: its target is the phase-zero mode
and its left line is (1.5).  It is not legitimate to obtain the physical
reverse return by changing only \(t\) to \(-t\), because the target mode and
its selected left line must also be changed.

The failure is already exact at quarter separation.  In a common exact
transverse-coordinate calculation, the period means of the persistent
sum-wake and zero-initial difference-wake returns to the phase-zero parent
are respectively

\[
 -{3\over8},\qquad -{3\over16},
 \tag{4.8}
\]

so its secular coefficient per period is

\[
 2\pi\left(-{3\over8}-{3\over16}\right)=-{9\pi\over8}.
 \tag{4.9}
\]

For the physical reverse return, with the phase-\(\pi/2\) right/left line
as target, the exact means are instead

\[
 +{3\over8},\qquad +{3\over16},
 \tag{4.10}
\]

and hence

\[
 2\pi\left(+{3\over8}+{3\over16}\right)=+{9\pi\over8}.
 \tag{4.11}
\]

Thus the pair-symmetric secular coefficient at this pair is exactly zero,
not \(-3\pi/4\).  The checker reconstructs both sum and difference means
using exact rational trigonometric moments.

There is also a coordinate-free reason this cancellation must occur at the
limiting rigid-rotation problem.  Let \(L\) denote the skew \(L^2\) linear
operator at \(\epsilon=0\), and expand a perturbation as

\[
 v=\alpha v_1+\alpha^2v_2+\alpha^3v_3+O(\alpha^4),
 \tag{4.12}
\]

with \(v_2(0)=v_3(0)=0\).  Using
\(\langle B(u,v),v\rangle=0\) and the skewness of \(L\), the order-four
energy identity is

\[
 {d\over ds}\left(\|v_2\|_2^2
       +2\operatorname{Re}\langle v_1,v_3\rangle\right)=0,
 \qquad
 2\operatorname{Re}\langle v_1,v_3\rangle=-\|v_2\|_2^2.
 \tag{4.13}
\]

Indeed,

\[
 v_2'=Lv_2-B(v_1,v_1),\qquad
 v_3'=Lv_3-B(v_1,v_2)-B(v_2,v_1).
\]

If \(b(u,v,w)=\langle B(u,v),w\rangle\), incompressibility gives
\(b(u,v,w)=-b(u,w,v)\).  Hence

\[
 {d\over ds}\|v_2\|_2^2
   =-2\operatorname{Re}b(v_1,v_1,v_2),
\]

while the derivative of the cross term is the opposite number:
the \(L\)-terms cancel by skewness,
\(b(v_2,v_1,v_1)=0\), and
\(b(v_1,v_2,v_1)=-b(v_1,v_1,v_2)\).

Therefore a bounded persistent quadratic wake cannot produce a nonzero
**aggregate retained-energy secular drift** at \(\epsilon=0\).  Individual
mode coordinates can and do have nonzero directed secular coefficients;
their energy-weighted sum must cancel.  This identity does not remove the
finite-time wake, and it does not extend without an error term to the
strained \(\epsilon>0\) cocycle.

More precisely, consider an isolated unit-normalized two-parent pair for
which the quadratic sum/difference wakes are bounded and there are no
other secular retained paths.  If \(K_{ij}\) is the physical directed
coefficient returning the shared wake to parent \(i\), then the coefficient
of \(|z_i|^2|z_j|^2\) in (4.13) is \(K_{ij}+K_{ji}\).  These coefficients
are real in the present real limiting-mode normalization.  Since the two
parent weights are arbitrary, (4.13) gives the exact antisymmetry

\[
 \boxed{K_{ji}=-K_{ij}.}
 \tag{4.14}
\]

This is the correct general reverse rule.  It is not the formal
substitution \(t\mapsto-t\) in (4.7), which leaves the wrong target line
fixed.

In particular, when parent \(1\) is the phase-zero target and parent \(2\)
has half-angle parameter \(t>0\), (4.7) and (4.14) give the physical reverse
coefficient for every fixed interior separation,

\[
 K_{21}(t)=-\mathfrak K(t)
 ={\pi t^3(t^2+2)(t^4+2t^3+3t^2+2t+4)\over(1+t^2)^5}.
\]

This is a general statement, not a quarter-pair coincidence.  For example,
at \(t=1\) it gives \(+9\pi/8\), whereas the invalid fixed-target
substitution gives \(\mathfrak K(-1)=+3\pi/8\).

The antisymmetry admits a small positive repair rather than a coercive
no-go.  Take three fixed-ring phases with half-angle parameters

\[
 (t_1,t_2,t_3)=\left(-3,-1,{1\over2}\right).
 \tag{4.15}
\]

Direct exact evaluation of the shared-wake normal form gives, in units of
\(\pi\),

\[
 K_{12}=-{1746\over3125},\qquad
 K_{13}={227409\over78125000},\qquad
 K_{23}=-{3861\over25000},
 \qquad K_{ji}=-K_{ij}.
 \tag{4.16}
\]

The strictly positive squared-amplitude vector

\[
 w=\left({3861\over25000},
          {227409\over78125000},
          {1746\over3125}\right)
 \tag{4.17}
\]

satisfies \(Kw=0\) exactly.  Thus the frozen limiting diagonal normal form
\(\dot z_i=z_i\sum_jK_{ij}|z_j|^2\) has a nontrivial positive balanced
three-ray state.  This does **not** cancel the nonzero quadratic wake, prove
that the balance survives finite \(\epsilon\), handle radial/normal
thickening or additive quartets, or construct a localized invariant graph.

The safe conclusions are now:

1. C155 and the directed kernel (4.7) remain exact for their stated target.
2. C156's quadratic fixed-ring null theorem (3.7) remains exact.
3. The former pair-symmetric coercivity formula and cubic
   phase-cancellation no-go are retracted; (4.14)--(4.17) give the exact
   antisymmetric replacement and one finite positive balance.  Different
   radii/charges, additive quartets, finite
   \(\epsilon\), slowly varying amplitudes, localization, and the neutral
   wake invariant graph all remain open.

Accordingly C156 is a quadratic fixed-ring obstruction plus a corrected
directed-cubic audit boundary, not a cubic coercivity theorem and not an
unforced one-cell stage theorem.
