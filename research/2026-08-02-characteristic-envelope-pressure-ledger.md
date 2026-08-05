# Characteristic-envelope pressure ledger

**Date:** 2026-08-02
**Status:** exact local algebra and scale ledger; SELF/ARITH-CHECKED
**Scope:** one two-wave colour with real scalar envelopes. This is not a
Navier--Stokes transition theorem.

## 1. Setup

Let \((e,v,t)\) be an orthonormal frame, let

\[
 w=W e,\qquad W>0,\qquad
 a=w-\delta t,\qquad b=w+\delta t,
\tag{1.1}
\]

and take the two carrier covectors

\[
 k=\Lambda t+\frac Q2 r,\qquad
 \ell=\Lambda t-\frac Q2 r,\qquad
 w\cdot r=1,\qquad
 \delta=\frac{Q}{2\Lambda}.
\tag{1.2}
\]

Then \(k\cdot a=\ell\cdot b=0\), \(k+\ell=K=2\Lambda t\),
and

\[
 (a\cdot\ell)b+(b\cdot k)a
 =-\frac{Q^2}{\Lambda}t=cK,
 \qquad c=-\frac{Q^2}{2\Lambda^2}=-2\delta^2.
\tag{1.3}
\]

For real scalar envelopes \(A,B\), put \(P=AB\). The common-high-sum
interaction between

\[
 Aa e^{ik\cdot x}\quad\hbox{and}\quad Bb e^{i\ell\cdot x}
\]

is

\[
 N_{\rm hi}=e^{iK\cdot x}\{F+i cKP\},
\qquad
F=A(a\cdot\nabla B)b+B(b\cdot\nabla A)a.
\tag{1.4}
\]

This formula includes the envelope derivatives that are invisible in a
frozen-symbol calculation.

## 2. Necessary and sufficient pressure condition

On a simply connected domain (with decay, or with no zero total Fourier
mode), \(N_{\rm hi}\) is a gradient if and only if

\[
 \boxed{\quad
 \nabla\times F=0,\qquad
 K\times(F-c\nabla P)=0.
 \quad}
\tag{2.1}
\]

Indeed,

\[
 \nabla\times N_{\rm hi}
 =e^{iK\cdot x}
 \{\nabla\times F+iK\times(F-c\nabla P)\}.
\tag{2.2}
\]

Both braces in (2.1) must vanish separately because the envelopes are
real. Thus the scalar condition

\[
 A(a\cdot\nabla B)+B(b\cdot\nabla A)=0
\tag{2.3}
\]

is not sufficient. It only removes one component of \(F\); it does not
account for the transverse derivative of the variable coefficient in the
nominal pressure term \(icKP\).

For complex envelopes the necessary-and-sufficient condition is only the
single complex equation obtained by setting the right side of (2.2) to
zero; its real and imaginary parts need not split as in (2.1). Variable
envelope phases are better treated as changes of the carrier phases and
polarizations, and are outside this constant-\((k,\ell,a,b)\) ledger.
On a torus, curl-free also permits a harmonic constant vector. The
equivalence with "pressure" therefore assumes that the envelope spectrum
does not contain the cancelling frequency \(-K\), as is automatic under
the intended low-envelope/high-carrier separation.

## 3. Exact simplification for divergence-free scalar packets

Suppose each enveloped carrier is itself exactly divergence-free:

\[
 a\cdot\nabla A=0,\qquad b\cdot\nabla B=0.
\tag{3.1}
\]

Writing subscripts for directional derivatives in the orthonormal frame,
one obtains

\[
 P_e=\frac{\delta}{W}(BA_t-AB_t),
\qquad
 F=2W^2P_e\,e-2\delta^2P_t\,t.
\tag{3.2}
\]

Subtract the explicit pressure \(cP e^{iK\cdot x}\). The exact
remaining source is

\[
 N_{\rm hi}-\nabla(cP e^{iK\cdot x})
 =e^{iK\cdot x}R,
\tag{3.3}
\]

\[
 \boxed{\quad
 R=2(W^2+\delta^2)P_e\,e+2\delta^2P_v\,v.
 \quad}
\tag{3.4}
\]

In particular, the Leray-projected high sum is exactly

\[
 \mathbb P N_{\rm hi}=\mathbb P(e^{iK\cdot x}R).
\tag{3.5}
\]

Since \(R\) is real and perpendicular to \(K\), (2.1) now gives the
sharp criterion

\[
 \boxed{\quad \mathbb P N_{\rm hi}=0
 \quad\Longleftrightarrow\quad
 P_e=P_v=0.\quad}
\tag{3.6}
\]

Thus the carrier product may depend only on the fast direction \(t\).
There is no nonzero compact lateral product with an exactly pressure-only
common high sum.

The failure is especially transparent for equal envelopes. Conditions
(3.1) for \(A=B=f\) imply \(f_e=f_t=0\), so \(f=f(v)\). Formula
(3.4) becomes

\[
 R=2\delta^2(f^2)_v v.
\tag{3.7}
\]

It vanishes only for a constant envelope. At \(H=-1\), the old
"envelope-dark" equation is just \(w\cdot\nabla f=0\); (3.7) is the
missing lateral pressure-gradient charge.

For example, in the \(H=-1\) frame take \(e=n\), \(v=m\), and
\(A=B=f(m\cdot x)\). Both carrier packets are divergence-free and (2.3)
holds, but

\[
 N_{\rm hi}
 =-i\frac{Q^2}{\Lambda}f^2\,t\,e^{i2\Lambda t\cdot x},
\qquad
 \nabla\times N_{\rm hi}
 =-i\frac{Q^2}{\Lambda}(f^2)'(m\times t)
 e^{i2\Lambda t\cdot x}.
\tag{3.8}
\]

This is a direct counterexample to sufficiency of (2.3).

For a pure lateral Fourier modulation
\(f^2(v)=\widehat P(\xi)e^{i\xi v}\), the exact Leray multiplier is

\[
 \left|\mathbb P_{2\Lambda t+\xi v}
 \left[-i\frac{Q^2}{\Lambda}t\,\widehat P(\xi)\right]\right|
 =\frac{Q^2}{\Lambda}
 \frac{|\xi|}{\sqrt{4\Lambda^2+\xi^2}}
 |\widehat P(\xi)|.
\tag{3.9}
\]

For \(|\xi|\ll\Lambda\), this is
\(Q^2|\xi|/(2\Lambda^2)\) to leading order. A lateral buffer of
width \(R_v\) therefore has absolute high-sum size

\[
 O\!\left(\frac{A_0^2Q^2}{\Lambda^2R_v}\right),
\tag{3.10}
\]

not zero. Relative to an \(O(A_0^2Q)\) low beat, it is
\(O(Q/(\Lambda^2R_v))\), hence at most
\(O((\Lambda R_v)^{-1})\) when \(Q\leq\Lambda\).

## 4. Classification of exact scalar-envelope solutions

Set \(\alpha=\delta/W\), \(x=e\cdot x\), and \(z=t\cdot x\).
The general local solutions of (3.1) are

\[
 A=A_0(s_+,v),\qquad B=B_0(s_-,v),
 \qquad s_+=z+\alpha x,\quad s_-=z-\alpha x.
\tag{4.1}
\]

If their high sum is exactly pressure and \(AB\ne0\) on a connected
open set, (3.6) and separation of the independent variables imply

\[
 A=C(v)e^{\kappa s_+},\qquad
 B=D(v)e^{\kappa s_-},\qquad C(v)D(v)=C_0.
\tag{4.2}
\]

Consequently \(AB=C_0e^{2\kappa z}\). Bounded or compactly supported
real envelopes force the nontrivial exact case to collapse: a nonzero
lateral cutoff is impossible, and a bounded two-sided longitudinal
profile has \(\kappa=0\). Reciprocal factors \(C,D\) do not localize
the product and make one carrier large where the other is small.

## 5. A surviving characteristic-crossing core--buffer cell

Although exact pressure and compact interaction are incompatible, there
is an exact divergence-free one-colour cell with a completely explicit
buffer charge. Let \(\chi,\eta\) be smooth cutoffs and set

\[
 A=\chi(s_+/L)\eta(v/R_v),\qquad
 B=\chi(s_-/L)\eta(v/R_v).
\tag{5.1}
\]

Each carrier is divergence-free and each self-advection vanishes exactly.
The product is localized to the intersection of two characteristic slabs.
Its derivatives are

\[
\begin{aligned}
 P_e&=\frac{\alpha}{L}\eta^2
   (\chi_+'\chi_- -\chi_+\chi_-'),\\
 P_t&=\frac1L\eta^2
   (\chi_+'\chi_- +\chi_+\chi_-'),\\
 P_v&=\frac2{R_v}\chi_+\chi_-\eta\eta'.
\end{aligned}
\tag{5.2}
\]

Equations (3.4)--(3.5) are an exact high-sum error formula. In particular,
using the \(L^2\) contraction of the Leray projector,

\[
 \|\mathbb P N_{\rm hi}\|_2
 \leq
 2(W^2+\delta^2)\|P_e\|_2
 +2\delta^2\|P_v\|_2.
\tag{5.3}
\]

Pointwise, relative to the \(O(Q)\) desired low interaction, the two
buffer scales are

\[
 \boxed{\quad
 \text{oblique endcaps: }O((\Lambda L)^{-1}),
 \qquad
 \text{common lateral collar: }
 O\!\left(\frac{Q}{\Lambda^2R_v}\right).
 \quad}
\tag{5.4}
\]

The low-difference interaction is also exact:

\[
 N_{\rm lo}
 =e^{i(k-\ell)\cdot x}
 \{F+i\,2QWP\,e\}.
\tag{5.5}
\]

Thus \(\mathbb P(e^{i(k-\ell)\cdot x}i2QWP e)\) is the localized
target force, while \(\mathbb P(e^{i(k-\ell)\cdot x}F)\) is an
explicit endcap error. From (3.2) and (5.2), its size relative to the
raw low interaction is again \(O((\Lambda L)^{-1})\), with smaller
\(O(Q/(\Lambda^2L))\) terms.

There is a geometric cost hidden in the first estimate. The intersection
of \(|s_+|\lesssim L\) and \(|s_-|\lesssim L\) has extent

\[
 |x|\lesssim \frac{L}{\alpha}
 =\frac{2\Lambda W L}{Q}.
\tag{5.6}
\]

If this physical extent is prescribed to be \(R_e\), then
\(L\asymp\alpha R_e\), and the first relative error in (5.4) becomes

\[
 O((Q R_e)^{-1}),
\tag{5.7}
\]

with no \(\Lambda\)-gain. The apparent
\(O((\Lambda L)^{-1})\) gain therefore buys a tube elongated by
\(\Lambda/Q\). This tradeoff cannot be removed within two real scalar
envelopes because the coefficient of \(P_e\) in (3.4) is strictly
positive.

## 6. A compact-tube variant using the exact Leray projection

There is a second useful geometry if exact scalar-envelope
incompressibility is relaxed before applying the Leray projector. For
equal envelopes, direct calculation from (1.4) gives

\[
 A=B=f,\qquad P=f^2,\qquad
 F=W^2P_e\,e-\delta^2P_t\,t.
\tag{6.1}
\]

Equations (2.1) then show that the raw common high sum is exactly pressure
if and only if \(P_e=P_v=0\). In particular, every real longitudinal
profile \(f=f(t)\), including a compact bump, has a pressure-only raw high
sum. It is not an admissible velocity packet by itself:

\[
 \nabla\cdot(fa e^{ik\cdot x})=-\delta f_t e^{ik\cdot x},
\qquad
 \nabla\cdot(fb e^{i\ell\cdot x})=+\delta f_t e^{i\ell\cdot x}.
\tag{6.2}
\]

This defect is small in the high-carrier regime and can be removed exactly
by projecting each packet. For an envelope Fourier mode \(\xi t\),

\[
 \left|(\mathbb P_{k+\xi t}-I)a\right|
 =\frac{\delta|\xi|}{|k+\xi t|},
\qquad
 \left|(\mathbb P_{\ell+\xi t}-I)b\right|
 =\frac{\delta|\xi|}{|\ell+\xi t|}.
\tag{6.3}
\]

Hence, for \(|\xi|\lesssim R^{-1}\ll\Lambda\), the solenoidal correction
has relative size

\[
 O\!\left(\frac{\delta}{\Lambda R}\right)
 =O\!\left(\frac{Q}{\Lambda^2R}\right).
\tag{6.4}
\]

The frequency split behind this statement is elementary and does not
require a formal symbol expansion. If \(c_+\) is defined by

\[
 \mathbb P(fa e^{ik\cdot x})
 =e^{ik\cdot x}(fa+c_+),
\]

then for every \(s,N\geq0\),

\[
\begin{aligned}
 \|c_+\|_{H^s}
 \lesssim{}&
 \frac{W}{\Lambda}\|\partial_e f\|_{H^s}
 +\frac{\delta}{\Lambda}\|\partial_t f\|_{H^s}\\
 &+C_N\Lambda^{-N}\|f\|_{H^{s+N}},
\end{aligned}
\tag{6.5}
\]

and similarly for \(c_-\). On \(|\xi|\leq\Lambda/2\), this follows from

\[
 |(\mathbb P_{k+\xi}-I)a|
 =\frac{|a\cdot\xi|}{|k+\xi|}
 \leq\frac{2}{\Lambda}
 (W|\xi_e|+\delta|\xi_t|).
\tag{6.6}
\]

The last term in (6.5) controls \(|\xi|>\Lambda/2\). For a Gevrey
cutoff it can be made super-polynomially small, although uniform constants
at derivative order \(M\) still have to be entered in the cascade ledger.

One derivative falling on a carrier in the quadratic interaction loses
one factor \(\Lambda\). The resulting high, low, and self-interaction
errors are therefore of absolute size

\[
 O\!\left(\frac{A_0^2Q}{\Lambda R}\right),
\tag{6.7}
\]

or \(O((\Lambda R)^{-1})\) relative to the desired
\(O(A_0^2Q)\) low beat. Formula (6.7) is a band-limited product estimate;
a compact Gevrey cutoff requires the corresponding pseudodifferential
estimate plus control of its high-frequency tail.

This suggests a genuinely compact one-colour core--buffer ansatz:

\[
 f(x)=g((e\cdot x)/L_e)\,
      \rho((t\cdot x)/R_t,(v\cdot x)/R_v),
\qquad
 u_\pm=\mathbb P(fa_\pm e^{ik_\pm\cdot x}).
\tag{6.8}
\]

Take \(g=1\) on a long axial core and put its transition in remote
endcaps. In the core, the compact cross-section \(\rho\) costs
\(O((\Lambda R_t)^{-1})\) through (6.2)--(6.7), while its \(v\)-cutoff
has the smaller raw high charge
\(O(Q/(\Lambda^2R_v))\). The axial endcaps have
\(P_e=O(L_e^{-1})\), so (6.1)--(2.1) give the unavoidable relative cost

\[
 O((Q L_e)^{-1}).
\tag{6.9}
\]

Unlike (5.1), (6.8) is finite-energy on \(\mathbb R^3\) and exactly
divergence-free by definition. What is not yet proved is a weighted
pseudodifferential estimate showing that its nonlocal solenoidal
correction and the endcap pressure tail obey (6.7)--(6.9) in the active
core through the required Gevrey orders.

## 7. What survives and what remains open

The strongest exact algebraic building block is the
characteristic-crossing cell (5.1):

* it is exactly divergence-free;
* each carrier is an exact shear with zero self-advection;
* its desired low beat and every high/endcap charge are given by exact
  formulas;
* the common lateral collar is suppressed by
  \(Q^2/(\Lambda^2R_v)\); and
* the more dangerous \(P_e\) charge can be moved into distant oblique
  endcaps, at the cost of a tube of length \(\Lambda L/Q\).

This is an interaction-localized cell, not yet a compact finite-energy
cell on \(\mathbb R^3\): each individual characteristic slab is
unbounded along its own velocity direction. It is directly meaningful
on a compatible periodic cell. A Euclidean construction still needs
remote closure plus a solenoidal correction. On a fixed torus, the
characteristic slopes must also satisfy the lattice periodicity and
packing constraints.

The compact projected tube (6.8) is the stronger PDE candidate: it removes
the infinite-slab defect and preserves the favorable
\(O((\Lambda R_t)^{-1})\) cross-sectional scale. Its estimates are not yet
an exact transition lemma because the global Leray, pressure, heat, and
Gevrey bounds remain to be supplied.

Three further caveats are load-bearing.

1. The relative estimates use a nondegenerate projected low beat, as at
   \(H=-1\). They are not uniform as the low polarization becomes pressure.
2. Vanishing of \(R\) in the core does not imply vanishing of
   \(\mathbb P(e^{iK\cdot x}R)\) there. The Leray pressure tail from the
   collar/endcaps must be estimated globally, including periodic images.
3. Viscous evolution does not preserve arbitrary characteristic cutoffs.
   The heat error and the remote-endcap part of the compact Leray
   correction have not been controlled by a transition estimate.

The next load-bearing estimate is therefore not another pointwise
cancellation identity. It is an all-order weighted
pseudodifferential/pressure/heat estimate for (6.8), showing that the
cross-sectional corrections (6.5)--(6.7) and the remote endcap source do
not erase the core BAS gain or the preseeded child. Formulas (3.4), (5.5),
and (6.5) give the exact sources and symbol bound for that estimate.

## 8. Mechanical verification

Run python checks/characteristic_envelope_pressure.py. The checker verifies
(1.3), (3.2)--(3.4), the characteristic invariants, the exact lateral
Fourier multiplier, the \(L\)-versus-physical-length tradeoff, and the
exact per-mode Leray correction in the compact projected-tube variant.
