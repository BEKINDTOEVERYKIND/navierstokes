# C140--C141: exact passive-scalar closure and cubic active return of the paired \(A_2\) gate

**Date:** 2026-08-05  
**Status:** exact fixed-plane Navier--Stokes reduction and exact Fourier
return coefficients; independently re-derived twice; localized
conversion/focus remains open  
**Checker:**
[checks/paired_gate_passive_return_c140_c141.py](../checks/paired_gate_passive_return_c140_c141.py)

## 0. Claim boundary

This note stays on the C114--C136 one-cell geometry.  It asks the precise
question left by C117 and C134: when the fixed paired gate creates an
unavoidable \(O(b^2)\) complement mode, does that mode contaminate the
active \(A_2\) coordinates at the same order, or only after another
interaction?

For the bare real four-mode gate the answer is exact:

\[
 \text{child}=O(b),\qquad
 \text{complement wake}=O(b^2),\qquad
 \text{wake-fed active return}=O(b^3).
 \tag{0.1}
\]

The cubic coefficient is nonzero.  At the scheduled
\(b_j=(j+1)^{-2}\), (0.1) matches the raw exponent pair reserved in BAFL:
\(O((j+1)^{-4})\) in the retained wake and
\(O((j+1)^{-6})\) before the active chart.

This statement refers only to the two **raw short-time exponents**.  It
does not establish the backward-weighted BAFL norm, which is measured
along a localized conversion/focus trajectory that has not been
constructed.

There is also a structural explanation.  The entire fixed \(A_2\) root
lattice lies in the plane perpendicular to \(N=(1,1,1)\).  The C116 gate
is an exact 2D3C Navier--Stokes solution: its planar velocity is an
autonomous heat-decaying two-dimensional flow, while its \(N\)-component
is a passive scalar.  Thus the raw gate does not feed its wake back into
the planar pump.

This positive statement is deliberately narrow.  Smooth periodic data in
the fixed-plane 2D3C class are globally regular; within the C116 split the
passive scalar cannot deplete the autonomous planar pump or supply the
scheduled three-dimensional volume focus.  This does not refute the
abstract C118/C119 normal form or its intended off-plane leaves
\(q_m\pm Kr_i\) when \(m\ne0\), which are outside the fixed root lattice.
The bare gate may still be a module before localization, an off-plane
leaf, or a rotated conversion exits this invariant class.  Controlling the
\(O(b^2)\) passive wake during that exit is one necessary conversion/focus
part of BAFL.

## 1. Conventions

Use the C114 roots

\[
 k_1=(1,-1,0),\qquad k_2=(0,1,-1),\qquad
 k_3=(-1,0,1),\qquad N=(1,1,1),
 \tag{1.1}
\]

and put

\[
 k_c=k_1+k_2=-k_3,\qquad
 e_1=k_1+k_c=(2,-1,-1),\qquad
 e_2=k_2+k_c=(1,1,-2).
 \tag{1.2}
\]

For Fourier dictionaries \(a,b\), let

\[
 {\cal B}(a,b)_k
 =-iP_k\sum_{p+q=k}(a_p\mathbin\cdot q)b_q.
 \tag{1.3}
\]

This is the ordered Navier--Stokes bilinear term in the sign convention of
the C114--C117 checker.  The unnormalised helices are

\[
 H_s(k)=N+is\,{k\times N\over\sqrt2}.
 \tag{1.4}
\]

The explicit C116 coefficients are

\[
 A_1=iH_+(k_1)+H_-(k_1),\qquad
 A_2=-iH_+(k_2)+H_-(k_2),
 \tag{1.5}
\]

with conjugate coefficients at \(-k_1,-k_2\).

Let

\[
 {\cal H}=\{\pm k_1,\pm k_2,\pm k_c\},
 \tag{1.6}
\]

and denote the fixed Fourier projections onto \({\cal H}\) and its
complement by \(P\) and \(Q\).

## 2. C140: the paired gate is exactly 2D3C

Every wavevector in

\[
 \Lambda_{A_2}=\mathbb Zk_1+\mathbb Zk_2
 \tag{2.1}
\]

is perpendicular to \(N\).  A real divergence-free field supported in
this lattice is independent of the \(N\)-coordinate and splits uniquely as

\[
 u(x,t)=v(x,t)+\phi(x,t)N,\qquad v\mathbin\cdot N=0.
 \tag{2.2}
\]

Substitution into unforced Navier--Stokes gives the exact 2D3C system

\[
 \begin{aligned}
  \partial_t v+\mathbb P_{2D}(v\mathbin\cdot\nabla v)
      &=\nu\Delta v,\\
  \partial_t\phi+v\mathbin\cdot\nabla\phi
      &=\nu\Delta\phi.
 \end{aligned}
 \tag{2.3}
\]

Indeed, \(N\mathbin\cdot\nabla=0\), so neither
\(\phi N\mathbin\cdot\nabla v\) nor
\(\phi N\mathbin\cdot\nabla(\phi N)\) contributes, while
\((v\mathbin\cdot\nabla\phi)N\) is purely \(N\)-directed.  In particular,
\(\phi\) has no feedback into \(v\).

For a general root coefficient

\[
 A_i=x_iH_+(k_i)+y_iH_-(k_i),
 \tag{2.4}
\]

the split is visible coefficientwise:

\[
 A_i=(x_i+y_i)N
      +i(x_i-y_i){k_i\times N\over\sqrt2}.
 \tag{2.5}
\]

For (1.5),

\[
 \begin{aligned}
 \widehat\phi_0(k_1)&=1+i,&
 \widehat v_0(k_1)&=-(1+i){k_1\times N\over\sqrt2},\\
 \widehat\phi_0(k_2)&=1-i,&
 \widehat v_0(k_2)&=(1-i){k_2\times N\over\sqrt2}.
 \end{aligned}
 \tag{2.6}
\]

The planar field \(v_0\) is divergence-free and supported on the single
shell \(|k|^2=2\).  Choose the orientation of the planar quarter-turn
\(J\) so that
\(\operatorname{curl}_{2D}(J\nabla\psi)=-\Delta\psi\).  Writing
\(v_0=J\nabla\psi\) then gives

\[
 -\Delta\psi=2\psi,\qquad
 \omega_0=2\psi,\qquad
 v_0\mathbin\cdot\nabla\omega_0
 =2J\nabla\psi\mathbin\cdot\nabla\psi=0.
 \tag{2.7}
\]

Hence its projected Euler nonlinearity vanishes and

\[
 \boxed{v(t)=e^{-2\nu t}v_0}
 \tag{2.8}
\]

is exact.  The full C116 evolution is therefore reduced to the linear
passive-scalar equation

\[
 \boxed{
 \partial_t\phi
 +e^{-2\nu t}v_0\mathbin\cdot\nabla\phi
 =\nu\Delta\phi.}
 \tag{2.9}
\]

At pump frequency \(K\), \(2\nu\) in (2.8)--(2.9) becomes
\(2\nu K^2\).  In the Euler case,

\[
 \phi(t)=e^{-t\,v_0\cdot\nabla}\phi_0.
 \tag{2.10}
\]

For positive viscosity, \(\Delta\) need not commute with
\(v_0\cdot\nabla\), so (2.10) is not asserted; (2.9) is the exact
nonautonomous equation.

This proves more than a finite Fourier cancellation but less than a stage
map.  The passive scalar generally develops an infinite planar Fourier
ladder.  It is \(L^2\)-nonexpansive when \(\nu\ge0\), but may mix in stronger
norms.  Moreover, smooth periodic data in the whole fixed lattice (2.1)
remain 2D3C and are globally regular.  The passive component cannot
deplete the autonomous planar pump, and the \(N\)-independent circuit
cannot by itself supply the scheduled three-dimensional volume focus.
The off-plane C118/C119 leaves \(q_m\pm Kr_i\), \(m\ne0\), are not part
of (2.1).

## 3. The full-versus-retained Navier--Stokes jet

Let \(u\) solve full unforced Navier--Stokes from (1.5), and let
\(\bar u\) solve the finite retained equation

\[
 \partial_t\bar u=\nu\Delta\bar u
                  +P{\cal B}(\bar u,\bar u),
 \qquad \bar u(0)=u(0).
 \tag{3.1}
\]

Write

\[
 u(t)=\sum_{m\ge0}c_mt^m,\qquad
 \bar u(t)=\sum_{m\ge0}\bar c_mt^m.
 \tag{3.2}
\]

The exact recurrence is

\[
 (m+1)c_{m+1}
 =\nu\Delta c_m+\sum_{a+b=m}{\cal B}(c_a,c_b),
 \tag{3.3}
\]

with \(P\) applied on the right for \(\bar c_{m+1}\).  Because
\(\Delta\) is Fourier diagonal, \(P\Delta Q=0\).  The first C116
interaction lies entirely in \(P\), so

\[
 Qc_0=Qc_1=0,\qquad Pc_m=\bar c_m\quad(0\le m\le2).
 \tag{3.4}
\]

The first complement Taylor coefficient is independent of viscosity:

\[
 \begin{aligned}
 Qc_2(e_1)&=(-9-9i)N,\\
 Qc_2(e_2)&=(-9+9i)N,
 \end{aligned}
 \tag{3.5}
\]

with conjugates at \(-e_1,-e_2\).  At the next order,

\[
 3P(c_3-\bar c_3)
 =P\left\{
 {\cal B}(c_0,Qc_2)+{\cal B}(Qc_2,c_0)
 \right\}.
 \tag{3.6}
\]

Exact convolution gives

\[
 \begin{aligned}
 P(c_3-\bar c_3)(k_c)&=18\sqrt2\,iN,\\
 P(c_3-\bar c_3)(-k_c)&=-18\sqrt2\,iN.
 \end{aligned}
 \tag{3.7}
\]

Thus

\[
 Qu(t)=t^2Qc_2+O(t^3),\qquad
 Pu(t)-\bar u(t)
 =t^3P(c_3-\bar c_3)+O(t^4).
 \tag{3.8}
\]

There are nonzero \(O(t^2)\) coefficients at the parent wavevectors too,
but they occur identically in the retained equation (3.1).  They are
passive-scalar redistribution, not leakage and not depletion of the
autonomous planar pump.  The cubic statement concerns the difference
\(Pu-\bar u\).

Heat preserves wavevector, and the recurrence (3.3) proves that the two
displayed leading coefficients are independent of viscosity.  No
all-order viscosity-independent Dyson coefficient is asserted.  The
fixed-projector conclusion does not extend automatically to moving
projectors or localized envelopes.

## 4. C141: phase tuning cannot cancel the cubic return

The cubic return is not special to (1.5).  Take the general pair (2.4),
impose Fourier reality at the negative roots, and require exact
first-difference cancellation:

\[
 D=x_1\overline{y_2}-y_1\overline{x_2}=0.
 \tag{4.1}
\]

The terminal coefficient is

\[
 T=y_1x_2-x_1y_2.
 \tag{4.2}
\]

If \(T\ne0\), then \(y_1y_2\ne0\).  Put \(r_i=x_i/y_i\).
Equation (4.1) gives

\[
 r_1=\overline{r_2},\qquad
 r_1=a+ic,\qquad r_2=a-ic,
 \tag{4.3}
\]

where \(a,c\) are real, and

\[
 T=-2ic\,y_1y_2.
 \tag{4.4}
\]

Thus terminal activity is equivalent to \(c\ne0\).

Let \(D_m=\partial_t^mu(0)\) for the Euler gate, and let \(W_2\) be
exactly the \(\{\pm e_1,\pm e_2\}\) portion of \(D_2\).  Put
\(d_i=x_i-y_i\).  Direct use of

\[
 (k_1\times N)\mathbin\cdot k_2=-3,\qquad
 (k_2\times N)\mathbin\cdot k_1=3
 \tag{4.5}
\]

gives

\[
 \begin{aligned}
 D_1(k_c)&=3\sqrt2\,T N,\\
 W_2(e_1)&=-9d_1T N,\\
 W_2(e_2)&=+9d_2T N.
 \end{aligned}
 \tag{4.6}
\]

The two reality paths \(e_1+(-k_1)=k_c\) and
\(e_2+(-k_2)=k_c\) add rather than cancel.  Their exact sum is

\[
 \boxed{
 D_{3,\mathrm{wake}}(k_c)
 =-{27\over\sqrt2}
   \left(|d_1|^2+|d_2|^2\right)T N.}
 \tag{4.7}
\]

Equivalently, using (4.3)--(4.4),

\[
 \boxed{
 D_{3,\mathrm{wake}}(k_c)
 =27\sqrt2\,i\,c\,y_1y_2
   \left(|y_1|^2+|y_2|^2\right)
   \left((a-1)^2+c^2\right)N.}
 \tag{4.8}
\]

Every factor on the right is nonzero when \(T\ne0\).  Therefore no choice
inside the four-mode, first-difference-cancelled family can remove the
cubic wake return while retaining the child.

The numerical constants in (4.6)--(4.8) use the \(K=1\), unnormalised
helicity, and bilinear-sign conventions of Section 1.  Replacing every
wavevector by \(Kk\) multiplies \(D_1,W_2,D_{3,\mathrm{wake}}\) by
\(K,K^2,K^3\), respectively.  Unit-normalised helices rescale the
amplitude prefactors, and reversing the evolution sign reverses the
corresponding odd-interaction sign.  The support ordering and
nonvanishing conclusion are invariant under these convention changes.

This is pointwise nonvanishing, not a uniform lower bound over a
degenerating polarization family.  For example, with
\(y_1=y_2=1,a=1,c=\varepsilon\), both the child and the ratio of return to
child degenerate as \(\varepsilon\to0\).  A uniform stage estimate needs a
fixed gate such as (1.5), or a compact normalized nondegenerate
polarization chart.

The checker verifies (3.3)--(3.7) at three exact viscosities, checks
(4.6)--(4.8) for general exact complex samples, and gives a rigorous
25-point interpolation certificate for (4.8) on the balanced
\(y_1=y_2=1\) family.  All arithmetic is in
\(\mathbb Q(\sqrt2)+i\mathbb Q(\sqrt2)\).

## 5. Consequence for the factorial stage

For the fixed nondegenerate C116 gate, use dimensionless interaction time
\(\tau\).  Equations (3.5)--(3.8) give

\[
 \|\text{child}\|\asymp\tau,\qquad
 \|\text{wake}\|\asymp\tau^2,\qquad
 \|\text{wake-fed active return}\|\asymp\tau^3
 \tag{5.1}
\]

for sufficiently small \(\tau\).  Stopping when the child reaches
\(b_j=n^{-2}\) therefore gives

\[
 \|\text{wake}\|=\Theta(n^{-4}),\qquad
 \|\text{raw active return}\|=\Theta(n^{-6}).
 \tag{5.2}
\]

The \(\Theta\)-constants are uniform for the normalized schedule: C127 has
\(\mu_j=\nu(j!)^{-2}\to0\), hence \(\mu_j\) stays in a compact interval,
and the fixed smooth 2D3C gate has uniform local \(C^4\)-in-time bounds on
that interval.  Together with the viscosity-independent nonzero leading
coefficients, Taylor's theorem gives a common short-time window.  This
uniformity concerns only the homogeneous normalized gate; it says nothing
about localization or conversion remainders.

The earlier worst-case reading of C134, in which the entire
\(O(n^{-4})\) wake immediately passes through an \(O(n^2)\) active
inverse, is too pessimistic for the homogeneous fixed projector.  The raw
gate therefore has exactly the two short-time powers reserved in the BAFL
ledger.  This is not yet a backward-weighted estimate.

One newly isolated problem appears when the construction exits (2.3).  If
\(P_j(t)\) is a moving active projector and \(w_j=Q_j(t)u_j\), then even a
linearized equation \(\partial_tu=A_j(t)u+\cdots\), with \(A_j\) the exact
Leray-projected linearization about the candidate trajectory, contains the
direct conversion term

\[
 \partial_t(P_ju)\supset
 \left(\dot P_jQ_j+P_jA_jQ_j\right)w_j.
 \tag{5.3}
\]

For fixed Fourier \(P\) and diagonal heat, \(\dot P=0\) and the first
nonzero \(PAQw\) path is the extra interaction in (3.6).  Localization,
off-plane routing, and the active focus make (5.3) nontrivial and can
promote an \(O(n^{-4})\) wake directly into the active chart.

Accordingly, the homogeneous short-gate Taylor grading is closed, and it
isolates the following necessary **localized conversion-exposure
estimate**:

\[
\boxed{
 \left\|
 \int_{I_j^{\rm conv/foc}}
 {\cal U}^{\rm act}_j(t_j^+,s)
 \left(\dot P_jQ_j+P_jA_jQ_j\right)w_j(s)\,ds
 \right\|
 \le Cb_j^3=Cn^{-6},}
 \tag{LCE}
\]

together with the already required \(O(b_j^2)\) terminal retained-wake
bound.  Here \({\cal U}^{\rm act}_j\) includes the backward focus weight;
\(t_j^+\) is the stage exit time.  The definition is meaningful only
after a candidate localization and focus map have been specified.

No such map or estimate is constructed in this note.  LCE is a necessary
named subproblem of BAFL, not a theorem and not an equivalent replacement
for full BAFL.  Once \(P_j,A_j\), and the candidate trajectory are fixed,
it is intended to include the linearized localization commutator,
pressure-mediated conversion, off-plane activation, and exposure to the
\(q_j^{3/2}\) focus.  It does not by itself control newly generated
conversion leakage, quadratic wake--wake feedback, closure of the
entrance/wake state classes, or construction of the depletion, conversion,
and focus modules, nor any residual source outside the already-created
wake \(w_j\).  Proving LCE, or finding one explicit term in (5.3)
that necessarily exceeds \(Cn^{-6}\), is the next exact test exposed by
the paired-gate calculation; the complete one-cell target remains BAFL.

## 6. Claim ledger

* **C140 (EXACT, fixed-plane):** the C116 gate is an exact 2D3C unforced
  Navier--Stokes flow: its monochromatic planar component heat-decays and
  its \(N\)-component solves a linear passive-scalar equation.  Relative
  to the fixed six-root projection, the first complement wake is
  quadratic and its first active return is cubic, with the exact leading
  Taylor coefficients (3.5), (3.7), independent of viscosity.
* **C141 (EXACT, bare four-mode family):** every nondegenerate
  first-difference-cancelled paired gate has the nonzero cubic wake-return
  coefficient (4.7)--(4.8).  Extra modes and later propagated
  cancellations are not excluded.
* **LCE (OPEN):** preserve the extra interaction factor while converting
  the passive child into the localized three-dimensional active focus.
  This is a necessary conversion/focus component of the one-cell BAFL
  theorem; full BAFL has the additional obligations listed above, and no
  stage map or blow-up claim is made.
