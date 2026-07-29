# Time reversal of the CDP complete inverse cascade: a heat-clock no-range theorem for a flat terminal force

Date: 2026-07-29

## Decision

This note audits the following proposed route to Clay alternative (D):

1. time-reverse the complete inverse Navier--Stokes cascade of
   Cheskidov--Dai--Palasek (CDP);
2. replace their double-exponential, infinite-energy ledger by a
   finite-energy cascade with polynomial internal carrier; and
3. cancel the reversal defect by an anti-diffusive corrector, followed by
   an all-order Borel correction that makes the remaining force
   \(C^\infty\)-flat at the terminal time.

There is a rigorous **no-range obstruction to importing the CDP
transition and fixed-point theorem in this way**.

* The exact time reverse \(v(t)=-u(T-t)\) has projected force
  \[
  {\cal N}_\nu(v)=-2\nu\Delta v.
  \]
  The published CDP branch has divergent \(L^2\) energy at its branching
  time.  Its reverse therefore cannot be a classical preterminal path
  driven by a bounded, let alone flat, force.
* On a carrier of frequency \(\Lambda\), over a stage of length \(\tau\),
  put
  \[
  \Theta=\nu\Lambda^2\tau .
  \]
  The exact same-terminal-trace correction that changes a reversed
  heat factor into a forward heat factor has relative amplitude
  \(2\sinh\Theta\).  A perturbative correction requires
  \(\Theta\to0\).
* CDP Proposition 5.1 is proved in the opposite range:
  \[
  N_{\rm hi}^{-2}\ll t_k,
  \qquad\text{hence}\qquad
  N_{\rm hi}^2t_k\longrightarrow\infty .
  \]
  Their concrete choice costs
  \(\exp(cN_{\rm lo}^{\,2b-4})\) under backward heat evolution.
  Proposition 6.2 is a forward semigroup estimate and supplies no
  terminal right inverse for this multiplier.
* For the proposed finite-energy polynomial ledger,
  \(\Theta_j\to0\).  That is the one range in which a Gevrey/Borel
  endpoint correction could be perturbative.  It is also the range in
  which the heat-saturated CDP transfer and the scale window used to
  prove it disappear.
* Renormalizing the high-wave amplitude does not recover the published
  theorem.  At \(\Theta\ll1\) it turns the cell into an ordinary
  inertial-stress transition whose estimates must be reproved.  Forcing
  the exact CDP scale separation into the polynomial ledger makes the
  viscous loss of a single order-one handoff diverge.

Thus there is no exact reduction to CDP Proposition 6.3.  The only
remaining version of this route is a **new theorem**: an all-order,
small-heat-clock, inertial transition with tame homological inverses and
super-algebraically matched endpoint jets.  The CDP rank-one stress
geometry may inform such a theorem, but neither their principal
transition estimate nor their fixed point proves it.

This is an import obstruction, not a proof that every possible forward
Navier--Stokes cascade or every possible construction for alternative
(D) is impossible.  No singular solution is claimed here.

---

## 1. What the primary paper proves

The source audited here is A. Cheskidov, M. Dai, and S. Palasek,
[*Instantaneous Type I blow-up and non-uniqueness of smooth solutions of
the Navier--Stokes
equations*](https://arxiv.org/abs/2511.09556v2), arXiv:2511.09556v2
(12 January 2026).  CDP normalize the viscosity to one.  The following
features are decisive.

1. **Orientation.**  Theorem 1.1 starts from a smooth value at \(T_*\)
   and constructs a branch that is classical for \(t>T_*\), with
   Type-I growth as \(t\downarrow T_*+\).  Energy is injected from
   infinite wavenumber and then transferred from high to low frequency.
   This is an instantaneous inverse cascade, not a forward terminal
   cascade.

2. **Energy divergence is necessary, not optional.**  Remark 1.4 states
   that their branch belongs to neither \(L^\infty_tL^2_x\) nor
   \(L^2_tH^1_x\) near \(T_*\).  Theorem 2.4 proves more generally that
   any instantaneous Type-I branch from a smooth value, in their
   setting, must satisfy
   \[
   \lim_{t\downarrow T_*+}\|u(t)\|_2=\infty .
   \tag{1.1}
   \]

3. **Shells and heat factors.**  Section 4.1 uses
   \[
   N_{j,k}\simeq m_*A^{b^k}
   \tag{1.2}
   \]
   (with a minor within-shell variation in dimension two).  The
   approximate shell in Section 4 has the form
   \[
   \bar v_k(t)
   =
   -\sum_jN_{j,k}e^{-N_{j,k}^2t}\Delta\psi_{j,k},
   \tag{1.3}
   \]
   and the exact daughter is generated forward:
   \[
   v_k(t)
   =
   -\int_0^t e^{(t-s)\Delta}
   {\mathbb P}\operatorname{div}
   \big(\bar v_{k+1}\otimes\bar v_{k+1}\big)(s)\,ds .
   \tag{1.4}
   \]

4. **The proved transition is heat saturated.**  Proposition 5.1
   chooses a separation time with
   \[
   N_{1,k+1}^{-2}\ll t_k\ll N_{J_d,k}^{-3},
   \tag{1.5}
   \]
   concretely \(t_k=N_{J_d,k}^{-4}\).  Hence, writing
   \(N_{\rm hi}=N_{1,k+1}\) and \(N_{\rm lo}=N_{J_d,k}\),
   \[
   N_{\rm hi}^2t_k\to\infty,
   \qquad
   N_{\rm hi}^2\gg N_{\rm lo}^3 .
   \tag{1.6}
   \]
   The factor \(1-e^{-2N_{\rm hi}^2t_k}\) is consequently saturated.

5. **The residual is small in a singular weighted norm, not flat.**
   Proposition 5.3 gives
   \[
   \partial_tv-\Delta v+
   {\mathbb P}\operatorname{div}(v\otimes v)
   ={\mathbb P}\operatorname{div}f
   \tag{1.7}
   \]
   with, schematically,
   \[
   \|\nabla^nf(t)\|_{C^\kappa}
   \lesssim_n\epsilon_0
   \big(t^{-1-n/2+\alpha}+1\big).
   \tag{1.8}
   \]

6. **The corrector is forward and initial-value.**  In Section 6,
   \(w(0)=0\), and Proposition 6.2 estimates \(S(t,t')\) only for
   \(0<t'\leq t\).  Proposition 6.3 is the contraction
   \[
   {\cal F}(w)(t)
   =
   -\int_0^tS(t,t')
   \big(w\otimes w+f+2U^{1/N_0}\mathbin\odot v\big)(t')\,dt'.
   \tag{1.9}
   \]
   It is not a terminal-value theorem.

Those statements are enough to decide the proposed import.  No
unpublished property of the CDP construction is used below.

---

## 2. The flat-force target and its energy consequence

On \(\mathbb T^3\), write the projected Navier--Stokes operator as

\[
{\cal N}_\nu(u)
:=
\partial_tu-\nu\Delta u+
{\mathbb P}\operatorname{div}(u\otimes u).
\tag{2.1}
\]

For a preterminal path \(u\in
C^\infty(\mathbb T^3\times[0,T))\), a sufficient condition for extending
\(f={\cal N}_\nu(u)\) by zero for \(t\geq T\) is

\[
\|\partial_t^qf(t)\|_{C^m}
\leq C_{m,q,N}(T-t)^N
\quad
\text{for every }m,q,N\geq0 .
\tag{2.2}
\]

This is the relevant form of terminal flatness for the periodic
alternative (D) in Fefferman's official
[*Existence and Smoothness of the Navier--Stokes
Equation*](https://www.claymath.org/wp-content/uploads/2022/06/navierstokes.pdf).
The zero extension is smooth and has the time decay required there.

There is an elementary necessary condition that immediately excludes the
unmodified reverse of CDP.

### Lemma 2.1 (bounded force implies finite energy and dissipation)

Suppose \(u\) is a smooth solution on \([0,T)\) of
\({\cal N}_\nu(u)=f\), with \(f\in L^1(0,T;L^2)\).  Then

\[
\sup_{t<T}\|u(t)\|_2
\leq
\|u(0)\|_2+\int_0^T\|f(s)\|_2\,ds,
\tag{2.3}
\]

and

\[
\nu\int_0^T\|\nabla u(s)\|_2^2\,ds<\infty .
\tag{2.4}
\]

#### Proof

The energy identity is

\[
\frac12\frac d{dt}\|u\|_2^2+
\nu\|\nabla u\|_2^2
=\langle f,u\rangle .
\tag{2.5}
\]

Dropping dissipation and dividing by \(\|u\|_2\), by regularization at
zeros if necessary, gives
\[
\frac d{dt}\|u\|_2\leq\|f\|_2.
\]
This proves (2.3).  Integrating (2.5), using (2.3), and applying
\(\int|\langle f,u\rangle|
\leq\sup\|u\|_2\int\|f\|_2\), proves (2.4). \(\square\)

A \(C^\infty\)-flat force is bounded and integrable.  In contrast, the
CDP branch obeys (1.1).  Time reversal preserves the instantaneous
\(L^2\) norm.  Therefore:

> The published CDP solution cannot itself be time-reversed into a
> preterminal classical path with flat force.  Any successful correction
> must cancel or replace enough of its leading \(L^2\) content to restore
> (2.3); the unchanged reversed solution cannot be the final path.

---

## 3. Exact reversal identity

The sign needed to preserve the Euler quadratic term is important.

### Lemma 3.1 (Navier--Stokes time reversal)

Let \(u(s)\) solve

\[
\partial_su-\nu\Delta u+
{\mathbb P}\operatorname{div}(u\otimes u)=0
\tag{3.1}
\]

on \(0<s<\delta\).  Set

\[
s=T-t,\qquad v(t)=-u(s).
\tag{3.2}
\]

Then

\[
\boxed{
{\cal N}_\nu(v)
=2\nu\Delta u(T-t)
=-2\nu\Delta v(t).
}
\tag{3.3}
\]

#### Proof

Since \(\partial_tv=\partial_su\),
\(\Delta v=-\Delta u\), and the quadratic term is even under
\(u\mapsto-u\),

\[
{\cal N}_\nu(v)
=\partial_su+\nu\Delta u+
{\mathbb P}\operatorname{div}(u\otimes u)
=2\nu\Delta u.
\]
\(\square\)

If no corrector is added, even convergence of (3.3) to zero in \(C^0\)
is impossible for a reversed branch with unbounded \(L^\infty\) norm.
Indeed, on the torus,

\[
\|v-\bar v\|_\infty
\leq C\|\Delta v\|_\infty,
\tag{3.4}
\]

and the mean is constant for the original unforced solution.

For \(U=v+z\), the exact corrected residual is

\[
\begin{aligned}
{\cal N}_\nu(v+z)
= {}&-2\nu\Delta v
+(\partial_t-\nu\Delta)z\\
&+{\mathbb P}\operatorname{div}
\big(v\otimes z+z\otimes v+z\otimes z\big).
\end{aligned}
\tag{3.5}
\]

Thus a correction must solve a terminal, backward-parabolic problem
unless the nonlinear terms furnish a new forward-oriented principal
balance.

---

## 4. The exact one-carrier anti-diffusive cost

The obstruction already appears on a divergence-free shear Fourier
mode, for which the quadratic term vanishes identically.

Let \(e_\lambda\) satisfy

\[
-\Delta e_\lambda=\lambda^2e_\lambda,\qquad
\operatorname{div}e_\lambda=0,\qquad
{\mathbb P}\operatorname{div}
(e_\lambda\otimes e_\lambda)=0.
\tag{4.1}
\]

For example, one may take
\(e_\lambda(x)=a\sin(k\cdot x)\) with \(a\cdot k=0\) and
\(|k|=\lambda\).

The forward heat solution issuing from \(Ae_\lambda\) at \(s=0\) is

\[
u(s)=Ae^{-\nu\lambda^2s}e_\lambda .
\tag{4.2}
\]

Its reversed velocity is

\[
v(t)=-Ae^{-\nu\lambda^2(T-t)}e_\lambda .
\tag{4.3}
\]

### Lemma 4.1 (same-trace correction)

Among fields in the one-mode subspace, the unique unforced forward
solution \(V\) with the same terminal trace \(V(T)=v(T)=-Ae_\lambda\)
is, in the variable \(s=T-t\),

\[
V(s)=-Ae^{+\nu\lambda^2s}e_\lambda .
\tag{4.4}
\]

Consequently the exact correction \(z=V-v\), with \(z(T)=0\), is

\[
\boxed{
z(s)=-2A\sinh(\nu\lambda^2s)e_\lambda .
}
\tag{4.5}
\]

At the other end of a stage of length \(\tau\),

\[
\boxed{
\frac{\|z(\tau)\|}{|A|\,\|e_\lambda\|}
=2\sinh\Theta,
\qquad
\Theta:=\nu\lambda^2\tau .
}
\tag{4.6}
\]

In particular,

\[
2\sinh\Theta
=
\begin{cases}
2\Theta+O(\Theta^3),&\Theta\to0,\\
\asymp1,&\Theta\asymp1,\\
\asymp e^\Theta,&\Theta\to\infty.
\end{cases}
\tag{4.7}
\]

#### Proof

Writing \(V(t)=q(t)e_\lambda\), the forward heat equation gives
\(q_t+\nu\lambda^2q=0\).  Since \(s=T-t\), this is
\(q_s=\nu\lambda^2q\).  The terminal value \(q(0)=-A\) yields
(4.4), and subtracting (4.3) gives (4.5). \(\square\)

This lemma is not a universal lower bound against a corrector using
cross-frequency nonlinear interactions.  It is an exact calculation for
the **anti-diffusive modewise correction proposed in the import**.
Using nonlinear cross-frequency cancellation instead would require a new
principal transition theorem.

The same calculation rules out leaving a small, merely nonzero residual.
Let \(q_v(s)=-Ae^{-\mu s}\), where \(\mu=\nu\lambda^2\), and let
\(q=q_v+z\), with \(z(0)=0\).  Define the forward heat defect in the
\(s\)-variable by

\[
r(s)=-q_s(s)+\mu q(s).
\tag{4.8}
\]

Variation of constants gives the exact identity

\[
\int_0^\tau e^{-\mu s}r(s)\,ds
=
-A(1-e^{-2\Theta})-e^{-\Theta}z(\tau).
\tag{4.9}
\]

Consequently, if \(|z(\tau)|\leq\varepsilon|A|\),

\[
\boxed{
\|r\|_{L^\infty(0,\tau)}
\geq
\mu|A|\,
\frac{\big(1-e^{-2\Theta}-\varepsilon e^{-\Theta}\big)_+}
     {1-e^{-\Theta}} .
}
\tag{4.10}
\]

For \(\Theta\to\infty\) and fixed \(\varepsilon\) (in particular for a
perturbative correction), the right side is
\((1-o(1))\mu|A|\).  Thus a correction that remains perturbative at the
other end of a heat-saturated stage leaves a defect of the same size as
the raw reversed-viscosity term.  Setting \(r=0\) in (4.9) recovers the
exponentially large correction (4.5).

---

## 5. The heat-clock no-range theorem

The preceding calculation gives a short exact obstruction.

### Proposition 5.1 (no direct CDP reversal/fixed-point import)

There is no sequence of CDP stages on which both of the following hold:

1. the principal transition is used in the scale range proved in CDP
   Proposition 5.1; and
2. the reversed high-shell heat factors are repaired, with their terminal
   traces fixed, by a perturbative modewise corrector.

#### Proof

CDP require

\[
N_{{\rm hi},k}^{-2}\ll t_k.
\tag{5.1}
\]

Therefore

\[
\Theta_k=N_{{\rm hi},k}^2t_k\longrightarrow\infty.
\tag{5.2}
\]

By Lemma 4.1, the relative correction at the start of that interval is
\(2\sinh\Theta_k\to\infty\), not \(o(1)\).  Conversely, a perturbative
correction requires \(\Theta_k\to0\).  The two ranges are disjoint.
\(\square\)

Allowing a small nonzero corrected force does not create an overlap:
(4.10) shows that a perturbative correction in the range (5.2) retains
an order-\(\nu N_{{\rm hi},k}^2A_k\) carrier defect.

For CDP's concrete choice

\[
t_k=N_{{\rm lo},k}^{-4},
\qquad
N_{{\rm hi},k}\simeq N_{{\rm lo},k}^{\,b},
\tag{5.3}
\]

with \(b>2\), the terminal inverse heat cost is

\[
\exp\!\big(cN_{{\rm hi},k}^2t_k\big)
=
\exp\!\big(cN_{{\rm lo},k}^{\,2b-4}\big).
\tag{5.4}
\]

This is superpolynomial and is absent from every estimate in CDP; their
fixed point allows only the mild power loss appearing in Proposition
6.2.

This is also visible directly in the orientation of the CDP semigroup
estimate.  Proposition 6.2 controls \(S(t,t')\) for \(t'\leq t\), with a
mild factor \((t/t')^\epsilon\).  A terminal problem asks for
\(S(t,T)\) with \(t<T\).  Even before the drift is included, its carrier
multiplier is

\[
e^{+\nu\lambda^2(T-t)}.
\tag{5.5}
\]

The contraction (1.9) therefore cannot be read backward.  Replacing its
forward Duhamel integral by a terminal integral is a different theorem,
not a change of variables in Proposition 6.3.

---

## 6. The finite-energy polynomial ledger lies in the opposite range

Consider the proposed outer scale, amplitude, normalized carrier, and
turnover time

\[
\ell_j=r^{-j},\qquad
K_j=(j+1)^A,\qquad
a_j=\ell_j^{-\gamma}K_j^\gamma,
\qquad
1<\gamma<\frac32,
\tag{6.1}
\]

\[
\Lambda_j=\frac{K_j}{\ell_j},
\qquad
\tau_j=\frac{\ell_j}{a_j}
=\ell_j^{1+\gamma}K_j^{-\gamma}.
\tag{6.2}
\]

The energy of a three-dimensional packet occupying volume
\(\asymp\ell_j^3\) is

\[
E_j\asymp a_j^2\ell_j^3
=\ell_j^{3-2\gamma}K_j^{2\gamma}.
\tag{6.3}
\]

Because \(\gamma<3/2\), \(\sum_jE_j<\infty\).

The carrier heat clock is

\[
\boxed{
\Theta_j
=\nu\Lambda_j^2\tau_j
=\nu\ell_j^{\gamma-1}K_j^{2-\gamma}.
}
\tag{6.4}
\]

Since \(\gamma>1\), the exponential decay of
\(\ell_j^{\gamma-1}\) beats the polynomial factor in \(K_j\), and

\[
\Theta_j\to0.
\tag{6.5}
\]

The carrier dissipation over one stage has the favorable ledger

\[
E_j\Theta_j
\asymp
\nu\ell_j^{2-\gamma}K_j^{\gamma+2},
\qquad
\sum_jE_j\Theta_j<\infty .
\tag{6.6}
\]

Thus finite energy and finite carrier dissipation are mutually
compatible in this small-clock range.  What fails is the CDP import.

First, the normalized heat integral is

\[
\nu\Lambda_j^2
\int_0^{\tau_j}e^{-2\nu\Lambda_j^2s}\,ds
=\frac{1-e^{-2\Theta_j}}2
=\Theta_j+O(\Theta_j^2).
\tag{6.7}
\]

It is unsaturated and tends to zero.

Second, identify the low envelope and high carrier in the exact CDP
window as

\[
N_{\rm lo}=\ell_j^{-1},
\qquad
N_{\rm hi}=K_j\ell_j^{-1}.
\tag{6.8}
\]

The necessary nonemptiness condition in (1.6) becomes

\[
N_{\rm hi}^2\gg N_{\rm lo}^3
\quad\Longleftrightarrow\quad
K_j^2\gg\ell_j^{-1}.
\tag{6.9}
\]

But

\[
K_j^2\ell_j=(j+1)^{2A}r^{-j}\longrightarrow0.
\tag{6.10}
\]

This is the opposite inequality.  Therefore no choice of the fixed
polynomial exponent \(A\) places this ledger in the range of CDP
Proposition 5.1.

The other natural identification fails as well.  If adjacent physical
cascade carriers are used,

\[
N_{\rm lo}=\Lambda_j,\qquad N_{\rm hi}=\Lambda_{j+1},
\tag{6.11}
\]

then \(\Lambda_{j+1}/\Lambda_j=r(1+O(j^{-1}))\), and hence

\[
\frac{N_{\rm hi}^2}{N_{\rm lo}^3}
=
\frac{r^2(1+O(j^{-1}))}{\Lambda_j}
\longrightarrow0,
\tag{6.12}
\]

again opposite to (1.6).

---

## 7. Amplitude renormalization does not reopen the range

It is natural to try to compensate for the small factor in (6.7) by
changing the high-wave amplitude.  The following stage calculation
isolates what that does.

Let a parent have amplitude \(a\), scale \(\ell\), turnover time
\(\tau=\ell/a\), and energy \(E\asymp a^2\ell^3\).  Let a high wave of
amplitude \(b\) and frequency \(\lambda\) have the heat profile
\(e^{-\nu\lambda^2t}\).  Contracting its covariance with a parent strain
of size \(a/\ell\) gives the stage-transfer scale

\[
\begin{aligned}
\Delta E
&\asymp
\frac a\ell\,b^2\ell^3
\int_0^\tau e^{-2\nu\lambda^2t}\,dt\\
&=
b^2\ell^3\,
\frac{1-e^{-2\Theta}}{2\Theta},
\qquad
\Theta=\nu\lambda^2\tau.
\end{aligned}
\tag{7.1}
\]

Demanding \(\Delta E\asymp E\) forces

\[
\boxed{
\frac{b^2}{a^2}
\asymp
\frac{2\Theta}{1-e^{-2\Theta}}.
}
\tag{7.2}
\]

There are two regimes.

* If \(\Theta\ll1\), then \(b/a=1+O(\Theta)\).  The transfer in (7.1)
  comes from maintaining an inertial covariance over the full turnover
  time.  Heat decay is perturbative.  This may be a viable design
  principle, but it is not the heat-saturated CDP module, and the
  principal error estimates must be proved again.
* If \(\Theta\gg1\), then \(b/a\asymp\Theta^{1/2}\).  The high packet is
  energetically inflated.  Its integrated viscous loss is
  \[
  \begin{aligned}
  {\cal D}_{\rm hi}
  &\asymp
  \nu\lambda^2b^2\ell^3
  \int_0^\tau e^{-2\nu\lambda^2t}\,dt\\
  &\asymp E\Theta .
  \end{aligned}
  \tag{7.3}
  \]
  This is the high heat-cost regime.  In the polynomial ledger forced
  into the CDP window, that cost diverges, as shown below.

The same conclusion follows without assuming a constant-amplitude wave.
If an order-\(E\) transfer is accomplished against strain
\(\lesssim a/\ell\), then necessarily

\[
\int_{\rm stage}\int_{\rm cell}|W|^2\,dx\,dt
\gtrsim E\tau.
\tag{7.4}
\]

For spectral localization at \(\lambda\), Bernstein/Poincare scaling
then gives

\[
\nu\int_{\rm stage}\|\nabla W\|_2^2\,dt
\gtrsim E\,\nu\lambda^2\tau
=E\Theta .
\tag{7.5}
\]

Spatial or temporal intermittency does not remove (7.4): reducing the
occupied measure requires the covariance amplitude to increase by the
reciprocal factor.

Finally, suppose an auxiliary normalized carrier \(H_j\) is inserted so
that

\[
N_{\rm hi}=\frac{H_j}{\ell_j}
\tag{7.6}
\]

does satisfy the exact CDP separation.  From (6.9),

\[
H_j^2\gg\ell_j^{-1}.
\tag{7.7}
\]

Its heat clock obeys

\[
\Theta_j^{H}
=
\nu\left(\frac{H_j}{\ell_j}\right)^2\tau_j
\gg
\nu\ell_j^{\gamma-2}K_j^{-\gamma}.
\tag{7.8}
\]

Combining (6.3), (7.5), and (7.8),

\[
{\cal D}^{H}_j
\gtrsim E_j\Theta_j^H
\gg
\nu\ell_j^{1-\gamma}K_j^\gamma
\longrightarrow\infty.
\tag{7.9}
\]

Thus even one late order-one handoff has divergent viscous cost in this
specific polynomial ledger if the carrier is inflated enough to enter
the proved CDP window.

---

## 8. What Borel correction can and cannot do

Small heat clock is favorable for an all-order correction, but relative
smallness is not flatness.  Before correction, the reversed viscous
defect on the \(j\)-th physical carrier has size

\[
\|\nabla^m(-2\nu\Delta v_j)\|_\infty
\asymp
\nu a_j\Lambda_j^{m+2}
=
\nu\ell_j^{-(\gamma+m+2)}
K_j^{\gamma+m+2}.
\tag{8.1}
\]

This diverges for every fixed \(m\).  Relative to the inertial scale
\(a_j^2/\ell_j\), however,

\[
\frac{\nu a_j\Lambda_j^2}{a_j^2/\ell_j}
=\nu\Lambda_j^2\tau_j
=\Theta_j\to0.
\tag{8.2}
\]

Hence a successful construction would have to cancel an asymptotically
small **relative** defect through an increasing number of orders until
the remaining **physical** defect is beyond all algebraic orders.

There is a useful conditional ledger calculation.

### Lemma 8.1 (conditional diagonal Borel budget)

Assume that a nonlinear homological construction exists through order
\(M\), and that after taking \(m\) spatial and \(q\) temporal
derivatives its stage-\(j\) remainder has the form

\[
\|\partial_t^qR_{j,M}\|_{C_x^m}
\leq
P_{m,q}(\ell_j^{-1},K_j,a_j,\tau_j^{-1})\,
C^M M!^\sigma K_j^{dM}\Theta_j^{M+1},
\tag{8.3}
\]

where \(P_{m,q}\) is a fixed monomial for fixed \(m,q\), and
\(C,\sigma,d\) are independent of \(j,M\).  Then choosing
\(M_j=\lfloor cj\rfloor\), for sufficiently small fixed \(c>0\), gives

\[
\|\partial_t^qR_{j,M_j}\|_{C_x^m}
\leq C_{m,q}e^{-c_{m,q}j^2}.
\tag{8.4}
\]

If the stages accumulate at \(T\) with
\(T-t_j\asymp\tau_j\), the summed remainder is \(C^\infty\)-flat at
\(T\).

#### Proof

From (6.4),

\[
\log\Theta_j
=-(\gamma-1)(\log r)j+O(\log j).
\tag{8.5}
\]

With \(M_j\asymp j\), the last factor in (8.3) therefore contributes
\(-cj^2+O(j\log j)\) to the logarithm.  Stirling's formula and the
polynomial carrier give

\[
\log(C^{M_j}M_j!^\sigma K_j^{dM_j})
=O(j\log j).
\tag{8.6}
\]

Every fixed physical prefactor \(P_{m,q}\) contributes only \(O(j)\) to
the logarithm.  This proves (8.4).

Moreover,

\[
\log\tau_j=-(1+\gamma)(\log r)j+O(\log j).
\tag{8.7}
\]

Thus \(e^{-c j^2}=O(\tau_j^N)\) for every fixed \(N\).  The same is true
after summing the Gaussian tail over all later stages. \(\square\)

Lemma 8.1 is a **budget, not an existence theorem**.  Its hypothesis is
precisely the missing analytic result.  In particular:

1. Borel's theorem realizes a prescribed compatible jet.  It does not
   solve the nonlinear Navier--Stokes homological equations that
   determine that jet.
2. Resonances, pressure inversion, localization commutators, and
   derivative losses must obey a uniform Gevrey-type bound such as
   (8.3).  CDP prove estimates to each fixed derivative order after
   choosing scale separation; they do not prove a diagonal
   \(M_j\to\infty\) estimate in the small-clock range.
3. In the CDP range \(\Theta_j\to\infty\), powers of \(\Theta_j\) do not
   produce a decreasing asymptotic series.  Borel summing the exact
   one-mode series simply recovers the exponentially large correction
   (4.5).
4. A cutoff between the uncorrected reverse and a terminal Borel layer
   introduces a seam term of size \(\chi' z\).  Since those seams also
   accumulate at \(T\), their mismatch must already satisfy a bound like
   (8.4).  Moving the mismatch to a thinner layer does not make it flat.

Therefore polynomial carriers make an all-order correction
**arithmetically conceivable**, but CDP's fixed point does not supply the
homological construction required by Lemma 8.1.

---

## 9. Why Proposition 6.3 cannot be imported after changing the ledger

The failure is not merely that CDP happened to choose a different
frequency sequence.

Their contraction uses all of the following simultaneously:

\[
\|v(t)\|_\infty\lesssim t^{-1/2},
\qquad
\|f\|_Y\lesssim\epsilon_0,
\tag{9.1}
\]

the logarithmic exposure estimate in their equation (5.19), the forward
semigroup estimate of Proposition 6.2, and the zero initial condition
\(w(0)=0\).  The small ball is

\[
\|w\|_X
=
\sup_t\left(
t^{(1-\alpha)/2}\|w(t)\|_\infty
+
t^{(2-\alpha)/2}\|\nabla w(t)\|_{C^\kappa}
\right)
<\delta .
\tag{9.2}
\]

Changing to the polynomial ledger affects every input:

* (1.5) has no admissible scale interval by (6.10);
* the heat integral used to identify the daughter is smaller by
  \(\Theta_j\);
* amplitude compensation changes the velocity and product estimates;
* the desired corrector has terminal rather than initial data;
* all-order flatness requires estimates with differentiation order
  increasing with \(j\), not merely a small fixed \(Y\)-norm; and
* the reversed linear carrier demands (5.5), which is absent from
  Proposition 6.2.

One could reuse the algebraic rank-one identity

\[
\sum_j a_{j,k+1}^2\theta_j\otimes\theta_j
=
2{\cal D}\sum_jN_{j,k}\psi_{j,k}+pI
\tag{9.3}
\]

as a stress factorization in a new construction.  But once the
time-integrated stress, shell amplitudes, heat clock, endpoint condition,
and error norms are changed, invoking Proposition 6.3 is circular: its
hypotheses are exactly what the new transition still has to prove.

---

## 10. The precise remaining theorem target

The audit leaves one sharply stated, genuinely new route:

> **Small-clock all-order inertial transition.**  Construct localized
> divergence-free stages with (6.1)--(6.6) such that:
>
> 1. nonlinear stress transfers an order-\(E_j\) packet to the next
>    smaller scale while \(\Theta_j\to0\);
> 2. the outgoing state and all pressure/velocity endpoint jets match the
>    incoming data for stage \(j+1\);
> 3. the homological inverse through order \(M\) obeys a uniform
>    Gevrey-tame estimate of the form (8.3);
> 4. the stage energies and dissipations obey
>    \(\sum E_j+\sum E_j\Theta_j<\infty\); and
> 5. the resulting physical residual satisfies (2.2).

If such a theorem were proved, Lemma 8.1 would convert its remainder
ledger into terminal flatness, and the usual uniqueness-before-\(T\)
argument would make it relevant to alternative (D).  CDP provide useful
stress geometry and an existence proof in the opposite, heat-saturated
orientation.  They do not provide this theorem.

---

## 11. Final verdict

\[
\boxed{
\begin{array}{c}
\text{CDP proved transition: }\Theta_k\to\infty,\\[2mm]
\text{perturbative terminal anti-heat correction: }\Theta_k\to0,\\[2mm]
\text{polynomial finite-energy carrier: }\Theta_j\to0.
\end{array}}
\tag{11.1}
\]

The first and second lines have no common parameter range.  The third
line agrees with the correction range but violates the scale window and
heat saturation on which the CDP principal estimate rests.  Inflating a
separate carrier enough to restore that window gives the divergent
single-stage dissipation (7.9).

Accordingly, time reversal plus an anti-diffusive/Borel corrector is **not
an exact reduction of the desired forward cascade to the 2026 CDP
fixed-point theorem**.  It either incurs an exponential backward-heat
cost, or moves to a small-clock inertial regime for which the needed
all-order transition theorem is presently absent.

No claim about resolution of the Navier--Stokes Millennium problem is
made.
