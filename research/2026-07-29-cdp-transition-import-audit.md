# Cheskidov--Dai--Palasek as a transition module: exact import boundary

Date: 2026-07-29

## Result and claim boundary

This note audits the principal construction and corrector in
Cheskidov--Dai--Palasek (CDP), arXiv:2511.09556v2, against the proposed
forward finite-time forced cascade.  It does **not** prove a Navier--Stokes
singularity.

There is one exact positive import and four exact obstructions to importing
the whole CDP module.

1. **The rank-one geometric recursion imports.**  Once a strictly positive
   covariance has been specified, the CDP geometric lemma factors it
   smoothly into fixed rank-one tensors.  Their velocity-potential
   implementation preserves divergence and localization.  This solves the
   principal covariance factorization, including an isotropically gauged
   version of the three-beat target.
2. **The heat-assisted dynamics do not operate at the polynomial carrier.**
   On our stage the carrier heat exposure tends to zero.  Making it order
   one forces an exponentially growing normalized carrier.
3. **The scale window in the proved CDP principal estimate is even
   stronger.**  It requires
   \[
     N_{\rm hi}^{-2}\ll t_k\ll N_{\rm lo}^{-3},
   \]
   hence \(N_{\rm hi}^2\gg N_{\rm lo}^3\).  Neither adjacent outer scales
   nor a polynomial internal carrier satisfy this.
4. **Time reversal cannot give a terminally flat force.**  The
   nonlinear-compatible reversal has residual \(2\nu\Delta u\).  On the
   torus, if that residual even tended to zero, elliptic inversion would
   rule out blowup.
5. **The CDP corrector is a forward initial-value corrector, not an endpoint
   inverse.**  Reversing it incurs the backward heat multiplier.  In the
   CDP scale window this multiplier grows faster than every power, while a
   forward corrector supplies no prescribed outgoing child.

There is also an independent work obstruction.  The exact CDP daughter
covariance is a symmetric derivative modulo pressure.  Its total
contraction with a constant trace-free parent strain is zero.  Thus it can
program a daughter source, but it cannot by itself both drain the affine
parent and produce the daughter.  A work-carrying boundary/wake component
must remain.

The verdict is consequently narrow:

> CDP can be imported as a **principal positive-stress programmer** after a
> high-frequency child already exists.  Its proved heat hierarchy and
> forward corrector cannot be used as an exact polynomial-carrier
> parent-to-child transition.  Such an import would require a new
> inertial-time, endpoint-controlled parametrix.

---

## 1. What CDP actually construct

CDP normalize viscosity to one and introduce shells \(N_{j,k}\).  In
dimension at least three all directions within a generation have the same
frequency scale, while successive generations satisfy, schematically,

\[
 N_{k+1}=N_k^b,\qquad b>1\ \hbox{large}.
\tag{1.1}
\]

Their high-frequency velocity potentials are

\[
 \psi_{j,k}
 =N_{j,k}^{-2}\,
   \varphi_k*
   \left(
     a_{j,k}\phi_{j,k}\theta_j
     \sin(N_{j,k}\eta_j\cdot x)
   \right),
\qquad \theta_j\cdot\eta_j=0.
\tag{1.2}
\]

The approximate shell is

\[
 \bar v_k
 =-\sum_j N_{j,k}e^{-N_{j,k}^2t}\Delta\psi_{j,k}.
\tag{1.3}
\]

Thus its leading amplitude is \(N_{j,k}\), and heat removes it on the time
scale \(N_{j,k}^{-2}\).  Its self-stress contains

\[
 N_{j,k}^2e^{-2N_{j,k}^2t}
 a_{j,k}^2\theta_j\otimes\theta_j.
\tag{1.4}
\]

The time integral of the scalar coefficient is

\[
 \int_0^t N^2e^{-2N^2s}\,ds
 =\frac{1-e^{-2N^2t}}2.
\tag{1.5}
\]

This is the core of the complete high-to-low transfer: once
\(N^2t\gg1\), the covariance has deposited an order-one low stress.

The exact low shell is then defined causally by

\[
 v_k(t)
 =-\int_0^t e^{(t-s)\Delta}
   \mathbb P\operatorname{div}
   (\bar v_{k+1}\otimes\bar v_{k+1})(s)\,ds.
\tag{1.6}
\]

This direction is important.  The already present high shell creates the
lower one as time increases.

---

## 2. The exact algebraic import

Let

\[
 {\cal D}V=2\nabla\mathbin{\odot}V-2(\operatorname{div}V)I.
\tag{2.1}
\]

CDP choose fixed rational directions \(\theta_j\) whose tensors
\(\theta_j\otimes\theta_j\) span the symmetric matrices and generate \(I\)
with positive coefficients.  The standard geometric lemma supplies smooth
functions \(\Gamma_j\) such that

\[
 S=\sum_{j=1}^6\Gamma_j(S)^2
       \theta_j\otimes\theta_j
\tag{2.2}
\]

for every \(S\) in a fixed neighborhood of \(I\) in
\(\operatorname{Sym}_3\).

Their recursion is

\[
 a_{j,k+1}
 =c^{-1/2}\chi_{k+1}
   \Gamma_j(I+cS_k),
\qquad
 S_k=2{\cal D}\sum_jN_{j,k}\psi_{j,k}.
\tag{2.3}
\]

On the support of \(S_k\), where \(\chi_{k+1}=1\), this gives the **exact**
identity

\[
 \boxed{
 \sum_j a_{j,k+1}^2\theta_j\otimes\theta_j
 =2{\cal D}\sum_jN_{j,k}\psi_{j,k}+pI .
 }
\tag{2.4}
\]

The pressure \(pI\) is harmless after Leray projection.  CDP also prove

\[
 \|\nabla^na_{j,k+1}\|_\infty
 \lesssim_n N_{{\rm lo},k}^{\,n},
\qquad
 \|\nabla^n\psi_{j,k}\|_\infty
 \lesssim_n N_{j,k}^{-2+n},
\tag{2.5}
\]

with nested support.

This imports into the proposed transition as follows.

### Algebraic import lemma

Let \(Q(x)\) be a smooth symmetric target covariance on one transition
region.  Suppose that for some \(\rho>0\)

\[
 \left\|\rho^{-1}Q-I\right\|_{L^\infty}<c_0,
\tag{2.6}
\]

where \(c_0\) is the radius in the geometric lemma.  Then, with the fixed
six CDP directions,

\[
 Q(x)=\sum_{j=1}^6 a_j(x)^2
        \theta_j\otimes\theta_j,
\qquad
 a_j=\rho^{1/2}\Gamma_j(\rho^{-1}Q).
\tag{2.7}
\]

The amplitudes depend smoothly on \(Q\), preserve its Gevrey class, and
remain positive.  If an Euler--Reynolds stress \(R\) is not positive,
choose

\[
 Q=\rho I-R,\qquad
 \rho>\|R\|_{\rm op}
\tag{2.8}
\]

with a strict margin.  The isotropic term disappears under Leray
projection.  Hence positivity is not the obstruction to the proposed
three-beat covariance.

The CDP potential (1.2) also provides a useful localized,
exactly-divergence-compatible way to implement the oscillatory factors.
What (2.7) does **not** give is an exact finite-frequency identity between
the full quadratic velocity and the desired child.  That equality holds
only for the principal averaged covariance; oscillatory, envelope,
pressure, and cross-family terms remain to be solved dynamically.

In particular, (2.7) factors a tensor field.  It does not prove that the
six amplitudes lie in the image of the three material doublets, that their
charge chains reach the required endpoint, or that the outgoing velocity
has the prescribed phase and circulation.

---

## 3. Polynomial carrier versus the CDP heat clock

Use the polynomial-carrier cascade ledger

\[
 \ell_j=r^{-j},\qquad
 K_j=(j+1)^A,\qquad
 a_j=\ell_j^{-\gamma}K_j^\gamma,
\qquad 1<\gamma<\frac32,
\tag{3.1}
\]

with physical carrier and stage time

\[
 \Lambda_j=\frac{K_j}{\ell_j},
\qquad
 \tau_j=\frac{\ell_j}{a_j}
 =\ell_j^{1+\gamma}K_j^{-\gamma}.
\tag{3.2}
\]

The heat exposure of that carrier during the whole stage is

\[
\boxed{
 \Theta_j
 :=\nu\Lambda_j^2\tau_j
 =\nu\ell_j^{\gamma-1}K_j^{2-\gamma}.
}
\tag{3.3}
\]

Since \(\gamma>1\), the exponential factor
\(\ell_j^{\gamma-1}=r^{-(\gamma-1)j}\) beats every polynomial.  Therefore

\[
 \Theta_j\longrightarrow0.
\tag{3.4}
\]

Equation (1.5) shows why this is incompatible with a faithful CDP module.
Their principal transfer saturates only when the high wave experiences an
order-one heat clock.  At (3.4),

\[
 1-e^{-2\Theta_j}=2\Theta_j+O(\Theta_j^2),
\tag{3.5}
\]

so the heat-assisted deposited covariance vanishes in normalized stage
variables.

If the same \(K_j\) is required to be both the Kelvin carrier and the CDP
heat carrier, the condition \(\Theta_j\ge c>0\) gives

\[
 \boxed{
 K_j\ge
 \left(\frac c\nu\right)^{1/(2-\gamma)}
 \ell_j^{-(\gamma-1)/(2-\gamma)}.
 }
\tag{3.6}
\]

The right side is exponential in \(j\), not polynomial.

One might introduce a separate normalized heat carrier \(H_j\), while
retaining the amplitude and stage time in (3.1)--(3.2).  Its physical
frequency is \(H_j/\ell_j\), and order-one heat exposure requires

\[
 \nu\left(\frac{H_j}{\ell_j}\right)^2\tau_j\ge c.
\tag{3.7}
\]

Equivalently,

\[
 \boxed{
 H_j\ge
 \left(\frac c\nu\right)^{1/2}
 \ell_j^{-(\gamma-1)/2}K_j^{\gamma/2}.
 }
\tag{3.8}
\]

This is still exponential.  Thus separating the CDP carrier from the
Kelvin carrier does not retain a polynomial-frequency transition.

This calculation does not prove that every inertial high-to-low mechanism
needs (3.6).  It proves that the **heat-decay mechanism actually used by
CDP** does.

---

## 4. The stronger frequency window in the CDP proof

CDP Proposition 5.1 introduces a separation time \(t_k\) satisfying

\[
 N_{{\rm hi},k}^{-2}
 \ll t_k
 \ll N_{{\rm lo},k}^{-3}.
\tag{4.1}
\]

They take \(t_k=N_{{\rm lo},k}^{-4}\) and
\(N_{{\rm hi},k}\simeq N_{{\rm lo},k}^{\,b}\) with \(b>2\) sufficiently
large.  Merely for the interval in (4.1) to be nonempty one needs

\[
 \boxed{
 N_{{\rm hi},k}^{\,2}\gg
 N_{{\rm lo},k}^{\,3}.
 }
\tag{4.2}
\]

There are two natural identifications with our cascade, and both fail.

### 4.1 Internal carrier over an outer envelope

Take

\[
 N_{\rm lo}=\ell_j^{-1},
\qquad
 N_{\rm hi}=K_j\ell_j^{-1}.
\tag{4.3}
\]

Then (4.2) becomes

\[
 \boxed{K_j^2\gg\ell_j^{-1}.}
\tag{4.4}
\]

No polynomial \(K_j=(j+1)^A\) satisfies (4.4) because
\(\ell_j^{-1}=r^j\).

### 4.2 Adjacent physical cascade carriers

Take

\[
 N_{\rm lo}=\Lambda_j,\qquad
 N_{\rm hi}=\Lambda_{j+1}.
\tag{4.5}
\]

But

\[
 \frac{\Lambda_{j+1}}{\Lambda_j}
 =r\left(1+\frac1{j+1}\right)^A
 =r(1+O(j^{-1})).
\tag{4.6}
\]

Consequently

\[
 \frac{\Lambda_{j+1}^2}{\Lambda_j^3}
 =\frac{r^2(1+O(j^{-1}))}{\Lambda_j}
 \longrightarrow0,
\tag{4.7}
\]

the opposite of (4.2).

The upper clock \(N_{\rm lo}^{-3}\) is tied to CDP's short-time error
analysis, not asserted here as a universal law of Navier--Stokes.
Therefore (4.4) is a rigorous no-go for importing their **proved
principal module**, not for every conceivable redesigned high-to-low
cell.  Removing it requires a new transition estimate on the inertial
turnover time.

### 4.3 The required inflation is fatal to the dissipation ledger

Suppose an auxiliary high carrier \(H_j/\ell_j\) is used to realize a
strictly positive covariance which effects an order-one normalized
transfer during one outer stage.  The invariant requirement is on the
**time-integrated** stress, not on its instantaneous size.  Since the
parent strain has size \(a_j/\ell_j\), changing an order
\(a_j^2\ell_j^3\) amount of energy requires, up to fixed efficiency
constants,

\[
 \boxed{
 \int_{\text{stage}}\int_{\text{cell}}|w_j|^2\,dx\,dt
 \gtrsim a_j^2\ell_j^3\tau_j.
 }
\tag{4.8}
\]

Indeed, contraction with a strain of size \(a_j/\ell_j\) bounds the total
work by that strain times the left side of (4.8), and
\(\tau_j=\ell_j/a_j\).  Neither spatial nor temporal intermittency evades
(4.8): reducing occupied space-time volume requires increasing the
amplitude by the inverse square root to preserve the same integrated
positive covariance.  Maintaining covariance of size \(a_j^2\) throughout
the stage is the stronger special case.

For a field spectrally localized at the carrier \(H_j/\ell_j\), its
viscous loss over the stage is therefore bounded below at the scaling
level by

\[
\begin{aligned}
 D_j^{\rm hi}
 &\gtrsim
 \nu\left(\frac{H_j}{\ell_j}\right)^2
 \int_{\text{stage}}\int_{\text{cell}}|w_j|^2\,dx\,dt\\
 &\gtrsim
 \nu\left(\frac{H_j}{\ell_j}\right)^2
 \tau_j\,a_j^2\ell_j^3\\
 &=E_j\,
 \frac{\nu H_j^2}{a_j\ell_j},
 \qquad
 E_j\asymp a_j^2\ell_j^3.
\end{aligned}
\tag{4.9}
\]

The CDP proof window, under the identification (4.3), demands

\[
 H_j^2\gg\ell_j^{-1}.
\tag{4.10}
\]

Since

\[
 a_j\ell_j=\ell_j^{1-\gamma}K_j^\gamma,
\tag{4.11}
\]

the relative heat loss then obeys

\[
 \frac{D_j^{\rm hi}}{E_j}
 \gg
 \nu\ell_j^{\gamma-2}K_j^{-\gamma}
 \longrightarrow\infty.
\tag{4.12}
\]

More strongly,

\[
 D_j^{\rm hi}
 \gg
 \nu\ell_j^{1-\gamma}K_j^\gamma
 \longrightarrow\infty.
\tag{4.13}
\]

Thus the particular frequency inflation needed for the proved CDP
principal estimate is fatal to the finite-dissipation ledger, not merely
an aesthetic departure from a polynomial carrier.  Avoiding this
conclusion requires abandoning (4.1), changing the positive covariance
target, or supplying comparably large external work; the latter cannot be
a terminally flat perturbation of the proposed cascade.

---

## 5. Time reversal cannot produce a flat Clay force

The following elementary lemma is independent of the details of CDP.

### Lemma 5.1 (time-reversal flat-force obstruction)

Let \(u\) be a smooth unforced Navier--Stokes solution on
\(\mathbb T^3\times(0,\delta]\),

\[
 \partial_su-\nu\Delta u+
 \mathbb P\operatorname{div}(u\otimes u)=0,
\tag{5.1}
\]

whose mean is bounded and whose \(L^\infty\) norm is unbounded along a
sequence \(s_n\downarrow0\).  Set

\[
 v(t)=-u(T-t),\qquad T-\delta<t<T.
\tag{5.2}
\]

Then the Navier--Stokes residual of \(v\) is

\[
 \boxed{
 F_v
 :=\partial_tv-\nu\Delta v+
 \mathbb P\operatorname{div}(v\otimes v)
 =2\nu\Delta u(T-t).
 }
\tag{5.3}
\]

In particular, \(F_v\) cannot tend to zero in \(C^0\) as \(t\uparrow T\);
hence it cannot extend \(C^\infty\)-flatly by zero.

#### Proof

The sign in (5.2) is the only one that preserves the sign of the Euler
quadratic term under time reversal.  Direct substitution of (5.1) gives
(5.3).

Assume \(\|F_v(t)\|_\infty\to0\).  Then
\(\|\Delta u(s)\|_\infty\to0\).  For a mean-zero field on the torus, the
periodic Green kernel gives

\[
 \|u(s)-\bar u(s)\|_\infty
 \le C\|\Delta u(s)\|_\infty.
\tag{5.4}
\]

The spatial mean \(\bar u(s)\) is conserved for unforced Navier--Stokes.
Thus \(u(s)\) remains bounded in \(L^\infty\), a contradiction. \(\square\)

CDP's branch is unbounded as \(s\downarrow0+\).  Therefore its
time-reoriented branch cannot be converted into a forward terminal
singularity by declaring the reversal error to be a smooth-flat external
force.

At a single polynomial-carrier stage the same defect has the physical
size

\[
 \|\nabla^mF_v\|_\infty
 \sim
 \nu a_j\Lambda_j^{m+2}
 =
 \nu\ell_j^{-(\gamma+m+2)}
 K_j^{\gamma+m+2},
\tag{5.5}
\]

before any new hierarchy is solved.  It grows exponentially in every
fixed physical \(C^m\) norm.  Smallness relative to the inertial term,
\(\Theta_j\to0\), is not terminal flatness.

---

## 6. Why the forward corrector does not repair the endpoint

After constructing the principal inverse cascade \(v\), CDP obtain

\[
 \partial_tv-\Delta v+
 \mathbb P\operatorname{div}(v\otimes v)
 =\mathbb P\operatorname{div}f,
\tag{6.1}
\]

where Proposition 5.3 gives, schematically,

\[
 \|\nabla^nf(t)\|_{C^\kappa}
 \lesssim_n
 \epsilon_0
 \left(t^{-1-n/2+\alpha}+1\right).
\tag{6.2}
\]

This residual is singular, not flat.  CDP do not use it as an external
force.  They solve a corrector \(w\) with

\[
\begin{aligned}
 \partial_tw-\Delta w
 &+\mathbb P\operatorname{div}
 \left(
 w\otimes w+2(U+v)\mathbin{\odot}w
 \right)\\
 &=-\mathbb P\operatorname{div}
 \left(f+2U\mathbin{\odot}v\right),
 \qquad
 w|_{t=0}=0.
\end{aligned}
\tag{6.3}
\]

Their Proposition 6.2 is a forward semigroup bound.  Its mild loss
\((t/t')^\epsilon\) is useful because \(0<t'\le t\).  It gives no right
inverse for a prescribed value of \(w\) at the far endpoint.

For the scalar heat equation, the obstruction is already visible on one
Fourier mode:

\[
 \widehat w_n(t)
 =e^{-\nu|n|^2(t-t')}\widehat w_n(t').
\tag{6.4}
\]

Solving backward multiplies by

\[
 e^{+\nu|n|^2(t-t')}.
\tag{6.5}
\]

The explicit high-shell factors in the CDP principal part are precisely
\(e^{-N_{\rm hi}^2t}\).  In their principal window,
\(N_{\rm hi}^2t_k\to\infty\).  Directly reversing those factors therefore
costs

\[
 \exp(cN_{\rm hi}^2t_k)\longrightarrow\infty.
\tag{6.6}
\]

For their concrete choice
\(t_k=N_{\rm lo}^{-4}\),
\(N_{\rm hi}=N_{\rm lo}^b\), this is

\[
 \exp\!\left(cN_{\rm lo}^{\,2b-4}\right),
\qquad b>2,
\tag{6.7}
\]

which is not a polynomial carrier loss and is not compatible with a
Gevrey-tame endpoint graph.  Proposition 6.2 is a one-sided forward
estimate and supplies no terminal inverse that cancels this explicit heat
cost.  This does not assert a lower bound for every conceivable
linearized transition operator; it rules out reversing the heat factors
and corrector that CDP actually prove.

Keeping the corrector forward avoids (6.7), but then its outgoing state is
whatever the initial-value problem produces.  CDP need no recurrent
endpoint.  Our construction needs the outgoing child, its phase-center
profile, pressure jet, and wake charges to match the next stage through
increasing order.  Cutting off the unmatched forward state introduces
seam terms of the same order as the state; (6.2) supplies no
\(e^{-cj^2}\) endpoint mismatch.

Thus the CDP corrector can be borrowed only **after** a correctly oriented
principal transition and its outgoing state have been specified.  It
cannot create that endpoint condition.

---

## 7. Positive covariance is not positive work

Let \(S\) be a constant symmetric trace-free strain.  The CDP covariance
identity (2.4) has the form

\[
 Q=2{\cal D}V+pI.
\tag{7.1}
\]

For a periodic \(V\), or a compactly supported \(V\) on
\(\mathbb R^3\),

\[
\begin{aligned}
 \int Q:S
 &=2\int{\cal D}V:S+\int p\,I:S\\
 &=4\int(\nabla\mathbin{\odot}V):S\\
 &=0.
\end{aligned}
\tag{7.2}
\]

The isotropic gauge does not contribute because \(\operatorname{tr}S=0\).

For unresolved high waves with covariance \(Q\), the exchange with the
parent is proportional to \(-\int Q:S\).  Therefore the pure CDP
daughter covariance performs zero net work against an exactly affine
parent.  Making \(Q\) strictly positive does not change that fact.

If the parent strain varies on scale \(\ell\) and the covariance is
localized on diameter \(\delta\), subtracting the constant value at the
center gives

\[
 \left|\int Q:S(x)\,dx\right|
 \le
 C\delta\|\nabla S\|_\infty\|Q\|_{L^1},
\tag{7.3}
\]

unless a boundary/wake part with nonzero first moment is retained.

This distinguishes two tasks that the proposed transition had combined:

* CDP's rank-one recursion can prescribe the daughter stress divergence;
* an additional work-carrying component must drain the affine parent.

The outgoing annular wake is therefore not merely a localization error.
It is required by the stress--work identity.

---

## 8. Consequence for the three-beat transition

The three-beat calculation supplies a rank-five principal chart for the
low child strain.  The CDP geometric lemma is compatible with that chart
after adding an isotropic gauge: it gives six positive rank-one rays and
smooth amplitudes.  This is useful and exact at the tensor-factorization
level.

It does not supply the missing dynamic statement for three reasons.

1. In dimension three CDP suppress cross-polarization interactions by
   separating pipe supports.  The three-beat child, in contrast, is made
   by controlled unequal-radius cross interactions and must retain the
   associated charge chains.
2. CDP evolve an already existing higher shell into a lower shell.  The
   proposed Zeno recurrence must causally amplify and hand off the next
   higher child before the reset can occur.
3. The CDP low stress is derivative-exact but work-neutral against the
   affine pump.  The three-beat transition requires simultaneous parent
   drain and child creation.

Accordingly, the smallest honest theorem target is now:

> **Inertial CDP-factorized transition theorem.**  Given the amplified
> three-beat carrier and a work-carrying wake source, construct on one
> turnover interval a localized material-phase flow whose principal
> covariance is factored by (2.7), whose full finite charge band is solved
> forward, and whose projected outgoing child and wake jets have a
> Gevrey-tame right inverse.  No estimate may use an order-one heat clock or
> the scale window (4.1).

That theorem would genuinely import the reusable CDP geometry while
replacing the three pieces that point in the wrong direction: the heat
clock, infinite-frequency reservoir, and initial-value-only corrector.

No GPU run is presently needed for this audit.  The decisive gates above
are exact inequalities.  A computation becomes relevant only after the
inertial finite-band transition operator has been written; it should then
measure the charged endpoint Jacobian and corrector norm as the carrier and
bandwidth grow.

---

## Primary sources

* A. Cheskidov, M. Dai, and S. Palasek,
  *Instantaneous Type I blow-up and non-uniqueness of smooth solutions of
  the Navier--Stokes equations*, arXiv:2511.09556v2:
  https://arxiv.org/abs/2511.09556
* C. L. Fefferman,
  *Existence and smoothness of the Navier--Stokes equation*, official Clay
  problem description:
  https://www.claymath.org/wp-content/uploads/2022/06/navierstokes.pdf
