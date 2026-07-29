# Near-identity analytic return: the extra time integral and the honest recurrence gate

Date: 2026-07-29

## Status and verdict

This note audits the strongest remaining version of the one-carrier
near-identity proposal.  It does **not** construct a Navier--Stokes
singularity.

There is one correction which changes the local endpoint argument.
The exact rank-five sideband chart is a chart for the **time derivative
of the charge-zero child velocity**.  If the charge-zero strain is absent
at entry and the partner amplitudes are initial-data controls \(c\), then

\[
 U_0(t)=t\,{\cal L}_Kc+O(t^2),
\qquad
 \log DX_0(T)
 =\frac{T^2}{2}\nabla{\cal L}_Kc+O(T^3).
\tag{0.1}
\]

Thus a logarithmic carrier deformation \(\sigma\) requires

\[
 |c|\asymp \frac{\sigma}{T^2},
\tag{0.2}
\]

not \(\sigma/T\).  A stage with \(T\asymp\sigma\) is therefore not a
small-control inverse theorem: its controls are \(O(\sigma^{-1})\).
The best generic analytic balance is \(T\asymp\sqrt{\sigma}\), for which
the controls are order one.

For one fixed small ratio \(r>1\), this is not a divergence in the stage
index.  With the shifted carrier

\[
 \ell_j=r^{-j},\qquad K_j=(j+j_0)^A,
\tag{0.3}
\]

the exact logarithmic handoff

\[
 \sigma_j
 =\log r+A\log\left(1+\frac1{j+j_0}\right)
\tag{0.4}
\]

lies in one compact interval once \(j_0\gg A/\log r\).  Consequently the
target size does not obstruct a finite-dimensional endpoint submersion
uniform in \(j\); the required uniform \(K\)-tame \(C^2\) bound is a
separate analytic gate.  What
fails is the claim that making \(r-1\) arbitrarily small automatically
makes the partner controls and analytic reset easy.

There are two honest ways forward.

1. The recurrent packet may **enter with the low strain already present**.
   Then its carrier deformation is first order,
   \(\log DX_0(T)=TS_{\rm in}+O(T^2)\), and a stage
   \(T\asymp\sigma_j\) is possible.  The missing theorem is a return map
   for the complete low-strain/high-charge packet, not birth of the strain
   from a pure high carrier.
2. A from-rest construction may use \(T\asymp\sqrt{\sigma_j}\), carry a
   low strain of size \(O(\sqrt{\sigma_j})\) between stages, and prove a
   two-level pipeline return.  Its physical duration is
   \(T\ell_j/a_j\), not \(\sigma_j\ell_j^2\).

The pressure audit supplies another exact boundary.  Cancelling all
pressure multipoles through degree \(M\) requires
\(M^2+2M-3\) scalar conditions.  A tight single-carrier packet cannot
perform even the degree-two cancellation: its covariance is asymptotically
transverse to the one fast direction and cannot be isotropic.  The
necessary controls live in a nonperturbative annular, multi-directional
wake.  Increasingly many moment cancellations also cannot remain in one
compact nonzero Gaussian packet class; the correction must escape to a
growing normalized halo.

The conclusion is therefore a sharper theorem target, not a solution:
the near-identity route survives only as an **entry-strain or two-level
return problem coupled to an outgoing pressure-control wake**.

The finite-\(K\) rank and Taylor statements used below are checked in
[`checks/linearized_shear_endpoint_series.py`](../checks/linearized_shear_endpoint_series.py).
The pressure statements are proved in
[`2026-07-29-pressure-multipole-control-audit.md`](2026-07-29-pressure-multipole-control-audit.md).

---

## 1. The exact shifted-carrier ledger

Put

\[
 \delta=\log r,\qquad n_j=j+j_0,\qquad
 k_j=\frac{K_j}{\ell_j}=r^j n_j^A.
\tag{1.1}
\]

The physical carrier ratio is

\[
 \frac{k_{j+1}}{k_j}
 =r\left(1+\frac1{n_j}\right)^A,
\tag{1.2}
\]

so the required logarithmic deformation is exactly (0.4).  The elementary
bounds

\[
 \frac1{n+1}\le
 \log\left(1+\frac1n\right)\le\frac1n
\tag{1.3}
\]

give

\[
 \delta+\frac{A}{n_j+1}
 \le \sigma_j
 \le \delta+\frac{A}{n_j}.
\tag{1.4}
\]

In particular, if

\[
 j_0\ge\frac{A}{\varepsilon\delta},
\tag{1.5}
\]

then

\[
 \delta\le\sigma_j\le(1+\varepsilon)\delta
\quad\hbox{for every }j\ge0.
\tag{1.6}
\]

This stronger choice, rather than merely
\(j_0\gtrsim A/\varepsilon_*\), is useful for analytic recurrence:
spatial rescaling restores analytic radius by the factor \(r=e^\delta\),
whereas the material phase must pay for the full \(\sigma_j\).

If the normalized stage duration is \(T_j\), the physical duration is

\[
 \boxed{\Delta t_j=T_j\frac{\ell_j}{a_j}.}
\tag{1.7}
\]

For the Kelvin/Gavrilov amplitude choice

\[
 a_j=k_j^\gamma
 =r^{\gamma j}n_j^{A\gamma},
\qquad 1<\gamma<\frac32,
\tag{1.8}
\]

this becomes

\[
 \Delta t_j
 =T_j\,r^{-(1+\gamma)j}n_j^{-A\gamma}.
\tag{1.9}
\]

Hence both choices

\[
 T_j=c_0\sigma_j
\quad\hbox{and}\quad
 T_j=c_0\sqrt{\sigma_j}
\tag{1.10}
\]

give a summable physical time.  No factor \(\ell_j^2\) belongs in this
inertial time calculation.

The same shift makes all large-carrier estimates uniform from the first
retained stage.  If

\[
 M_j\lesssim\frac{n_j^2}{\log n_j},
\qquad A>4,
\tag{1.11}
\]

then

\[
 \max\left\{
 \sup_j\frac{M_j}{K_j},
 \sup_j\frac{M_j^2}{K_j},
 \sup_j K_j^{-1}
 \right\}
\longrightarrow0
\quad\hbox{as }j_0\to\infty.
\tag{1.12}
\]

Thus aliases, mixed viscous derivatives, and the finite-\(K\) correction
to the child chart can all be made uniformly small.  The endpoint time
power, not the polynomial carrier, is the issue below.

---

## 2. What the rank-five chart differentiates

Let the real base carrier be the shear with Fourier modes

\[
 \pm Ke_3,\qquad a=e_1,
\tag{2.1}
\]

and let the three partner pairs start at

\[
 r_\alpha=q_\alpha-Ke_3,
\qquad
\begin{aligned}
q_1&=(-45,-36,20),\\
q_2&=(-4,-5,9),\\
q_3&=(1,1,1).
\end{aligned}
\tag{2.2}
\]

For a transverse partner amplitude \(b\in r_\alpha^\perp\), the exact
matched child coefficient is

\[
 {\cal L}_{q_\alpha,K}b
 =
 P_{q_\alpha}\left[
 (e_1\cdot q_\alpha)b+(b\cdot q_\alpha)e_1
 \right].
\tag{2.3}
\]

The six corresponding symmetric-gradient columns have rank five for
every positive integer \(K\).  In uniformly normalized transverse bases,
their least nonzero singular value has a positive lower bound for all
sufficiently large \(K\), because the matrices converge to a rank-five
limit.

Let \(c\in\mathbb R^6\) be the partner coordinates and suppose that the
charge-zero child is zero initially.  The charge-zero Fourier equation
at \(t=0\) gives

\[
 \partial_t U_0(0)
 ={\cal L}_Kc.
\tag{2.4}
\]

The complete linearized propagation about the real shear can revisit the
low child through

\[
 q-Ke_3\longleftrightarrow q
\longleftrightarrow q+Ke_3\longleftrightarrow\cdots.
\tag{2.5}
\]

Exact rational Taylor calculation shows that every low-child coefficient
through order \(21\) is \(K^0\), and that the combined rank-five endpoint
minor remains nonzero at the tested carriers
\(K=128,256,512,1024\).  For the time-power argument only (2.4) is needed:
analytic dependence gives

\[
 U_0(t;c)-U_0(t;0)
 =t\,{\cal L}_Kc
 +O(t^2|c|+t|c|^2)
\tag{2.6}
\]

in any fixed lower analytic radius on which the flow exists uniformly.

Let \(X_0\) be the flow of the complete charge-zero field.  It obeys

\[
 \partial_tDX_0
 =(\nabla U_0)(t,X_0)DX_0.
\tag{2.7}
\]

If the entry low gradient vanishes at the selected core point, integrating
(2.6) once more gives

\[
 \boxed{
 \Pi_{\mathrm{Sym}_0}\log DX_0(T;c)
 =
 \frac{T^2}{2}{\cal A}_Kc+
 O(T^3|c|+T^2|c|^2),
 }
\tag{2.8}
\]

where \({\cal A}_K\) is the rank-five symmetric child chart.  Harmless
Fourier and real-cosine normalizations change only fixed constants.

Equation (2.8) is the missing integration in a direct use of the
instantaneous rank calculation.  The chart controls charge-zero
**acceleration**.  Material deformation is its time integral.

---

## 3. Quantitative local endpoint theorem and its true scale

The following statement is rigorous for every fixed finite Galerkin band.
For the full profile PDE it is conditional on the uniform analytic
Cauchy theorem described in Section 5.

### Proposition 3.1: fixed-\(K\) from-rest finite-band submersion

Fix an integer \(K\) outside the three exceptional child values and retain a finite charge/slow
band which contains (2.1)--(2.2).  For all sufficiently small \(T>0\),
define

\[
 {\cal E}_{K,T}(c)
 =
 \Pi_{\mathrm{Sym}_0}\log DX_0(T;c).
\tag{3.1}
\]

There are constants \(T_*(K),R_*(K),C_*(K)>0\) such that

\[
 D{\cal E}_{K,T}(0)
 =\frac{T^2}{2}{\cal A}_K+O(T^3),
\tag{3.2}
\]

the five-dimensional range has a right inverse of norm at most

\[
 \frac{C_*(K)}{T^2},
\tag{3.3}
\]

and every target \(Y\in\mathrm{Sym}_0^3\) satisfying

\[
 |Y|\le R_*(K)T^2
\tag{3.4}
\]

is attained by some \(|c|\le C_*(K)|Y|/T^2\).

#### Proof

For fixed \(K\), the finite Galerkin vector field and its flow are analytic in the initial
coordinates.  Equation (2.4), followed by (2.7), proves (3.2).
The exact finite-\(K\) rank theorem gives a positive lower singular value.
After dividing (3.1) by \(T^2/2\), its derivative is
\({\cal A}_K+O(T)\) and its second derivative is bounded on a fixed
control ball.  The quantitative inverse-function theorem gives
(3.3)--(3.4). \(\square\)

Uniformity as \(K\to\infty\) additionally requires a \(K\)-uniform
\(C^2\) estimate for the phase-quotiented endpoint map.  The exact
first derivative is uniformly onto, the linear propagator audit removes
the long \(K\)-chain, and the rational Taylor checker finds \(K^0\)
low-output coefficients through order \(21\).  These facts are strong
evidence for the required bound, but they are not an all-order nonlinear
proof.  Therefore every use of Proposition 3.1 uniformly along
\(K_j\to\infty\) below is conditional on that tame endpoint estimate.

For the target

\[
 Y_j=\sigma_jG_*,
\tag{3.5}
\]

where \(G_*\) is the required trace-free logarithmic deformation,
Proposition 3.1 requires

\[
 T_j^2\gtrsim\sigma_j,
\qquad
 |c_j|\asymp\frac{\sigma_j}{T_j^2}.
\tag{3.6}
\]

Consequently:

* \(T_j\asymp\sigma_j\) gives
  \(|c_j|\asymp\sigma_j^{-1}\);
* \(T_j\asymp\sqrt{\sigma_j}\) gives
  \(|c_j|\asymp1\);
* a fixed \(T_0\) gives \(|c_j|=O(\sigma_j)\), but does not by itself
  reset the analytic radius after infinitely many stages.

Because (1.6) keeps \(\sigma_j\) in one fixed compact interval, every one
of these control sizes is uniform in \(j\) after \(\delta>0\) is fixed.
The point is that only the square-root choice remains uniformly
order-one as the near-identity parameter is itself made smaller.

---

## 4. The analytic-radius reset inequality

Suppose a Cauchy estimate for the selected analytic profile path gives

\[
 \rho_{\rm out}
 \ge \rho_{\rm in}-\Lambda(c)T.
\tag{4.1}
\]

After stage \(j\), changing from the old normalized variable
\(y=x/\ell_j\) to the next one
\(z=x/\ell_{j+1}=r y\) multiplies the available strip radius by \(r\).
A packet ball with one fixed radius \(\rho_*\) is therefore invariant
under this estimate only if

\[
 r\big(\rho_*-\Lambda(c_j)T_j\big)\ge\rho_*,
\tag{4.2}
\]

or

\[
 \boxed{
 \rho_*
 \ge
 \frac{r}{r-1}\Lambda(c_j)T_j
 \sim\frac{\Lambda(c_j)T_j}{\delta}.
 }
\tag{4.3}
\]

It is not enough that each individual stage leaves some positive radius.
Without (4.2), repeating the same loss infinitely often exhausts the
analytic class.

For a generic quasilinear analytic estimate one expects, at minimum,

\[
 \Lambda(c)\lesssim C(1+|c|)
\tag{4.4}
\]

as a sufficient loss rate.  Combining (3.6) and (4.3) makes the quantity
to minimize

\[
 T+\frac{\sigma}{T}.
\tag{4.5}
\]

Its minimum is

\[
 2\sqrt{\sigma}
\quad\hbox{at}\quad
 T=\sqrt{\sigma}.
\tag{4.6}
\]

Since \(\sigma_j\asymp\delta\), the generic reset estimate therefore
requires

\[
 \boxed{\rho_*\gtrsim\delta^{-1/2}.}
\tag{4.7}
\]

The alternatives are worse:

\[
\begin{array}{c|c|c}
T&|c|&\hbox{radius required by (4.3)--(4.4)}\\ \hline
c_0\delta& O((c_0^2\delta)^{-1})&
O((c_0\delta)^{-1})\\
C\sqrt\delta&O(1)&O(\delta^{-1/2})\\
T_0&O(\delta)&O(T_0/\delta).
\end{array}
\tag{4.8}
\]

This is a gate for the **generic analytic Cauchy proof**, not a theorem
that the selected trajectory actually loses radius at that rate.
A structurally stable invariant class could do better.  But Gaussian and
periodic heat-kernel profiles do not make (4.7) free: their usual strip
norms grow like

\[
 \|e^{-|\cdot|^2}\|_{\rho}
 \asymp e^{C\rho^2}.
\tag{4.9}
\]

Likewise, for the periodic heat kernel

\[
 G_\tau(x)
 =\sum_{n\in\mathbb Z^3}e^{-\tau|n|^2}e^{in\cdot x},
\]

the weighted Fourier norm contains

\[
 \sup_n e^{\rho|n|-\tau|n|^2}
 =e^{\rho^2/(4\tau)+O(1)}.
\tag{4.10}
\]

At (4.7), this is \(e^{C/\delta}\), and that large norm feeds back into
the ordinary Cauchy constant \(\Lambda\).  Thus a tight Gaussian plus a
generic CK theorem does not close merely by taking \(r\) closer to one.
One needs a profile-specific symmetrizer, an entry-strain recurrence, or
a packet whose normalized spatial width grows at least like
\(\delta^{-1/2}\).  The last option is strongly nonlocal relative to the
thin ratio \(r-1\).

The polynomial carrier term in (0.4) does not change this conclusion if
(1.5) is imposed.  Without (1.5), the ratio
\(\sigma_j/\delta\) is large in the initial tail and (4.3) is even less
favorable.

---

## 5. Entry strain changes the time power

Suppose instead that the recurrent packet enters with a charge-zero
gradient

\[
 \nabla U_0(0)=S_{\rm in}\ne0.
\tag{5.1}
\]

Then

\[
 \Pi_{\mathrm{Sym}_0}\log DX_0(T)
 =
 TS_{\rm in}
 +\frac{T^2}{2}{\cal A}_Kc+
 O(T^2\|S_{\rm in}\|^2+T^3|c|).
\tag{5.2}
\]

If \(S_{\rm in}=G_*/c_0+O(\sigma_j)\) and
\(T_j=c_0\sigma_j\), the first term realizes
\(\sigma_jG_*\).  The sideband chart may then correct the endpoint
**low strain** itself:

\[
 \nabla U_0(T)-\nabla U_0(0)
 =T{\cal A}_Kc+O(T^2),
\tag{5.3}
\]

which is first order in \(T\).  Order-\(\sigma_j\) strain errors can be
removed with \(c=O(1)\).

This is the mathematically consistent short-stage architecture:

1. an incoming low strain performs the immediate carrier deformation;
2. the one-carrier sidebands adjust and replenish the low strain for the
   following stage;
3. the full low-strain/high-charge state, rather than a pure high carrier,
   belongs to the recurrent packet class.

It is not supplied by the instantaneous rank theorem.  A proof must
construct the incoming strain at the previous stage and show that its
five coordinates return with a uniform right inverse.  The alternative
\(T_j\asymp\sqrt{\sigma_j}\) suggests a two-level pipeline in which both
the carried strain and its one-step increment are \(O(\sqrt{\sigma_j})\).

---

## 6. Why material transport does not create sequential pulses

Putting every charged packet in the common material coordinate is exactly
what removes the \(O(K)\) low--high transport.  It also has a consequence
for the proposed sequencing.

Let \(X\) be the complete charge-zero flow and \(y=X^{-1}(t,x)\).
Every one-carrier packet has the form

\[
 A_{h,n}(t,y)
 \exp\big(i(Kh\,w+\xi_n)\cdot y\big).
\tag{6.1}
\]

At geometric-optics order all such packets have the same group velocity:
the low material velocity.  Their centers are stationary relative to one
another in \(y\).  Thus initially separated parent and partner packets do
not pass through one another sequentially merely because the material
phase is used.

Sequential activation requires an additional mechanism:

* a genuinely different low transport for different packets;
* an internal instability/dormancy clock;
* prescribed time envelopes whose defects are cancelled by a WKB
  hierarchy; or
* an external control, whose size must enter the terminal-force ledger.

None follows from the one-phase change of variables.  Fortunately,
sequencing is not needed for the local rank theorem.  With controls in a
small fixed ball, simultaneous partner--partner interactions contribute
to the \(O(|c|^2)\) remainder in (2.8), and the quantitative
inverse-function theorem absorbs them.  What sequencing was supposed to
solve is the global terminal wake, not the five-dimensional derivative.

---

## 7. Pressure moments force an outgoing halo

Let a localized velocity have stress \(T=u\otimes u\), pressure source

\[
 f=\partial_i\partial_jT_{ij},
\qquad -\Delta p=f.
\tag{7.1}
\]

For a homogeneous harmonic polynomial \(H\) of degree \(m\),

\[
 {\mathfrak m}_H(u)
 =
 \int Hf
 =
 \int (u\otimes u):D^2H.
\tag{7.2}
\]

Only these harmonic moments enter the Newtonian far field.  Since
\(\dim{\cal H}_m=2m+1\), the number through degree \(M\) is

\[
 \boxed{
 D_M=\sum_{m=2}^M(2m+1)
 =M^2+2M-3.
 }
\tag{7.3}
\]

This is smaller than the number of all monomials, but it still diverges.
A fixed list of six scalar partner amplitudes cannot cancel it.  Allowing
analytic envelope coefficients supplies enough formal dimensions, and
keeping all active corrections on the same fast direction preserves the
exact charged--charged \(K\)-cancellation.

Positivity of the physical stress gives a stronger obstruction.  Degree
two is dark precisely when

\[
 \int u\otimes u=cI.
\tag{7.4}
\]

For a one-carrier field with slow band \(M\ll K\), incompressibility gives

\[
 \|w\cdot u\|_2
 \lesssim\frac{M}{K}\|u\|_2.
\tag{7.5}
\]

Equation (7.4) would require one third of the total covariance in the
\(w\) direction, contradicting (7.5) for \(K/M\) large unless \(u=0\).
Therefore the active one-carrier family cannot cancel even its quadrupole
by adding more envelopes with the same fast direction.  The pressure
corrector needs an order-one longitudinal component: an annular low wake
or at least one additional fast direction.

The companion pressure audit constructs, for every finite \(M\), a
positive annular stress/velocity right inverse with \(O(M^2)\) controls.
It can cancel moments through

\[
 M_{{\rm p},j}\asymp (j+j_*)^2
\tag{7.6}
\]

while retaining order-one normalized energy.  At a fixed larger radial
ratio this gives an \(e^{-cj^2}\) pressure tail, and \(A>4\) still absorbs
the Gevrey-two coefficient count.

There is an exact compactness reason this annular piece cannot be hidden
inside one tight nonzero Gaussian packet.

### Proposition 7.1: increasing moment darkness is non-tight

Let \(f_j\) satisfy a uniform Gaussian bound

\[
 |f_j(x)|\le Ce^{-c|x|^2}
\tag{7.7}
\]

and enough uniform local derivative bounds to be precompact in
\(L^1(\mathbb R^3)\).  Suppose that every harmonic moment through degree
\(M_j\) vanishes, where \(M_j\to\infty\).  Then every \(L^1\) subsequential
limit \(f\) has

\[
 \int Hf=0
\quad\hbox{for every harmonic polynomial }H.
\tag{7.8}
\]

If, more strongly, all polynomial moments vanish (as would be required
for a compactly supported pressure rather than only a dark harmonic
exterior), then \(f=0\).

#### Proof

The Gaussian bound permits passage to the limit in every fixed moment.
For all polynomial moments, the Fourier transform of \(f\) extends to an
entire function and every derivative at the origin is zero.  Analytic
unique continuation gives \(\widehat f\equiv0\). \(\square\)

For pressure darkness alone, (7.8) says that the Newtonian potential of
the limit has no exterior multipoles; it does not force an arbitrary
signed \(f\) to vanish because Laplacians of localized functions are
harmonically dark.  Physical covariance positivity, (7.4)--(7.5), is the
additional obstruction for the one-carrier velocity class.

The practical conclusion is unambiguous: the moment controls must form a
separate halo/wake whose normalized location or internal complexity grows
with \(j\).  A pressure-dark analytic core plus a non-tight annular
corrector is consistent.  One uniformly tight one-carrier Gaussian packet
is not.

---

## 8. Terminal charged modes cannot simply be erased

The finite-band profile evolution is an analytic ODE.  Its time-\(T\)
map is locally invertible.  The charge-zero subspace is invariant.
Consequently, if an exact finite-band unforced trajectory lands with
every nonzero charge equal to zero, backward uniqueness implies that it
had no nonzero charge at entry.

This elementary observation rules out interpreting smooth time windows
as exact disappearance of the parent and partner modes.  In the full
viscous problem an analogous claim requires a backward-uniqueness
theorem, but the finite-band obstruction already identifies the correct
design principle:

> every nonzero terminal charge must recur as part of the next packet,
> be transported into a dynamically valid outgoing wake, or be removed
> by a residual whose size is honestly counted.

The existing Gavrilov lattice is an exact stationary Euler storage wake,
but it does not accept an arbitrary charged terminal state.  Landing on
its stationary manifold is an infinite-dimensional endpoint condition.
The five-dimensional strain submersion does not solve it.

Near identity helps only if the **renormalized complete state** has a
trapping or invariant-graph estimate.  A bound of the form

\[
 \|Z_{j+1}\|_{\cal A}
 \le(1+C\sigma_j)\|Z_j\|_{\cal A}+C\sigma_j
\tag{8.1}
\]

is insufficient when \(\sigma_j\to\delta>0\): its iteration grows
exponentially in \(j\).  The polynomial-carrier pressure ledger tolerates
uniform or polynomial deterioration of normalized constants, but not
this exponential growth.  One needs a genuine contraction after
renormalization, an exact packet reset, or a conserved/coercive analytic
quantity for the selected path.

---

## 9. Revised theorem blueprint

The surviving prize-level target can now be stated without hiding the
time integral.

> **Entry-strain analytic return theorem.**  Fix one sufficiently small
> \(\delta=\log r>0\), choose \(A>4\), and take
> \(j_0\ge C A/\delta\).  Construct a nonempty analytic packet class
> \({\cal P}\) containing:
>
> * an order-one charge-zero strain at entry;
> * one common material carrier and its partner/envelope state;
> * an annular multi-directional pressure-control wake; and
> * all terminal charged modes required for recurrence.
>
> For every \(\sigma\in[\delta,(1+\varepsilon)\delta]\) and every
> \(K\ge K_0\), the normalized evolution on
> \(T=c_0\sigma\) must admit controls which:
>
> 1. realize the exact carrier and amplitude ratios
>    \(e^\sigma,e^{\gamma\sigma}\);
> 2. return all five low-strain coordinates to \({\cal P}\), with a
>    uniformly bounded endpoint right inverse;
> 3. map the complete charged state into the recurrent component of
>    \({\cal P}\), rather than deleting it;
> 4. cancel pressure moments through
>    \(M_{{\rm p},j}\asymp j^2\) using the annular wake, while keeping its
>    extra fast directions dynamically separated from the core;
> 5. satisfy the analytic-radius reset inequality after rescaling;
> 6. have a renormalized trapping or graph-transform constant at most
>    polynomial in \(j\); and
> 7. leave only an \(e^{-cj^2}\) physical residual, including endpoint
>    viscous jets and the two cutoff seams.

An alternative theorem may start with zero low strain, but then it must
use \(T\gtrsim\sqrt{\sigma}\) and return a carried strain of size
\(O(\sqrt{\sigma})\).  The physical time is still summable by (1.9).

The exact achievements available now are:

* the shifted ledger (1.1)--(1.12);
* the finite-\(K\) rank-five acceleration chart;
* the all-generation single-carrier \(K\)-null identity;
* the finite-band \(K\)-tame linear propagator;
* the variable-phase pressure derivative count; and
* the finite-dimensional annular pressure-moment right inverse.

The missing result is not another carrier optimization.  It is the
complete entry-strain/charged-wake return theorem above.  A GPU endpoint
calculation can falsify its finite-band transversality, but cannot repair
the \(T^2\) time power, the analytic reset inequality, or the need for a
dynamically generated pressure-control wake.
