# The infinite two-harmonic charge ladder: a bounded forward semigroup and a mandatory terminal wake

**Date:** 2026-08-03
**Status:** exact frozen-ladder estimate, Dyson tail, short-pulse
source/wake ledger, and formal total-interaction majorant; nonlinear
material endpoint theorem remains open
**Scope:** the infinite hierarchy forced by C95.  This note does not
construct a Navier--Stokes singularity.

## 1. Verdict

The extreme-charge leakage in C95 is a finite-support obstruction, not a
forward-semigroup obstruction.

For every fixed slow Fourier sector \(q\), linearization about the real
two-harmonic one-phase bath is a nearest/next-nearest charge operator.  On
the nonzero charges its matrix edges are bounded independently of the
carrier \(K\) and independently of the extreme charge.  Consequently it
generates a uniformly bounded propagator on every analytic or finite-order
Gevrey charge space

\[
 X_{\rho,\sigma,s}^{p}
 =
 \left\{Z:
 \left\|
 \langle h\rangle^s e^{\rho |h|^{1/\sigma}} Z_h
 \right\|_{\ell_h^p}<\infty\right\},
 \qquad \rho\ge0,\quad \sigma\geq1,\quad s\in\mathbb R.
\tag{1.1}
\]

More precisely,

\[
 \|U_q(t,s)\|_{X_{\rho,\sigma,s}^{p}\to
 X_{\rho,\sigma,s}^{p}}
 \le
 \exp\left(\int_s^t\Lambda_q(\tau)\,d\tau\right),
\tag{1.2}
\]

where \(\Lambda_q\) depends on \(q\), the four bath amplitudes, and the
weight, but not on \(K\), a charge cutoff, or the WKB depth \(M\).
Put \(x_q(s,t)=\int_s^t\Lambda_q(\tau)\,d\tau\).  Starting from finite
charge support, the part which requires more than \(N\) total bath edges
has the factorial Dyson bound

\[
 \left\|1_{\{|h|>H+2N\}}U_q(t,s)Z\right\|
 \le
 e^{x_q(s,t)}
 \frac{x_q(s,t)^{N+1}}{(N+1)!}\|Z\|.
\tag{1.3}
\]

Thus the infinite ladder is instantly populated, as C95 requires, but its
extreme tail is factorially small on a short pulse.  There is no
charge-edge amplification analogous to \(K^N\).

Sequential partners enter by Duhamel.  Their \(L^1_t\) size, rather than
their possibly large \(L^\infty_t\) pulse height, controls the linear
wake.  The rank-five zero-charge response persists with an \(O(T)\)
correction.  At the same time, the invertible positive-second-harmonic
edge in C95 gives an order-one charged terminal wake whenever the
short-pulse low response is order one.  Carrying that wake is therefore
the natural endpoint; a bare terminal state is not obtained from the
partner-only short-pulse chart.

There are three important boundaries.

1. The estimate is for the charged-to-charged block.  The full Euler
   linearization also has a charge-zero-to-charged \(O(K)\)
   phase-translation block.  A uniform triangular model requires moving
   that block into a material-phase gauge and controlling the remaining
   low-flow blocks.  This note does not derive that coupled gauge.  The
   physical velocity reconstruction has the already known polynomial
   \(K\) conditioning.
2. Forward viscosity is harmless, but there is no bounded all-charge
   two-ended inverse on finite-radius Gevrey spaces.  Backward heat
   multiplies charge \(h\) by \(e^{\theta T h^2+O(|h|)}\), which is
   unbounded between any two such spaces.  A construction must retain or
   damp the terminal charged wake, control it by a more restricted
   endpoint map, or work on a finite band.
3. The constant in (1.2) is proportional to the fixed slow frequency
   \(|q|\).  It is uniform on the three first-order C95 sectors, not on
   the direct sum of all nonlinear slow sectors.  The latter costs one
   slow derivative and is only covered here by the formal
   interaction-order ledger.

At total interaction order \(r\), with the bath edges themselves counted,
the existing \(Cr^2\) one-phase coefficient ledger still closes a
\(C^r(r!)^2\) formal majorant through every finite \(M\).  This statement
does not yet give the same estimate for a control-derivative hierarchy in
which infinitely many bath edges have first been resummed.

The exact formulas and scalar inequalities are checked in
[the charge-ladder checker](../checks/two_harmonic_charge_ladder.py).

## 2. Exact nonzero-charge operator

Use the C95 coordinates

\[
 w=e_1,\qquad
 B(\theta,t)=
 \sum_{m\in{\cal M}}D_m(t)e^{im\theta},
 \qquad
 {\cal M}=\{-2,-1,1,2\},
\tag{2.1}
\]

with \(D_m\cdot w=0\) and the reality condition
\(D_{-m}=\overline{D_m}\).  The sine convention changes harmless scalar
phases only.  Fix a nonzero slow frequency \(q\), assume the lifted modes
below are nonaliased, and write

\[
 p_h=q+hKw,\qquad Z_h\cdot p_h=0.
\tag{2.2}
\]

The symmetrized Euler interaction of \(D_m e^{imK w\cdot x}\) with
\(Z_h e^{ip_h\cdot x}\), suppressing the scalar Fourier factor, is

\[
 {\mathbb P}_{p_{h+m}}
 \left[
 (D_m\cdot p_h)Z_h+(Z_h\cdot mKw)D_m
 \right].
\tag{2.3}
\]

For \(h\ne0\), transversality gives

\[
 KZ_h\cdot w=-\frac{Z_h\cdot q}{h}.
\tag{2.4}
\]

Since \(D_m\cdot w=0\), (2.3) becomes the exact identity

\[
 \boxed{
 T_{m,h}^{q,K}Z_h
 =
 {\mathbb P}_{p_{h+m}}
 \left[
 (D_m\cdot q)Z_h
 -\frac{m}{h}(Z_h\cdot q)D_m
 \right].}
\tag{2.5}
\]

In particular,

\[
 \|T_{m,h}^{q,K}\|
 \le
 (1+|m|)|q|\,|D_m|
 \le 3|q|\,|D_m|,
 \qquad h\ne0.
\tag{2.6}
\]

There is no \(K\).  The ratio term improves rather than worsens as
\(|h|\to\infty\).  Define the charged block by deleting charge-zero input
and output:

\[
 ({\cal L}^{cc}_{q,K}Z)_j
 =
 -i
 \sum_{\substack{m\in{\cal M}\\j-m\ne0}}
 1_{\{j\ne0\}}T_{m,j-m}^{q,K}Z_{j-m}.
\tag{2.7}
\]

The omitted charged-to-zero row is retained as the bounded observation

\[
 {\cal O}_{q,K}Z
 =
 -i\sum_{m\in{\cal M}}T_{m,-m}^{q,K}Z_{-m}.
\tag{2.8}
\]

Only four input charges occur in (2.8), so it has the same \(K\)-uniform
bound.

## 3. Weighted convolution estimate

Put \(a=1/\sigma\in(0,1]\) and

\[
 W_h=\langle h\rangle^s e^{\rho|h|^a}.
\tag{3.1}
\]

For \(|m|\le2\),

\[
 |h+m|^a\le |h|^a+|m|^a,
\qquad
\langle h+m\rangle^s
\le (1+|m|)^{|s|}\langle h\rangle^s.
\tag{3.2}
\]

Therefore the shift \(S_mZ=(Z_{j-m})_j\) obeys, for every
\(1\le p\le\infty\),

\[
 \|S_mZ\|_{X_{\rho,\sigma,s}^p}
 \le
 (1+|m|)^{|s|}
 e^{\rho|m|^{1/\sigma}}
 \|Z\|_{X_{\rho,\sigma,s}^p}.
\tag{3.3}
\]

Combining (2.6) and (3.3) gives

\[
 \boxed{
 \|{\cal L}^{cc}_{q,K}(t)\|
 \le\Lambda_q(t):=
 |q|\sum_{m\in{\cal M}}
 (1+|m|)^{|s|+1}
 e^{\rho|m|^{1/\sigma}}|D_m(t)|.}
\tag{3.4}
\]

This proves (1.2) by the Dyson series, including for time-dependent bath
amplitudes.  The same edge-by-edge triangle estimate applies after every
finite Galerkin truncation, so the bound is uniform over the charge
cutoff.  The reverse inviscid propagator has the same estimate with the
time integral reversed.

For \(1\le p<\infty\), and on the corresponding weighted \(c_0\) space
when using a sup norm, forward diagonal heat can be included by bounded
perturbation of its contractive \(C_0\)-semigroup.  Unequal damping of the
two bath harmonics merely makes \(D_m(t)\) time dependent and decreases the
integral in (3.4).  On raw weighted \(\ell^\infty\), the heat family is not
strongly continuous at zero; fixed Galerkin estimates remain valid, but no
\(C_0\)-semigroup assertion is made.

### Factorial extreme tail

Every application of \({\cal L}^{cc}\) changes charge by at most two.
If \(Z\) is supported on \(|h|\le H\), the first \(N\) Dyson terms remain
inside \(|h|\le H+2N\).  Ordered time integration gives

\[
 \|U^{(n)}(t,s)\|
 \le\frac{(\Lambda_*T)^n}{n!},
\qquad
\Lambda_*=\sup_{[s,t]}\Lambda_q,\quad T=t-s.
\tag{3.5}
\]

Consequently,

\[
\begin{aligned}
 \left\|1_{\{|h|>H+2N\}}U(t,s)Z\right\|
 &\le
 \sum_{n=N+1}^{\infty}
 \frac{(\Lambda_*T)^n}{n!}\|Z\|\\
 &\le
 e^{\Lambda_*T}
 \frac{(\Lambda_*T)^{N+1}}{(N+1)!}\|Z\|.
\end{aligned}
\tag{3.6}
\]

C95's injective extreme edge shows that successive formal edge
coefficients generically populate the tail.  Equation (3.6) proves that
its existence alone causes no analytic or Gevrey loss.

### The scheduled charge cutoff is already flat

Let the retained charge radius be \(R>H\) and put

\[
 N_R=\left\lfloor\frac{R-H}{2}\right\rfloor.
\tag{3.7}
\]

Using \(n!\ge(n/e)^n\), (3.6) gives, for
\(x=\Lambda_*T\),

\[
 \left\|1_{\{|h|>R\}}U(t,s)Z\right\|
 \le
 e^x
 \left(\frac{ex}{N_R+1}\right)^{N_R+1}\|Z\|.
\tag{3.8}
\]

Thus a bounded short-pulse clock \(x=O(1)\) produces

\[
 \log\frac{\|1_{\{|h|>R\}}UZ\|}{\|Z\|}
 \le -\frac R2\log R+O(R).
\tag{3.9}
\]

At the repository's planned depth

\[
 R_j=M_j\sim\kappa\frac{j^2}{\log j},
\tag{3.10}
\]

this becomes

\[
 \boxed{
 \|1_{\{|h|>M_j\}}U_jZ\|
 \le
 \exp\big(-\kappa(1-o(1))j^2\big)\|Z\|.}
\tag{3.11}
\]

Hence the **linear bath-ladder truncation** itself meets the required
\(e^{-cj^2}\) flatness scale.  This conclusion uses fixed first-order
slow sectors and a bounded pulse clock.  It does not estimate nonlinear
slow-sector proliferation or a long window for which
\(\Lambda_*T\) grows comparably to \(M_j\).

## 4. Sequential partner pulses

Let \(P_j(t)\) be the prescribed charge-\(-1\) partner in slow sector
\(q_j\), with its reality conjugate, and assume the three time supports
are disjoint.  For the endpoint-rank statements, use fixed pulse shapes
so that

\[
 \|P_j\|_{L^1_t}\le C_P\left|\int P_j(t)\,dt\right|
\tag{4.0}
\]

on the selected two-dimensional polarization coordinates.  Put
\(P=\sum_jP_j\).  Split the perturbation into the prescribed partner and
the generated charged wake \(Z^c\).  At the linearized level,

\[
\begin{aligned}
 \partial_tZ^c
 &={\cal L}^{cc}Z^c+\Pi_c{\cal L}P,\\
 \partial_tY_0
 &={\cal O}(Z^c+P).
\end{aligned}
\tag{4.1}
\]

Here \(Y_0\) collects the desired low children.  Equation (4.1) is a
**conditional triangular model**: assume the dangerous low-to-bath phase
translation has been moved into a material-phase variable and the
remaining low-flow blocks are uniformly bounded.  Deriving that gauge
for the coupled localized flow is outside this note.  Under this
assumption, those bounded blocks contribute to the \(O(T)\) correction
below.  The ungauged physical-velocity system is not claimed to obey
(4.1).  Duhamel and (3.4) give

\[
 \|Z^c(T)\|_X
 \le
 e^{\Lambda_*T}
 \left(
 \|Z^c(0)\|_X+
 \Lambda_*\|P\|_{L^1_tX}
 \right).
\tag{4.2}
\]

Thus a pulse of height \(T^{-1}\) and width \(T\) has order-one linear
wake, not an order-\(K\) or order-\(T^{-1}\) wake.  Moreover,

\[
 Y_0(T)
 =
 \int_0^T{\cal O}(t)P(t)\,dt
 {}+O\left(
 T e^{\Lambda_*T}\Lambda_*^2
 \|P\|_{L^1_tX}
 \right).
\tag{4.3}
\]

For a frozen bath the leading term is exactly the six-column rank-five
C95 chart.  It remains an \(O(T)\) perturbation of that chart when the
bath has a uniform \(O(T)\) variation across the pulse.  Hence its
endpoint minor persists uniformly for a sufficiently short pulse under
that hypothesis.

Disjoint supports remove simultaneous partner--partner products.  They
do not remove interaction of a later partner with an earlier retained
wake.  That term is quadratic in the partner controls at the zero-control
Jacobian and belongs to the higher-order source hierarchy below.  Treating
an order-one partner path as part of the background instead gives the
one-slow-derivative fixed-star estimate of the earlier single-carrier
note, not the order-zero bath estimate (3.4).  This distinction is
load-bearing.

### A charged wake is visible already at first order

For a charge-\(-1\) partner, the positive second bath harmonic sends the
same integrated polarization to charge \(+1\):

\[
 Z^c_{+1}(T)
 =
 \int_0^T T_{2,-1}^{q,K}(t)P(t)\,dt
 {}+O\left(
 T e^{\Lambda_*T}\Lambda_*^2
 \|P\|_{L^1_tX}
 \right).
\tag{4.4}
\]

C95 proves that \(T_{2,-1}^{q_j,K}\) is invertible for all three checked
directions and every integer \(K\ge62\).  Its smallest singular value is
uniformly positive in the large-\(K\) regime: the transverse planes and
the matrix converge to a nonsingular limit, and the remaining bounded
set of integer \(K\)'s is finite.  For a frozen bath, (4.4) factors as
the displayed map applied to \(\int P\).  The same conclusion holds for
a time-dependent bath on fixed pulse shapes if
\(T_{2,-1}(t)=T_{2,-1}(t_*)+O(T)\) across the pulse, with a uniform
modulus.  Under that hypothesis, for sufficiently short pulses satisfying
(4.0), a nonzero leading low response from a partner-only source leaves
a nonzero charge-\(+1\) endpoint wake of comparable control order.

This is not a no-go against extra infinite-dimensional controls.  It says
that the six-coordinate partner chart by itself is not a bare-endpoint
chart.  Its correct codomain contains the charged wake.

## 5. Formal Gevrey-two interaction majorant

The uniform bath propagator shows that extreme-charge edges carry no
hidden \(K\)-factor.  In a formal total-interaction expansion which counts
each such edge, the exact one-phase incompressibility identity bounds an
order-\(r\) quadratic source by \(Cr^2\); this is the same count used in
the forward multiphase parametrix.  The relevant exact convolution is

\[
 r^2\sum_{p=1}^{r-1}
 \frac{(p!)^2((r-p)!)^2}{(r!)^2}
 =
 r^2\sum_{p=1}^{r-1}\binom rp^{-2}
 \le3.
\tag{5.1}
\]

For completeness, the two endpoint terms contribute exactly \(2\) when
\(r\ge3\) (at \(r=2\) they coincide and contribute \(1\)).
For \(r\ge4\), every interior binomial coefficient is at least
\(\binom r2\), so the interior contribution is at most

\[
 \frac{4(r-3)}{(r-1)^2}\le1.
\tag{5.2}
\]

Let \(\Gamma\) denote the time-integrated normalized bilinear coefficient.
For bounded profiles \(\Gamma\le CT\).  For a fixed-\(L^1\),
height-\(T^{-1}\) partner pulse, the partner--wake part contributes
\(O(\|P\|_{L^1})\), still independently of \(K\) and \(M\); exact partner
self-interaction in the frozen plane-wave model and disjoint partner
supports remove the apparent \(\int|P|^2\) term.  Suppose the
total-interaction coefficients obey

\[
 a_1\le A,\qquad
 a_r
 \le
 E\Gamma\,r^2
 \sum_{p=1}^{r-1}a_pa_{r-p},
 \qquad
 E=e^{\Lambda_*T}.
\tag{5.3}
\]

Then induction using (5.1) gives

\[
 \boxed{
 a_r\le
 A B^{r-1}(r!)^2,
 \qquad
 B\ge3E\Gamma A.}
\tag{5.4}
\]

The same estimate holds for the terminal coefficient \(a_r(T)\), because
the endpoint wake is the Duhamel value of the same forward equation.
Sequential partner--wake products enter the convolution in (5.3); they
do not introduce a carrier factor when all their constituent interaction
edges are counted.

Equation (5.4) is an all-finite-order Gevrey-two majorant, uniform through
any requested depth \(M\).  Its precise scope matters.  It counts total
interaction trees, with each bath edge included in the interaction order,
and uses the previously derived \(Cr^2\) coefficient bound at that order.
It is not automatically a majorant for derivatives in the partner
controls after the bath semigroup has been fully resummed: those
coefficients already have infinite charge support.  Proving that stronger
statement requires a weighted bilinear estimate on the resummed tails.
Nor is (5.4) a nonlinear Cauchy theorem for an untruncated material
profile.  A selected-path proof still has to control slow derivatives,
variable pressure, localization, and the control-dependent flow without
repeated radius loss.

## 6. No bounded all-charge viscous endpoint inverse

Let

\[
 \theta=\nu K^2>0.
\tag{6.1}
\]

On high charge,

\[
 \nu|q+hKw|^2T
 =
 \theta T h^2+O(|h|).
\tag{6.2}
\]

Forward heat is contractive.  Already for the free heat block, a
two-ended map accepting an arbitrary complete terminal charged state
would require the inverse multiplier

\[
 e^{\nu|q+hKw|^2T}.
\tag{6.3}
\]

Even if the input and output Gevrey radii differ by an arbitrary finite
amount \(\delta>0\), testing (6.3) on a single charge gives the ratio

\[
 \exp\left(
 \theta T h^2-\delta |h|^{1/\sigma}+O(|h|)
 \right)\longrightarrow\infty.
\tag{6.4}
\]

Thus backward heat is unbounded between every pair of ordinary
finite-radius Gevrey charge spaces.  A bounded finite-shift bath
perturbation does not repair this: the positive-time parabolic evolution
remains a compact operator on the weighted \(\ell^p\) sequence space for
\(1\le p<\infty\), and on the corresponding weighted \(c_0\) sup space,
so it cannot have a bounded inverse in infinite dimension.  On raw
\(\ell^\infty\), the heat family is not strongly continuous at zero; the
single-charge multiplier test (6.4) still rules out the bounded inverse,
but no \(C_0\)-semigroup assertion is made there.  No analytic weighting
supplies a bounded all-charge two-ended inverse.  This does not exclude a
specially constrained terminal control map.  On a finite band
\(|h|\le M\), the familiar cost is instead

\[
 \exp(C\theta TM^2).
\tag{6.5}
\]

This cleanly separates the two conclusions:

* extreme-charge leakage is harmless for a retained **forward** wake;
* there is no bounded general all-charge terminal-erasure inverse in the
  viscous problem.

## 7. What has and has not closed

The C95 hierarchy now has a complete frozen linear answer.

1. The two-harmonic bath generates a bounded convolution operator on
   analytic and Gevrey charge spaces, uniformly in \(K\) and charge depth.
2. The infinite extreme tail has a factorial short-pulse estimate.
3. Sequential partner sources preserve the rank-five low Jacobian and
   produce a controlled but mandatory charged terminal wake.
4. The formal total-interaction hierarchy retains the \(C^r(r!)^2\)
   majorant through every finite \(M\); the fully resummed nonlinear
   control-derivative hierarchy remains open.
5. A bounded viscous all-charge two-ended inverse on ordinary
   finite-radius Gevrey spaces is excluded; a restricted controlled
   endpoint or finite-band construction is not.

The remaining theorem is not a charge-edge estimate.  It is a
control-dependent localized material-phase endpoint map on a state space
which includes both the zero-charge global wake and the charged wake.
It must preserve a positive slow/Gevrey radius, retain the rank-five
minor, and meet the C90 stage-to-stage carry gate.  The present result
removes extreme-charge leakage itself from that list of possible fatal
losses.
