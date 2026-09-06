# C192: the certified zero-drift return expands by more than \(3000\)

**Date:** 2026-08-28

**Status:** outward-rounded strengthening of the C159 principal cocycle and
the C185/C189 operator-norm consequence; explicit shorter deficit and
preparation clocks; no essential-spectral-radius, finite-band, viscous,
nonlinear, UVSR, or singularity claim

**Checker:**
[checks/strong_zero_drift_gain_c192.py](../checks/strong_zero_drift_gain_c192.py)

## 0. Result and claim boundary

Let \(M\) be the one-period transverse Kelvin monodromy on C152's
zero-drift orbit and C159's returning \(m\ne0\) covector.  In C159's
returning orthogonal frame, put

\[
                         w=(1,3/20)^T.
\tag{0.1}
\]

The same outward-certified phase/covector tube used by C159 in fact gives

\[
 \boxed{Mw>3000w\quad\hbox{componentwise},\qquad
        \rho(M)>3000>e^8.}
\tag{0.2}
\]

No unstable amplitude column is integrated.  The checker subdivides the
certified path into 2048 cells, bounds the four entries of the cooperative
generator from below on each cell, and propagates a positive lower vector.
This avoids the interval wrapping that prevented C153 from certifying its
floating monodromy.

Using only the C189-approved abstract spectral inclusion, (0.2) upgrades
C185's robust operator-norm form to

\[
 \boxed{\|G_{rT}\|_{L^2\to L^2}\ge3000^r>e^{8r}
        \qquad(r\in\mathbb N).}
\tag{0.3}
\]

No statement about \(r_{\rm ess}\) is made.  Equation (0.3) is still an
inviscid, unrestricted operator-norm lower bound.  It does not by itself
identify a finite-frequency retained vector or prove fixed-energy
concentration.  Thus it does not reverse C191's normalized-endpoint
verdict.  Its concrete contribution is that the formerly long
finite-frequency window is no longer excluded at the power level.

## 1. Certified generator on the C159 tube

C159 writes the transverse amplitude equation in unit-period time as

\[
                              z'=B(s)z,\qquad 0\le s\le1.
\tag{1.1}
\]

Its checker certifies the phase variables and scalar covector coordinate
on 32 degree-40 Taylor panels, including all residual, panel-jump,
parameter, and initial-algebraic errors.  The resulting global tube is

\[
 \|\Delta(\cos a,\sin a,\cos b,\sin b)\|_\infty<2\times10^{-6},
 \qquad |\Delta\gamma|<8\times10^{-4}.
\tag{1.2}
\]

C192 reruns that proof and reevaluates the raw entries of \(B\) on 64
subcells of every Taylor panel.  Thus the cell length is

\[
                              \Delta={1\over2048}.
\tag{1.3}
\]

The C192 adversarial pass found an inherited implementation defect in
C159's shared interval primitives: Python Decimal unary negation and
absolute value can use the ambient 28-digit context. The primitives now
use exact copy-negate and copy-absolute operations instead. C159's full
tube certificate and C192's product were both rerun after that repair.

If \(L_j\) is the matrix of directed lower endpoints on cell \(j\), the
complete computation gives the global rounded bounds

\[
 (L_j)_{11}>-5,\qquad (L_j)_{12}>40,\qquad
 (L_j)_{21}>{97\over100},\qquad (L_j)_{22}>-{9\over2}.
\tag{1.4}
\]

In particular every actual \(B(s)\) and every \(L_j\) is Metzler.

## 2. Cellwise positive comparison

For each cell define, with a directed positive safety charge,

\[
 \alpha_j=\max\{0,-(L_j)_{11},-(L_j)_{22}\}+10^{-60},
 \qquad C_j=L_j+\alpha_jI.
\tag{2.1}
\]

Then \(C_j\ge0\) entrywise, and the certificate gives

\[
                         0\le\alpha_j<5,\qquad
                         \alpha_j\Delta<{5\over2048}<{1\over400}.
\tag{2.2}
\]

Because all powers of \(C_j\) are entrywise nonnegative,

\[
\begin{aligned}
 e^{L_j\Delta}
 &=e^{-\alpha_j\Delta}e^{C_j\Delta}\\
 &\ge(1-\alpha_j\Delta)
       \sum_{d=0}^{6}{(C_j\Delta)^d\over d!}
       =:Q_j
\end{aligned}
\tag{2.3}
\]

entrywise.  This uses only \(e^{-x}\ge1-x\) for \(x\ge0\) and positivity
of the omitted exponential-series terms.

Let \(\Phi_j(t,s)\) be the actual two-parameter cell propagator. Since
\(B(s)-L_j\ge0\), the positive Duhamel identity

\[
 \Phi_j(\Delta,0)-e^{L_j\Delta}
 =\int_0^\Delta
   \Phi_j(\Delta,\sigma)\{B(\sigma)-L_j\}
   e^{L_j\sigma}\,d\sigma
\tag{2.4}
\]

is entrywise nonnegative.  Therefore

\[
                         \Phi_j(\Delta,0)y\ge Q_jy
                         \qquad(y\ge0).
\tag{2.5}
\]

All arithmetic in the product \(Q_{2047}\cdots Q_0w\) is directed
downward.  The checker obtains

\[
 Q_{2047}\cdots Q_0w>
 \begin{pmatrix}
 3137.3029581527775142338\\
 503.6727017585062824359
 \end{pmatrix}.
\tag{2.6}
\]

The second component of \(3000w\) is \(450\).  Equations (2.5)--(2.6)
prove (0.2).  The returning frame is orthogonal and has the same column
lengths at each period endpoint, so the componentwise comparison also
gives physical Euclidean fiber-norm enlargement.  Positivity of \(M\)
then yields

\[
                         M^rw>3000^rw.
\tag{2.7}
\]

Finally the checker proves \(e^8<3000\) by an exact rational Taylor upper
sum with a geometric tail after order 32.  No floating logarithm enters
the load-bearing comparison.

## 3. Infinite-dimensional operator consequence

C185/C189 permits the abstract-level Shvydkoy spectral inclusion for this
velocity-form bicharacteristic-amplitude system.  From (2.7), the maximal
principal Lyapunov exponent satisfies

\[
                         \mu_{\max}>{\log3000\over T}>{8\over T}.
\tag{3.1}
\]

Applying exactly the already-approved operator-norm step gives

\[
 r(G_{rT})\ge3000^r,
 \qquad
 \|G_{rT}\|_{L^2\to L^2}\ge r(G_{rT})\ge3000^r,
\tag{3.2}
\]

which proves (0.3).  This neither uses nor establishes the citation-held
essential-spectral-radius identity.

## 4. The C191 clocks shorten by a factor forty

Retain C188's C180-preserving schedule

\[
 q=n^8,\qquad b=n^{-5/2},\qquad H=n^{51/2},\qquad n=j+1\ge2.
\tag{4.1}
\]

For the raw C182 deficit define

\[
 R_\Delta^{(192)}=\left\lceil{3\over8}\log n\right\rceil
                  =\left\lceil{3\over64}\log q\right\rceil.
\tag{4.2}
\]

Then

\[
 \boxed{
 \|G_{R_\Delta^{(192)}T}\|_{2\to2}
 \ge e^{8R_\Delta^{(192)}}
 \ge n^3=q^{3/8}.}
\tag{4.3}
\]

Using C159's \(3<T<76/25\), the required inertial action obeys

\[
 \boxed{
 {9\over64}\log q
 <TR_\Delta^{(192)}
 <{57\over400}\log q+{76\over25}.}
\tag{4.4}
\]

C191's upper coefficient was \(57/10\); the new coefficient \(57/400\)
is exactly forty times smaller.

If the same floor is also assigned to the formal preparation factor, put

\[
 R_*^{(192)}=\left\lceil{57\over16}\log n\right\rceil
             =\left\lceil{57\over128}\log q\right\rceil.
\tag{4.5}
\]

Then

\[
 \boxed{
 \|G_{R_*^{(192)}T}\|_{2\to2}
 \ge e^{8R_*^{(192)}}
 \ge n^{57/2}=Hq^{3/8},}
\tag{4.6}
\]

and

\[
 \boxed{
 {171\over128}\log q
 <TR_*^{(192)}
 <{1083\over800}\log q+{76\over25}.}
\tag{4.7}
\]

Again the logarithmic coefficient is exactly forty times smaller than
C191's \(1083/20\).

## 5. Explicit scalar feasibility with the shorter clock

C188 proves \(\log n\le6n^{1/14}\).  Hence

\[
                       R_*^{(192)}<23n^{1/14}.
\tag{5.1}
\]

Replacing C191's bound \(856n^{1/14}\) by (5.1), the declared C176/C188
nonviscous collar satisfies

\[
 \boxed{
 {C_{\rm col}(R_*^{(192)})^{7/2}q^{-1}\over b^3}
 <C_{\rm col}23^{7/2}n^{-1/4}.}
\tag{5.2}
\]

The pump retains at least \(99/100\) of its amplitude throughout this
clock under the explicit condition

\[
 \boxed{
 {152\over25}\nu(j!)^{-3/2}
 \left({57\over16}\log n+1\right)\le{1\over100}.}
\tag{5.3}
\]

Indeed, with \(\mu_j=\nu(j!)^{-3/2}\), condition (5.3) gives
\(x=2\mu_jTR_*^{(192)}\le1/100\).  The normalized physical time \(s_*\)
and actual stage time \(t_*=(j!)^{-35/2}s_*\) therefore obey

\[
 TR_*^{(192)}\le s_*\le{100\over99}TR_*^{(192)},
 \qquad e^{-2\mu_js_*}\ge{99\over100},
\tag{5.4}
\]

\[
 (j!)^{-35/2}TR_*^{(192)}
 \le t_*\le
 {100\over99}(j!)^{-35/2}TR_*^{(192)}.
\tag{5.5}
\]

The separate normalized viscous-collar line becomes

\[
 \boxed{
 C_{\nu{\rm col}}\nu23^{3/2}(j!)^{-3/2}n^{437/28}.}
\tag{5.6}
\]

Equations (5.2)--(5.6) remain conditional scalar ledgers.  They do not
prove finite-band viscous persistence or actual lower stage coverage.

## 6. Exact finite-frequency threshold exposed by C192

The shorter raw clock changes the power test for a future explicit
parametrix. From (4.4), suppose a retained packet at carrier \(q\) has a
remainder, measured in entrance units, bounded by

\[
 C_0(1+t)^Pq^{-1/2}e^{\Gamma t}
\tag{6.1}
\]

through \(t=TR_\Delta^{(192)}\), where \(C_0,P,\Gamma\) are independent
of \(q\). Relative to the required \(q^{3/8}\)
gain, (4.4) gives

\[
 C_0e^{76\Gamma/25}(1+t)^P
 q^{-7/8+(57/400)\Gamma}.
\tag{6.2}
\]

Thus the power closes whenever

\[
                         \boxed{\Gamma<{350\over57}}.
\tag{6.3}
\]

If an order-\(-1\) pseudodifferential remainder gives \(q^{-1}\) instead,
the corresponding exact threshold is

\[
                         \boxed{\Gamma<{550\over57}}.
\tag{6.4}
\]

These are conditional arithmetic implications, not estimates for the
actual remainder.  C176 makes the flow jets polynomial on this orbit, and
the value \(6\) lies strictly below \(350/57\), but no landed theorem yet
bounds every symbol/pressure remainder with exponential rate at most six.
That is now the finite-frequency constant that must be proved; no new
named gate is introduced.

## 7. What C192 changes and what remains open

C192 changes three facts quantitatively:

1. the certified one-period principal multiplier is \(>3000\), not merely
   \(>e^{1/5}\);
2. the raw \(q^{3/8}\) action coefficient falls from \(57/10\) to
   \(57/400\); and
3. a first-order \(q^{-1/2}\) packet remainder now has a nonempty explicit
   survival window, \(\Gamma<350/57\).

It does not prove that the unrestricted C192 operator-norm vector is the
C180/C182 packet, that a finite band follows the ray, that viscosity and
projection retain it, that the final \(L^\infty/L^2\) quotient increases,
or that C125, depletion, wake, UVSR, terminal singular-center tracking, or
a Navier--Stokes singularity closes.  The C185 citation remains restricted
to operator norm until the independent Shvydkoy paper-body check lands.

## 8. Verification boundary

The checker:

1. reruns C159's complete outward phase/covector certificate;
2. reevaluates raw \(B_{ij}\) intervals on 2048 cells;
3. proves the positive comparison (2.3) with directed rounding through
   degree six;
4. certifies both component ratios above \(3000\);
5. proves \(e^8<3000\) with rational arithmetic; and
6. verifies the exact clock, factor-forty, collar, heat, and conditional
   finite-frequency arithmetic in (4.2)--(6.4).

It does not check the external spectral-inclusion theorem, construct a
finite-frequency packet, certify the conditional remainder exponent in
(6.1), prove viscosity, or establish a nonlinear stage.
