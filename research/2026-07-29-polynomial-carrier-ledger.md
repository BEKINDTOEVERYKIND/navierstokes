# Polynomial-carrier refinement: flat residuals without an exponential seed gap

Date: 2026-07-29

## Claim boundary

This note does **not** construct the localized return cell.  It removes one
unnecessary difficulty from the conditional flat-force program.

The internal geometric-optics carrier need not grow exponentially with the
outer cascade.  A polynomial carrier

\[
 K_j=j^A
\]

combined with a much deeper, Gevrey-controlled truncation

\[
 M_j\asymp \frac{j^2}{\log j}
\]

still makes the normalized residual \(e^{-c j^2}\).  At the same time,

\[
 \frac{K_{j+1}}{K_j}=1+\frac Aj+O(j^{-2}),
\]

so the extra carrier handoff becomes a slowly varying modulation rather than
a second fixed-ratio cascade.

This is conditional on an all-order packet construction whose coefficients
obey a Gevrey bound.  It does not prove that the coefficient equations are
solvable.  In particular, the carrier cannot be a stationary Euler wave
with a stationary viscous corrector: the linearized Euler energy identity
puts \(\Delta W\) outside that static range.  Slow energy decay and
time-dependent endpoint jets must be part of the hierarchy.

## 1. Modified physical ledger

Let

\[
 \ell_j=r^{-j},\qquad K_j=j^A,\qquad
 a_j=\ell_j^{-\gamma}K_j^\gamma,
 \qquad 1<\gamma<\frac32.
\]

The amplitude is chosen this way because the exact affine Kelvin amplifier
obeys \(a\propto k^\gamma\), and the physical carrier is

\[
 k_j=\frac{K_j}{\ell_j}.
\]

The required one-stage carrier compression is therefore

\[
 q_j:=\frac{k_{j+1}}{k_j}
 =r\left(\frac{j+1}{j}\right)^A
 =r\left(1+\frac Aj+O(j^{-2})\right).
\]

An affine Kelvin interval with compression \(q_j\) gives exactly

\[
 \frac{a_{j+1}}{a_j}=q_j^\gamma.
\]

Thus the polynomial factor is compatible with the exact amplifier law.  It
does require the localized endpoint map to tolerate a nonautonomous
\(O(j^{-1})\) modulation.

The stage time is

\[
 \tau_j=\frac{\ell_j}{a_j}
 =\ell_j^{1+\gamma}K_j^{-\gamma}.
\]

The standard blow-up and energy quantities become

\[
\begin{aligned}
 a_j\ell_j
 &=\ell_j^{1-\gamma}K_j^\gamma\longrightarrow\infty,\\
 E_j
 &\asymp a_j^2\ell_j^3
 =\ell_j^{3-2\gamma}K_j^{2\gamma},\\
 \operatorname{Re}_j
 &=\frac{a_j\ell_j}{\nu}
 =\nu^{-1}\ell_j^{1-\gamma}K_j^\gamma
 \longrightarrow\infty.
\end{aligned}
\]

Polynomial factors do not affect summability against the exponential
\(\ell_j\)-powers.  Hence

\[
 \sum_j\tau_j<\infty,\qquad
 \sum_jE_j<\infty
\]

in the same sharp window \(1<\gamma<3/2\).

The normalized viscosity is

\[
 \varepsilon_j=\frac{\nu}{a_j\ell_j}
 =\nu\ell_j^{\gamma-1}K_j^{-\gamma}.
\]

At the carrier, the heat action during one normalized stage is

\[
 \theta_j=\varepsilon_jK_j^2
 =\nu\ell_j^{\gamma-1}K_j^{2-\gamma}.
\]

It is exponentially small in \(j\), despite the polynomial carrier.  Even
if an order-\(M_j\) construction creates harmonics \(nK_j\) for
\(n\leq C M_j\), the largest heat parameter is bounded by
\(C M_j^2\theta_j\to0\).

The carrier-dominated stage dissipation is

\[
 D_j\asymp
 \nu(a_jk_j)^2\ell_j^3\tau_j
 =\nu\ell_j^{2-\gamma}K_j^{\gamma+2}.
\]

It is summable, and

\[
 \frac{D_j}{E_j}
 \asymp\theta_j\longrightarrow0.
\]

## 2. Gevrey truncation gives a flat residual

Assume the order-\(m\) WKB/viscous coefficients satisfy, in each fixed
packet-plus-wake seminorm,

\[
 \|V_m\|\leq C^{m+1}(m!)^\sigma
\]

for some finite Gevrey order \(\sigma\geq1\).  A typical WKB remainder after
order \(M\) is then bounded by

\[
 C^M(M!)^\sigma K_j^{-M}.
\]

Choose an integer \(A>2\sigma\) and

\[
 M_j=\left\lfloor\eta\frac{j^2}{\log(e+j)}\right\rfloor
\]

with fixed \(\eta>0\).  Stirling's formula gives

\[
\begin{aligned}
 \log\!\left(C^{M_j}(M_j!)^\sigma K_j^{-M_j}\right)
 &\leq
 M_j\bigl(\log C+\sigma\log M_j-A\log j\bigr)\\
 &=
 -\eta(A-2\sigma)j^2+o(j^2).
\end{aligned}
\]

Thus, after decreasing the constant,

\[
 C^{M_j}(M_j!)^\sigma K_j^{-M_j}
 \leq e^{-c j^2}.
\]

The viscous expansion is easier.  Since

\[
 \log\theta_j
 =-(\gamma-1)(\log r)j+O(\log j),
\]

one has

\[
 C^{M_j}(M_j!)^\sigma\theta_j^{M_j}
 \leq
 \exp\!\left(-c\frac{j^3}{\log j}\right).
\]

Mixed monomials in \(K_j^{-1}\) and \(\theta_j\) obey the same
\(e^{-c j^2}\) bound, provided the bivariate coefficient estimates are
Gevrey-tame uniformly in the stage.

## 3. Terminal flatness survives physical rescaling

A normalized residual \(R_j\) produces physical force

\[
 f_j=\frac{a_j^2}{\ell_j}R_j.
\]

For fixed spatial and temporal derivative orders, conversion to physical
variables costs powers of

\[
 \ell_j^{-1},\qquad
 \tau_j^{-1},\qquad
 K_j.
\]

Their logarithm is only \(O(j+\log j)\).  Therefore

\[
 \|\partial_x^m\partial_t^n f_j\|_\infty
 \leq \exp(O_{m,n}(j+\log j)-c j^2).
\]

Moreover,

\[
 T-t_j\asymp\tau_j
 =r^{-(1+\gamma)j}j^{-A\gamma}.
\]

It follows that, for every fixed \(m,n,N\),

\[
 \|\partial_x^m\partial_t^n f_j\|_\infty
 =O((T-t_j)^N).
\]

Thus the residual and any equally accurate seam force extend
\(C^\infty\)-flatly by zero at \(T\).

## 4. What this changes

With an exponential carrier \(K_j=\ell_j^{-\kappa}\), the next stage asks
for a fixed additional internal frequency ratio.  The exact Kelvin
amplifier naturally preserves a fixed internal carrier under exact scaled
return, so the origin of that extra seed is a serious endpoint problem.

For \(K_j=j^A\), the extra factor is only

\[
 \frac{K_{j+1}}{K_j}=1+O(j^{-1}).
\]

The amplitude law acquires only the subexponential factor \(j^{A\gamma}\),
which leaves every energy, time, Reynolds, and force-flatness inequality
unchanged at the exponent level.

The corrected target is therefore a slowly nonautonomous invariant graph
or endpoint submersion, rather than a stage-independent fixed point with a
separately generated exponential carrier.  The \(O(j^{-1})\) modulation
must be incorporated **exactly** (or with an independently flat error) in
each one-step map.  A merely finite Gevrey expansion in \(j^{-1}\) is not
flat: its optimal remainder is generally only exponential in \(j\).

This does **not** make the handoff automatic.  The phase-resolved return map
must still:

1. compress the selected carrier by
   \(q_j=r(1+A/j+O(j^{-2}))\);
2. land on the next localized parent/seed state while retaining the
   work-carrying wake;
3. solve the WKB and viscous coefficient equations through order
   \(M_j\asymp j^2/\log j\);
4. satisfy Gevrey bounds with no repeated derivative loss and no small
   divisors growing faster than the chosen \(j^A\) budget.

The last requirement is now the decisive analytic gate.  A GPU
optimization of a low-order carrier cannot establish it.

There are two further non-scaling compatibilities.

* The exact energy identity requires
  \(E'=-\nu\|\nabla u\|_2^2+O((T-t)^\infty)\); the outgoing wake must retain
  the remaining energy while allowing the summable viscous loss.
* A bare steady Euler packet cannot occupy a seam collar, because its
  Navier--Stokes residual is \(-\nu\Delta U\), which is not flat.  The
  matched state must contain the full time-dependent viscous jet or an
  active carrier bath passed to the next wake.

## 5. Sharpened conditional theorem

For the forced Clay alternative, it would suffice to construct a localized,
phase-resolved Kelvin--Reynolds packet-plus-wake stage whose nonautonomous
return maps admit:

* a uniformly tame invariant sequence or graph evaluated at the exact
  compression
  \(q_j=r(K_{j+1}/K_j)\);
* an asymptotic expansion in the genuinely small parameters
  \(K_j^{-1}\) and \(\theta_j\);
* coefficient bounds \(C^m(m!)^\sigma\) for some finite \(\sigma\);
* uniform solvability through
  \(M_j\asymp j^2/\log j\);
* exact matching of the retained wake, or a wake coupling already bounded
  by \(e^{-c j^2}\).

Under those hypotheses, the ledger above converts the truncated equation
and seam errors into a smooth periodic force satisfying the decay
conditions in Clay alternative (D).  The missing result remains the
localized return theorem itself, not the scalar scaling.

## Primary sources

* C. L. Fefferman, *Existence and smoothness of the Navier--Stokes
  equation*, official Clay problem description:
  https://www.claymath.org/wp-content/uploads/2022/06/navierstokes.pdf
* A. D. D. Craik and W. O. Criminale, *Evolution of wavelike disturbances
  in shear flows: a class of exact solutions of the Navier--Stokes
  equations*, Proc. R. Soc. Lond. A 406 (1986), 13--26:
  https://doi.org/10.1098/rspa.1986.0061
* S. Daneri and L. Székelyhidi Jr., *Non-uniqueness and h-principle for
  Hölder-continuous weak solutions of the Euler equations*:
  https://arxiv.org/abs/1603.09714
