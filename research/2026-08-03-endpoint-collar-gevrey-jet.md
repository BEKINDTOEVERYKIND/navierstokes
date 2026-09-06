# Gevrey-2 majorant for the exact endpoint-collar jet

**Date:** 2026-08-03

**Status:** exact bidegree and factorial ledger; differential-polynomial
Gevrey majorant self-derived.  A spatially localized transition realizing
the prescribed jet remains open.

**Scope:** the growing-order caveat in C87.  This note estimates the formal
forward Navier--Stokes jet of the exact endpoint collar.  It does not splice
the collar to the active packet or prove convergence of the complete
cascade.

## 1. Outcome

The apparent time-Gevrey loss from the Laplacian is paid by the small
normalized viscosity.  If a collar is expanded through time order \(M\),
the relevant gate is
\[
                         C\mu M^2<1,                            \tag{1.1}
\]
where \(\mu=\operatorname{Re}^{-1}\).  Under (1.1), all endpoint jets
through \(M\) fit a single Gevrey-2 majorant
\[
                         C^{n+1}(n!)^2.                         \tag{1.2}
\]

For the cascade schedule, \(\mu_j\) is exponentially small while
\(M_j\asymp j^2/\log j\), so (1.1) holds with room to spare.  In fact the
existing carrier heat gate \(\mu_j(K_jM_j)^2\to0\) is strictly stronger.

Thus the exact Navier--Stokes collar does not introduce a new scalar
factorial obstruction.  The remaining endpoint problem is geometric and
operator-valued: the three-phase inverse must realize this prescribed jet
while retaining positivity, phase transport, and the global wake.

## 2. Viscosity-degree decomposition

Use the dimensionless projected equation
\[
 \partial_\tau U+\mathbb P(U\cdot\nabla U)=\mu\Delta U,
 \qquad U(0)=V,
 \qquad \mathbb P(V\cdot\nabla V)=0.                           \tag{2.1}
\]
Write
\[
 U_n=\partial_\tau^nU(0)
     =\sum_{q=0}^n\mu^qU_n^{[q]}.                              \tag{2.2}
\]
Then
\[
 U_0^{[0]}=V,
 \qquad U_n^{[0]}=0\quad(n\ge1),                              \tag{2.3}
\]
and the exact recurrence is
\[
\begin{aligned}
 U_{n+1}^{[q]}
 &=\Delta U_n^{[q-1]}\\
 &\quad-\mathbb P\sum_{a=0}^n{n\choose a}
       \sum_{q_1+q_2=q}
       U_a^{[q_1]}\cdot\nabla U_{n-a}^{[q_2]},                \tag{2.4}
\end{aligned}
\]
with the first term omitted when \(q=0\).

## 3. Exact differential bidegrees

Every differential monomial in \(U_n^{[q]}\), for \(1\le q\le n\), has
\[
 \boxed{
 \text{total spatial derivative order }D=n+q,
 \qquad
 \text{number of seed factors }L=n-q+1.}                       \tag{3.1}
\]
This follows by induction.

* A Laplacian sends
  \((n,q,D,L)\) to \((n+1,q+1,D+2,L)\), preserving (3.1).
* A quadratic product combines orders \(a,n-a\), adds one derivative, and
  sends \(q_1,q_2\) to \(q=q_1+q_2\).  Its derivative order is
  \[
  (a+q_1)+(n-a+q_2)+1=(n+1)+q,
  \]
  and its seed degree is
  \[
  (a-q_1+1)+(n-a-q_2+1)=(n+1)-q+1.
  \]

The pure heat branch is \(q=n\), \(D=2n\), \(L=1\).  The branch with only
one viscous insertion is \(q=1\), \(D=n+1\), and contains the repeated
linearized-Euler transport of \(\Delta V\).

## 4. The factorial gate

Assume the seed has a Gevrey-2 spatial bound in a fixed Sobolev algebra,
\[
 \|\nabla^mV\|_{H^s}
 \le A B^m(m!)^2,
 \qquad s>3/2.                                                  \tag{4.1}
\]
Products, the Leray projector, and the binomial time differentiation in
(2.4) preserve the same factorial class.  The exact time convolution is
\[
 \sum_{a=0}^n{n\choose a}(a!)^2((n-a)!)^2
 \le3(n!)^2.                                                    \tag{4.2}
\]

The only new cost is the \(q\) extra spatial derivatives in (3.1).
For \(1\le q\le n\),
\[
 {((n+q)!)^2\over(n!)^2}
 =\prod_{r=1}^q(n+r)^2
 \le(2n)^{2q}=(4n^2)^q.                                       \tag{4.3}
\]
An induction on (2.4), using (4.2) for each quadratic split and absorbing
the finite tensor contractions into \(C^{n+q}\), therefore gives
\[
 \|U_n^{[q]}\|_{H^s}
 \le C^{n+1}(n!)^2(Cn^2)^q.                                   \tag{4.4}
\]
Here \(C\) depends on the fixed seed and Gevrey norm, but not on
\(n,q,\mu\).

After restoring \(\mu^q\),
\[
\begin{aligned}
 \|U_n\|_{H^s}
 &\le C^{n+1}(n!)^2
       \sum_{q=1}^n(C\mu n^2)^q\\
 &\le C^{n+1}(n!)^2
       {C\mu n^2\over1-C\mu n^2},                             \tag{4.5}
\end{aligned}
\]
whenever \(C\mu n^2<1\).  In particular, if
\[
                         C\mu M^2\le{1\over2},                 \tag{4.6}
\]
then every \(1\le n\le M\) satisfies (1.2), with an additional small
factor \(O(\mu n^2)\) relative to the seed scale.

The pure heat branch illustrates why the small factor matters.  Spatial
Gevrey-2 alone gives
\[
 \|\mu^n\Delta^nV\|
 \lesssim C^n\mu^n((2n)!)^2
 \lesssim C^n(n!)^2(C\mu n^2)^n,                              \tag{4.7}
\]
so its apparent time-Gevrey-4 growth collapses back into the
Gevrey-2 budget under (4.6).

## 5. Cascade comparison

For the polynomial-carrier schedule,
\[
 \mu_j={\nu\over a_j\ell_j}
       =e^{-c_\mu j}\operatorname{poly}(j),
 \qquad
 M_j\asymp{j^2\over\log j}.                                  \tag{5.1}
\]
Hence
\[
                         \mu_jM_j^2\longrightarrow0.           \tag{5.2}
\]
The already-required largest-harmonic heat gate is
\[
                         \mu_j(K_jM_j)^2\longrightarrow0.      \tag{5.3}
\]
Since \(K_j\ge1\), (5.3) implies (5.2).  Matching the exact collar jet
therefore requires no stronger exponent window than the existing viscous
WKB hierarchy.

## 6. Remaining theorem

This ledger supplies a growing-order target with the correct factorial
class.  It does not construct the endpoint inverse.  One still has to:

1. solve the nonzero material charges with terminal data prescribed by
   \(U_n\) rather than zero;
2. absorb the zero charge through the C86 global stress/wake equation;
3. keep the leading covariance uniformly positive during the splice; and
4. prove the resulting endpoint derivative has the robust affine-capture
   rank.

The benefit is that viscosity no longer appears as an unspecified
all-order jet loss.  Its exact bidegree, factorial price, and sufficient
smallness gate are (3.1), (4.3), and (4.6).
