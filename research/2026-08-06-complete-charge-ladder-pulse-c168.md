# C168: the complete charge ladder retains a fixed active pulse gain

**Date:** 2026-08-06
**Status:** independently audited exact complete-ladder prescribed-pump
linear theorem and nonzero-radial-pulse robustness; no unforced or
nonlinear Navier--Stokes stage
**Checker:**
[checks/complete_charge_ladder_pulse_c168.py](../checks/complete_charge_ladder_pulse_c168.py)

## 0. Claim boundary

C166's source/first-daughter space is not invariant under its radial
cosine gate. This note removes that truncation. On the same horizontal
sheet, the radial gate is an explicit bounded nearest-neighbour
(Jacobi-type) operator on the **entire** integer normal-charge ladder.
The subsequent tangential gate is an exact square-zero operator on that
whole ladder. The radial operator is not asserted to be self-adjoint in
the physical coefficient \(\ell^2\) metric.

Two positive statements survive.

1. With zero radial action, a unit tangential pulse has the exact
   reality-completed perturbation point gain

   \[
      G_{\rm point}=\sqrt2,
   \]

   and half-lattice perturbation-coefficient energy \(3/2\).
2. This is not an artifact of setting the first pulse to zero. A radial
   pulse of normalized action \(1/100\), followed by the same unit
   tangential pulse, includes every repeated charge path and still has a
   rigorously certified perturbation point gain greater than \(13/10\).
   Its half-lattice perturbation-coefficient energy is less than \(3\),
   independently of every cascade parameter.

The mechanism is an active prescribed-pump stretch, not a conservative
charge-star focus. It gives only a fixed constant, uses a switched pump,
and lives in an \(x_t\)-independent 2D3C class. In the full nonlinear
2D3C equations the total tangential velocity is a passive scalar and
obeys the viscous maximum principle. Thus the calculation does not evade
the autonomous nonlinear resource problem.

In particular, no unforced/heat-compatible switching, \(q^2\) source
sheet, \(q\)-way normalization, localization, pressure estimate, BAFL
bound, or one-cell Navier--Stokes stage is claimed.

## 1. One horizontal sheet and all normal charges

Use the C164--C166 frame and write

\[
 \xi=Ae_r\cdot x,\qquad
 \theta={A\over\sqrt3}e_z\cdot x,\qquad
 k_m=A\left(e_r+{m\over\sqrt3}e_z\right),
 \quad m\in\mathbb Z.                              \tag{1.1}
\]

An orthonormal basis of \(k_m^\perp\) is

\[
 \rho_m={m e_r-\sqrt3e_z\over\sqrt{m^2+3}},
 \qquad e_t.                                       \tag{1.2}
\]

The complex positive-horizontal sheet and its reality completion are

\[
 V[a,b](\xi,\theta)
 =e^{i\xi}\sum_{m\in\mathbb Z}
       (a_m\rho_m+b_m e_t)e^{im\theta},
 \qquad v=2\mathop{\rm Re}V.                       \tag{1.3}
\]

The C166 source is \(m=1\), so \(\rho_1=e_\sigma\). Its daughter
orientations have the opposite sign:

\[
 \rho_2=-e_{+,\perp},\qquad \rho_0=-e_{-,\perp}.   \tag{1.4}
\]

Let the real gate have signed Fourier modes

\[
 \left(\pm {A\over\sqrt3}e_z,E\right),             \tag{1.5}
\]

each with coefficient \(E\). Its spatial velocity is
\(2\cos\theta\,E\), exactly the signed-pair normalization of C166. For
the symmetric projected Euler symbol

\[
 {\cal S}(k,a;q,E)
 =P_{k+q}\{(a\cdot q)E+(E\cdot k)a\},              \tag{1.6}
\]

normalized time is gate amplitude times \(At\).

## 2. Exact full radial nearest-neighbour operator

First take \(E=e_r\). For a source \((k_m,\rho_m)\) and signed gate
\(q=\sigma Ae_z/\sqrt3\), \(\sigma\in\{+1,-1\}\), the unprojected vector
in (1.6), divided by \(A\), is

\[
 { (m-\sigma)e_r-\sqrt3e_z\over\sqrt{m^2+3}}.      \tag{2.1}
\]

Its coefficient along \(\rho_{m+\sigma}\) is

\[
 \boxed{
 j_{m+\sigma,m}
 ={m^2+2\over
   \sqrt{m^2+3}\sqrt{(m+\sigma)^2+3}}.}            \tag{2.2}
\]

For a tangential source, the corresponding coefficient is one.
Consequently the complete radial operator on
\(\ell^2(\mathbb Z)\oplus\ell^2(\mathbb Z)\) is

\[
 H_r=\begin{pmatrix}J&0\\0&L\end{pmatrix},         \tag{2.3}
\]

where

\[
 \begin{aligned}
 (Ja)_n={}&{(n-1)^2+2\over
  \sqrt{(n-1)^2+3}\sqrt{n^2+3}}a_{n-1}\\
 &+{(n+1)^2+2\over
  \sqrt{(n+1)^2+3}\sqrt{n^2+3}}a_{n+1},\\
 (Lb)_n={}&b_{n-1}+b_{n+1}.                        \tag{2.4}
 \end{aligned}
\]

The term ``Jacobi-type'' means tridiagonal here, not symmetric in the
physical coefficient metric. Indeed,

\[
 {j_{m+1,m}\over j_{m,m+1}}
 ={m^2+2\over (m+1)^2+2}.                         \tag{2.4a}
\]

The detailed-balance weight \(w_m=m^2+2\) symmetrizes \(J\), but it is
unbounded and is not uniformly equivalent to coefficient \(\ell^2\).
No physical energy conservation is inferred from that formal weighted
symmetry.

Starting from the C166 source, the exact radial endpoint is

\[
 a(\tau)=e^{-i\tau J}\delta_1,\qquad b(\tau)=0.    \tag{2.5}
\]

This contains every walk on the integer ladder, including repeated
returns. The exponential is a bounded-operator exponential on both
\(\ell^1(\mathbb Z)\) and \(\ell^2(\mathbb Z)\); the bounds needed below
are proved in Section 4. The normalization against C166 is exact:

\[
 \begin{array}{c|cccc}
 \text{edge}&1\to2&1\to0&2\to1&0\to1\\ \hline
 j&\dfrac{3\sqrt7}{14}&\dfrac{\sqrt3}{2}
   &\dfrac{3\sqrt7}{7}&\dfrac{\sqrt3}{3}\\[1mm]
 \text{C166 coordinate}&-a_+&-a_-&-c_+&-c_-
 \end{array}.                                      \tag{2.6}
\]

The entries in the first row of (2.6) are positive in the \(\rho_m\)
basis. The C166 coordinates are their negatives because both of C166's
daughter radial vectors are \(-\rho_m\). Thus the signs in the second
row come solely from (1.4); no one-sided gate or doubled coefficient has
entered.

## 3. Exact full tangential square-zero pulse

Switch the prescribed gate polarization to \(E=e_t\). A tangential
source is dark, while a radial-plane source obeys

\[
 A^{-1}{\cal S}
 \left(k_m,\rho_m;{\sigma A\over\sqrt3}e_z,e_t\right)
 =-{\sigma\over\sqrt{m^2+3}}e_t.                  \tag{3.1}
\]

Hence

\[
 H_t=\begin{pmatrix}0&0\\K&0\end{pmatrix},
 \qquad
 (Ka)_n=-{a_{n-1}\over\sqrt{(n-1)^2+3}}
        +{a_{n+1}\over\sqrt{(n+1)^2+3}},           \tag{3.2}
\]

and, on the full ladder,

\[
 H_t^2=0,\qquad
 e^{-isH_t}(a,b)=(a,b-isKa).                       \tag{3.3}
\]

With the Fourier evolution convention
\(\partial_t(a,b)=-iA H_\bullet(a,b)\), the complete
radial-then-tangential endpoint is

\[
 \boxed{
 a=e^{-i\tau J}\delta_1,\qquad b=-isK a.}          \tag{3.4}
\]

For \(\tau=0\),

\[
 K\delta_1={1\over2}\delta_0-{1\over2}\delta_2.    \tag{3.5}
\]

At \(\theta=\pi/2,\ \xi=-\pi/2\), the unit pulse \(s=1\) gives

\[
 V[\delta_1,-iK\delta_1]=\rho_1-e_t,               \tag{3.6}
\]

which is real. The reality-completed perturbation has size \(2\sqrt2\)
there. The initial reality-completed source perturbation has supremum
\(2\), proving

\[
 \boxed{G_{\rm point}=\sqrt2.}                    \tag{3.7}
\]

The half-lattice perturbation-coefficient energy is
\[
 1+{1\over4}+{1\over4}={3\over2}.                 \tag{3.8}
\]

## 4. Nonzero radial pulse: an all-path certificate

For adjacent integers \(n=m\pm1\), (2.2) satisfies

\[
 j_{n,m}<\sqrt{m^2+3\over n^2+3}<{3\over2}.       \tag{4.1}
\]

The last inequality follows for both signs from

\[
 9\big((m\pm1)^2+3\big)-4(m^2+3)
 =5\left(m\pm{9\over5}\right)^2+{39\over5}>0.     \tag{4.2}
\]

Thus the row and column sums are bounded by \(3\), and Schur's bound
gives

\[
 \|J\|_{\ell^1\to\ell^1}\le3,\qquad
 \|J\|_{\ell^\infty\to\ell^\infty}\le3,\qquad
 \|J\|_{\ell^2\to\ell^2}\le3.                  \tag{4.3}
\]

Every entry of \(K\) is at most \(1/\sqrt3<3/5\), and each row and
column has two entries, so

\[
 \|K\|_{\ell^1\to\ell^1},
 \|K\|_{\ell^2\to\ell^2}<{6\over5}.               \tag{4.4}
\]

Put \(\tau=1/100,\ s=1\), and

\[
 h=e^{-i\tau J}\delta_1-\delta_1.                 \tag{4.5}
\]

The complete bounded-operator exponential series, not a finite section,
gives

\[
 \|h\|_1\le e^{3/100}-1
 < {1\over1-3/100}-1
 ={3\over97}<{1\over32}.                          \tag{4.6}
\]

This geometric majorant sums all positive-length walks, including every
repeated charge. In particular it gives the literal off-source tail

\[
 \sum_{m\ne1}|a_m|\le\|h\|_1<{1\over32}.          \tag{4.7}
\]

At the point in (3.6), the complex endpoint differs from
\(\rho_1-e_t\) by at most

\[
 \|h\|_1+\|Kh\|_1
 <\left(1+{6\over5}\right){1\over32}
 ={11\over160}.                                    \tag{4.8}
\]

Taking real parts cannot enlarge the error. Since

\[
 \sqrt2>{7\over5},\qquad
 {7\over5}-{11\over160}={213\over160}>{13\over10}, \tag{4.9}
\]

the reality-completed perturbation has size greater than \(13/5\).
Relative to the initial perturbation supremum \(2\),

\[
 \boxed{G_{\rm point}>{13\over10}}                \tag{4.10}
\]

for the genuine nonzero-radial circuit. Its normalized action is
\(101/100\). Moreover,

\[
 \begin{aligned}
 \|a\|_2&<e^{3/100}<{33\over32},\\
 \|a\|_2^2+\|Ka\|_2^2
 &<\left(1+{36\over25}\right)\left({33\over32}\right)^2
 ={66429\over25600}<3.                            \tag{4.11}
 \end{aligned}
\]

The full reality-completed perturbation-coefficient energy is therefore
less than \(6\).

## 5. Why this does not supply the one-cell theorem

The gain is powered by a prescribed tangential gate. It is an orthogonal
velocity component added by a square-zero linear response, not the
\(q\)-way, approximately energy-preserving concentration required by
C161. It remains a fixed constant and may draw energy from the pump.

There is also a nonlinear structural warning. Every field here is
independent of \(x_t\). For an exact 2D3C Navier--Stokes solution write
\(u=u_{rz}+w e_t\). Then \(u_{rz}\) solves the two-dimensional
Navier--Stokes equations and

\[
 \partial_t w+u_{rz}\cdot\nabla w=\nu\Delta w,
 \qquad
 \|w(t)\|_\infty\le\|w(0)\|_\infty.               \tag{5.1}
\]

The linear term in (3.1) is the derivative of this passive-scalar
transport with respect to an in-plane perturbation of a prescribed
tangential background. It does not make the total nonlinear tangential
field an autonomous amplifier. Any stage use must leave this closed
2D3C setting or find a reservoir exchange compatible with the maximum
principle.

The remaining obligations are unchanged: unforced decaying-pump
implementation, the \(q^2\)-source/\(q\)-daughter normalization, a
coherent \(\sqrt q\) point focus, localization and pressure control, and
the BAFL endpoint estimate. C168 proves only that C166's fixed active
gain is not destroyed by the complete integer charge ladder.
