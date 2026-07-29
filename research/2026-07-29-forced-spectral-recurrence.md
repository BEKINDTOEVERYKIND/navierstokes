# Forced spectral recurrence: a surviving prize-level route

Date: 2026-07-29

## Status

No solution of the Millennium problem is claimed here.

This note records a kill-or-advance result:

1. The direct combination of Palasek's smooth dormant forcing with the
   hyperbolically transported vortex layers of
   Córdoba–Martínez-Zoroa–Zheng (CMZ) cannot reach the ordinary Laplacian.
   The obstruction is a principal heat cost, not a fixed-order error.
2. A bounded-ratio cascade driven by a shape-preserving spectral instability
   survives the energy, viscosity, time, smooth-data, smooth-force, and
   critical-norm ledgers.
3. That surviving architecture reduces the problem to one nonlinear
   renormalized return-map theorem. This is a genuine existence problem, not
   another parameter optimization.

## 1. Why a smooth force is prize-relevant

Fefferman's official problem statement asks for any one of four alternatives.
Alternatives (C) and (D) permit a prescribed smooth force. In the periodic
case (D), the data are spatially smooth and periodic and every space-time
derivative of the force must decay faster than any power of time.

Thus a true forced Navier–Stokes blowup with terminally flat forcing would
resolve the problem in the breakdown direction. The periodic target avoids
the extra spatial-tail conditions of the whole-space formulation.

Source:
[official Clay problem description](https://www.claymath.org/wp-content/uploads/2022/06/navierstokes.pdf).

## 2. The superexponential shell window is formally real

Palasek's 2026 forced Obukhov model uses

\[
N_k=N_0^{b^k},\qquad A_k=N_k^\beta
\]

and a three-dimensional intermittency parameter \(a\le 5/2\). For the full
Laplacian, the formal window is

\[
1<b<\frac a2,\qquad 2b<\beta<a\le\frac52. \tag{2.1}
\]

Indeed, the parent strain is \(A_{k-1}=N_k^{\beta/b}\), so

\[
\frac{\nu N_k^2}{A_{k-1}}
=\nu N_k^{2-\beta/b}\longrightarrow0
\]

exactly when \(\beta>2b\). The unfavorable high-high feedback has the shell
size

\[
N_k^{-2(a-\beta)(b-1)}\longrightarrow0
\]

when \(\beta<a\).

The dormant seed and its viscosity-cancelling force carry the factor

\[
\exp\!\left[-c\frac{A_{k-1}}{A_{k-2}}\right]
=\exp(-cN_k^\theta),\qquad
\theta=\frac{\beta(b-1)}{b^2}>0, \tag{2.2}
\]

which defeats every polynomial spatial or temporal derivative. This is why
Palasek's force is smooth and flat even though it cancels a large relative
heat rate.

Source:
[Palasek, *Finite-time blow-up in an elementary model of the 3D
Navier–Stokes equations*](https://arxiv.org/abs/2605.13827).

## 3. Why the direct CMZ retrofit fails

CMZ construct a forced blowup for fractional dissipation
\(|\nabla|^d\) only when

\[
d<\frac{22-8\sqrt7}{9}\approx0.0926655.
\]

Their published ansatz locks the layer frequency and vorticity amplitude:

\[
M_k=N^{R^k},\qquad |\omega_k|\sim M_k^q,\qquad
q=\sqrt{\frac{2d}{7}},\quad R=\sqrt{\frac{2}{7d}}.
\]

Because \(q/R=d\), the parent strain acting on layer \(k\) is exactly
\(M_k^d\). At \(d=2\), this locking would force \(R=1/\sqrt7<1\), so it
cannot even define an increasing frequency ladder.

One might try to decouple the scales using (2.1) and remove fixed-order
errors by a WKB hierarchy of increasing order. That does reopen the formal
exponent window, but it does not fix the principal deformation.

The Palasek seed must gain at least

\[
K_k\gtrsim
\exp\!\left(c\frac{A_{k-1}}{A_{k-2}}\right)
=\exp(cN_k^\theta). \tag{3.1}
\]

For the CMZ hyperbolic pullback, the largest support width and a cutoff
frequency scale as

\[
\frac{K_k}{L_k},\qquad K_kL_k.
\]

Remaining inside the affine core of the parent requires

\[
\frac{K_k}{L_k}\lesssim N_{k-1}^{-1}. \tag{3.2}
\]

Keeping the deformed cutoff below the nominal carrier requires

\[
K_kL_k\lesssim N_k. \tag{3.3}
\]

Together, (3.2)–(3.3) imply

\[
K_k^2N_{k-1}\lesssim N_k, \tag{3.4}
\]

which permits only polynomial \(K_k\), contradicting (3.1).

Even if carrier dominance is abandoned, the ordinary Laplacian sees the
principal active frequency \(K_kL_k\). From (3.2),
\(L_k\gtrsim K_kN_{k-1}\), so heat domination by the parent would require

\[
(K_kL_k)^2\ll A_{k-1}
\quad\Longrightarrow\quad
K_k^4N_{k-1}^2\ll N_{k-1}^{\beta}. \tag{3.5}
\]

Again \(K_k\) must be polynomial. Leaving the dormancy force switched on
does not help: its absolute size then contains

\[
|\Delta\omega_k|
\gtrsim (K_kL_k)^2\frac{A_k}{K_k}
=A_kK_kL_k^2,
\]

which is not flat.

This is a principal-symbol obstruction. Increasing WKB order only reduces
lower-order residuals and cannot remove (3.5).

Source:
[CMZ, *Finite time blow-up for the hypodissipative Navier Stokes
equations...*](https://arxiv.org/abs/2407.06776).

## 4. What changes with a spectral amplifier

A material hyperbolic amplifier gains amplitude by deforming its wavevector.
A genuine unstable eigenmode or Floquet mode can instead gain amplitude while
returning to the same spatial profile after every normalized cycle.

Albritton and Ożański rigorously construct helical vortex-column ring modes
whose inviscid growth rate remains \(O(1)\) as their azimuthal and axial
wavenumbers tend to infinity at fixed ratio. This supplies a real
scale-uniform amplification mechanism, although not the compact recurrent
cell needed below.

Source:
[Albritton–Ożański, *Linear and nonlinear instability of vortex
columns*](https://arxiv.org/abs/2310.20674v3).

The distinction is decisive: if the active physical frequency remains
\(N_j\), the heat exponent during \(M_j\) normalized cycles is proportional
to

\[
M_j\,\frac{\nu N_j^2}{S_j},
\]

not \(M_j\nu K_j^2N_j^2/S_j\).

## 5. A bounded-ratio dormant spectral ladder

Let

\[
a_j=a_0r^{-j},\qquad
S_j=S_0\mu^j,\qquad
\eta_j=\eta_0q^{-j}, \tag{5.1}
\]

where \(a_j\) is core radius, \(S_j\) is strain, and
\(\eta_j=a_j/R_j\) is toroidal slenderness. Use a smooth preloaded seed

\[
\delta_j=\exp(-cj^2). \tag{5.2}
\]

It is \(C^\infty\) in space because every spatial derivative costs only an
exponential in \(j\), while (5.2) is Gaussian in \(j\). A dormancy force
may cancel waiting heat while the seed is tiny. If the activation times
satisfy

\[
T-t_j\asymp j^A\mu^{-j},
\]

then \(j\asymp|\log(T-t)|\), and

\[
e^{-cj^2+Cj}
=O((T-t)^K)
\]

for every \(K\). The dormancy force and all of its derivatives therefore
extend smoothly and flatly through \(T\).

If the unstable exponent per normalized cycle is \(\sigma>0\), reaching a
fixed endpoint section costs only

\[
M_j\sim\frac{c}{\sigma}j^2. \tag{5.3}
\]

The complete toroidal ledger, including the increasing number of core
lengths \(1/\eta_j\), is

\[
\begin{aligned}
\frac{E_{j+1}}{E_j}
 &=q\frac{\mu^2}{r^5},&
\frac{D_{j+1}}{D_j}
 &=q\frac{\mu}{r^3},\\
\frac{\tau_{j+1}}{\tau_j}
 &\sim\frac1\mu,&
\frac{\|u_{j+1}\|_3}{\|u_j\|_3}
 &=q^{1/3}\frac{\mu}{r^2},\\
\frac{\mathrm{Re}_{j+1}}{\mathrm{Re}_j}
 &=\frac{\mu}{r^2},&
\frac{R_{j+1}}{R_j}
 &=\frac qr.
\end{aligned} \tag{5.4}
\]

There is a strict nonempty window

\[
r^2<\mu<r^{5/2},\qquad
1<q<\frac{r^5}{\mu^2}. \tag{5.5}
\]

In this window:

- energy and total dissipation are summable, even after the polynomial
  residence factor \(M_j\);
- the total physical time is finite;
- major radius and volume shrink;
- Reynolds number and \(L^3\) grow;
- accumulated waiting-heat and finite-time viscous errors are summable.

An exact dyadic design point is

\[
r=32,\qquad \mu=2048,\qquad q=2,
\]

for which

\[
\frac{E_{j+1}}{E_j}=\frac14,\quad
\frac{D_{j+1}}{D_j}=\frac18,\quad
\frac{\mathrm{Re}_{j+1}}{\mathrm{Re}_j}=2,\quad
\frac{R_{j+1}}{R_j}=\frac1{16}.
\]

This is an algebraic pass, not a construction.

## 6. The exact remaining theorem

Let \(\Phi_\chi^t\) denote the normalized Navier–Stokes flow at local
Reynolds parameter \(\chi\), and let \(\mathcal S_{r,\mu,q}\) rescale
frequency, strain, slenderness, translation, rotation, and phase.

The target is a compact profile \(P\), an unstable direction \(H\), a
terminal section \(\Sigma\), and a trapping neighborhood \(\mathcal U\)
such that:

1. **Shape-preserving gain.** The linearized flow about \(P\) has an
   unstable eigenline or Floquet line \(H\) with exponent
   \(\sigma_\chi\to\sigma_\infty>0\), while all active frequencies remain
   polynomial in the cell frequency during an arbitrarily long residence.
2. **Nonlinear return.** For every sufficiently small
   \(\delta>0\), the orbit from \(P+\delta H\) reaches \(\Sigma\), and
   \[
   \mathcal S_{r,\mu,q}
   \Phi_\chi^{T(\delta)}(P+\delta H)\in\mathcal U.
   \]
   The rescaled output contains the next parent and a nonzero projection
   onto its unstable direction.
3. **Uniform trapping.** After quotienting the geometric symmetries, the
   return maps satisfy a graph-transform estimate uniformly as
   \(\chi\to\infty\).
4. **Energy and sidebands.** Parent drain, wakes, pressure tails, helicity,
   and every sideband obey the ratios in (5.4), with strict margin.
5. **Exact or all-order closure.** Any residual assigned to the external
   force is terminally flat in every \(C^m\) norm. A fixed algebraic error
   per generation is not enough.

This is the **renormalized unstable-manifold return theorem**. Proving it
would convert the ledger into a Clay alternative-(D) breakdown
construction. Failing it for a broad enough compact profile class would
also be valuable because it would remove the last currently
exponent-compatible forced recurrence identified here.

The newest complete high-to-low cascade construction shows that
high-high stress programming is not the wholly unknown part; however, its
solutions inject energy instantaneously from infinite wavenumber and are
not Leray-class finite-time blowups. It does not solve the Clay problem.

Source:
[Cheskidov–Dai–Palasek, *Instantaneous Type I blow-up and non-uniqueness
of smooth solutions of the Navier–Stokes
equations*](https://arxiv.org/abs/2511.09556).

## 7. Computational decision

Do not run another broad DNS or optimize a scalar gain. The next useful
computation must approximate the renormalized return map in Section 6 and
measure all of the following simultaneously:

- unstable projection at the daughter;
- normalized profile distance after rescaling;
- old-parent drain;
- total off-dictionary energy;
- frequency inflation during the long residence;
- convergence under resolution and box-size refinement.

A GPU run is justified only after specifying a compact candidate cell and a
fixed output section. Without those, a larger simulation would be
incremental evidence rather than a test of the missing theorem.
