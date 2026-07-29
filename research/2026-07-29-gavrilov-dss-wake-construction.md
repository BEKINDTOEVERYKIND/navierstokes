# An exact discretely self-similar Euler wake from disjoint Gavrilov bubbles

Date: 2026-07-29

## Result and claim boundary

The stationary outer-wake existence problem identified in
[`2026-07-29-log-periodic-wake-equation.md`](2026-07-29-log-periodic-wake-equation.md)
has an explicit solution.

Starting from Gavrilov's nonzero \(C^\infty\), compactly supported steady
Euler flow, put one bubble strictly inside a fundamental spherical annulus
and repeat it by exact dilations.  The copies can have disjoint supports,
so their sum is again an exact stationary Euler flow.  For any
\(\gamma\in\mathbb R\) it obeys

\[
 U(rx)=r^{-\gamma}U(x).
\]

Consequently it gives a smooth nonconstant periodic solution of the
cylinder system on
\(\mathbb S^2\times\mathbb T_{\log r}\).  A disjoint reflected pair makes
the construction nonaxisymmetric about the dilation centre and cancels
total helicity exactly in every log cell.

For \(1<\gamma<3/2\), the one-sided inward lattice has finite energy,
increasing local Reynolds number, and summable turnover times.  It is
therefore an admissible **leading Euler storage wake** for the proposed
cascade.

This does **not** construct a Navier--Stokes singularity.  The bubbles are
steady for Euler, not for viscous Navier--Stokes.  The remaining
prize-level problem is a time-dependent annular transition that creates
the next bubble, transports the active carrier and centre harmonics, and
cancels viscosity to all required orders with a smooth terminal force.

## 1. Compact steady seed

Gavrilov proves that there is a nonzero pair

\[
 W\in C^\infty_c(\mathbb R^3;\mathbb R^3),
 \qquad P\in C^\infty(\mathbb R^3),
\]

satisfying

\[
 (W\cdot\nabla)W+\nabla P=0,
 \qquad \operatorname{div}W=0,
\tag{1.1}
\]

with the velocity supported in an arbitrarily small neighbourhood of a
circle.  The pressure can be made compactly supported as well.  Indeed,
\(\nabla P=0\) off the velocity support, so \(P\) is constant on the
unbounded component and may be normalized to vanish there.  More
specifically, in Gavrilov's localization
\[
 d\widetilde P=\chi(P)^2\,dP
\]
with \(\chi\) supported in a small positive pressure interval.  After the
exterior constant is subtracted, \(\widetilde P\) is supported in the
corresponding outer solid torus, which shrinks to the seed circle with the
velocity support.

The primary source is A. V. Gavrilov,
[*A steady Euler flow with compact
support*](https://arxiv.org/abs/1810.08020).  Its final modulation by a
smooth function of pressure is what makes the seed \(C^\infty\)-flat at
the boundary of its support.

Fix any \(r>1\).  Because the joint support of velocity and normalized
pressure may be made arbitrarily thin and the Euler equations are invariant
under rigid motions and scaling, choose the seed so that, for some
\(0<a<b\),

\[
 \operatorname{supp}(W,P)
 \subset\{x:a<|x|<b\},
 \qquad
 \frac ba<r.
\tag{1.2}
\]

The strict inequality is the disjointness margin.

## 2. Bi-infinite exact wake

For \(j\in\mathbb Z\), define

\[
\begin{aligned}
 W_j(x)&=r^{\gamma j}W(r^jx),\\
 P_j(x)&=r^{2\gamma j}P(r^jx).
\end{aligned}
\tag{2.1}
\]

Spatial scaling gives

\[
 (W_j\cdot\nabla)W_j+\nabla P_j
 =
 r^{(2\gamma+1)j}
 \big((W\cdot\nabla)W+\nabla P\big)(r^jx)=0,
\tag{2.2}
\]

and \(\operatorname{div}W_j=0\).

The radial support of \(W_j\) is contained in

\[
 a r^{-j}<|x|<b r^{-j}.
\tag{2.3}
\]

Condition \(b/a<r\) makes adjacent intervals in (2.3) disjoint; all
nonadjacent ones are then disjoint as well.  Hence

\[
 U(x)=\sum_{j\in\mathbb Z}W_j(x),
 \qquad
 \Pi(x)=\sum_{j\in\mathbb Z}P_j(x)
\tag{2.4}
\]

are locally finite sums on \(\mathbb R^3\setminus\{0\}\).  They are smooth
there.  At every point either all terms vanish or exactly one seed solves
(1.1), so there are no cross interactions:

\[
 (U\cdot\nabla)U+\nabla\Pi=0,
 \qquad
 \operatorname{div}U=0.
\tag{2.5}
\]

Changing index in (2.4) gives the exact discrete homogeneity

\[
 U(rx)=r^{-\gamma}U(x),
 \qquad
 \Pi(rx)=r^{-2\gamma}\Pi(x).
\tag{2.6}
\]

With \(s=\log|x|\) and \(\omega=x/|x|\), set

\[
 V(s,\omega)=e^{\gamma s}U(e^s\omega),
 \qquad
 q(s,\omega)=e^{2\gamma s}\Pi(e^s\omega).
\tag{2.7}
\]

Then \(V,q\) are \(L=\log r\)-periodic and solve the exact stationary
cylinder equations.  This is the nonperturbative periodic orbit that the
earlier wake audit left open.

## 3. Fluxes, nonaxisymmetry, and helicity

There is a spherical gap between every two supports.  On a sphere in such
a gap, \(U=\Pi=0\).  Each stationary conservation flux is independent of
radius, so all mass, momentum, angular-momentum, energy, and helicity
fluxes vanish.  The construction therefore satisfies the zero-flux
identities from the cylinder audit automatically.  It is a storage wake,
not the active energy-transfer region.

Gavrilov's original seed is axisymmetric about its own axis.  Translate a
sufficiently small seed away from the dilation centre before imposing
(1.2), and choose its intrinsic symmetry axis not parallel to the ray from
the dilation centre to the translated torus.  The scaled torus centres lie
on that ray while their intrinsic axes lie on distinct parallel lines.
Consequently no one-parameter rotation group about the dilation centre
preserves the lattice.  A generic finite collection of disjoint translated
seeds gives the same conclusion without relying on this geometric
description.

Exact helicity cancellation can be imposed without creating interactions.
For compatibility with a later componentwise compact anti-divergence, first
choose the Gavrilov pressure modulation so that the individual seed has
zero angular momentum,
\[
 \int_{\mathbb R^3}x\times W(x)\,dx=0.
\]
The sign-changing modulation that achieves this is constructed in
[`2026-07-29-gavrilov-active-transition-ledger.md`](2026-07-29-gavrilov-active-transition-ledger.md).
This extra choice is not needed for the stationary Euler equation, but it
prevents a future transition stress from being forced to bridge two
disjoint reflected supports.

Let \(R\) be an orientation-reversing orthogonal map, choose the first
seed's support disjoint from its reflected image, and put

\[
 W^{R}(x)=R\,W(R^Tx),
 \qquad
 P^{R}(x)=P(R^Tx).
\tag{3.1}
\]

For \(\det R=-1\),

\[
 \operatorname{curl}W^{R}(x)
 =-R(\operatorname{curl}W)(R^Tx),
\tag{3.2}
\]

so

\[
 W^{R}\cdot\operatorname{curl}W^{R}
 =-(W\cdot\operatorname{curl}W)(R^Tx).
\tag{3.3}
\]

Choose the joint velocity-pressure support disjoint from its reflected
image.  The disjoint sum \(W+W^R\), with pressure \(P+P^R\), is again a
compact steady Euler seed and has zero total helicity.  Repeating this
paired seed gives exact cancellation in every logarithmic cell.  Generic
placement within the two reflection half-spaces retains nonaxisymmetry.
Other finite-group symmetries can be imposed by the same disjoint-orbit
construction.

## 4. One-sided finite-energy terminal wake

For a physical wake localized inside an outer scale, use only \(j\ge0\):

\[
 U_+(x)=\sum_{j=0}^{\infty}W_j(x),
 \qquad
 \Pi_+(x)=\sum_{j=0}^{\infty}P_j(x).
\tag{4.1}
\]

The pair is supported in a bounded ball, is smooth away from the
accumulation point \(x=0\), and solves stationary Euler on the punctured
space.  Its cell quantities are exact:

\[
\begin{aligned}
 \|W_j\|_2^2
 &=r^{(2\gamma-3)j}\|W\|_2^2,\\
 \|W_j\|_3^3
 &=r^{(3\gamma-3)j}\|W\|_3^3,\\
 \|\nabla W_j\|_2^2
 &=r^{(2\gamma-1)j}\|\nabla W\|_2^2,\\
 \operatorname{Re}_j
 &\asymp \nu^{-1}r^{(\gamma-1)j},\\
 \tau_j
 &\asymp r^{-(1+\gamma)j}.
\end{aligned}
\tag{4.2}
\]

Thus

\[
 1<\gamma<\frac32
\tag{4.3}
\]

simultaneously gives finite total kinetic energy, increasing Reynolds
number, divergent \(L^3\) strength, and summable turnover times.  This is
exactly the sharp scalar window found independently in the cascade
ledger.

Disjointness also permits an exact polynomial **amplitude** modulation
without creating an Euler error.  For any positive sequence \(K_j\),
replace (2.1) on the one-sided lattice by

\[
 W_j^{K}(x)=r^{\gamma j}K_j^\gamma W(r^jx),
 \qquad
 P_j^{K}(x)=r^{2\gamma j}K_j^{2\gamma}P(r^jx).
\tag{4.4}
\]

Each pair is still an amplitude-scaled steady Euler solution, and cross
terms still vanish.  Taking \(K_j=(j+1)^A\) gives the same polynomial
amplitude factor as in the flat-force ledger; the additional polynomial
factors do not disturb energy or time summability.  Exact discrete
homogeneity is then replaced by an exact nonautonomous bubble lattice.

This observation does **not** realize the ledger's internal carrier
\(k_j=K_j/\ell_j\).  The bubble's spatial derivative scale remains
\(\ell_j^{-1}\), not \(K_j/\ell_j\), and its viscous cost therefore lacks
the carrier factor \(K_j^2\).  An oscillatory microstructure or a
Gevrey-tame carrier corrector is still needed for that part of the proposed
Navier--Stokes stage.

There is, however, an exact **packed-bubble carrier**.  Let
\(\ell_j=r^{-j}\), let \(K_j\) be a positive integer, and put
\[
 \delta_j=c_0\frac{\ell_j}{K_j},
 \qquad
 a_j=\ell_j^{-\gamma}K_j^\gamma .
\tag{4.5}
\]
A fixed subannulus of diameter comparable to \(\ell_j\) contains
\(N_j\asymp K_j^3\) disjoint balls of radius comparable to \(\delta_j\).
Translate a spatially \(\delta_j\)-scaled copy of the joint
velocity-pressure support into each ball and set
\[
\begin{aligned}
 {\mathcal W}_{j,m}(x)
 &=a_j W\!\left(\frac{x-x_{j,m}}{\delta_j}\right),\\
 {\mathcal P}_{j,m}(x)
 &=a_j^2 P\!\left(\frac{x-x_{j,m}}{\delta_j}\right).
\end{aligned}
\tag{4.6}
\]
Every pair in (4.6) solves stationary Euler, all cross terms vanish, and
different annuli remain disjoint.  If
\({\mathcal W}_j=\sum_{m=1}^{N_j}{\mathcal W}_{j,m}\), then
\[
\begin{aligned}
 \|{\mathcal W}_j\|_2^2
 &\asymp \ell_j^{3-2\gamma}K_j^{2\gamma},\\
 \|{\mathcal W}_j\|_3^3
 &\asymp \ell_j^{3-3\gamma}K_j^{3\gamma},\\
 \|\nabla{\mathcal W}_j\|_2^2
 &\asymp \ell_j^{1-2\gamma}K_j^{2\gamma+2},\\
 k_j&\asymp K_j/\ell_j,\\
 \operatorname{Re}_{j,\mathrm{carrier}}
 &\asymp \nu^{-1}\ell_j^{1-\gamma}K_j^{\gamma-1}.
\end{aligned}
\tag{4.7}
\]
Thus the leading Euler carrier bath and all scalar derivative counts are
exact, not asymptotic.  Over the ledger stage time
\[
 \tau_j=\ell_j^{1+\gamma}K_j^{-\gamma},
\]
its nominal viscous loss is
\[
 \nu\tau_j\|\nabla{\mathcal W}_j\|_2^2
 \asymp
 \nu\ell_j^{2-\gamma}K_j^{\gamma+2},
\tag{4.8}
\]
exactly the polynomial-carrier scaling.  For
\(K_j=(j+1)^A\), the exponential factors in \(\ell_j\) still give finite
energy and divergent carrier-scale Reynolds number throughout
\(1<\gamma<3/2\).

The bath is not an active cascade: its bubbles do not interact and each
has zero external flux.  Formulae (4.5)--(4.8) remove the *static leading
Euler carrier* obstruction, but not the time-dependent creation,
low-mode-strain coupling, viscous endpoint jets, or all-order seam
matching.

For the DSS amplitudes,
\[
 \|P_j\|_1=r^{(2\gamma-3)j}\|P\|_1.
\]
The same geometric factor times a polynomial occurs for
\(K_j=(j+1)^A\), so the pressure remains summable in local \(L^1\).  Hence
the partial sums
converge in
\[
 U_+^{(N)}\longrightarrow U_+\quad\hbox{in }L^2,
 \qquad
 U_+^{(N)}\otimes U_+^{(N)}+\Pi_+^{(N)}I
 \longrightarrow U_+\otimes U_++\Pi_+I\quad\hbox{in }L^1.
\]
For every compactly supported smooth vector test function \(\varphi\),
each finite partial sum obeys
\[
 \int_{\mathbb R^3}
 \big(U_+^{(N)}\otimes U_+^{(N)}+\Pi_+^{(N)}I\big):
 \nabla\varphi\,dx=0.
\]
The \(L^1\) stress convergence passes this identity to the limit.  The
\(L^2\) convergence likewise passes
\(\int U_+^{(N)}\cdot\nabla\psi\,dx=0\) to the limit for every scalar test
\(\psi\).  Thus no point source is hidden at the accumulation point:
\((U_+,\Pi_+)\) is a finite-energy weak stationary Euler solution on all
of \(\mathbb R^3\).  It is unbounded at \(x=0\) and is not a smooth
Navier--Stokes state there.  At each preterminal time, a proposed cascade
would contain only finitely many active cells, hence a smooth velocity.

## 5. What remains

The explicit wake removes the need to solve a global cylinder Liouville or
Floquet-bifurcation problem.  It also supplies an infinite-dimensional
endpoint state: each disjoint bubble can retain outgoing energy, phase
centre data, and finite-moment corrections without interacting with the
other stored cells at principal Euler order.

Three hard requirements remain.

1. **Viscous endpoint jets.**  A bare steady bubble has residual
   \(-\nu\Delta W_j\).  Across the shrinking tail these residuals are not
   smooth at the terminal time.  The active stage must incorporate their
   forward heat/Navier--Stokes jets or an oscillatory Reynolds bath.
2. **Active annular transition.**  The zero-flux bubbles cannot transfer
   energy between scales.  A time-dependent transition must drain the
   parent, create the smaller child, and export the complementary data
   into the disjoint wake.
3. **All-order tame closure.**  The transition must solve the finite-band
   Kelvin--heat centre block and transverse corrector equations through
   \(M_j\asymp j^2/\log j\), with Gevrey bounds and
   \(e^{-c j^2}\) residual/seam errors.

The remaining prize-level target is therefore a **local active transition
theorem**, not existence of a stationary wake.
