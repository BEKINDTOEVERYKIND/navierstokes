# Counterfactual success portfolio: what a genuine result would probably look like

**Date:** 2026-08-30

**Status:** strategy and pre-registered computational discriminators, grounded
in C188 and C192--C196 and the landed Stage-5 data; not a theorem or a
claimed numerical candidate

## 0. Why this note exists

The repository has repeatedly described the residue of one failed
implementation as “the one realistic path.” That wording is withdrawn.
The local gates are not an exhaustive partition of the mathematics.

Conditional on this program eventually producing a Clay-level result, the
most likely retrospective story is:

1. a coherent renormalized object is discovered, probably numerically;
2. the computation exposes a symmetry or a coercive bulk plus a finite bad
   block;
3. the full residual and finite block are interval-certified; and
4. the existing exact pump, scaling, pressure, and tail estimates are used
   as proof components.

It is substantially less likely that all current hand-built modules simply
snap together after enough additional local estimates.

## 1. Route A: the \(A_2\) same-witness cascade

C192--C195 make the operator mechanism real: the exact \(A_2\) orbit has a
large hyperbolic polarization cocycle and a quantitative shrinking-tube
bundle. C193 supplies the fixed-energy polarization swap, C194 supplies a
pressure-resolved first-order local bridge, and C196 supplies exact real
periodic solenoidal endpoint profiles.

C196 also gives the strict discriminator. A fixed-energy lower bound
\(c_0q^{3/2}\) requires at least \(c_0^2q^3\) dynamically retained Fourier
modes. A three-coordinate tube of physical half-width \(cq^{3/4}\) has at
most \(64c^3q^{9/4}\) modes once \(cq^{3/4}\ge\sqrt3/2\); a projective
tube in \(0<|k|\le bq\) has at most \(64b^3c^2q^{5/2}\) once its three
box half-widths exceed \(\sqrt3/2\); and C180's one-sided retained slab
has at most \(64\delta^3q^3/J^4\) under the same unit-cube condition.
Within this retained-band formulation, the three
creditable escape mechanisms are:

* a fixed-aperture hyperbolic/band theorem;
* an honest anisotropic child of volume \(J^4q^{-3}\), whose forced focus is
  \(q^{3/2}/J^2\); or
* a terminal nonlinear converter that genuinely creates the missing band.

Another central-ray gain estimate does not address this discriminator.

The next creditable \(A_2\) computation is one complete exact periodic
linearized-Euler witness with a full endpoint band and an **absolute** point
lower bound after subtracting C194. A relative concentration quotient alone
receives no stage credit.

## 2. Route B: a direct nonautonomous Navier--Stokes invariant graph

This route bypasses module assembly. Let \(S_T^\mu\) be dimensionless
unforced Navier--Stokes evolution and define the complete exit chart

\[
  ({\cal C}v)(y)=g^{-1}Q^Tv(a+q^{-1}Qy),\qquad
  \mu'=\frac qg\mu.                                      \tag{2.1}
\]

Search for a curve of full states, including all wake,

\[
 \boxed{{\cal F}(\mu)
   ={\cal C}S_{T(\mu)}^\mu X(\mu)-X((q/g)\mu)=0.}          \tag{2.2}
\]

For the axisymmetric parity-pinned search below, restrict this generic chart
to its symmetry stabilizer: \(a=0\) in the pinned coordinates and
\(Qe_z=e_z\) (the first experiment takes \(Q=I\)). Under this restriction
the axis and parity class are preserved and the scalings in (2.4) are exact.

A single autonomous profile is the wrong object: C188 proves that the
viscosity coordinate moves, and the Stage-5 data show that a large fraction
of the endpoint remains outside the selected child shell.

The most economical discovery class is axisymmetric-with-swirl on
\(\mathbb R^3\), which retains genuine three-dimensional stretching while
reducing the PDE to two spatial variables. With

\[
 G=ru^\theta,\qquad \eta=\omega^\theta/r,
\]

use the parity

\[
 G(r,-z)=G(r,z),\qquad \eta(r,-z)=-\eta(r,z),             \tag{2.3}
\]

and free-space decay. This pins the terminal center to the axis and retains
both swirl and poloidal velocity. The exact exit scalings are

\[
 G'(r,z)=\frac qgG_T(r/q,z/q),\qquad
 \eta'(r,z)=\frac1{gq^2}\eta_T(r/q,z/q).                 \tag{2.4}
\]

At the axis the search enforces the standard regularity conditions
\(G=O(r^2)\), \(\partial_r\eta=0\), and the corresponding even regular
streamfunction condition. Discovery may use a padded free-space elliptic
solve, but certification must reach spatial infinity rather than hide a
radial-wall residual.

The first bounded experiment is pre-registered as follows:

* \(q=6/5\), \(7/6<\log g/\log q<1.48\);
* three consecutive returns, not one selected-shell transfer;
* a \(128\times256\) axisymmetric grid followed by \(192\times384\);
* 12--20 independent symmetry-preserving restarts; and
* full-field endpoint residual in

  \[
  \|v\|_{{\cal X}_\mu}^2
   =\|v\|_{H^{1/2}}^2+\mu\|v\|_{H^{3/2}}^2
                    +\mu^2\|v\|_{H^{5/2}}^2.            \tag{2.5}
  \]

A candidate receives further work only if its relative residual is below
\(10^{-5}\), doubled-grid/half-step replay is below \(10^{-4}\), top-third
spectral and outer-collar residual fractions are each below \(10^{-7}\),
all parameters stay in the interiors of their boxes, and both swirl and
poloidal energies stay nonzero. Boundary-hitting, \(T\to0\), cutoff pile-up,
or a one-component collapse is a failed search.

The landed data justify exactly one such run. Seed B gives a resolved
one-step \(Q_3^{\rm oct}>1\). At the \(N=320\) \(Q_3\) peak, \(56.23\%\)
lies in the child shell, \(34.61\%\) in the parent band, and \(9.15\%\)
elsewhere. Separately, the best stored corrected shape-return seed has
fidelity about \(0.109\), and the best stored second-octave candidate is far
below one.
This is positive evidence for one-step strength and negative evidence for a
child-only fixed profile. It does not test (2.2).

The old optimizer also had two implementation defects. This checkpoint
repairs both in engines/return_map_opt.py: the seed-dependent time step is
no longer detached from automatic differentiation, and replayed
concentration and \(L^3\) quantities now use Euclidean vector speed rather
than componentwise magnitude. Score, field, normalization, and shift are
also snapshotted atomically at the best pre-step iterate across restarts.
The engine still optimizes selected shells
and is not represented as an implementation of (2.2). Its anti-cheat caps
are translated conservatively by \(1/\sqrt3\) and \(1/3\). The historical
independent audit script remains unchanged for reproducibility, and stored
results have not been rerun under the corrected engine diagnostic; they
receive no new numerical credit here.

## 3. Route C: a non-outgoing Euler-dominant Type-II profile

Let \(s=T-t\) and

\[
 u(x,t)=s^{-\alpha}U(x/s^\beta,\tau),\qquad
 \alpha+\beta=1,\qquad \tau=-\log s.                    \tag{3.1}
\]

Uniform finite total energy gives \(\beta\ge2/5\), while viscosity is lower
order only for \(\beta<1/2\). Chae--Wolf exclude the endpoint
\(\beta=2/5\) for their Type-I atomic/DSS class; the proposed search avoids
that endpoint and uses the strict strip

\[
 \boxed{\frac25<\beta<\frac12,\qquad
        1<\Gamma:=\frac\alpha\beta<\frac32.}             \tag{3.2}
\]

This is the same numerical corridor C188 found from the stage ledger, now
derived from a genuinely different similarity mechanism. Recent primary
work proves that finite-energy Euler self-similarity requires
\(\beta\ge2/5\), and that a smooth axisymmetric or locally outgoing global
profile cannot live in the strict sub-parabolic strip. Therefore a viable
profile must be fully three-dimensional, nonaxisymmetric, and non-outgoing;
a rotating or discretely self-similar relative equilibrium is a natural
first test.

Fix

\[
 \beta=\frac9{20},\qquad \alpha=\frac{11}{20},\qquad
 \Gamma=\frac{11}{9}.                                   \tag{3.3}
\]

For rotation generator \(J\), search \((U,P,\kappa)\) satisfying

\[
 {\mathbb P}\left[
 \kappa\{JU-(Jy\cdot\nabla)U\}
 +\frac{11}{20}U+\frac9{20}(y\cdot\nabla)U
 +(U\cdot\nabla)U\right]=0.                              \tag{3.4}
\]

The viscous coefficient then decays exactly as

\[
                              \mu(\tau)=\nu e^{-\tau/10}. \tag{3.5}
\]

The large-similarity-radius law used for **inner-profile discovery** is
\(|y|^{-11/9}\), with a logarithmic angular twist compensating the time
rotation.  This is not a proposed physical far field: the corresponding
global profile is not in \(L^2(\mathbb R^3)\).  Interpreted only as the
terminal inner asymptotic \(u_T(x)\sim |x|^{-\Gamma}\), it has finite local
energy, while

\[
 r^{-1}\int_{B_r}|u_T|^2\sim r^{2-2\Gamma}\longrightarrow\infty
 \qquad(r\downarrow0),                                  \tag{3.6}
\]

so it is not removed by the scalar epsilon-regularity heuristic. Its local
terminal \(L^3\) norm diverges, so critical-norm regularity criteria are not
being assumed away.  None of these inner asymptotics supplies finite-energy
initial data or an exact Navier--Stokes solution.

An **inner-profile discovery candidate** requires a full weighted residual
below \(10^{-10}\), the correct tail and twist, and the exact energy-flux
identity

\[
 \frac{\displaystyle\int_{\partial B_R}
   \left[\frac\beta2|U|^2y\cdot n+
   \left(\frac{|U|^2}{2}+P\right)U\cdot n\right]\,dS}
 {\displaystyle\int_{B_R}|U|^2\,dy}
 =\frac{5\beta}{2}-1=\frac18                            \tag{3.7}
\]

at three radii. It must also exhibit a resolution-stable finite bad
spectrum separated from a coercive tail. A floating residual minimum
without those tests is merely a failed discovery search.

Promotion beyond an inner-profile candidate has separate load-bearing
requirements.  One must certify a time-dependent cutoff or outer matching
that produces smooth finite-energy data; compute and control the complete
matching-annulus defect; and include the pressure/Biot--Savart nonlocal
coupling between the inner core and outer flow with explicit constants.  In
the unforced route the induced defect must be cancelled by the actual
dynamics, not relabeled as a force.  The matched solution must also satisfy
the terminal singular-center local-energy and critical-norm obligations.
Until all of these are proved, even an interval-certified weighted residual
for (3.4) is discovery evidence only and is **not progress toward a Clay
singular solution**.

The direct constraints are
[Constantin--Ignatova--Vicol](https://arxiv.org/abs/2602.17570) and
[Chae--Wolf](https://arxiv.org/abs/1706.02020).
[Pineau--Vicol](https://arxiv.org/abs/2607.09619) concerns the distinct
parabolic rotated problem and is cited only as contextual methodology, not
as an exclusion of (3.4).

## 4. Route D: exploit the allowed smooth force at all orders

The Clay formulation permits a smooth admissible force. In the periodic
alternative targeted here, the defect must extend past the terminal time
and be terminally flat so that the required long-time derivative decay can
be arranged. One concrete version worth testing retains the complete
charged and zero-charge wake and solves one all-order Gevrey endpoint
equation; the old finite-ladder shortcut does not do this.

The bounded discriminator is not a new module. Form increasing finite
sections of the full retained-wake endpoint derivative, enclose their
smallest singular values, and prove an analytic tail. Backward-heat growth
\(e^{cK^2}\) or a persistent adjoint compatibility kills the route.
Factorial/Gevrey growth compatible with \(C^M(M!)^2\) warrants a
Nash--Moser/Borel attack. No GPU profile search should precede this
singular-value test.

## 5. Portfolio verdict

The active portfolio is:

1. finish one strict \(A_2\) fixed-aperture or log-decorated endpoint test;
2. run the nonautonomous full-wake search (2.2);
3. independently search the nonaxisymmetric Type-II equation (3.4); and
4. run one retained-wake tame-inverse spectrum test for the forced route.

Accumulating more class-specific no-gos is not a plausible route to a
global regularity theorem. Bare \(L^2\) cocycle growth, another relative
concentration multiplier, a selected-shell optimizer, or finite-dimensional
numerics without a tail certificate receive no prize-level credit.
