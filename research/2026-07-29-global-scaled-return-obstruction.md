# Global scaled-return obstruction for a periodic Bloch amplifier

Date: 2026-07-29

## Decision

The exact infinite-chain Bloch instability is a genuine **local
amplifier**, but a globally periodic Bloch module cannot be the recurrent
Millennium cell.

There are three independent reasons.

1. On a fixed torus, a clean isotropic scale return changes helicity while
   preserving the normalized \(L^2\) energy.  Hence no inviscid finite-time
   orbit can connect a nonzero-helicity pump profile to a strictly finer
   scaled copy.
2. Even when helicity is zero, a clean finer copy lies in a proper Fourier
   sublattice.  Euler reversibility and invariance of that sublattice forbid
   entering it at a finite time.  The same obstruction holds for a
   uniformly strong, fixed-clock sequence whose normalized viscosity and
   force tend to zero.
3. A volume-filling periodic copy cannot beat viscosity at unbounded
   frequency with bounded energy.  Exact Navier--Stokes torus scaling
   multiplies the energy by the square of the scale ratio.

The current charge-lattice causal fork does not test a return.  Its
`diagonal_source_off` branch retains the old pump, and its diagonal lattice
has rank two.  A pass can only show continued growth of the already known
sum-chain instability about the **old** pump.

The surviving prize route is therefore a physically localized,
non-Beltrami (or helicity-balanced) heteroclinic with a leading infinite
tail.  It must shrink its occupied volume, drain the old pump, and emerge
as a new parent capable of amplifying a fresh transverse seed.  No existing
GPU run tests those properties.

No Navier--Stokes singularity is proved here.

## 1. Active-stage normalization

For a lattice-compatible active ansatz that can be placed on a fixed
normalized torus, write

\[
 u(t,x)=A_jv_j(\tau,y),\qquad
 y=N_jx,\qquad
 \tau=A_jN_j(t-t_j).
\]

The forced Navier--Stokes equation becomes

\[
 \partial_\tau v_j+B(v_j,v_j)
 =\varepsilon_j\Delta v_j+g_j,
 \qquad
 \varepsilon_j={\nu N_j\over A_j},
 \qquad
 g_j(\tau,y)
 ={ \mathbb P f(t_j+\tau/(A_jN_j),\,y/N_j)
    \over A_j^2N_j}.
\tag{1.1}
\]

A recurrent inviscid cell on a fixed normalized clock would have
\(\varepsilon_j\to0\) and \(g_j\to0\), together with enough strong
compactness to pass to an Euler stage.

The common-domain and normalized-force hypotheses must be verified by a
physical construction.  The substitution \(y=N_jx\) does not put an
arbitrary periodic solution on one fixed torus, and smoothness of the
physical force alone does not imply \(g_j\to0\) in the rescaled norm.

The local Bloch amplifier is compatible with this limit.  A single
unidirectional or Beltrami Fourier pump \(P_N\) satisfies

\[
 U_N(t)=A e^{-\nu N^2(t-t_0)}P_N
\tag{1.2}
\]

exactly in unforced Navier--Stokes.  On \(O((AN)^{-1})\) active times its
relative heat drift is \(O(\nu N/A)\), so a high-Reynolds stage shadows the
positive viscous/Euler Bloch eigenmode without an active pump force.
Vasudevan's continued-fraction theorem supplies the positive viscous
eigenvalue for type-I unidirectional chains.  This validates amplification;
it supplies no scale-changing endpoint.

## 2. Energy--helicity obstruction to a clean torus return

Let \(P\) be a real smooth periodic profile and let

\[
 {\cal D}_{r,O,\phi}P(x)
 =O\,P(rO^Tx+\phi),
\tag{2.1}
\]

where \(r\ge2\) is an integer and \(O\) is an orthogonal lattice symmetry.
The factor \(O\) rotates the velocity components together with the
coordinates, so incompressibility is preserved.  Curl then transforms
with the additional orientation factor \(\det O\).
With

\[
 E(P)={1\over2}\int_{\mathbb T^3}|P|^2,\qquad
 H(P)=\int_{\mathbb T^3}P\cdot\operatorname{curl}P,
\]

the covering map in (2.1) gives

\[
 E(c{\cal D}_{r,O,\phi}P)=c^2E(P),
\qquad
 H(c{\cal D}_{r,O,\phi}P)
 =(\det O)c^2rH(P).
\tag{2.2}
\]

Suppose a smooth Euler orbit has

\[
 v(0)=P,\qquad
 v(T)=c{\cal D}_{r,O,\phi}P.
\tag{2.3}
\]

Energy conservation gives \(c^2=1\).  If \(H(P)\ne0\), helicity
conservation then gives

\[
 1=(\det O)r,
\]

which is impossible for \(r>1\).  Thus

\[
\boxed{\text{no clean finite-time scaled return of a helical pump}.}
\tag{2.4}
\]

This applies in particular to a circular Beltrami pump.  Adding a
super-algebraically small daughter doublet in a topology that controls
helicity (for example \(H^{1/2}\), and hence \(H^s\) with \(s>1/2\)) does
not alter the contradiction: the limiting helicity is still that of the
pump.  Merely making the daughter small in \(L^2\) is insufficient, since a
vanishing-energy high-frequency tail can carry order-one helicity.

There is a quantitative wake version.  With the energy convention above,
start from a positive Beltrami pump of wavenumber \(N\), energy \(E\), and
helicity \(2NE\).  Suppose a child at \(rN\), with relative helicity sign
\(\sigma\in\{+1,-1\}\), carries energy \(\theta E\), while a
spectrally disjoint (equivalently here, helicity-orthogonal) wake carries
\((1-\theta)E\).  Helicity conservation forces

\[
 H_{\rm wake}=2(1-\sigma r\theta)NE.
\]

If the wake is supported below wavenumber \(K_{\rm wake}\), the spectral
realizability bound

\[
 |H_{\rm wake}|\le2K_{\rm wake}E_{\rm wake}
\]

implies

\[
 \boxed{
 {K_{\rm wake}\over N}
 \ge {|1-\sigma r\theta|\over1-\theta}.}
\tag{2.5}
\]

Consequently a same-chirality clean energy return
\(\theta\uparrow1\) requires an unbounded-frequency opposite-helicity
wake.  Choosing \(\theta=1/r\) avoids the helicity defect only by leaving
the order-one energy fraction \(1-1/r\) behind.  A recurrent Beltrami
handoff must therefore keep an order-one wake or lose compactness in a
helicity-controlling topology to a still finer helicity-carrying tail.

Viscosity can change helicity, but it does not evade (2.4) in a fixed-clock
stage limit for which \(\varepsilon_jT_j\to0\) and the relevant derivative
norms stay uniform.  An order-one viscous helicity change must instead use
either a nonvanishing viscous clock or a loss of those uniform bounds,
possibly through a transient high-frequency excursion.  That excursion
need not remain as a terminal wake, so this observation alone does not
exclude a noncompact intermediate mechanism.

## 3. Proper-sublattice entry is impossible on a finite Euler clock

Let \(\Lambda\subset\mathbb Z^3\) be the additive Fourier lattice of the
normalized parent and put

\[
 Y_{r,O}
 =\{v:\operatorname{supp}\widehat v
       \subset rO\Lambda\}.
\tag{3.1}
\]

Because \(rO\Lambda\) is an additive subgroup, \(Y_{r,O}\) is invariant
under Euler.  Smooth Euler evolution is locally unique forward and
backward.  Therefore

> If a smooth Euler solution belongs to \(Y_{r,O}\) at one finite time,
> it belongs to \(Y_{r,O}\) throughout its smooth lifespan.

Indeed, evolve the terminal datum backward inside the invariant subspace
and use uniqueness.  A primitive parent \(P\notin Y_{r,O}\) therefore
cannot reach a clean scaled endpoint in \(Y_{r,O}\) at finite time, even
when \(H(P)=0\).

The useful asymptotic version is as follows.  Fix \(T<\infty\) and
\(s>5/2\).  Suppose solutions of (1.1) satisfy

\[
\begin{aligned}
 &\sup_j\|v_j\|_{C([0,T];H^s)}<\infty,\\
 &\varepsilon_j\to0,\qquad
   g_j\to0\quad\hbox{in }C([0,T];H^{s-1}),\\
 &v_j(0)\to P\quad\hbox{in }H^{s-1},\qquad
   \operatorname{dist}_{H^{s-1}}(v_j(T),Y_{r,O})\to0.
\end{aligned}
\tag{3.2}
\]

The equation gives a uniform \(\partial_t v_j\) bound in \(H^{s-2}\).
Rellich compactness and interpolation therefore give, after passing to a
subsequence, strong convergence in \(C([0,T];H^{s-1})\) to an \(H^s\)
classical Euler solution \(v\).  Since \(P\ne0\), its terminal limit is
nonzero by Euler uniqueness.  Its endpoint lies in \(Y_{r,O}\), so its
initial datum must also lie there.  This contradicts a primitive
\(P\notin Y_{r,O}\).

Thus an asymptotic scale reset must lose at least one of:

- finite normalized residence time;
- a vanishing normalized force/viscosity defect;
- strong compactness;
- a vanishing off-daughter wake; or
- a nondegenerate fixed parent shape.

An infinite-time heteroclinic remains logically possible.  For a helical
profile, however, strong convergence at both ends still contradicts
energy and helicity.  Its only escape is again a nonvanishing or
noncompact wake.

## 4. Global periodic scaling cannot survive viscosity

For a volume-filling periodic module

\[
 U_j(x)=A_jP(N_jx),
\tag{4.1}
\]

the covering map preserves the normalized volume, so

\[
 E(U_j)=A_j^2E(P).
\tag{4.2}
\]

If \(U_j\) is the whole state, an energetically dominant state, or an
\(L^2\)-orthogonal stage projection, the forced energy inequality bounds
\(E(U_j)\) on every common finite physical time interval for smooth
forcing.  Hence \(A_j\) is bounded in that setting.  Linearized strain is
\(O(A_jN_j)\), whereas viscosity on a fixed-shape, comparable-frequency
Bloch class is \(O(\nu N_j^2)\).  More explicitly,

\[
 {1\over2}{d\over dt}\|w\|_2^2
 +\nu\|\nabla w\|_2^2
 +\langle w,(w\cdot\nabla)U_j\rangle=0
\]

gives, when the class has minimum frequency \(cN_j\),

\[
 {1\over2}{d\over dt}\|w\|_2^2
 +\big(c^2\nu N_j^2-C A_jN_j\big)\|w\|_2^2
\le0.
\tag{4.3}
\]

At sufficiently large \(N_j\) such a class is strictly damped.  Estimate
(4.3) does not cover a perturbation class containing a fixed or other
low-frequency central mode; the broader global-stage obstruction is the
dimensionless ratio \(\nu N_j/A_j\to\infty\) when \(A_j\) is bounded.

The exact unforced torus Navier--Stokes scaling makes the same point more
starkly:

\[
 u_r(t,x)=r\,u(r^2t,rx)
\tag{4.4}
\]

solves the same-viscosity equation, but

\[
 E(u_r)=r^2E(u)
\tag{4.5}
\]

because \(x\mapsto rx\) tiles the fixed torus rather than localizing one
copy.  With forcing, the corresponding force is
\(f_r(t,x)=r^3f(r^2t,rx)\).  Repeated copies cannot occur as successive
states of one bounded-energy finite-time trajectory.  This does not
forbid (4.4) as a map between separate solutions with separately scaled
initial data.

For a zero-helicity generalized shear one can preserve energy by omitting
the amplitude factor \(r\), but then the effective Reynolds number drops
by \(r^{-1}\) at every return.  The fixed-shape,
comparable-frequency Bloch instability eventually disappears.  Therefore
changing from a Beltrami pump to a zero-helicity global shear removes
(2.4), but not this global fixed-profile prize obstruction.

This is exactly where physical localization matters.  On
\(\mathbb R^3\), the same Navier--Stokes scaling has
\(\|u_r\|_2^2=r^{-1}\|u\|_2^2\), because it shrinks one copy instead of
tiling the domain.  Thin rings or compact intermittent blocks can exploit
that volume gain; a global charge lattice cannot.

## 5. Why steady pump forcing does not rescue the iteration

Holding a scale-\(N_j\) Laplacian eigenpump steady costs

\[
 f_j=-\nu\Delta U_j,
\qquad
 \|\nabla^m f_j\|_\infty
\asymp\nu A_jN_j^{m+2}.
\tag{5.1}
\]

An assembled terminal hierarchy of such forces cannot be \(C^\infty\) if
the strains \(A_jN_j\) become unbounded.  Although the **normalized** force
in (1.1) is only \(O(\varepsilon_j)\), the physical high-frequency force
in (5.1) is not terminally flat.

One low-frequency pump can be maintained by a smooth fixed force, but one
fixed Bloch geometry has only finitely many type-I central charges at the
pump scale.  Its remote sidebands are not a sequence of new scale-\(N_j\)
pumps.

More generally, a fixed \(C^\infty\) periodic force cannot repair the
order-one energy--helicity mismatch of a mean-zero, fixed-shape global
scaled block at terminal scales.  Uniformly near the terminal time its
Fourier coefficients satisfy

\[
 |\widehat f(k,t)|\le C_m(1+|k|)^{-m}
 \qquad\hbox{for every }m.
\tag{5.2}
\]

For a block whose spectrum stays in a uniformly controlled band near
\(N\), on one natural clock \(\Delta t\asymp(AN)^{-1}\), its
force-induced energy and helicity changes are bounded schematically by

\[
 \Delta E_f
 \lesssim A|\widehat f(N)|\Delta t
 \lesssim {|\widehat f(N)|\over N},
\qquad
 \Delta H_f
 \lesssim AN|\widehat f(N)|\Delta t
 \lesssim |\widehat f(N)|.
\tag{5.3}
\]

Both are \(O(N^{-\infty})\).  This estimate is not asserted for an
arbitrary localized broad-spectrum state with appreciable low modes; there
one must project the force onto the actual high-frequency block.  Under
the stated fixed-band, fixed-clock assumptions, a Clay-admissible smooth
force can keep a super-algebraically tiny dormant seed alive, but it
cannot supply the order-one invariant change needed by a clean
high-frequency return.

Thus the steady-forced unstable manifold is useful as a local spectral
model.  An admissible terminal cascade must switch to the unforced
heat-decaying parent (1.2) on active intervals, or make the pump residual
flat by a genuinely nonlinear localized construction.

## 6. The current diagonal fork is not a recurrence gate

In charge coordinates \(k=B(m,a,b)\), the retained diagonal set is

\[
 D=\{(m,a,b):a=b\}
   =\mathbb Z(1,0,0)+\mathbb Z(0,1,1).
\tag{6.1}
\]

Its physical lattice is

\[
 B D=\mathbb Zp+\mathbb Z(q_1+q_2),
\tag{6.2}
\]

which has rank two.  It is an exact 2D3C invariant system and is globally
regular.  More importantly, the `diagonal_source_off` fork keeps both:

- the old pump \(\mathbb Zp\);
- the populated sum chain \(q_1+q_2+\mathbb Zp\).

Hence its `child_internal_work = N(D,D)` includes the old-pump--child
interaction.  Positive growth after deleting the original \(q_1,q_2\)
chains is exactly what the already certified positive sum-chain eigenvalue
predicts.  It does not show:

1. drain or disappearance of the old pump;
2. concentration into a stationary shorter-scale parent;
3. growth after the old pump is removed; or
4. amplification of a fresh transverse seed by the alleged new parent.

A causal fork can distinguish a slaved \(2\lambda\) response from an
independent \(\sigma\)-mode about the old pump.  It cannot establish a
parent-to-parent return.

## 7. Surviving theorem gate

The theorem above closes clean, primitive, volume-filling, fixed-profile
returns under uniform \(H^s\) compactness and a vanishing normalized
defect.  It does not rule out low-frequency transient mediation, a
long-clock viscous mechanism, or a tail that is \(L^2\)-small but
noncompact in a helicity-controlling topology.  A prize-relevant positive
lemma outside the closed class must instead construct a localized orbit
with all of the following:

1. **localized volume scaling:** a minor scale \(\ell_j\) and occupied
   volume \(V_j\) for which \(A_j^2V_j\) stays bounded while
   \(A_j\ell_j/\nu\to\infty\);
2. **old-parent drain:** after the active interval, the old parent is a
   summable wake rather than the source of the measured child growth;
3. **new-parent test:** with the old parent removed, the endpoint block
   amplifies a newly inserted pair of noncoplanar, super-algebraically
   small dormant seeds;
4. **helicity closure:** either the parent class has zero net helicity, or
   the opposite-helicity wake required by (2.5) is included in the leading
   recurrent profile and its dissipation ledger;
5. **leading infinite tail:** localization and all Fourier sidebands are
   part of the main profile, not an off-dictionary error tending to zero;
6. **terminally flat residual (a design requirement):** every active-stage
   force other than the dormant seed control must be
   \(O(N_j^{-\infty})\) in all space-time seminorms.  The Fourier estimate
   above does not by itself construct or prove this property.

This is naturally a thin-ring/compact-column or other intermittent
heteroclinic theorem, not a charge-lattice theorem.

No present GPU run implements this gate.  If a concrete localized endpoint
is produced, the first decisive numerical fork should:

- delete the old parent, not merely the original transverse sources;
- retain only the proposed localized child parent;
- insert two fresh transverse seeds in the child's own scaled geometry;
- require converged positive gain for both seeds for at least one child
  clock; and
- verify that curvature, envelope, pressure, and off-profile residuals
  decrease under the proposed physical localization scaling.

Until such an ansatz exists, additional resolution of the present diagonal
fork is an incremental confirmation of a known instability, not a
Millennium-level gate.

## Primary references

- S. Vasudevan, *Instability of unidirectional flows for the 2D
  Navier--Stokes equations and related alpha-models*,
  [arXiv:2011.02244](https://arxiv.org/abs/2011.02244).
- H. R. Dullin, Y. Latushkin, R. Marangell, S. Vasudevan, and
  J. Worthington, *Instability of unidirectional flows for the 2D
  alpha-Euler equations*,
  [arXiv:1901.01367](https://arxiv.org/abs/1901.01367).
- N. Kishimoto and T. Yoneda, *Characterization of three-dimensional
  Euler flows supported on finitely many Fourier modes*,
  [arXiv:2110.08039](https://arxiv.org/abs/2110.08039).
- Z. Lin and C. Zeng, *Unstable manifolds of Euler equations*,
  [arXiv:1112.4525](https://arxiv.org/abs/1112.4525).
