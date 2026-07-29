# The one-phase Kelvin cell inverse: an infinite cokernel and an exact forward heat block

## Status

This note audits the operator-theoretic gap in the proposed one-carrier
Gevrey construction.  It does **not** prove a Navier--Stokes singularity.

There are two exact conclusions.

1. A one-phase shear carrier, even with arbitrary time-dependent amplitude,
   scalar phase shift, and two transverse polarizations, has an
   infinite-dimensional carrier-chain kernel and adjoint kernel.  Therefore
   its stationary fast operator has no right inverse modulo only finitely
   many resonances.  The time-dependent operator with both complete
   endpoints prescribed has the corresponding infinite family of moment
   conditions.  Finite-dimensional amplitude and phase modulation cannot
   remove them.
2. In affine Kelvin coordinates that same infinite block is not mysterious:
   it is an exactly invariant one-dimensional heat system.  It propagates
   all phase harmonics explicitly, including the viscous correction to the
   endpoint carrier.  On the polynomial-carrier truncation
   \(|m|\leq M_j\), its forward and backward multipliers are uniformly tame
   because \(\Theta_jM_j^2\to0\).

Thus the viable target is narrower than the earlier “cell inverse modulo
finite resonances” formulation:

> solve the fast equations as a **forward initial-value problem**, retain
> the full one-phase centre profile as part of the wake state, impose only
> finitely many projected endpoint conditions, and invert heat only on the
> stage-dependent finite harmonic band.  The discarded tail must be put
> into the \(C^\infty\)-flat force.

An exact two-ended inverse on a fixed Gevrey space is ruled out already by
the heat multiplier.  A finite-band, forward invariant graph is not ruled
out.

---

## 1. Frozen fast operator

Work on the mean-zero, divergence-free subspace of
\(\mathbb T^3\).  Let \(n\in\mathbb Z^3\) be primitive, and let

\[
 W(y,t)=Z(n\cdot y,t),\qquad n\cdot Z(\theta,t)=0.
 \tag{1.1}
\]

This includes a linearly polarized carrier

\[
 W=a(t)p\sin(n\cdot y+\delta(t)),\qquad p\cdot n=0,
 \tag{1.2}
\]

and a common-phase two-polarization carrier.  Every field (1.1) is
nonlinearly dark:

\[
 (W\cdot\nabla_y)W=0.
 \tag{1.3}
\]

The frozen spatial linearized Euler operator is

\[
 L_Wz
 =
 \mathbb P_y\big((W\cdot\nabla_y)z+(z\cdot\nabla_y)W\big).
 \tag{1.4}
\]

Its \(L^2\) adjoint on divergence-free fields is

\[
 L_W^*\psi
 =
 \mathbb P_y\big(-(W\cdot\nabla_y)\psi+(\nabla_yW)^T\psi\big).
 \tag{1.5}
\]

Define the one-phase subspace

\[
 X_n
 :=
 \left\{
 z(y)=H(n\cdot y):
 n\cdot H=0,\ \int_{\mathbb T}H\,d\theta=0
 \right\}.
 \tag{1.6}
\]

It contains two independent transverse components for every nonzero phase
harmonic.

### Proposition 1.1: exact zero block

For every carrier (1.1),

\[
 L_WX_n=\{0\},
 \qquad
 X_n\subseteq\ker L_W^*.
 \tag{1.7}
\]

More strongly, if \(\Pi_n\) denotes Fourier projection onto the modes
\(\{mn:m\in\mathbb Z\setminus\{0\}\}\), then

\[
 L_W\Pi_n=0,
 \qquad
 \Pi_nL_W=0.
 \tag{1.8}
\]

#### Proof

If \(z=H(n\cdot y)\in X_n\), then

\[
 (W\cdot\nabla)z=(W\cdot n)H'=0,
 \qquad
 (z\cdot\nabla)W=(z\cdot n)Z'=0,
\]

which proves the first part of (1.7).

For the other half of (1.8), use Fourier variables.  The carrier has
support only at wave numbers \(\ell n\).  The coefficient of \(L_Wz\) at
an output \(mn\) can therefore involve only inputs \((m-\ell)n\).
Incompressibility gives
\(\widehat z((m-\ell)n)\cdot n=0\), while every carrier coefficient is
orthogonal to \(n\).  Both the transport and stretching coefficients
vanish term by term.  Thus \(\Pi_nL_W=0\).

Since \(\Pi_n\) is the \(L^2\)-orthogonal projection onto \(X_n\), for
every \(\psi\in X_n\) and every divergence-free \(z\),

\[
 \langle L_Wz,\psi\rangle
 =
 \langle \Pi_nL_Wz,\psi\rangle
 =0.
 \tag{1.9}
\]

Hence \(X_n\subseteq\ker L_W^*\). \(\square\)

For the sinusoidal carrier, the chain structure is completely explicit.
Writing

\[
 W_{\sigma n}=\frac{\sigma a}{2i}p,\qquad \sigma=\pm1,
\]

one has

\[
 \widehat{L_Wz}(k)
 =
 i\mathbb P_k\sum_{\sigma=\pm1}
 \left[
   (W_{\sigma n}\cdot k)\widehat z(k-\sigma n)
   +\sigma\big(\widehat z(k-\sigma n)\cdot n\big)W_{\sigma n}
 \right].
 \tag{1.10}
\]

When \(k=mn\), both coefficients in every summand vanish by
incompressibility.  Equation (1.10) is the velocity-form counterpart of
the standard Fourier-chain decomposition of the Euler linearization about
a periodic shear.

### Corollary 1.2: no finite-resonance stationary inverse

A necessary condition for

\[
 L_Wz=g
 \tag{1.11}
\]

is

\[
 \Pi_ng=0.
 \tag{1.12}
\]

The obstruction space \(X_n\) is infinite-dimensional.  In particular,
the carrier viscosity

\[
 \Delta_yW\in X_n
 \tag{1.13}
\]

cannot be produced by a stationary corrector.  This recovers and
strengthens the energy pairing obstruction: it gives one compatibility
condition for every transverse phase profile, not only one scalar energy
condition.

The conclusion is unchanged by adding a second polarization with the same
phase.  A common-phase “multiwave” can realize a rank-two covariance, but
it still has the entire zero block \(X_n\).

### 1.1 The simplest non-collinear exact two-wave also has an infinite cokernel

There is a natural attempt to escape the common-phase block.  On
\(\mathbb T^2\subset\mathbb T^3\), take

\[
 W_{A,B}(y)
 =
 \big(B\sin y_2,\ A\sin y_1,\ 0\big).
 \tag{1.14}
\]

With the convention
\(\nabla^\perp\psi=(\partial_2\psi,-\partial_1\psi)\), this is generated by

\[
 \psi=A\cos y_1-B\cos y_2,
 \qquad
 W_{A,B}=\nabla^\perp\psi,
 \qquad
 \omega=\partial_1W_2-\partial_2W_1=\psi.
 \tag{1.15}
\]

Hence \(W_{A,B}\cdot\nabla\omega=0\), so (1.14) is an exact stationary
Euler flow for arbitrary \(A,B\).  Equivalently, the cross interaction of
the two shears is a pressure gradient.  This is the smallest
non-collinear exact multiwave one would try.

It still does not have a finite-codimension stationary inverse.  In
vorticity variables the two-dimensional Euler equation preserves every
Casimir

\[
 {\cal C}_F(\omega)=\int_{\mathbb T^2}F(\omega)\,dy.
 \tag{1.16}
\]

Let \({\cal L}_\omega\) be the vorticity linearization at the equilibrium
(1.15).  Differentiating the identity
\(\langle F'(\omega),X(\omega)\rangle=0\), where \(X\) is the Euler vector
field and \(X(\omega)=0\), gives

\[
 {\cal L}_\omega^*F'(\omega)=0
 \tag{1.17}
\]

for every smooth \(F\).  Since \(\omega\) is nonconstant, these supply
infinitely many independent adjoint nullvectors.  Thus changing from one
phase to the simplest exact orthogonal two-wave merely replaces the
carrier-chain cokernel by the two-dimensional Casimir cokernel.

This obstruction remains valid for the operator on \(\mathbb T^3\), not
only for an artificially restricted equation.  The base has vertical
Fourier wave number zero, so its linearization preserves every vertical
Fourier sector.  Inside the zero sector, the vertical velocity is passive
and does not feed the horizontal two-dimensional block.  The horizontal
2D block is therefore reducing, and the adjoint nullvectors (1.17) extend
by zero to the full 3D operator.

The exact cancellation in (1.14) also fails under the anisotropic Kelvin
strain used by the cascade.  If the two wave numbers are \(k\) and \(m\),
their cross interaction is a gradient only when \(k=m\); the desired
strain evolves them at rates \(\alpha\) and \(\beta\), with
\(\alpha\ne\beta\).  A successful non-collinear cell must therefore be
genuinely three-dimensional and time-dependent, or retain its generated
cross harmonics in the wake.

---

## 2. Time dependence: the cokernel becomes endpoint data

Let \(W(t)\) have the form (1.1), with arbitrary smooth time dependence,
and consider

\[
 \mathcal D_Wz
 :=
 \partial_tz+L_{W(t)}z
 =
 g.
 \tag{2.1}
\]

Projection by \(\Pi_n\) and (1.8) give the exact equation

\[
 \partial_t\Pi_nz=\Pi_ng.
 \tag{2.2}
\]

Consequently,

\[
 \Pi_nz(T)-\Pi_nz(0)
 =
 \int_0^T\Pi_ng(t)\,dt.
 \tag{2.3}
\]

Neither amplitude modulation \(a(t)\) nor a scalar phase shift
\(\delta(t)\) changes this identity: \(X_n\) depends only on the carrier
line, not on its amplitude or origin.

On a fixed periodic fast torus, a smoothly varying integer direction
\(n(t)\) is necessarily constant.  The continuously changing physical
wave vector of a Kelvin mode is instead represented in material fast
coordinates, where its integer carrier line is fixed.  Thus ordinary
Kelvin wave-vector evolution does not evade (2.3).  A genuinely
non-collinear phase conversion can evade this particular fixed-line
identity, but it is then a multiwave transition and must retain the cross
harmonics it creates.

### Proposition 2.1: two-ended no-go

If the complete phase profiles at \(t=0\) and \(t=T\) are both fixed, then
(2.3) is an infinite-dimensional family of compatibility conditions.
In particular, the operator

\[
 \mathcal D_W:
 \{z:z(0)=z(T)=0\}\longrightarrow g
 \tag{2.4}
\]

has infinite-dimensional cokernel.

Suppose \(d\) scalar modulation parameters are allowed to alter \(g\).
Their contributions to the right side of (2.3) span a subspace of \(X_n\)
of dimension at most \(d\).  They cannot cancel a generic defect in
\(X_n\).  Thus a finite number of phase, amplitude, or timing controls
does not turn (2.4) into a finite-codimension inverse.

This is not an obstruction to the forward problem.  With only \(z(0)\)
prescribed, (2.1) is the ordinary linearized Euler initial-value equation.
Standard energy estimates give, for finite \(T\),

\[
 \|z(t)\|_{H^s}
 \leq
 C_{s,T,W}
 \left(
 \|z(0)\|_{H^s}
 +\int_0^t\|g(\tau)\|_{H^s}\,d\tau
 \right),
 \tag{2.5}
\]

with no derivative loss on \(z\).  The distinction between a stationary
or two-ended inverse and a forward propagator is therefore structural, not
terminological.

For a general time-dependent Euler carrier, the adjoint evolution gives
the analogous endpoint identities.  The one-phase case is special only
in making an infinite family of them stationary and explicit.

---

## 3. Exact affine Kelvin--heat subsystem

The centre block can be solved more strongly than by linearization.  Let

\[
 S=\operatorname{diag}(-\alpha,-\beta,\alpha+\beta),
 \qquad \alpha,\beta>0,
 \tag{3.1}
\]

and

\[
 U(x)=Sx,\qquad
 \theta=k(t)x_1,\qquad
 k'=\alpha k.
 \tag{3.2}
\]

For arbitrary periodic profiles \(h(t,\theta)\) and \(q(t,\theta)\), set

\[
 u(x,t)=Sx+e_2h(t,\theta)+e_3q(t,\theta).
 \tag{3.3}
\]

The oscillatory part is divergence-free and its self-interaction vanishes
exactly.  Substitution into Navier--Stokes, with the affine pressure
\(-x\cdot S^2x/2\), gives the closed system

\[
 \begin{aligned}
 \partial_th-\beta h
 &=\nu k(t)^2\partial_\theta^2h+g_2,\\
 \partial_tq+(\alpha+\beta)q
 &=\nu k(t)^2\partial_\theta^2q+g_3.
 \end{aligned}
 \tag{3.4}
\]

Thus (3.3)--(3.4) is an exact nonlinear invariant subsystem, not merely a
formal WKB equation.

Let

\[
 H(s,t):=\int_s^t k(\tau)^2\,d\tau.
 \tag{3.5}
\]

The \(m\)-th Fourier harmonics satisfy

\[
 \begin{aligned}
 h_m(t)
 &=
 e^{\beta(t-s)-\nu m^2H(s,t)}h_m(s)
 +\int_s^t
 e^{\beta(t-\tau)-\nu m^2H(\tau,t)}g_{2,m}(\tau)\,d\tau,\\
 q_m(t)
 &=
 e^{-(\alpha+\beta)(t-s)-\nu m^2H(s,t)}q_m(s)
 +\int_s^t
 e^{-(\alpha+\beta)(t-\tau)-\nu m^2H(\tau,t)}g_{3,m}(\tau)\,d\tau.
 \end{aligned}
 \tag{3.6}
\]

For \(h=A(t)\sin\theta\), (3.4) reduces to

\[
 A'=(\beta-\nu k^2)A,
 \tag{3.7}
\]

the exact viscous Kelvin amplitude equation.

If \(T=(\log r)/\alpha\), then

\[
 k(T)=rk(0),\qquad
 H(0,T)=\frac{k(0)^2(r^2-1)}{2\alpha},
 \tag{3.8}
\]

and, with \(\gamma=\beta/\alpha\),

\[
 \frac{A(T)}{A(0)}
 =
 r^\gamma
 \exp\left(
 -\frac{\nu k(0)^2(r^2-1)}{2\alpha}
 \right).
 \tag{3.9}
\]

Equation (3.9) is the correct endpoint carrier law.  Dropping the
exponential and calling it an error leaves an algebraic-in-stage residual,
not a \(C^\infty\)-flat one.  The heat factor must be included in the
renormalized amplitude sequence exactly.

Because the stage heat actions are summable in the proposed cascade, their
cumulative product converges to a nonzero constant.  Thus (3.9) changes
the multiplicative normalization but does not, by itself, kill the Kelvin
amplification.

### Corollary 3.1: the selected carrier endpoint is an exact submersion

Regard \(h_1\) as a complex coefficient, so that its two real components
encode amplitude and phase.  From (3.6),

\[
 \frac{\partial h_1(T)}{\partial h_1(0)}
 =
 \exp\big(\beta T-\nu H(0,T)\big)\ne0.
 \tag{3.10}
\]

Therefore, after the forcing contribution in (3.6) is fixed, any desired
outgoing first-harmonic coefficient has a unique incoming coefficient.
The inverse norm is

\[
 \exp\big(-\beta T+\nu H(0,T)\big),
 \tag{3.11}
\]

which is uniformly bounded in the small-heat cascade.  Viscosity does not
obstruct the selected endpoint carrier; it only changes its exact return
multiplier.  The obstruction concerns prescribing the *entire* centre
profile at both ends while refusing to carry its remaining harmonics as
wake data.

---

## 4. Full Gevrey endpoint inversion is impossible; finite-band inversion survives

Set

\[
 \Theta:=\nu H(0,T)>0.
 \tag{4.1}
\]

On the phase centre block, forward heat multiplies harmonic \(m\) by
\(e^{-\Theta m^2}\).  Its formal backward inverse multiplies by
\(e^{\Theta m^2}\).

For every ordinary Gevrey class of finite order \(\sigma\geq1\), whose
Fourier weight is of the form

\[
 e^{\rho|m|^{1/\sigma}},
 \tag{4.2}
\]

the backward multiplier \(e^{\Theta m^2}\) is unbounded even if one allows
a fixed loss of Gevrey radius.  Therefore the exact viscous endpoint map is
not onto any fixed ordinary Gevrey space.  A full-profile, two-ended
Gevrey inverse is impossible for every \(\nu>0\).

On the truncated band \(|m|\leq M\), however, the inverse norm is exactly
bounded by

\[
 e^{\Theta M^2}.
 \tag{4.3}
\]

For the polynomial-carrier ledger,

\[
 \ell_j=r^{-j},\qquad
 K_j=j^A,\qquad
 M_j\asymp\frac{j^2}{\log j},
 \tag{4.4}
\]

the Kelvin-compatible stage time and physical carrier are

\[
 \tau_j\asymp\ell_j^{1+\gamma}K_j^{-\gamma},
 \qquad
 k_j\asymp\frac{K_j}{\ell_j}.
\]

Up to the bounded factor caused by the within-stage compression,

\[
 \Theta_j
 \asymp \nu\tau_jk_j^2
 =
 \nu\ell_j^{\gamma-1}K_j^{2-\gamma}.
 \tag{4.5}
\]

Hence

\[
 \Theta_jM_j^2
 \lesssim
 r^{-(\gamma-1)j}
 \frac{j^{A(2-\gamma)+4}}{(\log j)^2}
 \longrightarrow0.
 \tag{4.6}
\]

Thus the centre-block heat inverse is uniformly tame on precisely the
finite harmonic band used by the flat-force truncation.  The parabolic
obstruction is fatal to an exact all-harmonic endpoint inverse, but not to
the stage-dependent finite-band construction.

This also explains why the exponentially small discarded WKB tail is not
an aesthetic convenience.  It is what permits the construction to avoid
backward inversion of the full heat semigroup.  The tail must remain in
the prescribed force, with the already-required
\(e^{-c j^2}\) flatness after physical rescaling.

---

## 5. Consequences for carrier creation and the return theorem

The calculation separates what one carrier can and cannot do.

### What is exact

* Affine compression changes \(k\) by the desired bounded factor.
* The growing transverse polarization changes amplitude by the Kelvin
  factor and the exact heat multiplier (3.9).
* Every collinear phase harmonic is propagated explicitly by (3.6).
* On \(|m|\leq M_j\), forward and backward propagation is uniformly tame
  by (4.6).

### What one phase cannot do

* The pure one-phase sector is nonlinearly dark.  By itself it cannot
  create a non-collinear child or a localized affine return strain.
* A common-phase second polarization changes the covariance but does not
  remove the infinite zero block.
* The simplest exact orthogonal two-wave is a two-dimensional steady Euler
  eigenflow and inherits the infinite Casimir cokernel (1.17).
* A stationary spatial corrector cannot absorb any nonzero defect in
  \(\Pi_ng\), including \(\Delta W\).
* Finite-dimensional modulation cannot enforce two complete endpoint
  profiles against a generic localized hierarchy.

The missing low-frequency return must therefore come from transverse
envelope dynamics, spatially separated colours, or a genuinely
non-collinear packet-plus-wake interaction.  Those mechanisms will
generically feed phase-centre harmonics.  The harmonics must be retained as
state variables, not declared to be finitely many resonances.

### 5.1 Compatibility with a disjoint Gavrilov-bubble wake

A scaled lattice of mutually disjoint compact steady Euler bubbles changes
the global wake problem in a favourable way, but it does not change the
local calculation above.

At a fixed stage, bubbles whose supports are disjoint from the active
transition have exactly zero quadratic interaction with it.  The
one-phase principal carrier inside the active region therefore still has
the zero block (1.8).  A slowly varying bubble background contributes
transport and stretching to the *time-dependent* amplitude equations; in
material coordinates this is part of the forward propagator, not a
uniform stationary inverse for \(L_W\).

Nor can one simply use an \(O(K^{-1})\) envelope or bubble perturbation as
a stationary spectral “splitting” of the normalized zero block.  The
resulting inverse costs \(O(K)\).  Unless the hierarchy is reorganized so
that the centre equation is solved at its natural next order as a
modulation equation, repeated use of that inverse produces exactly the
carrier loss excluded by the Gevrey ledger.  Thus the slow bubble terms
may govern the centre, but they must do so by propagation rather than by
pretending that \(L_W\) has become a uniformly Fredholm spatial operator.

The disjoint bubble lattice may supply exactly the missing global
freedom: the outgoing state need not return to one smooth finite-template
carrier, and an infinite collection of wake-bubble parameters can in
principle receive energy, flux, and off-child endpoint data.  What remains
unproved is the active transition map that converts the propagated
phase-centre profile into

1. the finitely projected smaller child bubble/carrier, and
2. admissible changes of the disjoint wake bubbles,

with a \(C^\infty\)-flat residual.  The present no-go applies only if that
infinite wake freedom is replaced by finitely many endpoint controls.

---

## 6. Revised prize-level theorem target

The following formulation survives all exact calculations in this note.

> **Forward finite-band Kelvin--Reynolds theorem.**  Construct a compact
> Gevrey Euler--Reynolds stage and a phase-resolved packet-plus-wake state
> for which:
>
> 1. the fast equations are solved forward, not by a stationary
>    \(L_W^{-1}\);
> 2. the entire one-phase centre profile through
>    \(|m|\leq M_j\) is included in the renormalized state;
> 3. the exact Kelvin--heat map (3.6), including its endpoint amplitude
>    factor, is part of the one-step return map;
> 4. only finitely many low/child projections are imposed at the outgoing
>    endpoint, while the complementary tail is carried by the wake;
> 5. transverse correctors obey Gevrey bounds with no repeated carrier
>    derivative loss;
> 6. harmonics above \(M_j\) and the physical seams have residual
>    \(e^{-c j^2}\), so their extension to the terminal time is
>    \(C^\infty\)-flat; and
> 7. the projected endpoint derivative is a submersion after the
>    infinite-dimensional centre block has been propagated rather than
>    inverted.

The next decisive algebraic test is a two-colour, non-collinear
linearization: compute whether its **projected** endpoint map onto the
localized affine child has full rank while all newly generated
carrier-chain components can be assigned to the outgoing wake.  A static
Fredholm inverse modulo finitely many modes is no longer a viable target.

---

## Primary references

* R. T. Craik and W. O. Criminale,
  [Evolution of wavelike disturbances in shear flows: a class of exact
  solutions of the Navier--Stokes
  equations](https://doi.org/10.1098/rspa.1986.0061).
* H. R. Dullin and J. Worthington,
  [Stability Theory of the 3-Dimensional Euler
  Equations](https://arxiv.org/abs/1903.09970), especially the Fourier-chain
  linearization about periodic shears.
