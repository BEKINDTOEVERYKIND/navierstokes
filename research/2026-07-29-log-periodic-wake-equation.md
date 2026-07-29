# The log-periodic Euler wake on \(\mathbb S^2\times\mathbb T\)

Date: 2026-07-29

## Claim boundary

This note derives the exact equation for a discretely homogeneous,
log-periodic *stationary outer wake*.  It does not construct such a wake and
does not construct an Euler or Navier--Stokes singularity.

The principal conclusions are:

1. The finite-energy/high-Reynolds window is exactly
   \[
   \frac25<c<\frac12,\qquad
   \gamma=\frac{1-c}{c}\in(1,3/2).
   \]
2. Stationary Euler becomes an autonomous first-order system on
   \(\mathbb S^2\times\mathbb T_L\), \(L=\log r\), with an invertible
   pressure operator throughout this open interval.
3. A nonconstant spatially log-periodic field cannot be an exact
   time-self-similar Euler profile of the corresponding exponent.  The
   similarity drift and stationary Euler terms have different radial
   weights; periodicity forces both to vanish separately.
4. Every stationary log-periodic wake in this window has zero net mass,
   momentum, angular-momentum, energy, and helicity flux through every
   sphere.  It can store a wake, but cannot itself carry a one-way cascade
   flux.  Annular transition regions must carry that flux.
5. Irrotational fields, all generalized Beltrami fields, purely radial
   fields, homogeneous tangential fields, separable rotating tangential
   fields, and a perturbative branch from zero are all excluded.  The
   surviving target is a genuinely three-dimensional, nonaxisymmetric,
   nonperturbative periodic shape of the cylinder system, or a Floquet
   bifurcation from a presently unknown genuinely three-dimensional
   homogeneous state.

## 1. Scaling and the physical window

Let \(\tau=T-t\), and use the Euler-dominated similarity scaling

\[
u(x,t)=\tau^{-(1-c)}U(y),\qquad
y=\frac{x}{\tau^c}.
\tag{1.1}
\]

The Euler profile equation is

\[
(1-c)U+c(y\cdot\nabla)U+(U\cdot\nabla)U+\nabla P=0,
\qquad \operatorname{div}U=0.
\tag{1.2}
\]

Set

\[
\gamma=\frac{1-c}{c}.
\tag{1.3}
\]

A velocity of size \(\ell^{-\gamma}\) on spatial scale \(\ell\) has

\[
\begin{aligned}
\text{energy per scale}&\asymp \ell^{3-2\gamma},\\
\text{local Reynolds number}&\asymp \nu^{-1}\ell^{1-\gamma},\\
\text{turnover time}&\asymp \ell^{1+\gamma}.
\end{aligned}
\tag{1.4}
\]

Thus nested scales \(\ell_j=r^{-j}\to0\) have summable energy precisely when
\(\gamma<3/2\), increasing Reynolds number precisely when \(\gamma>1\), and
summable turnover times for every \(\gamma>-1\).  Equations (1.3)--(1.4)
give

\[
1<\gamma<\frac32
\quad\Longleftrightarrow\quad
\frac25<c<\frac12.
\tag{1.5}
\]

In (1.1), viscosity is lower order by the factor
\(\tau^{1-2c}\to0\).  This is the reason the leading wake equation is Euler,
not stationary Navier--Stokes.

For a blowup embedding, \(U(y)\sim |y|^{-\gamma}\) is the *far-field*
behavior \( |y|\to\infty\).  At fixed physical \(x\ne0\), the similarity
amplitude cancels because \(c\gamma=1-c\).  For the log-periodic ansatz
(2.1), however, what remains is

\[
u(x,t)=|x|^{-\gamma}
V\!\left(\log|x|-c\log(T-t),\frac{x}{|x|}\right).
\tag{1.6}
\]

Thus a nonconstant \(V\) has a perpetually drifting log phase and no
pointwise terminal limit.  This is a second way to see why the object below
can only be a separately embedded stationary wake, not an exact similarity
profile.  Its physical amplitude \(|x|^{-\gamma}\) nevertheless has finite
local energy near \(x=0\) when \(\gamma<3/2\).  Equivalently, the divergent
far-field profile energy is exactly balanced by the similarity prefactor
and the growing profile cutoff.

## 2. Cylinder variables

Write

\[
y=\varrho\omega,\qquad
\varrho=e^s,\qquad
\omega\in\mathbb S^2,
\]

and let \(L=\log r\).  The discretely homogeneous ansatz is

\[
U(y)=e^{-\gamma s}V(s,\omega),\qquad
V(s+L,\omega)=V(s,\omega).
\tag{2.1}
\]

Split \(V\) into radial and tangential parts:

\[
V=v(s,\omega)\omega+w(s,\omega),\qquad
w\cdot\omega=0.
\tag{2.2}
\]

The stationary Euler pressure has the natural form

\[
P(y)=e^{-2\gamma s}q(s,\omega),\qquad
q(s+L,\omega)=q(s,\omega).
\tag{2.3}
\]

All differential operators below on \(\mathbb S^2\) are for the unit round
metric.  Direct spherical-coordinate calculation gives

\[
\boxed{
v_s+(2-\gamma)v+\operatorname{div}_{\!S}w=0.
}
\tag{2.4}
\]

For later use,

\[
\nabla P=e^{-(2\gamma+1)s}
\left((q_s-2\gamma q)\omega+\nabla_{\!S}q\right).
\tag{2.5}
\]

The convective acceleration is

\[
(U\cdot\nabla)U=e^{-(2\gamma+1)s}
\left(\mathcal N_r\omega+\mathcal N_T\right),
\tag{2.6}
\]

where

\[
\begin{aligned}
\mathcal N_r
&=v(v_s-\gamma v)+w\cdot\nabla_{\!S}v-|w|^2,\\
\mathcal N_T
&=v w_s+(1-\gamma)vw+\nabla^{S}_{w}w.
\end{aligned}
\tag{2.7}
\]

Hence stationary Euler on the punctured space is exactly the autonomous
cylinder system

\[
\boxed{
\begin{aligned}
v_s+(2-\gamma)v+\operatorname{div}_{\!S}w&=0,\\
v(v_s-\gamma v)+w\cdot\nabla_{\!S}v-|w|^2
+q_s-2\gamma q&=0,\\
v w_s+(1-\gamma)vw+\nabla^{S}_{w}w+\nabla_{\!S}q&=0.
\end{aligned}
}
\tag{2.8}
\]

When \(v,w,q\) are independent of \(s\), (2.8) is exactly the sphere
system studied by Shvydkoy for homogeneous stationary Euler flows.

The vorticity is

\[
\operatorname{curl}U=e^{-(\gamma+1)s}
\left[
\zeta\,\omega+
\omega\times
\left(w_s+(1-\gamma)w-\nabla_{\!S}v\right)
\right],
\tag{2.9}
\]

where \(\zeta=\operatorname{curl}_{S}w\).  Formula (2.9) makes helicity and
generalized Beltrami constraints explicit on the cylinder.

## 3. Pressure is uniquely determined in the open window

Taking divergence of stationary Euler gives

\[
-\Delta P=\partial_iU_j\,\partial_jU_i.
\tag{3.1}
\]

For (2.3),

\[
\Delta P=e^{-(2\gamma+2)s}\mathcal A_\gamma q,
\]

with

\[
\boxed{
\mathcal A_\gamma
=
\partial_s^2+(1-4\gamma)\partial_s+\Delta_S
+2\gamma(2\gamma-1).
}
\tag{3.2}
\]

On the Fourier--spherical mode
\[
e^{i\kappa_ns}Y_{\ell m}(\omega),\qquad
\kappa_n=\frac{2\pi n}{L},
\]
the multiplier is

\[
(i\kappa_n-2\gamma)(i\kappa_n+1-2\gamma)-\ell(\ell+1).
\tag{3.3}
\]

If it vanishes, its imaginary part first gives \(n=0\), since
\(1-4\gamma\ne0\).  The real part then requires

\[
2\gamma=\ell+1.
\tag{3.4}
\]

There is no integer \(\ell\) satisfying (3.4) for
\(1<\gamma<3/2\); the two neighboring resonances occur at the excluded
endpoints \(\gamma=1\) and \(\gamma=3/2\).  Therefore

> \(\mathcal A_\gamma\) is invertible on smooth
> \(L\)-periodic cylinder functions for every
> \(1<\gamma<3/2\).

Consequently \(q\) is a definite quadratic, nonlocal functional of
\((v,w)\).  There is no free log-periodic harmonic pressure mode in this
window.

## 4. Spatial log-periodicity cannot be an exact similarity profile

For (2.1),

\[
(1-c)U+c(y\cdot\nabla)U
=c e^{-\gamma s}V_s,
\tag{4.1}
\]

because \(1-c=c\gamma\).  Taking divergence first gives a rigorous pressure
decomposition without assuming the pressure scaling.  The Poisson source
in (3.1) has the form
\(e^{-(2\gamma+2)s}F(s,\omega)\), with \(F\) periodic.  By Section 3 there is
a unique periodic \(q_0\) for which
\[
P_0=e^{-2\gamma s}q_0
\]
solves the pressure Poisson equation.  Every other pressure solving that
equation is \(P=P_0+h\), where \(\Delta h=0\).  Define
\[
A=e^{-\gamma s}V_s,\qquad
R=(U\cdot\nabla)U+\nabla P_0
  =e^{-(2\gamma+1)s}R_0(s,\omega),
\tag{4.2}
\]
where \(R_0\) is periodic.  The exact similarity momentum equation is
\[
cA+R+\nabla h=0.
\tag{4.3}
\]

Taking curl removes the arbitrary harmonic pressure:
\[
c\,\operatorname{curl}A+\operatorname{curl}R=0.
\tag{4.4}
\]
The two terms have respectively the forms
\[
e^{-(\gamma+1)s}C_1(s,\omega),\qquad
e^{-(2\gamma+2)s}C_2(s,\omega),
\]
with periodic \(C_1,C_2\).  Evaluate (4.4) at \(s+L\), compare with (4.4) at
\(s\), and use \(\gamma>-1\).  The two radial multipliers are distinct, so
\[
\operatorname{curl}A=\operatorname{curl}R=0.
\tag{4.5}
\]

The \(s\)-derivative of (2.4) gives \(\operatorname{div}A=0\).  Thus \(A\)
is an irrotational, divergence-free discretely homogeneous field of degree
\(-\gamma\).  The Fourier--spherical calculation in Section 7.1 excludes
such a field for \(1<\gamma<3/2\), so
\[
A=0,\qquad V_s=0.
\tag{4.6}
\]

The choice of \(P_0\) and (3.1) also gives \(\operatorname{div}R=0\).
Together with (4.5), \(R\) is an irrotational, divergence-free discretely
homogeneous field of degree \(-(2\gamma+1)\).  Its harmonic potential has
degree \(-2\gamma\); the only possible periodic spherical resonance is
\(2\gamma=\ell+1\), precisely (3.4).  Section 3 excludes it in the open
window, hence \(R=0\).  Equation (4.3) finally gives \(\nabla h=0\).
Therefore the remaining equation is exactly stationary Euler with its
natural pressure:
\[
\mathcal E_\gamma(V,q_0)=0,\qquad V_s=0.
\tag{4.7}
\]

Thus:

> An exact similarity profile of the form (2.1) is necessarily continuously
> homogeneous.  Spatial log-periodicity alone does not evade homogeneous
> profile rigidity.

This is different from the usual discretely self-similar Euler ansatz in
the literature, which is periodic in *similarity time*.  A nonconstant
solution of (2.8) can only serve here as a stationary annular wake, with
separate time-dependent transition regions.

## 5. Exact flux constraints

Let

\[
b=q+\frac12(v^2+|w|^2)
\tag{5.1}
\]

be the cylinder Bernoulli function.  From stationary Euler,

\[
\boxed{
v(b_s-2\gamma b)+w\cdot\nabla_{\!S}b=0.
}
\tag{5.2}
\]

The physical Bernoulli function is \(e^{-2\gamma s}b\).  Let
\(n=x/|x|\), \(\Omega=\operatorname{curl}U\),
\(B=P+\frac12|U|^2\), and \(T=U\otimes U+PI\).  The five sphere fluxes used
below are
\[
\begin{aligned}
\mathsf M(\varrho)&=\int_{S_\varrho}U\cdot n,\\
\boldsymbol{\mathsf P}(\varrho)&=\int_{S_\varrho}Tn,\\
\boldsymbol{\mathsf J}(\varrho)&=\int_{S_\varrho}x\times(Tn),\\
\mathsf E(\varrho)&=\int_{S_\varrho}B\,U\cdot n,\\
\mathsf H(\varrho)&=\int_{S_\varrho}
\left[(U\cdot\Omega)U+
\left(P-\frac12|U|^2\right)\Omega\right]\cdot n .
\end{aligned}
\tag{5.3}
\]
Surface measure is understood in these formulas.  The local stationary
conservation laws make each flux independent of \(\varrho\).  Discrete
homogeneity also multiplies it by a power after \(s\mapsto s+L\).  A
constant flux with a nonunit multiplier must vanish.  Explicitly:

\[
\begin{array}{c|c|c}
\text{flux}&\text{radial scaling}&\text{zero unless}\\ \hline
\text{mass}&\varrho^{\,2-\gamma}&\gamma=2\\
\text{momentum}&\varrho^{\,2-2\gamma}&\gamma=1\\
\text{angular momentum}&\varrho^{\,3-2\gamma}&\gamma=3/2\\
\text{energy}&\varrho^{\,2-3\gamma}&\gamma=2/3\\
\text{helicity}&\varrho^{\,1-3\gamma}&\gamma=1/3 .
\end{array}
\tag{5.4}
\]

None of the exceptional exponents lies in the open window.  In particular,
for every \(s\),

\[
\int_{\mathbb S^2}v(s,\omega)\,d\omega=0,
\qquad
\int_{\mathbb S^2}b(s,\omega)v(s,\omega)\,d\omega=0.
\tag{5.5}
\]

The second identity is zero radial energy flux.  The local helicity
conservation law used above is

\[
\partial_t(U\cdot\Omega)+
\operatorname{div}\left[
(U\cdot\Omega)U+
\left(P-\frac12|U|^2\right)\Omega
\right]=0,
\qquad \Omega=\operatorname{curl}U.
\tag{5.6}
\]

Thus an exact stationary log-periodic wake cannot be the active one-way
energy cascade.  It can only be a zero-net-flux storage state.  The active
annular transition must carry the energy and helicity exchanges.

There is also a useful pressure identity.  Integrating the radial equation
in (2.8) over one cylinder period and using (2.4) gives

\[
\boxed{
2\gamma\int q
=2(1-\gamma)\int v^2-\int|w|^2.
}
\tag{5.7}
\]

All integrals in (5.7)--(5.8) are over
\(\mathbb T_L\times\mathbb S^2\).  Hence every nonzero solution in the
window has strictly negative mean pressure.  Meanwhile,

\[
\boxed{
\int b
=
\frac{2-\gamma}{2\gamma}\int v^2
+
\frac{\gamma-1}{2\gamma}\int|w|^2
>0.
}
\tag{5.8}
\]

Therefore a nonzero solution has \(q<0\) somewhere and \(b>0\) somewhere.
Any ansatz with nonnegative pressure or nonpositive Bernoulli function is
excluded.  In particular, \(b\equiv0\) forces \(V\equiv0\).  Here
``generalized Beltrami'' means
\(\operatorname{curl}U=m(x)U\), with \(m\) allowed to vary and change sign.
Such a stationary Euler field has spatially constant physical Bernoulli
function, and hence \(b\equiv0\) after removal of the irrelevant constant
pressure.  Thus (5.8) excludes every nonzero generalized Beltrami wake in
this window.

Finally, the helicity contained in a logarithmic shell scales like
\(\varrho^{2-2\gamma}\).  Since \(\gamma>1\), nested shells toward the
origin have increasing unsigned helicity.  A finite- or zero-helicity
embedding therefore needs exact bihelical cancellation in each log cell
(for example an orientation-reversing symmetry).

## 6. Streamline topology forced by Bernoulli

On the cylinder, particle paths obey

\[
\dot s=v,\qquad \dot\omega=w.
\tag{6.1}
\]

Equation (5.2) gives along such a path

\[
\frac{d}{dt}b=2\gamma v b
=2\gamma\dot s\,b.
\tag{6.2}
\]

If a streamline closes on the quotient cylinder after winding \(k\ne0\)
times around the \(s\)-circle, periodicity gives

\[
b_{\rm final}=b_{\rm initial}
=e^{2\gamma kL}b_{\rm initial},
\]

and hence

\[
b=0
\quad\text{on every log-radially winding streamline.}
\tag{6.3}
\]

Thus a wake with a large family of radial-throughput streamlines is forced
toward its zero-Bernoulli set.  Formula (5.8) rules out the extreme case in
which this set fills the entire cylinder.  A viable wake must therefore
have substantial non-winding/recirculating streamline regions; a globally
generalized-Beltrami escape is unavailable.

## 7. Tests of explicit classes

### 7.1 Irrotational fields

Because \(\mathbb R^3\setminus\{0\}\) is simply connected, an irrotational
and divergence-free field is \(U=\nabla\Phi\) with \(\Delta\Phi=0\).
Discrete homogeneity gives

\[
\Phi=e^{(1-\gamma)s}\phi(s,\omega),
\qquad \phi(s+L,\omega)=\phi(s,\omega).
\tag{7.1}
\]

Put \(\lambda=1-\gamma\in(-1/2,0)\).  Harmonicity becomes

\[
\phi_{ss}+(2\lambda+1)\phi_s+\Delta_S\phi
+\lambda(\lambda+1)\phi=0.
\tag{7.2}
\]

For a Fourier--spherical mode, (7.2) requires

\[
(\lambda+i\kappa_n)(\lambda+i\kappa_n+1)=\ell(\ell+1).
\tag{7.3}
\]

Since \(2\lambda+1=3-2\gamma>0\), the imaginary part forces \(n=0\).
The real equation then requires
\(\lambda=\ell\) or \(\lambda=-\ell-1\), impossible for
\(\lambda\in(-1/2,0)\).  Therefore

\[
U=0.
\tag{7.4}
\]

There is no irrotational log-periodic wake in the target window.

### 7.2 Beltrami fields

Suppose more generally that
\[
\operatorname{curl}U=m(x)U.
\tag{7.5}
\]
Then \(U\times\operatorname{curl}U=0\), so the stationary Euler identity
\[
U\times\operatorname{curl}U
=\nabla\left(P+\frac12|U|^2\right)
\]
shows that the physical Bernoulli function is constant on the connected
punctured space.  Subtracting this constant from \(P\) gives
\[
q=-\frac12(v^2+|w|^2),\qquad b=0.
\]
The positive identity (5.8) forces \(v=w=0\).  This argument permits
variable and sign-changing \(m\), and therefore excludes all generalized
Beltrami and bihelical-parallel variants.

For completeness, discrete homogeneity itself forces
\[
m(\varrho,\omega)=\varrho^{-1}\widetilde m(s,\omega)
\]
wherever \(U\ne0\), with \(\widetilde m\) periodic.  In particular a
constant-\(\mu\ne0\) Beltrami equation is already incompatible with the
radial weights; \(\mu=0\) is the irrotational case.

### 7.3 Purely radial fields

If \(w=0\), incompressibility alone gives

\[
v_s+(2-\gamma)v=0.
\]

Since \(\gamma\ne2\), periodicity forces \(v=0\).  No purely radial
log-periodic wake exists.

### 7.4 Tangential fields

For \(v=0\), the cylinder equations reduce to

\[
\operatorname{div}_{S}w=0,\qquad
\nabla^S_w w+\nabla_Sq=0,\qquad
q_s-2\gamma q=|w|^2.
\tag{7.6}
\]

If \(w,q\) are independent of \(s\), Shvydkoy's tangential rigidity theorem
applies: there are no \(C^1(\mathbb S^2)\) nonzero tangential homogeneous
solutions for \(\gamma>-1\).  In particular there are none in the target
window.

The simplest genuinely log-periodic test is a rotating latitude field

\[
w=a(s)\sin\theta\,e_\varphi .
\tag{7.7}
\]

The tangential equation gives

\[
q=\frac12a(s)^2\sin^2\theta+C(s).
\]

The radial equation then yields

\[
a'=(1+\gamma)a,\qquad C'=2\gamma C.
\tag{7.8}
\]

No nonzero solution of (7.8) is periodic.

More generally, for a separable family
\[
w=a(s)w_0(\omega),\qquad
q=a(s)^2q_0(\omega)+C(s),
\]
where \((w_0,q_0)\) is a steady Euler flow on \(\mathbb S^2\), periodicity
either fails exponentially or forces
\[
|w_0|^2+2\gamma q_0=\text{constant}.
\]
After absorbing the constant into \(C\), this is exactly the homogeneous
tangential relation excluded by Shvydkoy.  Thus no nonzero separable
tangential branch survives.

Equation (7.6) still leaves open a nonseparable periodic path through the
infinite-dimensional set of steady Euler flows on \(\mathbb S^2\).  Such a
path would have to be found nonperturbatively; it is not generated by a
simple rotation or amplitude modulation.  Once found, its overall velocity
amplitude could of course be rescaled by Euler's quadratic homogeneity.

### 7.5 Bifurcation from zero

After pressure is eliminated using the invertible operator
\(\mathcal A_\gamma\), stationary Euler on the cylinder is homogeneous
quadratic in \(V\).  If

\[
V_\varepsilon=\varepsilon V_1+o(\varepsilon),
\]

then the pressure is \(q_\varepsilon=\varepsilon^2q_2+o(\varepsilon^2)\),
and the order-\(\varepsilon^2\) equation says that \(V_1\) itself solves the
full stationary cylinder equation.

Therefore zero has no useful linear bifurcation theory:

> A small branch from zero exists only if a nonzero exact cylinder solution
> is already available.  Amplitude smallness gives no perturbative lever
> because Euler is quadratically homogeneous.

## 8. The only conventional Floquet bifurcation target

Let \((v_0,w_0,q_0)\) be a smooth \(s\)-independent homogeneous solution.
For a nonzero log Fourier mode
\[
(\phi,X,\pi)e^{i\kappa s},\qquad \kappa=2\pi n/L,
\]
the linearized cylinder equations are

\[
\begin{aligned}
(i\kappa+2-\gamma)\phi+\operatorname{div}_SX&=0,\\
(i\kappa-2\gamma)v_0\phi
+X\cdot\nabla_Sv_0+w_0\cdot\nabla_S\phi
-2w_0\cdot X+(i\kappa-2\gamma)\pi&=0,\\
i\kappa v_0X
+(1-\gamma)(\phi w_0+v_0X)
+\nabla^S_Xw_0+\nabla^S_{w_0}X+\nabla_S\pi&=0.
\end{aligned}
\tag{8.1}
\]

A conventional log-periodic bifurcation requires a nontrivial kernel of
(8.1) for some \(n\ne0\), followed by an adjoint transversality condition
as \(L\) or \(\gamma\) varies.

The available simple bases do not work:

* irrotational bases occur only at the discrete integer homogeneities
  identified by Shvydkoy, not at \(1<\gamma<3/2\);
* tangential bases are excluded for \(\gamma>-1\);
* Shvydkoy excludes all smooth axisymmetric homogeneous solutions for
  \(0<\gamma<2\).

Thus any base for (8.1) must itself be a genuinely three-dimensional,
nonaxisymmetric homogeneous Euler state.  No such explicit smooth base is
currently available in this interval.  In fact, Shvydkoy conjectures that
for every homogeneity \(\gamma>-1\) the only \(C^1\) homogeneous solutions
are the irrotational solutions at the permitted integer exponents.  If that
conjecture is correct, the Floquet base route is empty throughout the open
window.  The other possibility is a nonperturbatively found periodic shape
of (2.8), not connected in shape to any equilibrium; its overall amplitude
is freely rescalable.

## 9. Physical annular wake and pressure matching

Suppose a nonzero periodic cylinder solution exists and is embedded as a
stationary physical wake behind a time-dependent active core.  Its repeated
log cells describe nested annuli around the singular point,

\[
\ell_{j+1}<|x-x_*|<\ell_j,\qquad \ell_j=r^{-j}.
\]

Then

\[
E_j\asymp \ell_j^{3-2\gamma},
\qquad
\operatorname{Re}_j\asymp\nu^{-1}\ell_j^{1-\gamma},
\qquad
\Delta t_j\asymp\ell_j^{1+\gamma}.
\tag{9.1}
\]

The energy and stage times sum, while the Reynolds numbers diverge.  In the
similarity coordinates of the active core, the already deposited wake
appears in successive outer cells \(|y|\sim r^j\).  This is the attractive
part of the construction, but Section 4 and (1.6) show that the deposition
cannot be replaced by one exact self-similar formula.

The pressure is not a free matching parameter.  Section 3 fixes
\[
P_j\asymp\ell_j^{-2\gamma}q
\]
through an invertible elliptic operator, and Section 5 forces zero flux in
the exactly periodic interior of every cell.  Consequently the annular
interfaces must do all of the following:

1. carry the nonzero energy transfer needed by the active cascade;
2. match the order-\(\ell_j^{-2\gamma}\) pressure without introducing an
   uncontrolled nonlocal multipole;
3. retain the mismatch as part of the outgoing wake rather than discard it
   as an external force.

The zero momentum and angular-momentum fluxes in (5.4) are necessary moment
compatibilities for pressure localization, but they do not by themselves
control the full pressure tail.

### Conditional slowly modulated WKB refinement

Exact discrete self-similarity is not required for an all-order
construction.  Let the internal carrier count on annulus \(j\) be

\[
K_j=j^A,\qquad
M_j=\left\lfloor\frac{j^2}{\log j}\right\rfloor .
\tag{9.2}
\]

Assume, conditionally, an order-\(M\) WKB construction with Gevrey bound

\[
\|\mathrm{coefficient}_M\|
\leq C^M(M!)^\sigma
\]
and residual bounded by
\[
C^{M_j}(M_j!)^\sigma K_j^{-M_j}.
\]
Stirling's formula gives

\[
\log\|\mathrm{residual}_j\|
\leq
-(A-2\sigma+o(1))j^2.
\tag{9.3}
\]

Hence \(A>2\sigma\) gives an \(e^{-c j^2}\) residual.  Multiplication by any
fixed number of spatial or temporal derivative costs only
\(e^{O(j)}j^{O(1)}\), so (9.3) is flat at the accumulation time.  Meanwhile

\[
\frac{K_{j+1}}{K_j}=1+\frac{A}{j}+O(j^{-2}).
\tag{9.4}
\]

This replaces an exponential carrier jump by a slowly nonautonomous,
asymptotically discretely self-similar modulation.  It does **not** produce
the leading Euler wake or the annular transition: it is useful only after a
base cylinder cell and Gevrey-tame corrector theory have been proved.

## 10. Decisive next theorem

The surviving long-shot can now be stated without ambiguity.

> **Cylinder-wake theorem target.**  For some
> \(\gamma\in(1,3/2)\) and \(L>0\), construct a smooth nonconstant
> \(L\)-periodic solution of (2.8) which is genuinely three-dimensional,
> nonaxisymmetric, has an orientation-reversing zero-helicity symmetry, and
> admits pressure and annular transition estimates uniform under scaling.

Such a solution must satisfy all zero-flux conditions in Section 5 and
therefore serves only as the wake state.  A second theorem must construct a
time-dependent annular transition whose endpoint is the next scaled wake
and whose nonzero energy transfer is confined to that transition.

The conventional route to the first theorem is to find a genuinely
three-dimensional homogeneous base and a simple nonzero-log-mode kernel of
(8.1).  If no such base or kernel exists, the alternative is a
nonperturbative periodic shape of (2.8).  A Liouville theorem excluding
both would decisively kill this outer-wake route.

No GPU calculation is warranted until either an explicit homogeneous base
for (8.1) or a finite-dimensional symmetry reduction of (2.8) is found.

## Primary sources

* R. Shvydkoy, *Homogeneous solutions to the 3D Euler system*:
  https://arxiv.org/abs/1510.03378
* D. Chae and R. Shvydkoy, *On formation of a locally self-similar collapse
  in the incompressible Euler equations*:
  https://arxiv.org/abs/1201.6009
* D. Chae and T.-P. Tsai, *On discretely self-similar solutions of the Euler
  equations*:
  https://arxiv.org/abs/1304.7414
