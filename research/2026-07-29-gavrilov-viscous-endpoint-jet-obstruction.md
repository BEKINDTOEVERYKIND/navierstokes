# A second-jet obstruction to locally dressed Gavrilov endpoints

Date: 2026-07-29

## Result and claim boundary

Consider an exact compact steady Euler bubble at a stage endpoint.  Allow
its amplitude, scale, and centre to vary, and allow an otherwise arbitrary
time-dependent divergence-free corrector, but require the bubble and all
of its corrector jets to remain in one of a collection of mutually
disjoint balls.

The zeroth Navier--Stokes residual can always be cancelled by the first
time derivative of the corrector.  The next equation cannot.  If the
residual and its first time derivative vanish at the endpoint, then every
individual endpoint bubble \(W\) must satisfy

\[
 {\mathsf G}(W)_{ij}
 :=
 \int_{\mathbb R^3}\partial_kW_i\,\partial_kW_j\,dx
 =c\,\delta_{ij}.
\tag{0.1}
\]

This is a pressure-quadrupole condition, not the scalar energy condition.
It is unaffected by amplitude, dilation, translation, or rotation.

A canonical fixed-profile thin Gavrilov bubble does not satisfy (0.1).
For the localization

\[
 W_\varepsilon=\chi(p/\varepsilon)u,
 \qquad 0\ne\chi\in C^\infty_c((0,\infty)),
\tag{0.2}
\]

around Gavrilov's seed circle, axial symmetry gives

\[
 {\mathsf G}_\varepsilon
 =
 \operatorname{diag}(G_{\perp,\varepsilon},
                     G_{\perp,\varepsilon},
                     G_{\parallel,\varepsilon}),
\]

and the explicit thin-torus limit is

\[
 G_{\parallel,\varepsilon}-G_{\perp,\varepsilon}
 =
 \pi^2\varepsilon\int_0^\infty\chi(s)^2\,ds
 +o(\varepsilon)>0.
\tag{0.3}
\]

Consequently a bare thin Gavrilov endpoint cannot be made
Navier--Stokes-flat even through the first derivative of the residual by
correctors confined to its disjoint bubble ball.  The obstruction occurs
at the second velocity jet.  The missing jet has a nonzero harmonic
pressure quadrupole and therefore an algebraic spatial tail.

This does **not** obstruct a transition carrying a nonlocal pressure/centre
wake, a nonzero endpoint oscillatory bath, or a global corrector joining
the nominally disjoint bubbles.  It proves that at least one of those
features is mandatory.  It does not construct or rule out a
Navier--Stokes singularity.

---

## 1. Exact residual of the modulated disjoint-bubble ansatz

Let

\[
 (W\cdot\nabla)W+\nabla P=0,\qquad
 \operatorname{div}W=0,
 \qquad W,P\in C^\infty_c(\mathbb R^3),
\tag{1.1}
\]

and, for fixed \(Q\in O(3)\), put

\[
 y=Q^T\frac{x-c(t)}{\ell(t)},\qquad
 U(t,x)=A(t)QW(y),\qquad
 \Pi(t,x)=A(t)^2P(y).
\tag{1.2}
\]

The Euler part cancels exactly.  With the pressure included, the full
Navier--Stokes residual is

\[
\boxed{
 \begin{aligned}
 {\cal R}(U,\Pi)
 &:=\partial_tU+(U\cdot\nabla)U+\nabla\Pi-\nu\Delta U\\
 &=A'QW
   -A\frac{\ell'}{\ell}Q(y\cdot\nabla W)
   -\frac{A}{\ell}Q\big((Q^Tc')\cdot\nabla W\big)
   -\frac{\nu A}{\ell^2}Q\Delta W .
 \end{aligned}}
\tag{1.3}
\]

For finitely many copies whose joint velocity--pressure supports are
disjoint, the residual is exactly the sum of (1.3); there are no quadratic
cross terms.

Now add an arbitrary local divergence-free corrector \(z\) and use the
projected notation

\[
 B(v,w):=\mathbb P((v\cdot\nabla)w).
\tag{1.4}
\]

Then

\[
 {\cal N}(U+z)
 =
 {\cal N}(U)
 +\partial_tz-\nu\Delta z
 +B(U,z)+B(z,U)+B(z,z).
\tag{1.5}
\]

Suppose that \(z(t_*)=0\).  Equation (1.5) shows that
\(\partial_tz(t_*)\) can always cancel (1.3).  This observation can make
the endpoint look unobstructed if only the residual value is checked.
It does not provide freedom in the total velocity jet: once
the projected residual vanishes,

\[
 \partial_t(U+z)(t_*)=\nu\Delta(U+z)(t_*)
\tag{1.6}
\]

because the endpoint state is steady Euler.  Thus all choices of
\(A',\ell',c'\) merely change how the fixed right side of (1.6) is split
between the parameter tangent and \(z_t\).

The scalar work identity at the endpoint is correspondingly

\[
 \frac{d}{dt}\frac12\|U+z\|_2^2\bigg|_{t=t_*}
 =
 -\nu\|\nabla U(t_*)\|_2^2.
\tag{1.7}
\]

Amplitude decay can satisfy this one scalar equation.  The obstruction
below is the trace-free part of a \(3\times3\) tensor and survives that
choice.

## 2. The formal Navier--Stokes jet is unique

Set \(t_*=0\) and write

\[
 u_n:=\partial_t^nu(0).
\]

If every time jet of the projected force vanishes at zero, differentiating
the equation gives the exact recursion

\[
\boxed{
 u_{n+1}
 =
 \nu\Delta u_n
 -
 \sum_{a=0}^n\binom na B(u_a,u_{n-a}).
}
\tag{2.1}
\]

There is no endpoint shooting parameter in (2.1).  A decomposition into
similarity parameters and correctors only decomposes this uniquely fixed
total jet.

For a steady Euler endpoint \(u_0=W\),

\[
 u_1=\nu\Delta W
\tag{2.2}
\]

and

\[
 u_2
 =
 \nu^2\Delta^2W
 -
 \nu\mathbb P\left(
  (W\cdot\nabla)\Delta W
  +(\Delta W\cdot\nabla)W
 \right).
\tag{2.3}
\]

Every jet in (2.1) is a well-defined spatially global element of

\[
 H^\infty_\sigma(\mathbb R^3)
 :=
 \bigcap_{s\ge0}H^s_\sigma(\mathbb R^3).
\]

Here the subscript \(\sigma\) denotes the solenoidal subspace.

For example, for \(s>3/2\),

\[
 \|u_{n+1}\|_{H^s}
 \le
 \nu\|u_n\|_{H^{s+2}}
 +C_s\sum_{a=0}^n\binom na
   \|u_a\|_{H^{s+1}}\|u_{n-a}\|_{H^{s+1}}.
\tag{2.4}
\]

The standard vector-valued Borel construction, with successively
shrinking time cutoffs, realizes any such sequence as the time jet of a
path in \(C^\infty_tH^\infty_\sigma\).  Substitution into the polynomial
map \({\cal N}\) and Taylor's theorem then make its residual flat at the
endpoint in every \(H^s\), hence every spatial \(C^m\), norm.  The
Borel-realized path is not asserted to solve Navier--Stokes; only its
residual has all endpoint jets equal to zero.

This global formal construction is not the desired disjoint construction.
The Leray projection in (2.3) generally creates algebraic spatial tails.
The next theorem identifies its first unavoidable tail for a thin
Gavrilov bubble.

There is a useful all-order necessary condition for keeping the recursion
local.  In unprojected variables the \(n\)-th pressure jet obeys

\[
 -\Delta p_n
 =
 \partial_i\partial_j(T_n)_{ij},
 \qquad
 T_n
 =
 \sum_{a=0}^n\binom na u_a\otimes u_{n-a}.
\tag{2.5}
\]

If all velocity jets through \(u_{n+1}\) are confined to one ball and the
\(n\)-th force jet vanishes, then \(\nabla p_n\) vanishes outside that
ball.  Consequently, for every harmonic polynomial \(q\),

\[
 \int_{\mathbb R^3}T_n:\nabla^2q\,dx=0.
\tag{2.6}
\]

These harmonic-moment identities are necessary, not asserted sufficient,
for a compact all-order jet.  Theorem 3.1 is the degree-two condition in
(2.6) at \(n=1\).

## 3. Local second-jet theorem

> **Theorem 3.1 (local second-jet obstruction).**
> Let \(B_1,\ldots,B_M\) be open balls with pairwise disjoint closures, and
> let
> \[
>  W=\sum_{m=1}^MW_m,\qquad
>  W_m\in C^\infty_c(B_m;\mathbb R^3),
> \]
> where every \(W_m\) is divergence-free and is a compactly supported
> steady Euler flow.  Let \(\nu>0\).
>
> Suppose there is a \(C^2_tC^\infty_x\) divergence-free path \(u(t)\), a
> \(C^1_tC^\infty_x\) pressure \(p(t)\), and a \(C^1_tC^\infty_x\) force
> \(f(t)\) such that, for \(t\) near zero,
> \[
> \begin{aligned}
>  &\operatorname{supp}u(t)\subset\bigcup_{m=1}^MB_m,\\
>  &u(0)=W,\\
>  &\partial_tu+(u\cdot\nabla)u+\nabla p-\nu\Delta u=f,\\
>  &f(0)=\partial_tf(0)=0.
> \end{aligned}
> \tag{3.1}
> \]
> Normalize the pressure by a time-dependent spatial constant on the
> unbounded component.  Then, for every \(m\),
> \[
>  {\mathsf G}(W_m)=c_mI
> \tag{3.2}
> \]
> for some \(c_m\ge0\).
>
> In particular, if one \(W_m\) has nonisotropic \({\mathsf G}(W_m)\),
> there is no such path.  This remains true if \(u\) is presented as a
> time-dependent amplitude/scale/translation path plus arbitrary
> divergence-free correctors supported in the same balls and vanishing in
> value at the endpoint.

### Proof

Since \(W\) is steady Euler and its components do not interact,
the zeroth equation in (3.1), after applying the Leray projection, gives

\[
 u_1:=\partial_tu(0)=\nu\Delta W.
\tag{3.3}
\]

Differentiate (3.1) once and write \(p_1=\partial_tp(0)\).  Since
\(\partial_tf(0)=0\),

\[
 u_2
 +(u_1\cdot\nabla)W
 +(W\cdot\nabla)u_1
 +\nabla p_1
 -\nu\Delta u_1
 =0.
\tag{3.4}
\]

Taking divergence gives

\[
 -\Delta p_1
 =
 \partial_i\partial_jS_{ij},
\qquad
 S
 =
 W\otimes u_1+u_1\otimes W.
\tag{3.5}
\]

All velocity terms in (3.4) vanish outside the union of the balls.
Because \(u_2\) also vanishes there, (3.4) gives
\(\nabla p_1=0\) there.  The complement of finitely many disjoint closed
balls is connected, so the normalization makes \(p_1=0\) outside their
union.  In a collar of every \(\partial B_m\), \(p_1\) is zero.  It
therefore splits into smooth compactly supported pieces

\[
 p_1=\sum_mp_{1,m},\qquad p_{1,m}\in C^\infty_c(B_m),
\tag{3.6}
\]

and (3.5) splits componentwise with

\[
 S_m
 =
 \nu\left(
 W_m\otimes\Delta W_m+\Delta W_m\otimes W_m
 \right).
\tag{3.7}
\]

Let \(H=H^T\) be any constant trace-free matrix and put

\[
 q_H(x)=\frac12x\cdot Hx.
\]

This is a harmonic polynomial.  Testing the \(m\)-th equation in (3.5)
against \(q_H\), and integrating twice by parts, gives

\[
 0
 =
 \int(-\Delta p_{1,m})q_H\,dx
 =
 \int S_m:H\,dx.
\tag{3.8}
\]

On the other hand,

\[
 \begin{aligned}
 \int(S_m)_{ij}\,dx
 &=
 \nu\int\left(
  (W_m)_i\Delta(W_m)_j
  +\Delta(W_m)_i(W_m)_j
 \right)dx\\
 &=
 -2\nu\int
 \partial_k(W_m)_i\,\partial_k(W_m)_j\,dx\\
 &=-2\nu{\mathsf G}(W_m)_{ij}.
 \end{aligned}
\tag{3.9}
\]

Thus \({\mathsf G}(W_m):H=0\) for every trace-free symmetric \(H\).
The orthogonal complement of the trace-free symmetric matrices is
\(\mathbb RI\), proving (3.2). \(\square\)

The same proof applies on \(\mathbb T^3\) when the disjoint balls lie in
coordinate charts.  After \(p_{1,m}\) is shown to vanish in a collar of
its ball, lift that compactly supported piece to \(\mathbb R^3\) and use
the same quadratic harmonic test.

### What the theorem does and does not use

The proof permits arbitrary time derivatives of every local corrector.
It does not assume that the corrector is stationary, finite-dimensional,
small, or chosen by a particular anti-divergence.  It only uses:

1. an exact steady Euler endpoint;
2. vanishing of the force through its first time derivative;
3. positive viscosity; and
4. confinement of the first two velocity jets to disjoint balls.

If only the total far-field quadrupole were tested, several rotated
bubbles could make \(\sum_m{\mathsf G}(W_m)\) isotropic.  The support
argument in (3.6) is stronger: under genuinely disjoint local correction,
each ball must pass the condition separately.

## 4. The thin Gavrilov bubble fails the tensor condition

The following computation uses only the local Taylor expansion in
Gavrilov's construction.

> **Lemma 4.1 (anisotropy of a thin Gavrilov localization).**
> Normalize the radius of Gavrilov's seed circle to one and use cylindrical
> coordinates \((\rho,\varphi,z)\).  Let \((u,p)\) be the analytic local
> flow in a neighborhood of the circle, and let
> \[
>  W_\varepsilon=\chi(p/\varepsilon)u,
>  \qquad 0\ne\chi\in C^\infty_c((0,\infty)).
> \tag{4.1}
> \]
> For all sufficiently small \(\varepsilon>0\), this is a smooth compact
> steady Euler flow.  Its gradient covariance is axisymmetric and obeys
> \[
>  G_{\parallel,\varepsilon}-G_{\perp,\varepsilon}
>  =
>  \pi^2\varepsilon\int_0^\infty\chi(s)^2\,ds
>  +o(\varepsilon).
> \tag{4.2}
> \]
> In particular, it is not a scalar matrix.

### Proof

Put \(r^2=(\rho-1)^2+z^2\).  Gavrilov's Taylor expansion and the definition
of his function \(H\) give

\[
 \alpha(\rho,z)=2r^2+O(r^3),\qquad
 p=\frac{\alpha}{4},\qquad
 H(\alpha)=4\alpha+O(\alpha^2).
\tag{4.3}
\]

The cylindrical components of formula (6) in Gavrilov's paper therefore
have the \(C^1\) expansions

\[
 u_\rho=z+O(r^2),\qquad
 u_z=-(\rho-1)+O(r^2),\qquad
 u_\varphi=\frac{r}{\sqrt2}+O(r^2).
\tag{4.4}
\]

Set

\[
 \eta_1=\frac{\rho-1}{\sqrt\varepsilon},\qquad
 \eta_2=\frac{z}{\sqrt\varepsilon},\qquad
 R=|\eta|.
\]

On the fixed annulus containing the support of \(\chi(R^2/2)\),
\(\varepsilon^{-1/2}W_\varepsilon\) converges in \(C^1\) to the field with
cylindrical components

\[
 V_\rho=\chi(R^2/2)\eta_2,\qquad
 V_z=-\chi(R^2/2)\eta_1,\qquad
 V_\varphi=\frac1{\sqrt2}\chi(R^2/2)R.
\tag{4.5}
\]

The \(\rho,z\) derivatives contribute at order \(\varepsilon\) to
\({\mathsf G}\), because
\(\rho\,d\rho\,dz\,d\varphi
 =\varepsilon(1+O(\sqrt\varepsilon))\,d\eta\,d\varphi\).
The azimuthal derivatives contribute only \(O(\varepsilon^2)\).

Write

\[
 g(R)=R\chi(R^2/2).
\]

In polar coordinates on the \(\eta\)-plane,

\[
 V_\rho=g(R)\sin\theta,\qquad
 V_z=-g(R)\cos\theta,\qquad
 V_\varphi=\frac{g(R)}{\sqrt2}.
\]

Let \(E_\rho,E_z,E_\varphi\) denote the two-dimensional Dirichlet
energies of these three scalar components.  Direct integration in
\(\theta\) gives

\[
 \begin{aligned}
 E_\rho=E_z
 &=\pi\int_0^\infty
   \left(g'(R)^2+\frac{g(R)^2}{R^2}\right)R\,dR,\\
 E_\varphi
 &=\pi\int_0^\infty g'(R)^2R\,dR,
 \end{aligned}
\tag{4.6}
\]

and hence

\[
 E_z-E_\varphi
 =
 \pi\int_0^\infty\frac{g(R)^2}{R}\,dR
 =
 \pi\int_0^\infty\chi(s)^2\,ds.
\tag{4.7}
\]

After the outer azimuthal integration, the axial component of
\({\mathsf G}/\varepsilon\) is \(2\pi E_z\).  Either Cartesian transverse
component is \(\pi(E_\rho+E_\varphi)\).  Their difference is therefore
the right side of (4.2).  The \(C^1\) convergence above controls the
remainder by \(o(\varepsilon)\). \(\square\)

For a scaled, translated, and rotated copy

\[
 W_{A,\ell,c,Q}(x)
 =
 AQW_\varepsilon\left(Q^T\frac{x-c}{\ell}\right),
\tag{4.8}
\]

one has exactly

\[
 {\mathsf G}(W_{A,\ell,c,Q})
 =
 A^2\ell\,Q{\mathsf G}(W_\varepsilon)Q^T.
\tag{4.9}
\]

Thus no similarity parameter changes whether \({\mathsf G}\) is scalar.
Theorem 3.1 therefore applies to every member of a disjoint Gavrilov
lattice built from this fixed-profile thin seed, including its
packed-bubble versions.

## 5. The unavoidable pressure/velocity tail

The obstruction has a concrete far-field form.  For one bubble, let

\[
 S=\nu(W\otimes\Delta W+\Delta W\otimes W).
\]

The first pressure derivative solves

\[
 -\Delta p_1=\partial_i\partial_jS_{ij}.
\tag{5.1}
\]

Its integrated tensor is

\[
 M_{ij}:=\int S_{ij}\,dx=-2\nu{\mathsf G}(W)_{ij}.
\tag{5.2}
\]

With \(\Gamma(x)=(4\pi|x|)^{-1}\), the Newtonian solution has expansion

\[
 \begin{aligned}
 p_1(x)
 &=
 \partial_i\partial_j\Gamma(x)M_{ij}
 +O(|x|^{-4})\\
 &=
 -\frac{\nu}{2\pi|x|^3}
 \big(3\widehat x\otimes\widehat x-I\big):
 {\mathsf G}(W)
 +O(|x|^{-4}).
 \end{aligned}
\tag{5.3}
\]

For nonisotropic \({\mathsf G}(W)\), the leading coefficient is not
identically zero.  Hence

\[
 \nabla p_1(x)=O(|x|^{-4})
\]

with a nonzero leading angular profile.  Outside the original support,
equation (3.4) reduces to

\[
 u_2=-\nabla p_1.
\tag{5.4}
\]

The second formal velocity jet is therefore smooth and finite-energy but
not compactly supported.  This is the first pressure/centre wake that an
all-order transition must retain.

Cancellation of the leading tensor after summing several remote bubbles
only removes the first term at spatial infinity.  It does not make the
pressure vanish in the gaps.  If the correctors are confined to disjoint
balls, Theorem 3.1 forces every component's entire exterior pressure field
to vanish and recovers the componentwise condition (0.1).

## 6. Consequences for Borel and endpoint constructions

The exact conclusions are:

1. **A bare bubble has only one local viscous jet.**  
   The first jet \(\nu\Delta W\) is compact.  For a thin Gavrilov bubble,
   the second jet already has the nonlocal tail (5.4).

2. **Finite-dimensional modulation is irrelevant to this obstruction.**  
   Amplitude, scale, and translation derivatives can be absorbed into the
   first corrector derivative, but they do not change the total recursion
   (2.1) or the tensor (4.9).

3. **Ordinary Borel summation cannot restore disjoint support.**  
   It can realize the global sequence (2.1), and the resulting residual is
   flat.  It cannot replace the forced nonlocal \(u_2\) by a compact jet.
   The incompatibility precedes all factorial/Gevrey estimates.

4. **Disjoint symmetry copies do not cure a local endpoint.**  
   Reflections or cubic rotations can cancel a total quadrupole, but
   componentwise local correction still fails.  Such cancellation becomes
   relevant only if a global pressure/corrector field is allowed to join
   the copies.

5. **An admissible endpoint must be enlarged.**  
   At least one of the following is necessary:
   * retain the pressure/velocity tail as wake or centre data;
   * keep a nonzero oscillatory viscous bath at the endpoint, so the
     endpoint is not a bare Gavrilov bubble;
   * allow a global corrector or stress bridge between bubbles; or
   * leave a nonflat endpoint force.

The first three options remain possible.  The theorem says that the
"exactly disjoint bare fixed-profile thin bubbles plus locally supported
all-order correctors" ansatz is not one of them.

For Gavrilov localizations outside the fixed-profile thin family (4.1),
the exact test is still (0.1), but this note does not assert its sign.

## Primary source

* A. V. Gavrilov,
  [*A steady Euler flow with compact support*](https://arxiv.org/abs/1810.08020),
  especially the Taylor expansion in Remark 4 and formulas (6) and
  \(H(\alpha)\) used in Section 4 above.
