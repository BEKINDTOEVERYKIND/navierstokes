# Flat-force inverse design by Borel or meromorphic ansatz: exact gates and two remaining theorem targets

Date: 2026-07-29

## Decision

This note tests whether Clay alternative (D) can be reached without first
constructing a recurrent geometric cascade: prescribe a periodic
divergence-free velocity \(u\), singular at \(T<\infty\), and arrange that

\[
 {\cal N}(u):=\partial_tu+
 {\mathbb P}\operatorname{div}(u\otimes u)-\nu\Delta u
\tag{0.1}
\]

is \(C^\infty\)-flat at \(T\).  Here \(\mathbb P\) is the periodic Leray
projection.

There is no such shortcut in any of the usual finite-dimensional exact
families.  Energy, dissipation, and a fixed-Fourier-mode estimate rule out:

* affine or polynomial fields after imposing the Clay spatial conditions;
* a finite Fourier or finite fixed-profile ansatz;
* shears, potential flows, Beltrami eigenfields, and the known
  nonlinear-dark generalized Beltrami plane-wave spaces;
* a finite temporal Laurent ansatz with fixed spatial coefficients; and
* hiding an affine singularity with either a fixed cutoff or a
  beyond-all-orders shrinking cutoff.

Ordinary Borel summation does not alter this conclusion.  A flat residual
has zero formal power series, so every algebraic order must already solve
the *unforced* hierarchy.  Borel summation can realize a compatible
hierarchy or encode exponentially small terms; it cannot repair a failed
leading balance.

Two precise high-risk targets remain.

1. **All-order Euler--Reynolds parametrix.**  Interpolate between compact
   steady Euler bubbles at the low-frequency level, realize the resulting
   Reynolds stress by a single-carrier Mikado/Beltrami WKB expansion, and
   make adjacent stage jets match super-algebraically.  An exact gluing
   identity below then produces a flat force.  The scalar frequency ledger
   for this route actually closes with a polynomial carrier.  What is not
   known is the required all-order resonant cell inverse and nondegenerate
   endpoint construction.
2. **Meromorphic pole collision.**  Find a genuinely three-dimensional
   rational/meromorphic invariant manifold on the complexified torus whose
   poles reach the real torus at \(T\), while real energy stays bounded,
   total dissipation stays finite, and (0.1) is flat in every real
   \(C^m\) norm.  A local simple-pole calculation shows that viscosity
   forces a generic leading pole hypersurface to be characteristic for the
   *complexified* Laplacian.  This is not a construction, but it is a
   concrete symbolic search problem distinct from numerical optimization.

No singular Navier--Stokes solution is claimed.

## 1. Exact inverse-design criterion

The following elementary proposition explains why a flat residual would be
enough.

> **Proposition 1 (flat-residual criterion).**  Let
> \(u\in C^\infty(\mathbb T^3\times[0,T))\) be divergence-free and let
> \(f={\cal N}(u)\).  Suppose that, for every \(m,q,N\geq0\),
> \[
>  \|\partial_t^q f(t)\|_{C^m}
>  \leq C_{m,q,N}(T-t)^N.                         \tag{1.1}
> \]
> If some \(C^m\) norm of \(u(t)\) is unbounded as \(t\uparrow T\), then
> extending \(f\) by zero for \(t\geq T\) gives admissible smooth periodic
> forcing for Clay alternative (D).

Indeed, (1.1) makes the zero extension \(C^\infty\).  It is compactly
supported in future time and therefore satisfies the decay condition in
the official Clay statement.  The gradient part removed by \(\mathbb P\)
defines the pressure.  If a global smooth solution with the same data and
force existed, classical uniqueness on every compact subinterval of
\([0,T)\) would make it equal to \(u\).  Smoothness on the compact slab
\(\mathbb T^3\times[0,T]\) would then bound every \(C^m\) norm, a
contradiction.

Thus inverse design is legitimate.  The difficulty is entirely in (1.1),
not in the last uniqueness step.

## 2. Three exact necessary conditions

Write

\[
 E(t)=\frac12\|u(t)\|_2^2,\qquad
 D(t)=\|\nabla u(t)\|_2^2.
\]

For every smooth inverse-designed path,

\[
 E'(t)+\nu D(t)=\langle f(t),u(t)\rangle.          \tag{2.1}
\]

Since a flat force is bounded on \([0,T]\),

\[
 \|u(t)\|_2
 \leq \|u(0)\|_2+\int_0^t\|f(\tau)\|_2\,d\tau,    \tag{2.2}
\]

and consequently

\[
 \sup_{t<T}E(t)<\infty,\qquad
 \int_0^T D(t)\,dt<\infty.                        \tag{2.3}
\]

Moreover, because \(f\) is flat,

\[
 E(t)-E(T-)
 =\nu\int_t^T D(\tau)\,d\tau+O_N((T-t)^N)
\tag{2.4}
\]

for every \(N\), after changing \(N\) by one.  Near the singular time the
force supplies no algebraic-order energy input.

There is also a useful modewise statement.  With the usual Fourier
normalization, for every fixed \(k\in\mathbb Z^3\),

\[
\begin{split}
 \partial_t\widehat u(k)
 +\nu|k|^2\widehat u(k)
 +i\mathbb P_k k_j\widehat{u_j u}(k)
 =\widehat f(k),                                  \tag{2.5}\\
 |\widehat{u_i u_j}(k)|
 \leq C\|u(t)\|_2^2.                              \tag{2.6}
\end{split}
\]

The second line is just Cauchy--Schwarz for the Fourier convolution.
Equations (2.2), (2.5), and (2.6) imply that every fixed Fourier
coefficient is Lipschitz up to \(T\).  For \(k=0\),

\[
 \partial_t\widehat u(0)=\widehat f(0),
\tag{2.7}
\]

so the mean differs from its terminal value by a flat function.

Therefore any flat-force singularity must be a genuine escape to
\(|k|\to\infty\).  No fixed mode, finite collection of modes, or divergent
mean can carry it.

## 3. A finite-dimensional no-go theorem

Let \(V\subset C^\infty_\sigma(\mathbb T^3)\) be any fixed
finite-dimensional space and suppose

\[
 u(t)=\sum_{r=1}^d a_r(t)\phi_r\in V.              \tag{3.1}
\]

The \(L^2\) bound (2.2) bounds all coordinates \(a_r\).  Projecting the
equation onto \(V\) gives a finite-dimensional ODE

\[
 G a'=L a+Q(a,a)+F_V(t),                           \tag{3.2}
\]

where \(G\) is the positive Gram matrix, \(L\) is linear, \(Q\) is
quadratic, and \(F_V\) is smooth through \(T\).  A bounded solution of this
smooth finite-dimensional ODE extends through \(T\).  Repeated
differentiation then gives a smooth extension.

This argument does not require \(V\) to be invariant under the full
nonlinearity; projection of an alleged solution is enough.  It rules out,
in one stroke:

* every fixed finite Fourier template;
* every finite sum of fixed spatial profiles with Laurent, rational, or
  Borel-summed time coefficients; and
* every finite-dimensional exact-solution manifold whose parametrization
  remains nondegenerate in \(L^2\).

The last qualification matters.  A rational spatial function whose pole
approaches the real torus leaves every compact subset of its parameter
domain while possibly remaining \(L^2\)-bounded.  That is an
infinite-Fourier-tail mechanism and is considered in Section 7.

## 4. Audit of the standard exact ansatzes

### 4.1 Affine and polynomial fields

On \(\mathbb R^3\), put \(u(x,t)=M(t)x\), with
\(\operatorname{tr}M=0\).  Since \(\Delta u=0\),

\[
 \partial_tu+(u\cdot\nabla)u=(M'+M^2)x.            \tag{4.1}
\]

Whenever the non-gradient part of \(M'+M^2\) vanishes, a quadratic
pressure makes this an exact Euler and Navier--Stokes solution.  Miller
constructed explicit matrix solutions with

\[
 |M(t)|\asymp (T-t)^{-1}                           \tag{4.2}
\]

and finite-time blow-up.  These examples are important local strain
models, but they grow linearly at spatial infinity and have infinite
energy.  They are explicitly outside the Millennium data class.

They also cannot be put on the torus: periodicity of \(Mx\) gives
\(Mn=0\) for every \(n\in\mathbb Z^3\), hence \(M=0\).  More generally, a
periodic polynomial is bounded and therefore constant.  A nonconstant
polynomial on \(\mathbb R^3\) is neither rapidly decaying nor finite
energy.

A naive cutoff does not cure the problem.  If a strain of size
\(A\asymp (T-t)^{-1}\) is retained in a ball of radius \(L\), then

\[
\begin{array}{c|c}
\text{quantity}&\text{characteristic size}\\ \hline
\text{energy}&A^2L^5\\
\text{cutoff convection}&A^2L\\
\text{cutoff viscosity}&\nu A/L.
\end{array}                                       \tag{4.3}
\]

The last two entries have product \(\nu A^3\), independent of \(L\).
Thus no choice of one cutoff scale makes both errors flat.  Balancing them
gives the parabolic boundary scale and an error of order
\(A^{3/2}\), rather than a flat error.  This is a scaling statement for a
generic cutoff, not a universal lower bound after all possible nonlinear
correctors.  Such correctors would have to solve a localized
backward-profile or wake problem, which is precisely the hard missing
construction.

### 4.2 Potential and Beltrami pressure loopholes

The identity

\[
 (u\cdot\nabla)u
 =\nabla\frac{|u|^2}{2}-u\times\operatorname{curl}u
\tag{4.4}
\]

makes potential and Beltrami fields attractive because pressure removes
their nonlinearity.

On \(\mathbb T^3\), a curl-free and divergence-free field has only its
zero Fourier mode and is spatially constant.  Its evolution is
\(u'(t)=f_0(t)\), so flat forcing cannot make it singular.  On
\(\mathbb R^3\), an \(L^2\) curl-free, divergence-free field is zero by
the same Fourier argument.

If \(W\) is a periodic Beltrami eigenfield,

\[
 \operatorname{curl}W=\lambda W,\qquad
 \Delta W=-\lambda^2W,
\]

then \(u=a(t)W\) has projected residual

\[
 f=(a'+\nu\lambda^2a)W.                            \tag{4.5}
\]

Variation of constants shows that a smooth, flat right-hand side makes
\(a\) extend smoothly through \(T\).  Nadirashvili's Liouville theorem
also says that a finite-energy Beltrami field on \(\mathbb R^3\) is zero.

Known superposed plane-wave generalized Beltrami families similarly make
the nonlinearity vanish or become a gradient; their forced
Navier--Stokes dynamics reduce to linear variation of constants.  They
are valuable exact test fields but cannot blow up under smooth flat
forcing.

### 4.3 Shears and lower-dimensional reductions

For

\[
 u=(U(x_2,x_3,t),0,0),                             \tag{4.6}
\]

one has \((u\cdot\nabla)u=0\), and the equation is the forced heat
equation

\[
 \partial_tU-\nu\Delta_{2,3}U=f_1.                 \tag{4.7}
\]

Smooth data and smooth forcing give a smooth solution through every finite
time.  Common-direction multi-shears remain nonlinearly dark.  The
standard 2D and 2D3C reductions are also globally regular and cannot
supply three-dimensional vortex stretching.

## 5. Why ordinary Borel summation is not the missing cancellation

There are three different statements that are sometimes conflated.

1. The ordinary \(C^\infty\) Borel lemma realizes an arbitrary collection
   of finite jets at \(T\).  Its output is smooth at \(T\), so it cannot by
   itself realize a singular velocity.
2. Borel summation of the Navier--Stokes time-Taylor series reconstructs
   the unique local analytic solution.  Costin and Tanveer obtain
   \[
   u(t)=u_0+\int_0^\infty e^{-p/t}U(p)\,dp,         \tag{5.1}
   \]
   and translate longer continuation into a growth problem for \(U(p)\)
   as \(p\to\infty\).  This is a useful equivalent formulation of the
   evolution, not freedom to choose a different terminal branch.
3. A transseries may contain beyond-all-orders terms such as
   \(e^{-c/(T-t)}\).  These can produce a nonzero smooth-flat residual,
   but they have zero formal power series.

The third observation gives a hard formal gate.  If

\[
 {\cal N}(u)=O((T-t)^N)\quad\text{for every }N,     \tag{5.2}
\]

then substitution of any valid power, Puiseux, or logarithmic asymptotic
expansion of \(u\) makes every coefficient in the *unforced* formal
Navier--Stokes hierarchy vanish.  An exponentially small Stokes term
cannot cancel a nonzero algebraic leading residual.  A Borel procedure
can sum an all-order hierarchy only after that hierarchy has been solved.

There is an even simpler analytic version: a force analytic at \(T\) and
flat there is identically zero in a neighborhood of \(T\).  A nonzero
flat-force construction must therefore use terminal nonanalyticity, but
its algebraic asymptotics are still unforced.

## 6. Superflat concentration fails for a single high-Re packet

A tempting transseries trick is to make the support shrink faster than
every power, hoping to hide a divergent amplitude.  Consider the concrete
isotropic packet

\[
 u_s(x)=a(s)U\!\left(\frac{x-x_*}{\ell(s)}\right),
 \qquad s=T-t,
\tag{6.1}
\]

where \(U\) is a fixed nonconstant compactly supported divergence-free
profile.  Then

\[
 \|u_s\|_2^2=a^2\ell^3\|U\|_2^2,\qquad
 \|\nabla u_s\|_2^2=a^2\ell\|\nabla U\|_2^2.
\tag{6.2}
\]

Let

\[
 a=e^{A\Phi(s)},\qquad
 \ell=e^{-B\Phi(s)},\qquad
 \frac{\Phi(s)}{\log(1/s)}\longrightarrow\infty,
\tag{6.3}
\]

with \(A,B>0\).  A high-Reynolds-number packet requires

\[
 a\ell/\nu\longrightarrow\infty,
\quad\text{hence}\quad A>B.                        \tag{6.4}
\]

But then

\[
 \|\nabla u_s\|_2^2
 \asymp e^{(2A-B)\Phi(s)}.                         \tag{6.5}
\]

Since \(2A-B>0\) and \(\Phi(s)/\log(1/s)\to\infty\), the right-hand side
eventually exceeds \(s^{-2}\).  It is not integrable at \(s=0\), contrary
to (2.3).

Thus a single packet cannot combine beyond-all-orders concentration with
an Euler-dominated high-Reynolds blow-up.  Choosing \(A\leq B\) moves the
packet into the viscosity-dominated regime, where the dominant heat term
must itself be cancelled to all orders.  If each spatial jet of the heat
term is merely made flat term by term, then every jet of \(u\) is flat as
well; a singular construction again needs a genuine PDE cancellation.

This does not exclude anisotropic, multistrand, or wake-carrying
concentration.  It does rule out the simplest exponential Borel cutoff
around an affine or fixed-profile core.

## 7. Laurent and rational ansatzes

There are two different Laurent ideas.

* A finite Laurent series in time with fixed smooth spatial coefficients
  lies in a fixed finite-dimensional space and is ruled out by Section 3.
  More directly, a nonzero negative-power leading coefficient makes the
  \(L^2\) energy diverge.
* A finite Laurent polynomial in
  \(z_j=e^{ix_j}\) is a finite Fourier field and is ruled out for the same
  reason.

A rational function of \(z=(e^{ix_1},e^{ix_2},e^{ix_3})\) is different:
while its poles remain off the unit torus it is smooth and has an infinite
Fourier tail.  A pole approaching the unit torus can collapse the
analyticity radius and move energy to \(|k|\to\infty\), exactly the escape
allowed by Section 2.

There is nevertheless a strong local viscous constraint on the most
plausible finite-pole case.  In the complexified spatial variables suppose,
at a regular point of a *simple* pole hypersurface \(q=0\),

\[
 u\sim a\,q^{-1},\qquad n=\nabla q.
\tag{7.1}
\]

The leading divergence equation gives

\[
 a\cdot n=0.                                      \tag{7.2}
\]

Freezing the leading residue and linearizing \(q\) at the pole,

\[
 \Delta(aq^{-1})
 =2(n\cdot n)a\,q^{-3}+\text{lower poles}.
\tag{7.3}
\]

The nominal most singular self-advection term is proportional to
\((a\cdot n)a\) and vanishes by (7.2).  A pressure pole has leading
gradient parallel to \(n\); it cannot cancel a nonzero tangential vector
\(a\).  Therefore a nonzero tangential leading residue requires

\[
 \boxed{\,n\cdot n=0\,}                            \tag{7.4}
\]

unless another pole component creates a same-order tangential
cancellation.  Higher-order poles have a more complicated hierarchy:
derivatives of their residues and cross interactions can precede the
viscous order, so (7.4) is not asserted for them without further
hypotheses.
Equation (7.4) is the characteristic cone of the *complexified* Laplacian.
There are no real characteristics, but complex null normals do exist.

For a periodic rational ansatz this condition has an exact algebraic form.
Write \(z_j=e^{ix_j}\) and
\(D_j=z_j\partial_{z_j}\).  If an irreducible Laurent polynomial \(q(z)\)
defines a generically smooth pole hypersurface, then (7.4) requires

\[
 q\ \big|\ \sum_{j=1}^3(D_jq)^2                  \tag{7.5}
\]

in the Laurent polynomial ring.  This follows because
\(\partial_{x_j}q=iD_jq\) and a polynomial vanishing on the irreducible
hypersurface \(q=0\) is divisible by \(q\).

Condition (7.5) already excludes every genuine binomial denominator.  If

\[
 q=c_1z^k+c_2z^\ell,\qquad k\ne\ell,
\]

then on \(q=0\), writing \(c_1z^k=A\) and \(c_2z^\ell=-A\), one has

\[
 \sum_j(D_jq)^2=|k-\ell|^2A^2\ne0.               \tag{7.6}
\]

Irreducible trinomials are excluded as well.  After removing one monomial,
write

\[
 q=c_0+c_1z^a+c_2z^b.                            \tag{7.7}
\]

If \(a,b\) are rationally dependent, \(q\) is a one-variable Laurent
polynomial in a single character and factors over \(\mathbb C\); its
irreducible factors are binomials.  If \(a,b\) are independent, the
character map \(z\mapsto(z^a,z^b)\) is surjective onto
\((\mathbb C^*)^2\).  On \(q=0\), set \(A=c_1z^a\) and
\(c_2z^b=-c_0-A\).  Then

\[
 \sum_j(D_jq)^2
 =\big((a-b)A-bc_0\big)\cdot
   \big((a-b)A-bc_0\big),                         \tag{7.8}
\]

where the dot is the complex bilinear extension of the Euclidean
integer-lattice form.  As \(A\) varies, its leading coefficient is
\(|a-b|^2>0\), so it cannot vanish identically.  At this stage of the
audit, the smallest unresolved irreducible characteristic denominator had
at least four Laurent monomials.  The subsequent note
[`2026-07-29-laurent-null-pole-no-go.md`](2026-07-29-laurent-null-pole-no-go.md)
closes the remaining finite-support search: no nonunit irreducible Laurent
polynomial with any number of monomials can satisfy (7.5).

This calculation does not prove that a pole collision is possible.
Rather, it turns the rational route into a sharp algebraic search:

1. find a nontrivial Laurent denominator satisfying (7.5), with conjugate
   pole surfaces disjoint from the real torus for \(t<T\);
2. satisfy divergence and the full next-order pole equations with genuine
   three-dimensional cross interactions;
3. make the real energy uniformly bounded and the real dissipation
   integrable;
4. force all uncancelled real-axis residual jets to be flat; and
5. avoid the shear/generalized-Beltrami reductions, which revert to a
   linear heat equation.

After that all-support Laurent no-go, the uncovered pole variants are
coupled divisors with same-order tangential cancellations, higher-order
pole hierarchies, or non-Laurent essential singularities.  None currently
has a closed residue system or an energy-dissipation ledger.

## 8. Exact all-order seam identity

Borel/Whitney gluing can be stated without hand-waving.  Let

\[
 Q(v)=\mathbb P\operatorname{div}(v\otimes v)
\]

and let \(a,b\) be two divergence-free paths on an overlap.  Put
\(d=b-a\) and

\[
 u=a+\chi(t)d,
\]

where \(0\leq\chi\leq1\).  Direct expansion gives the exact identity

\[
\boxed{
 {\cal N}(u)
 =(1-\chi){\cal N}(a)+\chi{\cal N}(b)
 +\chi'd-\chi(1-\chi)Q(d).
}
\tag{8.1}
\]

Suppose seams \(t_j\uparrow T\) have widths \(\delta_j\) and
\(\chi_j^{(q)}=O(\delta_j^{-q})\).  If, on every overlap, the raw residuals
are flat in \(s_j=T-t_j\), and for every \(m,q,L,N\)

\[
 \|\partial_t^q d_j\|_{C^{m+1}}
 \leq C_{m,q,L,N}\,s_j^N\delta_j^L,               \tag{8.2}
\]

then (8.1), Leibniz's rule, and choosing \(L\) larger than the number of
cutoff derivatives show that the glued residual is flat in every
space-time \(C^m\) seminorm.

This is the rigorous content of the proposed Borel seam method.  It also
identifies exactly what Borel's theorem does **not** provide: (8.2).
Ordinary smoothing can preserve super-algebraic endpoint compatibility,
but it cannot turn an order-one endpoint mismatch into compatibility.

## 9. Surviving theorem targets and priority

The broad inverse-design question can now be replaced by either of two
falsifiable statements.

### Target A: flat stage-shadowing theorem

Construct a sequence of bounded-energy, finite-dissipation,
divergence-free stage paths with:

1. a high norm diverging as \(t_j\uparrow T\);
2. raw projected residuals flat on the tail;
3. the all-order overlap estimate (8.2); and
4. no reduction to a fixed finite-dimensional or nonlinear-dark class.

Then Proposition 1 and (8.1) prove Clay alternative (D).  For the current
Kelvin--Reynolds program, (8.2) is the exact analytic form of the
increasing-order endpoint return requirement.

### Target B: coupled or higher-order complex poles

The single reduced finite-Laurent simple-pole target is now closed by the
all-support no-go cited after (7.8).  A separate pole attack would first
have to derive a new leading system in which coupled divisors cancel the
tangential viscous pole, or in which a higher-order residue hierarchy
precedes the null-normal condition.  Without such a system there is no
meaningful symbolic or GPU search.

## 10. Highest-priority bypass: Gavrilov bubbles plus one-carrier Euler--Reynolds WKB

There is a more concrete version of Target A which could, if its missing
parametrix theorem were proved, bypass an exact Kelvin state-to-state
orbit.

### 10.1 What works at principal order

Gavrilov constructs a nonzero \(C^\infty\) compactly supported steady
Euler field \(U\).  The construction first produces an analytic local
field and then multiplies it by an arbitrary smooth function of its
pressure.  Choosing that final function to be a compact
Gevrey-\(\sigma\) bump, with any \(\sigma>1\), gives the natural candidate
for a compact Gevrey bubble.  Analytic compact support is impossible, so
\(\sigma>1\) is essential.

Scaled and translated copies

\[
 U_{a,\ell,x_0}(x)=a\,U((x-x_0)/\ell)             \tag{10.1}
\]

are again steady Euler fields.  A smooth path \(v(t)\) between a parent
bubble and a smaller child plus disjoint wake bubbles can be made constant
in total momentum.  Its Euler defect

\[
 g=\partial_tv+\mathbb P\operatorname{div}(v\otimes v) \tag{10.2}
\]

can formally be written as

\[
 g=-\mathbb P\operatorname{div}R_0.               \tag{10.3}
\]

For a compactly supported symmetric \(R_0\) on \(\mathbb R^3\), this
requires the exact moment conditions

\[
 \int g\,dx=0,\qquad \int x\times g\,dx=0,         \tag{10.3a}
\]

obtained by integration by parts.  Thus a localized bubble interpolation
must conserve total momentum and angular momentum, or export the defect
into a nonlocal wake.  Symmetric bubble pairs can plausibly enforce these
finite conditions; an isotropic gauge cannot change them.

Adding \(\rho(t,x)I\) changes neither (10.3) nor the work against a
divergence-free field.  Hence, subject to constructing a suitably
localized symmetric anti-divergence, one can make
\(R=R_0+\rho I\) positive in the interpolation region.  The Mikado
geometric lemma then supplies leading oscillations \(W\) with

\[
 \langle W\rangle_y=0,\qquad
 \langle W\otimes W\rangle_y=R.                  \tag{10.4}
\]

Thus neither positivity nor the principal averaged equation kills the
proposal.  This is exactly the setting in which Daneri--Székelyhidi
introduced Mikado flows to recover arbitrary Reynolds stresses.

### 10.2 Three necessary corrections to the naive proposal

First, **constant total energy is incompatible with viscosity and flat
forcing**.  If the proposed smooth approximate velocity has
\(E'(t)=0\), then (2.1) gives

\[
 \langle f,u\rangle=\nu\|\nabla u\|_2^2.          \tag{10.5}
\]

The left side is flat because \(f\) is flat and \(\|u\|_2\) is bounded.
The right side is not flat for a nontrivial high-frequency carrier.
Therefore the isotropic gauge must prescribe a *viscously decreasing*
energy profile,

\[
 E'(t)=-\nu\|\nabla u\|_2^2+O((T-t)^\infty),      \tag{10.6}
\]

not a constant one.  In the cascade window the stage loss is summable, so
this correction is compatible with the scalar energy ledger.

Second, a pure steady Euler bubble cannot occupy an endpoint collar.  For
such a bubble,

\[
 \mathbb P\operatorname{div}(U\otimes U)=0,
\qquad
 {\cal N}(U)=-\nu\Delta U.                        \tag{10.7}
\]

For (10.1), its \(C^m\) size is of order

\[
 \nu a\ell^{-m-2},                                \tag{10.8}
\]

which diverges rather than becoming flat along the proposed shrinking
bubbles.  At a seam one must therefore do at least one of the following:

* match the full forward heat/Navier--Stokes jet beginning with
  \(\partial_tU=\nu\Delta U\);
* retain a nonzero high-frequency stress satisfying
  \(\mathbb P\operatorname{div}R=\nu\Delta U\); or
* put (10.8) into the external force, which is forbidden by flatness.

The second option means that the endpoint is not a bare Gavrilov bubble:
it includes an active oscillatory bath which must pass into the wake.  The
first option loses exact stationarity and compact support immediately,
although the deformation is perturbative on a high-Reynolds stage.

Third, viscosity has an exact cell solvability condition.  Let \(W\) be a
stationary fast Euler field and let

\[
 L_W z=
 \mathbb P_y\operatorname{div}_y(W\otimes z+z\otimes W)
\tag{10.9}
\]

be its stationary linearization.  Differentiating the Euler energy
identity at \(W\) gives

\[
 \langle L_Wz,W\rangle_y=0                       \tag{10.10}
\]

for every \(z\).  But

\[
 \langle-\Delta_yW,W\rangle_y
 =\|\nabla_yW\|_2^2>0.                            \tag{10.11}
\]

Consequently the static first viscous equation

\[
 L_Wz=\Delta_yW                                  \tag{10.12}
\]

has no solution.  Viscosity cannot be removed by a stationary
single-carrier corrector while keeping its energy fixed.  A viable WKB
hierarchy must include slow amplitude decay or energy transfer from the
low field, and (10.10) becomes precisely the order-by-order energy
compatibility (10.6).

### 10.3 Why standard convex integration has fatal frequency inflation

A usual convex-integration step cancels the current Reynolds stress by a
new wave whose amplitude is the square root of that stress.  Its transport
and cross-interaction errors become the next stress, which is cancelled at
a fresh, substantially larger carrier.  The published Euler and
Navier--Stokes schemes use a hierarchy

\[
 \lambda_{q+1}\gg\lambda_q,                       \tag{10.13}
\]

not an all-order expansion at one carrier.

If even a fixed factor \(C>1\) is spent at each of \(M_j\) correction
levels, then

\[
 \lambda_{\mathrm{final}}\geq C^{M_j}K_j.         \tag{10.14}
\]

For the truncation order in Section 10.4,
\(M_j\asymp j^2/\log j\), this makes
\(\log\lambda_{\mathrm{final}}\gtrsim j^2/\log j\).
The normalized viscosity decays only exponentially in \(j\), so

\[
 \varepsilon_j\lambda_{\mathrm{final}}^2
 \longrightarrow\infty.                          \tag{10.15}
\]

Intermittency can make a Reynolds stress small in integral norms, as in
the Buckmaster--Vicol construction, but it does not hide
\(\nu\Delta u\) from the \(C^\infty\) supremum seminorms required of the
Clay force.

Therefore no existing multi-carrier convex-integration iteration gives the
needed flat residual.  The proposal lives or dies on a genuinely different
object: a nonlinear WKB/Newton parametrix whose \(M\)-th truncation uses
only harmonics \(O(MK)\), rather than \(M\) nested carrier scales.

### 10.4 Polynomial carriers pass the complete scalar ledger

The earlier exponential-carrier choice is unnecessary.  Let the physical
stage scale be

\[
 \ell_j=r^{-j},\qquad
 a_j=\ell_j^{-\gamma},\qquad 1<\gamma<3/2,
\]

so the normalized viscosity is

\[
 \varepsilon_j=\frac{\nu}{a_j\ell_j}
 \asymp e^{-c_\varepsilon j}.
\tag{10.16}
\]

Choose an internal carrier and truncation order

\[
 K_j=j^A,\qquad
 M_j=\left\lfloor\frac{\eta j^2}{\log j}\right\rfloor.
\tag{10.17}
\]

Here \(K_j\) has two conceptually different possible roles.

* In the generic Euler--Reynolds construction it is only the approximation
  carrier.  The low bubble amplitude remains
  \(a_j=\ell_j^{-\gamma}\), and the principal wave amplitude is fixed by
  \(\sqrt{R_j}\), not by a Kelvin law \(K_j^\gamma\).
* If the fast field is required to be created by the exact Kelvin
  amplifier, its amplitude-frequency law adds the slowly varying factor
  \(K_j^\gamma=j^{A\gamma}\).  Thus the corresponding physical amplitude
  is schematically
  \(\ell_j^{-\gamma}K_j^\gamma\).  This polynomial factor is
  subexponential in \(j\): it does not change the exponential scale window,
  active-energy summability, or terminal flatness, and it makes normalized
  viscosity smaller.

The second interpretation is not obtained for free from an
Euler--Reynolds interpolation.  Creating the next carrier with the correct
phase, polarization, and amplitude, and making it recur at the next
endpoint, remains part of the missing cell inverse in Section 10.5.

Assume the order-\(M\) WKB coefficients obey a Gevrey bound

\[
 \|U_M\|\leq C^{M+1}(M!)^\sigma                  \tag{10.18}
\]

and that harmonic support grows at most linearly with \(M\).  The pure
geometric-optics remainder is then bounded schematically by

\[
 C^{M_j}(M_j!)^\sigma K_j^{-M_j}.
\]

Stirling's formula gives

\[
\begin{split}
 \log\!\left(C^{M_j}(M_j!)^\sigma K_j^{-M_j}\right)
 &=
 \sigma M_j\log M_j-A M_j\log j+O(M_j)\\
 &=
 -(A-2\sigma)\eta j^2+o(j^2).                    \tag{10.19}
\end{split}
\]

Thus any

\[
 \boxed{A>2\sigma}                                \tag{10.20}
\]

gives an \(e^{-cj^2}\) remainder.  Fixed physical space-time derivatives
cost only \(e^{O(j)}\) from stage rescaling and polynomial powers of
\(K_jM_j\), so (10.19) still implies flatness.  Indeed the remaining
physical time is comparable to
\(\ell_j^{1+\gamma}=e^{-c_Tj}\); hence \(e^{-cj^2}\) is smaller than every
power of \(T-t\).

The carrier heat parameter is

\[
 \theta_j=\varepsilon_jK_j^2
 =e^{-c_\varepsilon j}j^{2A}\longrightarrow0.    \tag{10.21}
\]

If the largest harmonic is \(O(M_jK_j)\), then even
\(\varepsilon_j(M_jK_j)^2\to0\).  A viscosity expansion truncated at
\(M_j\) has an error far smaller than (10.19).  Hence viscosity causes no
*scalar* frequency obstruction for a true one-carrier bivariate
\((K^{-1},\theta)\) expansion.

Compact cutoffs are also compatible: standard compact
Gevrey-\(\sigma\) functions exist for every \(\sigma>1\), and their
factorial derivative cost is already included in (10.18).  Gavrilov's
analytic local field followed by a pressure cutoff is structurally suited
to this choice, although a uniform Gevrey estimate for the complete
scaled/interpolated construction still has to be written down.

Finally,

\[
 \frac{K_{j+1}}{K_j}
 =1+\frac{A}{j}+O(j^{-2}).                        \tag{10.22}
\]

The internal carrier handoff is therefore a slow modulation.  The physical
carrier \(K_j/\ell_j\) still changes by the bounded factor
\(r(1+A/j+O(j^{-2}))\).  Polynomial or other subexponential corrections to
the bubble amplitude do not alter the exponent window; they must be
included in the endpoint state rather than discarded into the force.

### 10.5 The exact missing theorem

The scalar estimates leave one sharply isolated question.

> **One-carrier Gevrey parametrix problem.**  Given a compact
> Gevrey-\(\sigma\) Euler--Reynolds path with positive stress and a
> viscously admissible energy profile, construct a bivariate WKB expansion
> through arbitrary order \(M\) such that:
>
> 1. its fast harmonics are bounded by \(CMK\);
> 2. its coefficients satisfy (10.18) uniformly under bubble scaling;
> 3. every mean, energy, momentum, and helicity resonance is absorbed by
>    the low path or the outgoing wake;
> 4. its endpoint jets either retain a nondegenerate wake carrier or solve
>    the full heat/Navier--Stokes jet, without division by a vanishing
>    stress amplitude; and
> 5. the resulting stage endpoints satisfy the seam estimate (8.2),
>    including creation and recurrence of the next carrier rather than
>    merely its averaged Reynolds stress.

This theorem is not supplied by the Mikado geometric lemma.  For a single
straight Mikado pipe, the fast linearized equation contains transport
along closed streamlines and has an explicit infinite-dimensional
cokernel.  In local periodic coordinates take

\[
 W=e_1\phi(y_2,y_3).
\]

The adjoint of (10.9) on divergence-free fields is

\[
 L_W^*\psi
 =\mathbb P_y\big(-(W\cdot\nabla)\psi+
                  (\nabla W)^T\psi\big).         \tag{10.23}
\]

For every smooth scalar \(F\), with \(F(0)=0\) if transverse localization
is desired,

\[
 \psi_F=e_1F(\phi)
\]

is divergence-free and satisfies \(L_W^*\psi_F=0\): the transport term
vanishes, while
\((\nabla W)^T\psi_F=\nabla G(\phi)\) for \(G'=F\), and is removed by
\(\mathbb P_y\).  Varying \(F\) gives infinitely many independent
compatibility conditions.  A finite collection of spatially disjoint
pipes inherits these local cokernels.

Beltrami cells remove the leading quadratic interaction, but their
linearized stationary Euler operator is still non-elliptic and already has
the energy cokernel (10.10).  It is presently unproved that the special
forcing generated by the WKB hierarchy satisfies all of these
compatibilities with a Gevrey-tame inverse.

Accordingly:

* the **naive** recipe (constant energy, bare steady endpoints, and an
  ordinary finite convex-integration stack) is closed;
* the **modified one-carrier Gevrey recipe is not closed by scaling** and
  is a serious prize-level theorem target;
* proving or disproving its resonant cell inverse should take priority over
  further GPU tuning.

## Primary references

* C. Fefferman, [Existence and Smoothness of the Navier--Stokes
  Equation](https://www.claymath.org/wp-content/uploads/2022/06/navierstokes.pdf).
* O. Costin and S. Tanveer,
  [Borel summability of Navier--Stokes equation in \(\mathbb R^3\) and
  small time existence](https://arxiv.org/abs/math/0612063).
* E. Miller,
  [Finite-time blowup for smooth solutions of the Navier--Stokes
  equations on the whole space with linear growth at
  infinity](https://arxiv.org/abs/2103.12237).
* N. Nadirashvili,
  [Liouville theorem for Beltrami
  flow](https://arxiv.org/abs/1403.1414).
* A. Prugger and J. D. M. Rademacher,
  [Explicit superposed and forced plane wave generalized Beltrami
  flows](https://arxiv.org/abs/2003.07824).
* A. V. Gavrilov,
  [A steady Euler flow with compact
  support](https://arxiv.org/abs/1810.08020).
* S. Daneri and L. Székelyhidi Jr.,
  [Non-uniqueness and \(h\)-principle for Hölder-continuous weak
  solutions of the Euler equations](https://arxiv.org/abs/1603.09714).
* T. Buckmaster and V. Vicol,
  [Nonuniqueness of weak solutions to the Navier--Stokes
  equation](https://arxiv.org/abs/1709.10033).
