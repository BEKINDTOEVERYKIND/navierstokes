# Can a forward corrector replace the two-endpoint Euler inverse?

Date: 2026-07-29

## Verdict

A generic forward Navier--Stokes stability argument does **not** replace the
missing return theorem.  It solves an initial-value problem, but it does not
make the outgoing packet land on the input set for the next stage.  In the
present geometric-optics window it has two additional difficulties:

1. the active carrier lies deliberately below the heat cutoff, so viscosity
   supplies asymptotically no contraction over one stage;
2. the crude linearized estimate costs
   \(\exp(CK_j)\), while the previously proposed order
   \(M_j\asymp j\) only produces \(\exp(-c j^2)\).

The Cheskidov--Dai--Palasek forward corrector avoids such a loss by a
structure-specific semigroup theorem and by putting each active shell on its
heat time.  Neither feature is automatic for a forward, Euler-dominated
Kelvin amplifier.

There is, however, a sharper route that avoids the impossible problem of
fixing both complete endpoints.  One can seek a parameter-dependent
**discrete invariant graph** for the renormalized one-step map.  Its
coefficient equations are

\[
 \big(\lambda_{pq}I-L\big)Z_{pq}=G_{pq},
\]

where \(L\) is the derivative of the inviscid return map and
\(\lambda_{pq}\in(0,1)\) is the stage-to-stage scaling of a carrier/viscous
monomial.  This is a resolvent problem for initial data, not a right inverse
for a spacetime operator with two fixed endpoints.  If these resolvents are
Gevrey-tame, truncation at order \(M_j\asymp j\) leaves \(e^{-c j^2}\)
endpoint and equation errors.

For the forced Clay alternative, that equation error can simply be declared
to be the force.  It is terminally flat after physical rescaling, so no
parabolic corrector is required.  If one additionally wants every stage to
be an exact unforced Navier--Stokes segment, then a tame forward propagator
is needed to correct the interior residual before terminally-flat seam
forces join the stages.

This does not construct the leading localized Euler return cell.  It
replaces the overstrong all-order theorem that had been placed after that
cell by a more precise, potentially weaker one.

## 1. Stage scaling and the two small parameters

Use

\[
 \ell_j=r^{-j},\qquad
 a_j=\ell_j^{-\gamma},\qquad
 1<\gamma<\frac32,
\]

and normalized variables

\[
 y=\frac{x-x_j}{\ell_j},\qquad
 s=\frac{t-t_j}{\tau_j},\qquad
 \tau_j=\frac{\ell_j}{a_j}=\ell_j^{1+\gamma}.
\]

Writing \(u=a_jV\), the projected Navier--Stokes equation becomes

\[
 \partial_sV+\mathbb P\operatorname{div}(V\otimes V)
 -\varepsilon_j\Delta V=R,
 \qquad
 \varepsilon_j=\frac{\nu}{a_j\ell_j}
 =\nu\ell_j^{\gamma-1}.
\]

Choose an internal carrier

\[
 K_j=\ell_j^{-\kappa}.
\]

The first localization parameter and the heat action on that carrier are

\[
 h_j=K_j^{-1}=\ell_j^\kappa,\qquad
 \theta_j=\varepsilon_jK_j^2
 =\nu\ell_j^{\gamma-1-2\kappa}.
\]

Thus the nonempty carrier window is exactly

\[
 \boxed{0<\kappa<\frac{\gamma-1}{2}.}
\]

Put

\[
 \sigma_0=\min\{\kappa,\gamma-1-2\kappa\}.
\]

Then \(\rho_j:=\max(h_j,\theta_j)\lesssim\ell_j^{\sigma_0}\).  The
balanced choice

\[
 \kappa=\frac{\gamma-1}{3}
\]

gives \(\sigma_0=(\gamma-1)/3\).

The stage-to-stage parameter ratios are

\[
 q_h=\frac{h_{j+1}}{h_j}=r^{-\kappa},
 \qquad
 q_\theta=\frac{\theta_{j+1}}{\theta_j}
 =r^{-(\gamma-1-2\kappa)}.
\]

These two contractions, rather than a fixed-\(\varepsilon\) return profile,
are the natural base dynamics for an all-order construction.

## 2. What a one-sided forward corrector actually gives

Suppose \(V_j^{\rm app}\) has residual \(R_j\).  A forward correction with
zero incoming data solves

\[
\begin{split}
 \partial_sw-\varepsilon_j\Delta w
 &+\mathbb P\operatorname{div}
 \big(V_j^{\rm app}\otimes w+w\otimes V_j^{\rm app}
       +w\otimes w\big)\\
 &=-R_j,\qquad w(0)=0.
\end{split}
\]

If the linearized evolution family \(\mathcal U_j(s,\sigma)\) obeys a
uniform tame estimate

\[
 \sup_{0\le \sigma\le s\le S}
 \|\mathcal U_j(s,\sigma)\|_{X_K^m\to X_K^m}
 \le C_{j,m},
\]

then a standard contraction estimate gives

\[
 \|w(S)\|_{X_K^m}
 \lesssim C_{j,m}\|R_j\|_{L^1_sY_K^m}.
\]

This is an initial-value theorem.  The endpoint \(w(S)\) is generally not
zero.  Consequently, it proves neither an exact return nor the matching
condition for the next stage.

There is nevertheless one functional-analytic advantage.  A forward
symmetric-hyperbolic energy estimate need not lose derivatives of the
unknown.  It asks for one more derivative of the known coefficient, but
maps a source in \(H^m\) to a solution in \(H^m\).  By contrast, repeatedly
using a solver that loses \(d\) carrier derivatives at every asymptotic
order would change

\[
 h_j^n
\quad\hbox{into}\quad
 h_j^nK_j^{dn}=K_j^{(d-1)n}.
\]

For any integer \(d\ge1\), that is not an asymptotic expansion.  A viable
scheme must therefore use semiclassical norms with no repeated carrier
loss, a Nash--Moser mechanism that restores it, or a special triangular
identity.

A single fixed loss is harmless after all-order truncation:
\(K_j^D e^{-c j^2}=e^{-c j^2+O(j)}\).  A loss \(K_j^{Dn}\) at the \(n\)-th
coefficient is not.

## 3. Why ordinary parabolic stability is insufficient

The carrier condition \(\theta_j=\varepsilon_jK_j^2\to0\) says precisely
that heat damping over a normalized stage is

\[
 e^{-S\varepsilon_jK_j^2}=1-o(1).
\]

Thus viscosity does not contract carrier sidebands, envelope errors, or
low-frequency wake errors.  It only strongly damps frequencies well above
\(\varepsilon_j^{-1/2}\), whereas

\[
 K_j\ll\varepsilon_j^{-1/2}
\]

throughout the admissible window.

There is also a serious high-carrier stability cost.  A rank-one stress of
order one is produced by carrier velocities of order one.  In ordinary
Sobolev estimates this gives

\[
 \|\nabla V_j^{\rm app}\|_\infty\sim K_j
\]

during an order-one normalized interval, and hence only the crude bound

\[
 C_{j,m}\lesssim e^{CK_j}.
\]

With \(K_j=r^{\kappa j}\),

\[
 e^{CK_j}e^{-c j^2}\longrightarrow\infty.
\]

Therefore the \(M_j\asymp j\) flatness ledger cannot support a generic
forward energy estimate.  The required estimate must exploit the exact
Kelvin polarization, frequency separation, spatial intermittency, or a
triangular cancellation strongly enough to make

\[
 \log C_{j,m}=o(j^2).
\]

A polynomial loss \(C_{j,m}\lesssim K_j^{D_m}\), or more generally
\(\log C_{j,m}=O_m(j\log j)\), is acceptable.

This is not a universal no-go theorem: the \(e^{CK_j}\) estimate is an upper
bound, not proof that the proposed cell has such an unstable mode.  It is a
fatal obstruction to the *unstructured* forward-corrector argument.  A
structure-specific semigroup theorem is indispensable.

## 4. Why summably small endpoint errors are not automatically enough

Let \(e_j\) be the normalized error in the incoming state and let \(L\) be
the derivative of the inviscid renormalized return map.  Forward correction
alone leads schematically to

\[
 e_{j+1}=Le_j+b_j,\qquad |b_j|\lesssim e^{-c j^2}.
\]

If an unstable multiplier of \(L\) has modulus \(\Lambda>1\), then an error
created at one finite stage grows like

\[
 \Lambda^{j-k}b_k
\]

and eventually leaves every fixed tube, however small \(b_k\) was.  If
\(L\) is power bounded, the superexponential \(b_j\)'s are summable and the
error stays small.  If the outgoing error is spatially expelled into a wake
whose child projection is summable, the same conclusion can hold for the
projected child even when the full state is not stable.

Thus one needs one of the following, none supplied by ordinary parabolic
well-posedness:

1. a power-bounded or contracting renormalized child map;
2. a finite-codimensional unstable bundle plus a shooting/stable-manifold
   theorem;
3. a one-way packet/wake geometry that makes old endpoint errors
   asymptotically invisible to the child;
4. the invariant-graph construction below.

The heat equation cannot provide item 1 at the active carrier because
\(\theta_j\to0\).

## 5. A discrete invariant-graph alternative

Let

\[
 F_{h,\theta}:X\longrightarrow X
\]

denote one complete, phase-resolved, renormalized packet-plus-wake stage.
It is important that \(X\) retain the oscillatory phases and the wake; an
averaged Reynolds tensor alone loses information and need not define an
invertible inviscid map.

Assume the leading inviscid cell gives

\[
 F_{0,0}(P)=P,
\qquad
 L=D_zF_{0,0}(P).
\]

Seek a graph

\[
 Z(h,\theta)
 =
 P+\sum_{p+q\ge1}h^p\theta^qZ_{pq}
\]

that obeys the nonautonomous invariance equation

\[
 \boxed{
 F_{h,\theta}\big(Z(h,\theta)\big)
 =
 Z(q_hh,q_\theta\theta).
 }
\]

At bidegree \((p,q)\), all lower-order terms are already known and the new
coefficient solves

\[
 \boxed{
 \big(\lambda_{pq}I-L\big)Z_{pq}=G_{pq},
 \qquad
 \lambda_{pq}=q_h^pq_\theta^q\in(0,1).
 }
\]

This equation chooses the incoming correction \(Z_{pq}\).  It does not ask
for a solution with both \(h(0)=0\) and \(h(S)=0\), so the
infinite-dimensional endpoint cokernel of the latter problem is absent.

Several useful spectral facts follow.

* Neutral modulation multipliers at \(1\) do not resonate with any
  \(\lambda_{pq}\), because \(p+q\ge1\) implies
  \(\lambda_{pq}<1\).
* If the full inviscid packet-plus-wake return map is a local
  diffeomorphism, then \(L^{-1}\) is bounded and all sufficiently small
  \(\lambda_{pq}\) lie in the resolvent set.  Only finitely many low
  bidegrees can cause a spectral obstruction.
* Point resonances can generically be tested or avoided by changing
  \(r,\gamma,\kappa\).  Continuous spectrum crossing the positive numbers
  \(\lambda_{pq}\), or resolvent norms growing too quickly, is a genuine
  obstruction.
* A resonance with a finite Jordan block may sometimes be handled by
  polyhomogeneous \(j^a h^p\theta^q\) terms.  Such polynomial factors do
  not spoil \(e^{-c j^2}\) flatness.  This does not cure a continuous
  resonant band.

The needed bound is a scale-tame resolvent estimate, for example

\[
 \|(\lambda_{pq}I-L)^{-1}g\|_{X^m}
 \le
 C^{p+q+m+1}
 ((p+q)!)^\sigma(m!)^\sigma
 \|g\|_{Y^{m+d}},
\]

with no carrier loss repeated at every bidegree.

Under compatible nonlinear estimates, this gives

\[
 \|Z_{pq}\|_{X^m}
 \lesssim
 C^{p+q+m}((p+q)!)^\sigma(m!)^\sigma.
\]

This is the discrete replacement for the previously requested
two-endpoint spacetime right inverse.

## 6. All-order size and endpoint accuracy

Truncate the bivariate graph and stage trajectory at total degree

\[
 M_j=\lfloor\eta j\rfloor
\]

for a sufficiently small fixed \(\eta>0\).  Since
\(\rho_j\lesssim r^{-\sigma_0j}\),

\[
\begin{split}
 \log\!\left(
   C^{M_j}(M_j!)^\sigma\rho_j^{M_j+1}
 \right)
 &\le
 -\eta\sigma_0(\log r)j^2
 +O(j\log j)\\
 &\le -c j^2.
\end{split}
\]

Consequently the normalized equation residual and graph mismatch are
\(O(e^{-c j^2})\).

For a forced construction this is already sufficient: define the force to
be the projected equation residual of the truncated trajectory.  No
stability estimate about the high-carrier flow is used.

If one instead insists on an exact unforced segment in the interior of each
stage, and the forward propagator has only polynomial carrier loss, the
one-sided corrector has the same endpoint accuracy:

\[
 \|w_j(S)\|_{X_K^m}
 \lesssim K_j^{D_m}e^{-c j^2}
 \le e^{-c'j^2}.
\]

For a fixed physical derivative order, conversion from semiclassical to
ordinary derivatives costs only another factor \(K_j^{O(1)}\), and hence
does not change superexponential flatness.

If the only available propagator bound is \(e^{CK_j}\), the exact-segment
conclusion is false.  It does not invalidate the forced-residual route.
Increasing \(M_j\) from \(O(j)\) to its Gevrey-optimal value gives
roughly

\[
 \exp\big(-c\rho_j^{-1/\sigma}\big).
\]

For compactly localized Gevrey order \(\sigma>1\) and balanced
\(\rho_j\sim K_j^{-1}\), this is
\(\exp(-cK_j^{1/\sigma})\), still weaker than \(e^{CK_j}\).
Thus analyticity with favorable constants or a tame structured propagator
would be needed to rescue that case.

## 7. Interior residuals, stage seams, and a terminally-flat force

For Clay's forced alternative it is not necessary either to correct the
tiny interior equation error or to make the tiny final graph error vanish
by an unforced endpoint condition.  Use the former as the interior force.
Join two consecutive divergence-free stage states over a fixed normalized
seam interval using a smooth cutoff that is flat at both ends.  If their
normalized difference is

\[
 d_j=O(e^{-c j^2})
\]

in the packet-plus-wake seminorms, then the normalized seam residual is
also \(O(K_j^De^{-c j^2})\) for each fixed seminorm.

In physical variables a normalized residual \(R_j\) gives

\[
 f_j=\frac{a_j^2}{\ell_j}R_j.
\]

After \(m\) spatial and \(n\) temporal derivatives, let
\(K_j^{D_{m,n}}\) bound the additional fixed carrier powers in the chosen
semiclassical-to-ordinary norm conversion.  Then

\[
 \|\partial_x^m\partial_t^nf_j\|_\infty
 \lesssim
 \frac{a_j^2}{\ell_j}
 \ell_j^{-m}\tau_j^{-n}K_j^{D_{m,n}}e^{-c j^2}
 =
 r^{[2\gamma+1+m+(1+\gamma)n+\kappa D_{m,n}]j}e^{-c j^2}.
\]

Since

\[
 T-t_j\asymp\ell_j^{1+\gamma}
 =r^{-(1+\gamma)j},
\]

the right side is \(O((T-t_j)^N)\) for every fixed \(m,n,N\).
The sums of the interior residual forces and seam forces therefore extend
\(C^\infty\)-flatly by zero at \(T\).

An algebraic endpoint mismatch is not enough.  If
\(d_j\sim\rho_j=\ell_j^{\sigma_0}\), then sufficiently high physical
derivatives of the seam force diverge.  The discrete graph must be solved
to increasing order before the flat seam shortcut is legitimate.

The seam construction must act on the complete retained state.  If only
the child is matched while an order-one pressure or wake discrepancy is
discarded, the force is not flat.  Alternatively, one needs a proved
one-way wake theorem showing that the discarded component is already
superexponentially invisible in all child pressure/Leray seminorms.

## 8. Exact comparison with the CDP corrector

Cheskidov--Dai--Palasek first construct a principal inverse cascade \(v\)
which satisfies

\[
 \partial_tv-\Delta v+\mathbb P\operatorname{div}(v\otimes v)
 =\mathbb P\operatorname{div}f.
\]

Their Proposition 5.3 bounds the residual schematically by

\[
 \|\nabla^nf(t)\|_{C^\kappa}
 \lesssim
 \epsilon_0\big(t^{-1-n/2+\alpha}+1\big).
\]

They then solve a forward corrector equation with zero data at the
branching time.  Proposition 6.2 gives a structure-specific linearized
semigroup estimate with only the mild loss

\[
 (t/t')^\epsilon.
\]

The key input is their logarithmic exposure estimate: the number of active
heat shells between \(t'\) and \(t\) is divided by a large frequency-gap
parameter, yielding a coefficient of order \((\log A)^{-1}\) in
\(\log(t/t')\).  Choosing \(\epsilon<\alpha\) makes the Duhamel integral
converge, and Proposition 6.3 closes a fixed point.

Three features do not transfer automatically.

1. **Heat orientation.**  Their \(k\)-th carrier is active on
   \(t\sim N_k^{-2}\) and decays as \(e^{-N_k^2t}\).  A forward terminal
   cascade has \(\varepsilon_jK_j^2\to0\), so its carrier does not receive
   comparable damping.
2. **Direction of time.**  Their solution evolves away from the branching
   time \(t=0\).  If \(s=T-t\) is remaining time in a terminal cascade, the
   analogous forward loss is \((s'/s)^\epsilon\).  An error created at one
   fixed earlier stage is amplified as \(s\downarrow0\); flatness of errors
   generated only at later stages does not erase it.
3. **Endpoint condition.**  Their corrector is an initial-value fixed point.
   It is not designed to hit a recurrent outgoing packet.  That missing
   condition is exactly what the invariant graph or a stable-manifold
   theorem must supply.

What can be borrowed is the methodology: define the principal part through
velocity potentials, prove a semigroup estimate from the actual
frequency/support geometry rather than from
\(\|\nabla v\|_\infty\), and put the residual in a weighted tensor norm
adapted to the Duhamel operator.

## 9. Sharper sufficient theorem

The all-order part of the prize route would close if the following theorem
were proved after constructing the leading localized
Kelvin--Reynolds cell.

There exist a phase-resolved packet-plus-wake scale of spaces \(X^m\), a
leading state \(P\), parameters \(r,\gamma,\kappa\) in the window above,
and one-step maps \(F_{h,\theta}\) such that:

1. \(F_{0,0}(P)=P\), with the required child amplification, strict parent
   drain, exact helicity symmetry, circulation provenance, and finite-energy
   wake;
2. \(F_{h,\theta}\) has a Gevrey-tame two-parameter expansion obtained from
   forward initial-value equations, with no repeated carrier derivative
   loss;
3. every
   \(\lambda_{pq}I-D F_{0,0}(P)\), \(p+q\ge1\), has a compatible
   Gevrey-tame inverse, modulo explicitly treated finite resonances;
4. the full wake is in the graph, or its pressure/Leray coupling into the
   next child is \(O(e^{-c j^2})\);
5. smooth flat seams can be inserted without violating the spatial domain,
   symmetry, energy, and circulation ledgers.

Then the invariant graph truncated at \(M_j\asymp j\), with its interior
equation residual and seam residual declared to be the force, produces a
smooth, terminally-flat forced cascade.  Combined with the already checked
\(1<\gamma<3/2\) energy/dissipation ledger, this would yield the forced Clay
alternative, provided the resulting solution is unbounded at \(T\).

If exact unforced stage interiors are desired, add:

6. the linearized forward Navier--Stokes propagator about the truncated
   stage has
   \(\log C_{j,m}=o(j^2)\), preferably a polynomial carrier bound.

This optional hypothesis lets one apply the one-sided corrector before the
flat seam.  It is not needed for Clay's forced alternative.

The leading localized return cell and items 3--5 are unproved.  No numerical
run can certify them, although a discretized one-step map could usefully
screen for a bad resolvent, a large unstable spectrum, or catastrophic
\(e^{cK}\) growth.

## Primary source

* C. L. Fefferman, *Existence and smoothness of the Navier--Stokes
  equation*, official Clay problem description:
  https://www.claymath.org/wp-content/uploads/2022/06/navierstokes.pdf
* A. Cheskidov, M. Dai, and S. Palasek,
  *Instantaneous Type I blow-up and non-uniqueness of smooth solutions of
  the Navier--Stokes equations*:
  https://arxiv.org/abs/2511.09556
