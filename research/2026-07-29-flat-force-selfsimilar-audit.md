# Flat-force self-similar inverse design: no-go phase diagram and one surviving cascade class

Date: 2026-07-29

## Question

Can one prescribe a divergence-free velocity \(u\), singular at a finite time
\(T\), and define

\[
 f=\partial_tu+\mathbb P\nabla\!\cdot(u\otimes u)-\nu\Delta u
\]

so that \(f\) extends smoothly and flatly through \(T\)? Here \(\mathbb P\)
is the Leray projection, so singular gradient terms are correctly assigned to
the pressure rather than to the force.

If this worked on \(\mathbb T^3\), with \(f\) extended by zero after \(T\),
it would target Clay alternative (D). The official statement permits a smooth
periodic force whose every space-time derivative decays faster than any power
of time. A compactly time-supported smooth force is admissible.

The conclusion of this audit is:

1. A fixed, localized self-similar or discretely self-similar profile cannot
   make this construction work. Flat forcing disappears in the blow-up limit,
   and all nondegenerate one-profile balances reduce to profile classes already
   ruled out by Liouville, endpoint-regularity, energy, or backward-heat
   arguments.
2. Small support, a flat amplitude, or rapid oscillation does not hide the
   residual from the \(C^\infty\) force condition. All spatial jets are tested.
3. One sharply defined class survives the ledgers: a globally coupled,
   non-precompact, annular-wake generalized Euler cascade with scale
   \(\ell_j=r^{-j}\), amplitude \(a_j=\ell_j^{-\gamma}\), and
   \(1<\gamma<3/2\), followed by an increasing-order viscous correction.
   It must be zero-net-helicity (or exactly bihelical), multistrand rather
   than a one-tube reset, and infinite-sideband rather than finite Fourier.
   This is formal compatibility, not an existence result. Its missing theorem
   is stated below.

## 1. Flat force is invisible in every blow-up tangent

Under the Navier--Stokes parabolic rescaling

\[
 u_r(y,\tau)=r\,u(x_*+ry,T+r^2\tau),\qquad
 f_r(y,\tau)=r^3 f(x_*+ry,T+r^2\tau),
\]

every smooth force tends to zero locally as \(r\downarrow0\). If \(f\) is flat
at \(T\), the convergence is faster than every power, in every local
space-time seminorm. Thus a Type-I self-similar tangent must solve the
*unforced* equation. The force can seed the flow before the singular regime,
but it cannot supply a principal self-similar balance.

In standard similarity variables

\[
 s=T-t=e^{-\tau},\qquad y=\frac{x-x_*}{\sqrt{s}},\qquad
 u=s^{-1/2}U(y,\tau),
\]

one has

\[
 f=s^{-3/2}\mathcal R(U),
\]

\[
 \mathcal R(U)=
 \partial_\tau U+\frac12(U+y\cdot\nabla U)
 +\mathbb P\nabla\!\cdot(U\otimes U)-\nu\Delta U.
\]

Moreover,

\[
 \partial_t^n\partial_x^m f
 =
 e^{(n+(3+m)/2)\tau}\,\mathcal D_{m,n}\mathcal R(U),
\]

where \(\mathcal D_{m,n}\) is a finite combination of \(\tau\)-derivatives,
\(y\)-derivatives, and the dilation field \(y\cdot\nabla\).
Flatness of \(f\) therefore requires, for every \(K,m,n\),

\[
 \|\mathcal D_{m,n}\mathcal R(U(\tau))\|_\infty
 \lesssim_{K,m,n}
 e^{-[K+n+(3+m)/2]\tau}.
\]

The rescaled residual must be superexponentially small in all relevant
seminorms. A nonzero periodic residual cannot pass this condition. A strict
DSS ansatz therefore has to be an exact unforced periodic orbit of the
rescaled equation.

## 2. The one-profile scaling phase diagram

Let

\[
 u(x,t)=s^{-\alpha}U\!\left(\frac{x-x_*}{s^\beta}\right)
\]

with a localized, smooth, divergence-free profile. The time, nonlinear, and
viscous terms have respective sizes

\[
 s^{-(\alpha+1)},\qquad
 s^{-(2\alpha+\beta)},\qquad
 s^{-(\alpha+2\beta)}.
\]

A nondegenerate flat residual needs at least two leading exponents to match.
There are only four cases.

### 2.1 Triple balance

\[
 \alpha=\beta=\frac12.
\]

The profile obeys the backward Leray equation. A localized profile belongs to
\(L^3(\mathbb R^3)\), and the Nečas--Růžička--Šverák theorem makes it
trivial. Later results rule out much broader Lorentz and Morrey profile
classes. Equivalently, a localized SS/DSS core has bounded critical
\(L^\infty_tL^3_x\), which is incompatible with blow-up.

### 2.2 Time--viscosity balance

\[
 \beta=\frac12,\qquad \alpha<\frac12.
\]

The leading profile satisfies

\[
 \alpha U+\frac12y\cdot\nabla U-\nu\Delta U=0.
\]

Fourier transformation gives, along every ray,

\[
 \widehat U(r,\theta)
 =C(\theta)r^{2\alpha-3}e^{\nu r^2}.
\]

No nonzero tempered, let alone localized, profile exists.

### 2.3 Nonlinearity--viscosity balance

\[
 \alpha=\beta>\frac12.
\]

The leading profile is a finite-energy stationary unforced Navier--Stokes
field:

\[
 \mathbb P\nabla\!\cdot(U\otimes U)-\nu\Delta U=0.
\]

Taking the \(L^2\) inner product with \(U\) yields
\(\nu\|\nabla U\|_2^2=0\), hence a localized profile is zero.

### 2.4 Euler-dominated time--nonlinearity balance

\[
 \alpha=1-\beta,\qquad \beta<\frac12.
\]

Bounded core energy requires

\[
 3\beta-2\alpha=5\beta-2\ge0,
\]

so \(2/5\le\beta<1/2\). The leading Euler profile obeys

\[
 \alpha U+\beta y\cdot\nabla U+
 \mathbb P\nabla\!\cdot(U\otimes U)=0.
\]

Its \(L^2\) identity is

\[
 \left(\alpha-\frac{3\beta}{2}\right)\|U\|_2^2=0.
\]

Thus a nonzero localized fixed profile forces

\[
 (\alpha,\beta)=\left(\frac35,\frac25\right).
\]

This is the energy-conserving, Type-I Euler scale. Chae--Wolf exclude the
corresponding finite-energy Type-I atomic concentration and, as a corollary,
energy-conserving Euler DSS blow-up. A strict SS profile is a special case.

Degenerate choices in which a formally dominant term vanishes (stationary
Euler, harmonic, or homogeneous profiles) merely expose the next algebraic
term. A finite polyhomogeneous expansion then ends with a non-flat residual.
Avoiding this requires an infinite correction scheme or a non-precompact
renormalized orbit.

## 3. The all-jet obstruction to naive inverse design

Take a localized bump

\[
 u=a(s)U\!\left(\frac{x-x_*}{\ell(s)}\right).
\]

Ignoring cancellation, spatial derivatives of the three projected residual
terms have sizes

\[
 \begin{aligned}
 \|\partial_x^m\partial_tu\|_\infty
 &\sim \ell^{-m}
 \left(|a_s|+\left|a\frac{\ell_s}{\ell}\right|\right),\\
 \|\partial_x^m\mathbb P\nabla\!\cdot(u\otimes u)\|_\infty
 &\sim a^2\ell^{-m-1},\\
 \|\partial_x^m\Delta u\|_\infty
 &\sim a\ell^{-m-2}.
 \end{aligned}
\]

If the viscous term is made flat *term by term* for every \(m\), then
\(a\ell^{-q}\) is flat for every fixed \(q\). All spatial derivatives of
\(u\) then extend smoothly through \(T\). The same conclusion holds for
\(u=aW(k(s)x)\): requiring \(ak^{m+2}\) to be flat for every \(m\) also
makes every jet of \(u\) flat.

Therefore shrinking support, superflat amplitude, and oscillation cannot by
themselves produce a singular velocity with a \(C^\infty\)-flat force.
They only work if the individually singular time, nonlinear, and viscous
terms cancel to arbitrarily high order. Intermittency helps integral norms,
but the Clay force condition tests supremum-size space-time jets.

## 4. The surviving wake-carrying cascade window

The fixed-profile energy identity can be avoided if the active child does not
carry all the parent's energy. Let

\[
 \ell_j=r^{-j},\qquad a_j=\ell_j^{-\gamma},\qquad r>1.
\]

Use an Euler turnover time

\[
 \Delta t_j\sim\frac{\ell_j}{a_j}=\ell_j^{1+\gamma}.
\]

If the tail of the stage durations is identified with \(s=T-t\), this
corresponds to generalized similarity exponents

\[
 \beta=\frac{1}{1+\gamma},\qquad
 \alpha=\frac{\gamma}{1+\gamma}=1-\beta.
\]

The exact compatible window is

\[
 \boxed{1<\gamma<\frac32.}
\]

Equivalently, \(2/5<\beta<1/2\) and \(1/2<\alpha<3/5\): the open
Euler-dominated sector from the phase diagram. Unlike a fixed profile, the
front sheds energy into wakes, so its active packet is not subject to the
fixed-profile \(L^2\) identity.

The ledgers are:

\[
 \begin{array}{lll}
 \text{terminal time:}
 &\sum_j\Delta t_j<\infty,\\[2mm]
 \text{Reynolds number:}
 &a_j\ell_j=\ell_j^{1-\gamma}\to\infty,\\[2mm]
 \text{critical norm:}
 &\|u_j\|_3\sim a_j\ell_j=\ell_j^{1-\gamma}\to\infty,\\[2mm]
 \text{active energy:}
 &E_j\sim a_j^2\ell_j^3=\ell_j^{3-2\gamma}\to0,\\[2mm]
 \text{stage dissipation:}
 &\nu\|\nabla u_j\|_2^2\Delta t_j
 \sim\nu\ell_j^{2-\gamma},\\[2mm]
 \text{relative viscous loss:}
 &\varepsilon_j:=
 \frac{\nu}{a_j\ell_j}
 =\nu\ell_j^{\gamma-1}\to0,\\[2mm]
 \text{Type-I gradient size:}
 &(T-t_j)\|\nabla u_j\|_\infty\sim1.
 \end{array}
\]

Moreover

\[
 \|\omega_j\|_\infty\Delta t_j\sim1,
\]

so the Beale--Kato--Majda integral accumulates one order-one contribution per
stage and diverges, even though the total energy and viscous dissipation
ledgers converge.

The child-to-parent energy ratio is

\[
 q_E=\frac{E_{j+1}}{E_j}=r^{2\gamma-3}<1.
\]

The missing energy must enter a controlled wake. Since \(E_j\) is geometric,
the total wake-energy ledger can remain finite. The active energy tends to
zero, so this mechanism does not form the fixed atomic energy profile ruled
out at \(\gamma=3/2\). Its \(L^3\) norm diverges, so endpoint regularity does
not apply.

This is not an exact Navier--Stokes DSS symmetry. It is a generalized Euler
renormalization with a stage-dependent small viscosity
\(\varepsilon_j\). The flow is necessarily non-precompact in the usual
Navier--Stokes similarity variables.

## 5. Helicity and circulation obstructions

These are leading constraints on any recurrent cell, not bookkeeping details.

### 5.1 Helicity

For an isotropic localized packet,

\[
 H_j=\int u_j\cdot\omega_j\,dx
 \sim a_j^2\ell_j^2
 =\ell_j^{2-2\gamma}.
\]

For \(\gamma>1\), its characteristic magnitude grows without bound, with
ratio

\[
 \frac{H_{j+1}}{H_j}=r^{2\gamma-2}>1.
\]

Smooth Euler evolution conserves total helicity. Since viscosity is only an
\(O(\varepsilon_j)\) perturbation per stage and the force is beyond all
orders, neither can support an order-one helicity multiplication at every
stage. A recurrent leading cell must therefore satisfy one of:

1. exactly zero net helicity at leading order and at every solvability order;
2. a bihelical daughter state with equally leading opposite-helicity
   components included in the recurrent state.

An old, low-frequency wake cannot absorb an unbounded opposite helicity.
If a wake carries the cancellation, it must contain an equally leading
high-frequency companion, which is the second option in substance.

### 5.2 Circulation and vortex topology

The characteristic circulation grows:

\[
 \Gamma_j\sim a_j\ell_j=\ell_j^{1-\gamma},\qquad
 \frac{\Gamma_{j+1}}{\Gamma_j}=r^{\gamma-1}>1.
\]

A single material vortex tube cannot recur with this scaling: Kelvin
circulation and Helmholtz flux conservation preserve its circulation in the
Euler limit. The surviving cell must therefore be demonstrably one of:

1. non-tubular and strain-dominated, so the measured child loop is not the
   advected parent loop;
2. an aggregator of an increasing number of same-sign circulation strands;
3. a mechanism with leading viscous reconnection.

The third option is in tension with
\(\varepsilon_j\to0\). An order-one reconnection would require a thinner
layer in which viscosity is leading, invalidating the regular perturbative
viscosity hypothesis below. Thus circulation aggregation/topology is an
additional theorem gate for the Euler-cell route.

## 6. The exact missing all-order cell theorem

The following is a sufficient *hypothetical* theorem. No such theorem is
currently known for the true three-dimensional Euler/Navier--Stokes
nonlinearity.

On stage \(j\), write

\[
 x=x_j+\ell_j\xi,\qquad
 t=t_j+\frac{\ell_j}{a_j}\sigma,\qquad
 u(x,t)=a_jV_j(\xi,\sigma).
\]

The projected equation becomes

\[
 \partial_\sigma V_j+
 \mathbb P\nabla_\xi\!\cdot(V_j\otimes V_j)
 -\varepsilon_j\Delta_\xi V_j
 =
 \frac{\ell_j}{a_j^2}f,
\qquad
 \varepsilon_j=\frac{\nu}{a_j\ell_j}.
\]

Thus an inviscid unit cell followed by a uniform expansion in
\(\varepsilon_j\) is exactly the correction problem that must be solved.

### Hypothesis A: an inviscid handoff cell

There are \(r>1\), \(1<\gamma<3/2\), a time \(S>0\), and a smooth,
localized, divergence-free Euler trajectory \(V_0(\sigma)\),
\(0\le\sigma\le S\), in a packet-plus-wake phase space such that:

1. the incoming state is a normalized parent packet plus an admissible
   background;
2. the outgoing state contains a daughter
   \(r^\gamma P(r(\,\cdot-\zeta))\), with the same polarization/phase class,
   plus an admissible wake;
3. the daughter carries energy fraction \(q_E=r^{2\gamma-3}<1\);
4. the total leading helicity is exactly zero, or the recurrent state is
   explicitly bihelical;
5. the circulation increase is realized by a proved non-tubular or
   multi-strand topology, not by relabeling a single material tube;
6. translations, rotations, dilation, phase, and other neutral directions
   are fixed by modulation conditions.

### Hypothesis B: tame right-invertibility

Let

\[
 \mathcal L_{V_0}h
 =
 \partial_\sigma h+
 \mathbb P\nabla\!\cdot(V_0\otimes h+h\otimes V_0)
\]

with the incoming/outgoing packet-and-wake boundary conditions. After
quotienting the neutral modes and imposing the energy, helicity, circulation,
and phase solvability conditions, \(\mathcal L_{V_0}\) has a right inverse
\(\mathcal G\) on a scale of localized Banach spaces:

\[
 \|\mathcal G g\|_{X^m}
 \le C_m\|g\|_{Y^{m+d}},
\qquad
 C_m\le C^{m+1}(m!)^\sigma,
\]

for fixed derivative loss \(d\) and Gevrey index \(\sigma\), uniformly under
the stage rescalings. The quadratic Euler map obeys compatible tame
estimates. The endpoint return map has the same bounds.

### Hypothesis C: all-order viscous solvability

For

\[
 \mathcal F_\varepsilon(V)
 =
 \partial_\sigma V+
 \mathbb P\nabla\!\cdot(V\otimes V)-\varepsilon\Delta V,
\]

there are coefficients \(V_n\) and modulation/wake corrections satisfying

\[
 \mathcal L_{V_0}V_n
 =
 \Delta V_{n-1}
 -
 \sum_{\substack{p+q=n\\p,q\ge1}}
 \mathbb P\nabla\!\cdot(V_p\otimes V_q)
 +\text{modulation terms},
\]

with the same renormalized endpoint condition and Gevrey bounds
\(\|V_n\|_{X^m}\lesssim C^{n+m}(n!)^\sigma(m!)^\sigma\).
The helicity and circulation solvability conditions must hold at every order.

At stage \(j\), truncate at

\[
 M_j=\lfloor\kappa j\rfloor.
\]

Since

\[
 \varepsilon_j=\nu r^{-(\gamma-1)j},
\]

the normalized remainder obeys

\[
 C^{M_j}(M_j!)^\sigma\varepsilon_j^{M_j+1}
 \le e^{-c j^2}
\]

for sufficiently small fixed \(\kappa>0\) and all large \(j\). Factorial
Gevrey growth is harmless because \(j^2\) dominates \(j\log j\).

### Hypothesis D: wake closure and nonlocal separation

The accumulated wake must satisfy one of two strong alternatives.

1. **Closed background state.** Its constant velocity, strain, pressure jet,
   and any other leading child interaction are part of the renormalized state
   and return under the cell map.
2. **Causal separation.** After removing constant advection by translating
   the child center, all child--wake velocity, strain, pressure, and Leray-tail
   interactions are \(O(e^{-cj^2})\) in the stage-normalized seminorms.

For the second alternative, ordinary support separation is insufficient
because pressure is nonlocal. One needs exact cancellation or increasing
moment cancellation of the wake stresses, with quantitative control of the
periodic/Euclidean pressure kernel. Schematically, for the unmodeled old wake
\(W_{<j}\) in the child region,

\[
 \frac{\|W_{<j}-W_{<j}(x_j)\|_\infty}{a_j}
 +
 \frac{\ell_j\|\nabla W_{<j}\|_\infty}{a_j}
 +
 \text{normalized pressure-tail seminorms}
 \lesssim e^{-cj^2}.
\]

Any prescribed low-mode wake strain used as the amplifier belongs in the
closed state, not in the error.

There is a further reason to prefer the first alternative.  In the natural
geometric scaling, the older wakes become a nested annular train in
child-rescaled variables.  The \(n\)-th older wake has characteristic volume
\(r^{3n}\), velocity \(r^{-\gamma n}\), and vorticity
\(r^{-(1+\gamma)n}\).  Consequently

\[
 \sum_{n\ge0}r^{3n-q(1+\gamma)n}
\]

diverges when \(q<3/(1+\gamma)\).  Retaining this global annular state
therefore escapes the low-\(q\) vorticity integrability assumptions in
asymptotic-DSS Liouville theorems.  Exiling the wake while keeping a compact
periodic core risks putting the construction back inside those no-go
hypotheses.  The most defensible survivor is thus a wake-inclusive global
renormalized block, not a passive-wake core.

### Hypothesis E: flat gluing

The stage endpoint mismatch and time-cutoff residual have physical size

\[
 \frac{a_j^2}{\ell_j}e^{-cj^2}.
\]

After \(m\) spatial and \(n\) temporal derivatives, the worst scale factors
are

\[
 \ell_j^{-m}\Delta t_j^{-n}
 =
 r^{[m+(1+\gamma)n]j}.
\]

Multiplying by any negative power of
\(s_j=T-t_j\asymp\ell_j^{1+\gamma}\) costs only another \(e^{Cj}\).
The \(e^{-cj^2}\) remainder therefore beats every required space-time jet and
every power of \(s_j\). The summed force extends \(C^\infty\)-flatly by zero
at \(T\).

Under Hypotheses A--E, the stage intervals are locally finite on
\([0,T)\), the total energy and dissipation ledgers converge, while
\(\|u(t_j)\|_\infty\sim a_j\to\infty\). Smooth-solution uniqueness before
\(T\) would then turn the constructed singular solution into a Clay-(D)
breakdown datum.

Again, this paragraph is conditional. Hypothesis A, the return-cell
right inverse in Hypothesis B, the helicity/circulation conditions, and the
wake-pressure closure are not presently proved.

## 7. An exact affine Kelvin amplifier inside the surviving window

There is an exact true-Euler mechanism for the active low-to-high
amplification half of Hypothesis A.  It is a particularly simple member of
the known Kelvin/Craik--Criminale class of exact waves on affine base flows,
specialized here to expose the cascade exponent and full-Laplacian ledger.
Let

\[
 S=\operatorname{diag}(-\alpha,-\beta,\alpha+\beta),
 \qquad
 \frac{\beta}{\alpha}=\gamma\in(1,3/2),
\]

and define

\[
 u(x,t)
 =
 Sx+A(t)e_2\sin(k(t)x_1),\qquad
 p(x,t)=-\frac12x\cdot S^2x.
\]

The field is divergence free. Since \(S\) is symmetric,
\((Sx\cdot\nabla)Sx=S^2x=-\nabla p\). For
\(w=Ae_2\sin(kx_1)\),

\[
 (w\cdot\nabla)w=0,\qquad
 (w\cdot\nabla)Sx=-\beta w,
\]

and

\[
 (\partial_t+Sx\cdot\nabla)w
 =
 A'e_2\sin(kx_1)
 +A(k'-\alpha k)x_1e_2\cos(kx_1).
\]

Consequently,

\[
 k'=\alpha k,\qquad A'=\beta A
\]

makes \(u\) an exact Euler solution, with

\[
 k(t)=k_0e^{\alpha t},\qquad
 A(t)=A_0e^{\beta t}
 =A_0\left(\frac{k(t)}{k_0}\right)^\gamma.
\]

This exactly realizes the required amplitude-frequency law with no wave
self-advection and no off-shell leakage. The shear wave has zero
self-helicity. Its growing Eulerian wavelength-scale circulation
\(A/k\sim k^{\gamma-1}\) does not contradict Kelvin's theorem because the
wavelength-sized Eulerian loop is not the advected parent loop. Thus the
mechanism gives a concrete version of the non-tubular/strain-dominated option
in the circulation ledger.

The full Laplacian is also exactly auditable. For Navier--Stokes the same
ansatz solves the equation when

\[
 k'=\alpha k,\qquad A'=(\beta-\nu k^2)A,
\]

so, over the time \(h=(\log r)/\alpha\) needed to send \(k_0\) to \(rk_0\),

\[
 \frac{A(h)}{A_0}
 =
 r^\gamma
 \exp\!\left[-\frac{\nu k_0^2(r^2-1)}{2\alpha}\right].
\]

For a stage with \(k_0\sim\ell_j^{-1}\) and
\(\alpha\sim a_j/\ell_j\), the damping exponent is
\(O(\nu/(a_j\ell_j))=O(\varepsilon_j)\). It can be compensated by a
stage modulation, and it confirms rather than merely estimates the viscous
ledger.

What this solves:

1. an exact true-Euler low-mode-strain amplifier with any
   \(\gamma\in(1,3/2)\);
2. exact polarization and phase closure for one shear family;
3. an explicit full-Laplacian correction during the active amplification;
4. a zero-self-helicity and non-material-loop circulation mechanism.

What it does not solve:

1. \(Sx\), its quadratic pressure, and the plane wave are nonperiodic,
   nonlocalized, and infinite-energy;
2. a material support transported by \(S\) is squeezed by
   \(r^{-1}\) and \(r^{-\gamma}\) in two directions but stretched by
   \(r^{1+\gamma}\) in the third. It is not an isotropically smaller child;
3. the stationary affine strain is an infinite energy reservoir. In a finite
   packet, wave growth has to drain a parent strain and leave a controlled
   wake;
4. cutting off \(Sx\) and the wave creates leading transport, pressure, and
   divergence corrections at the boundary;
5. one growing shear polarization has
   \(w\cdot\nabla w\equiv0\), hence no Reynolds backreaction at all; a
   common-direction shear family remains nonlinear-dark and cannot prescribe
   the next strain tensor;
6. it creates no smaller localized strain packet and no recurrent reset.

The sharpened candidate for Hypothesis A is therefore a localized
**Kelvin-amplify / rank-one-reset** cell: a finite family of separated Kelvin
waves is amplified by an approximately affine parent strain; their
high-high stress drains that strain and creates the next smaller strain
packet; the old strain becomes a wake. The affine calculation makes the
amplification exact. Localization, energy transfer, finite-rank compatibility,
pressure closure, and the return map remain the prize-level part.

The qualifier “finite family” refers only to the carrier orientations.  Once
their amplitudes are localized, divergence corrections and nonlinear
interactions generate an infinite sideband hierarchy.  That hierarchy is
unavoidable: Kishimoto--Yoneda rule out a nonstationary exact finite-Fourier
Euler switch on \(\mathbb T^3\).  The common-direction shear family is
therefore a hard no-go for reset, while a finitely oriented but
infinite-sideband localized stress remains open.

There is also a nonempty full-Laplacian geometric-optics window.  If an
internal carrier number is chosen as

\[
 K_j=\ell_j^{-\kappa},\qquad N_j=K_j/\ell_j,
\]

then the first localization error is \(K_j^{-1}=\ell_j^\kappa\), whereas
heat relative to the parent strain is

\[
 \varepsilon_jK_j^2
 =\ell_j^{\gamma-1-2\kappa}.
\]

Thus

\[
 0<\kappa<\frac{\gamma-1}{2}
\]

keeps both effects perturbative.  The balanced choice
\(\kappa=(\gamma-1)/3\) makes them both
\(O(\ell_j^{(\gamma-1)/3})\).  Heat scaling therefore does not kill an
all-order localized Kelvin/rank-one expansion; construction and uniform
solvability do.

## 8. Comparison with the Cheskidov--Dai--Palasek inverse cascade

The 2025 Cheskidov--Dai--Palasek construction is the closest exact
Navier--Stokes technology to parts of Hypotheses A--E, but its time orientation
and energy source are the opposite of a forward Clay breakdown.

Their solution agrees with a classical background through a time \(T_*\),
has a smooth spatial value at \(T_*\), and blows up as \(t\downarrow T_*\)
*from the right*. For \(t>T_*\), an infinite hierarchy already present at
arbitrarily high wavenumber transfers energy high-to-low. They prove that the
inverse flux from infinite frequency is infinite and explicitly note that the
new branch lies in neither \(L^\infty_tL^2_x\) nor \(L^2_tH^1_x\) near
\(T_*\). A classical branch with the same data continues to exist; the result
is weak nonuniqueness, not failure of global smooth existence.

This cannot be turned into a Clay-(D) example by reversing time:

1. Navier--Stokes time reversal changes the forward heat operator into an
   anti-diffusive one.
2. Their robust nonlinear direction is high-high-to-low forward in time.
   A forward terminal singularity needs the causal production of a higher
   child from lower modes.
3. Their high-frequency reservoir is present at the branching time and
   injects unbounded energy. In the wake-cascade program, every high-frequency
   child has to be generated by the unique smooth evolution before \(T\), with
   a uniform energy ledger.

Nevertheless, several pieces are directly informative.

### Technology that can transfer

- Their rank-one geometric identity prescribes a low-frequency symmetric
  stress as
  \[
  \sum_j a_{j,k+1}^2\theta_j\otimes\theta_j
  =
  2\mathcal D\sum_jN_{j,k}\psi_{j,k}+pI.
  \]
  This is an exact high-high-to-low stress realization. It can inform the
  construction of a prescribed wake/background stress once the high child is
  already available.
- Their oscillatory profiles use velocity potentials before applying the
  Laplacian and Leray projection. This avoids losing localization immediately
  to a nonlocal projection and is relevant to Hypothesis D.
- In dimension at least three, disjoint periodic pipe supports, shrinking
  intersection sets of volume \(O(2^{-k})\), and separated heat time scales
  suppress cross-polarization and cross-shell errors. These are useful models
  for wake geometry and pressure-tail estimates.
- The individual linearly polarized waves
  \(\theta\sin(N\eta\cdot x)\), \(\theta\perp\eta\), have zero self-helicity.
  Their basis is therefore compatible with the zero-helicity side of the
  ledger, although cross-block helicity still has to be controlled.
- Once a principal inverse cascade is known, their forward linearized
  semigroup estimates and contraction argument absorb a small residual into
  an exact solution. This is evidence that a *forward-parabolic* corrector can
  close around a correctly oriented principal construction.

### Technology that does not supply the missing theorem

- The recursive amplitude identity creates the *lower* shell from an already
  existing higher shell. It does not provide Hypothesis A's causal
  parent-to-higher-daughter return.
- Their corrector is a forward initial-value problem with zero data at the
  branching time. It is not a two-endpoint, modulated right inverse for a
  recurrent cell and gives no control of the return-map cokernel in
  Hypothesis B.
- Viscosity is incorporated exactly through factors such as
  \(e^{-N_{j,k}^2t}\). They do not construct a uniform all-order expansion in
  the small stage viscosity \(\varepsilon_j\), as required in Hypothesis C.
- Their principal residual tensor is small only in a singular weighted norm;
  schematically its \(n\)-th spatial derivative is bounded by
  \(\epsilon_0(t^{-1-n/2+\alpha}+1)\). It is then canceled by the weak
  corrector. It is not a \(C^\infty\)-flat external force and cannot be used
  as Hypothesis E's gluing remainder.
- Periodic pipe localization and volume smallness control integral and
  Hölder norms, but do not by themselves give the superalgebraic
  child--old-wake pressure separation required at every physical jet.
- The oscillatory pipes are not a recurrent material vortex tube, which is
  helpful conceptually, but the construction does not prove the increasing
  circulation-strand aggregation demanded by
  \(\Gamma_{j+1}/\Gamma_j=r^{\gamma-1}>1\).

The most realistic reuse is therefore asymmetric: borrow their geometric
stress decomposition, potential localization, support separation, and
forward corrector *after* a genuine low-to-high Euler return cell has been
found. Reversing their completed solution does not provide that cell.

## 9. Relation to strict SS/DSS literature

The no-go conclusion is intentionally limited to the hypotheses actually
covered by known results:

- Nečas--Růžička--Šverák exclude backward Leray profiles in
  \(L^3(\mathbb R^3)\).
- Tsai and later work cover larger local-energy, Lorentz, Marcinkiewicz, and
  Morrey profile classes. In particular, the 2022 Morrey-space theorem of
  Jiu--Wang--Wei summarizes the extensions.
- Hou--Li exclude locally self-similar blow-up when the rescaled profile
  converges in \(L^p\), \(3<p<\infty\).
- Chae excludes locally asymptotically DSS Navier--Stokes blow-up for a
  time-periodic profile in
  \(C^1_s(L^3_y\cap C^2_y)\).
- Chae--Wolf exclude energy-conserving Euler DSS through their Type-I atomic
  energy-concentration theorem.
- Chae--Wolf's Navier--Stokes DSS removal theorem only removes general DSS
  scaling factors near one; it should not be misquoted as excluding every
  large-ratio DSS profile without integrability hypotheses.

Thus a wildly nondecaying or non-Morrey profile can remain outside a named
Liouville theorem. But embedding such a profile into Clay (C) or (D) while
keeping initial data and every force jet rapidly decaying creates an outer
matching residual. Unless that residual is canceled to all orders, it is not
flat. The wake-carrying cascade above states the concrete noncompact
alternative rather than relying on an uncontrolled tail.

## Primary sources

- Charles Fefferman, official Clay problem statement:
  https://www.claymath.org/wp-content/uploads/2022/06/navierstokes.pdf
- Nečas, Růžička, Šverák, *On Leray's self-similar solutions of the
  Navier--Stokes equations*:
  https://doi.org/10.1007/BF02551584
- Escauriaza, Seregin, Šverák, *\(L_{3,\infty}\)-solutions ... and backward
  uniqueness*:
  https://www.mathnet.ru/eng/rm609
- Hou and Li, *Nonexistence of local self-similar blow-up*:
  https://arxiv.org/abs/math/0603126
- Chae, *Remarks on asymptotically discretely self-similar solutions*:
  https://arxiv.org/abs/1306.0305
- Chae and Wolf, *On the Liouville type theorems for self-similar solutions*:
  https://arxiv.org/abs/1609.06962
- Chae and Wolf, *Removing discretely self-similar singularities*:
  https://arxiv.org/abs/1610.09464
- Guevara and Phuc, Marcinkiewicz/Morrey profile theorem:
  https://doi.org/10.1137/16M110099X
- Jiu, Wang, Wei, Morrey profile theorem:
  https://arxiv.org/abs/2006.15776
- Chae and Wolf, *Energy concentrations and Type I blow-up for the 3D Euler
  equations*:
  https://arxiv.org/abs/1706.02020
- Fabijonas and Holm, *Multi-frequency Craik--Criminale solutions of the
  Navier--Stokes equations*:
  https://arxiv.org/abs/nlin/0304049
- Cheskidov, Dai, and Palasek, *Instantaneous Type I blow-up and
  non-uniqueness of smooth solutions of the Navier--Stokes equations*:
  https://arxiv.org/abs/2511.09556
- Kishimoto and Yoneda, *Characterization of three-dimensional Euler flows
  supported on finitely many Fourier modes*:
  https://arxiv.org/abs/2110.08039
