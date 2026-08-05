# Reconstructed C107--C113: transported-writer fork, Clebsch transport, on-shell pumps, and Weber bookkeeping

Date: 2026-08-05
Status: **reconstructed SELF / exact algebra, conditional PDE boundary**
Checker: [`checks/transported_writer_reconcile_c107_c113.py`](../checks/transported_writer_reconcile_c107_c113.py)

## 0. Scope and provenance

The 2026-08-03 registry reserved C107--C136 for results reported during a
long in-session derivation, but the corresponding files were not preserved.
This note reconstructs the conservative content of C107--C113 from

* `audit/ATTEMPT-2026-08-03-transported-writer.md` (C102--C106),
* `audit/ATTEMPT-2026-08-03b-wake-slaving.md` (C137--C139), and
* exact Clebsch, Beltrami-shell, and Weber identities.

The labels below are therefore **reconstructed labels**, not a claim that
the lost prose has been recovered verbatim.  Nothing here proves a one-cell
Navier--Stokes stage map, a cascade, or blow-up.

Write

\[
 B(v,w):=(v\cdot\nabla)w,
 \qquad J(x)=\gamma (r\cdot x)h,
 \qquad r\cdot h=0,
\]

and let \(\mathbb P\) be the Leray projector.

## 1. C107 -- the wake-slaving/spatial-export fork is genuine

The exact rank-one identity behind C137 is

\[
 B(V,J)=\gamma(V\cdot r)h. \tag{1.1}
\]

It proves **directional protection**: at first order, advection of the
rank-one jet cannot create an \(r\)- or \(t\)-component.  It does not by
itself prove **profile slaving**.  To be a bounded pointwise multiple of the
affine jet, one needs

\[
 V\cdot r=c(t,x)(r\cdot x)+\text{controlled remainder},
 \qquad \|c\|_{L^\infty}<\infty. \tag{1.2}
\]

For a genuine affine **coefficient shift**, rather than merely an
in-direction profile, (c) must additionally be spatially constant on the
active core up to a remainder controlled in the chosen chart norm.

The distinction is exact, not semantic.  The divergence-free constant field
\(V=r\) gives \(B(V,J)=\gamma h\), which is nonzero on the plane
\(r\cdot x=0\), where \(J=0\).  Thus (1.1) alone cannot imply a bounded
pointwise multiplier relation \(B(V,J)=cJ\).

For the raw wakes in C103, with an orthonormal frame \((r,h,t)\),

\[
 W_t=w\,t,
 \qquad W_\zeta=\phi(h-Hr),
\]

one has exactly

\[
 B(W_t,J)=0,
 \qquad B(W_\zeta,J)=-H\gamma\phi\,h. \tag{1.3}
\]

The first channel is dark.  The second is in-direction, but it is an affine
coefficient shift only if the transported profile \(\phi\) has the required
\(s=r\cdot x\) factor with spatially constant quotient on the core (or if
the quotient's variation and the non-affine part are harmless in the chosen
chart norm).  Leray completion is nonlocal, and after conversion the next
pump is no longer the rank-one field \(J\), so neither issue is settled by
(1.1).

Consequently the rigorous continuation has two, and only two currently
identified, branches:

1. **Wake slaving.** Construct a weighted invariant graph on which the
   non-affine/off-chart component is quadratic, schematically
   \(
   \|(I-\Pi_{\rm chart})W\|_X\le C\delta^2\|J\|_X
   \), while the chart component changes the target parameters by
   \(O(\delta)\).  C139 checks the frozen first-order direction count, but
   not this nonlinear graph.
2. **Spatial export.** Move the collar and conversion wake out of the next
   active cell, impose the necessary zero charge/moments, and prove a
   stage-uniform Leray/heat tail estimate before the child becomes a full
   three-direction pump.

This is the precise C107 fork.  The auditor's first-order darkness identities
make the slaving branch plausible; they do not close it.  An export estimate
also requires actual separation, moment cancellation, periodic-image
control, and nonlinear persistence.  Those are hypotheses, not conclusions.

## 2. C108 -- exact Clebsch--Piola transport

Let \(X(t,a)\) be a smooth incompressible flow map,
\(F=\nabla_aX\), so \(\det F=1\).  If two scalars are transported,

\[
 \alpha(t,X(t,a))=\alpha_0(a),
 \qquad
 \beta(t,X(t,a))=\beta_0(a),
\]

then

\[
 \nabla_x\alpha\circ X=F^{-T}\nabla_a\alpha_0,
 \qquad
 \nabla_x\beta\circ X=F^{-T}\nabla_a\beta_0.
\]

The cofactor identity gives

\[
 (F^{-T}p)\times(F^{-T}q)
   ={1\over\det F}F(p\times q).
\]

Hence the Clebsch two-form

\[
 \omega=\nabla\alpha\times\nabla\beta
\]

obeys the exact Cauchy transport law

\[
 \omega(t,X(t,a))=F(t,a)\omega_0(a). \tag{2.1}
\]

This supplies an exact Lagrangian design language for a transported packet:
the desired terminal vorticity can be pulled back through \(F^{-1}\) and
encoded in initial scalar level sets.  It also preserves divergence without
an envelope-by-envelope Leray correction.

Boundary: (2.1) is a vorticity/two-form statement.  The associated velocity
is recovered nonlocally, and arbitrary transported Clebsch scalars do not by
themselves solve the full Euler or Navier--Stokes momentum equation.

## 3. C109 -- one global Clebsch pair cannot carry a helical pump

Suppose a velocity has a single-valued global Clebsch representation

\[
 u=\nabla\varphi+\alpha\nabla\beta,
 \qquad
 \omega=\nabla\times u=\nabla\alpha\times\nabla\beta.
\]

Then

\[
 u\cdot\omega
 =\nabla\varphi\cdot(\nabla\alpha\times\nabla\beta)
 =\nabla\cdot\bigl(\varphi\,\nabla\alpha\times\nabla\beta\bigr). \tag{3.1}
\]

For periodic single-valued potentials, or for decaying potentials with no
boundary flux, (3.1) implies

\[
 \int u\cdot\omega\,dx=0. \tag{3.2}
\]

A nonzero homochiral Beltrami pump has
\(\omega=\sigma K u\), and therefore helicity
\(\sigma K\|u\|_2^2\ne0\).  It cannot be represented globally by one such
Clebsch pair.  Multiple charts/pairs, multivalued potentials, or an explicit
helical remainder are necessary; their overlaps and cross-interactions must
remain in the stage ledger.  Thus Clebsch coordinates simplify transport
but are not an exact one-pair replacement for the pump.

## 4. C110 -- an on-shell homochiral pump is an exact unforced NS orbit

Let \(u_0\) be divergence-free and satisfy

\[
 \nabla\times u_0=\sigma K u_0,
 \qquad \sigma\in\{-1,1\},\quad K>0.
\]

Equivalently on the torus, the zero Fourier mode vanishes and every nonzero
Fourier mode of \(u_0\) lies on the single shell \(|k|=K\) with helicity
\(\sigma\).  Then

\[
 (u_0\cdot\nabla)u_0
 =\nabla {|u_0|^2\over2}-u_0\times(\nabla\times u_0)
 =\nabla {|u_0|^2\over2}, \tag{4.1}
\]

and \(\Delta u_0=-K^2u_0\).  Therefore

\[
 u(t)=e^{-\nu K^2t}u_0,
 \qquad
 p(t)=-e^{-2\nu K^2t}{|u_0|^2\over2} \tag{4.2}
\]

is an exact, unforced Navier--Stokes solution.  No finite-mode truncation is
being used in (4.2).

## 5. C111 -- exact Euler Weber identity

For a smooth Euler solution \(u\) and its flow map \(X\), set
\(F=\nabla_aX\).  Since

\[
 \partial_tF=(\nabla u\circ X)F,
 \qquad
 (\partial_t+u\cdot\nabla)u=-\nabla p,
\]

the pulled-back velocity one-form obeys

\[
 {d\over dt}\bigl(F^Tu\circ X\bigr)
 =F^T\nabla_x\left({|u|^2\over2}-p\right)\circ X
 =\nabla_a\left[\left({|u|^2\over2}-p\right)\circ X\right]. \tag{5.1}
\]

Integration yields the Weber formula

\[
 F(t,a)^Tu(t,X(t,a))
 =u_0(a)+\nabla_a\Phi(t,a), \tag{5.2}
\]

where

\[
 \Phi(t,a)=\int_0^t
 \left({|u|^2\over2}-p\right)(s,X(s,a))\,ds.
\]

Taking the exterior derivative removes \(\Phi\) and recovers the Cauchy
vorticity law.  This is an exact way to design a terminal Euler one-form in
Lagrangian coordinates; it is not a controllability theorem, because
\(F\), \(u\), and \(\Phi\) are coupled.

## 6. C112 -- viscosity is an exact Weber volume defect

For Navier--Stokes, the same computation gives

\[
 {d\over dt}\bigl(F^Tu\circ X\bigr)
 =\nabla_a\left[\left({|u|^2\over2}-p\right)\circ X\right]
  +\nu F^T(\Delta u)\circ X. \tag{6.1}
\]

Thus

\[
 F^Tu\circ X
 =u_0+\nabla_a\Phi
  +\nu\int_0^tF(s)^T(\Delta u)(s,X(s,a))\,ds. \tag{6.2}
\]

The last term is not generally a gradient and cannot be removed by pressure.
Any Weber/Clebsch terminal construction for the one-cell stage must retain
this impulse together with the pressure and heat wakes.  Smallness follows
only from a separately proved high-Reynolds estimate; (6.2) itself supplies
no small parameter.

## 7. C113 -- exact dark tangent and the named activation obstruction

Let \(E_{K,\sigma}\) be the homochiral curl eigenspace in C110 and

\[
 Q(u)=\mathbb P[(u\cdot\nabla)u].
\]

Equation (4.1) says \(Q\equiv0\) on the whole linear space
\(E_{K,\sigma}\).  Consequently, for every \(u,v\in E_{K,\sigma}\),

\[
 DQ(u)[v]
 =\mathbb P\{(u\cdot\nabla)v+(v\cdot\nabla)u\}=0. \tag{7.1}
\]

So a same-shell, same-helicity sideband is exactly dark: it cannot provide
just-in-time quadratic activation of another member of the same eigenspace.
An active terminal interaction must leave that manifold -- for example by
using opposite helicity, a different shell, or a localized envelope.  The
next note checks the smallest opposite-helicity choice on the integer
\(A_2\) hexagon and shows that reality creates an off-shell output of the
same leading size.

The exact named obstruction after C107--C113 is therefore **mode leakage at
activation**.  Clebsch--Piola and Weber identities organize transport, and
the homochiral shell gives an exact unforced decaying pump, but none proves
that the active sideband returns to the next pump while all nonterminal modes
stay within the summable wake allowance.

## 8. Claim ledger

| ID | Reconstructed claim | Status |
|---|---|---|
| C107 | Rank-one protection is directional, not profile slaving; the nonlinear continuation is a wake-slaving or spatial-export fork. | Identity and counterexample EXACT; either PDE branch OPEN |
| C108 | Advected Clebsch gradients satisfy the exact Piola/Cauchy law (2.1). | EXACT |
| C109 | One global single-valued Clebsch pair has zero integrated helicity and cannot carry a nonzero homochiral Beltrami pump. | EXACT under stated boundary hypotheses |
| C110 | A homochiral single-shell Beltrami field decays by (4.2) as an exact unforced NS solution. | EXACT |
| C111 | Euler pulled-back velocity satisfies the Weber identity (5.2). | EXACT |
| C112 | Navier--Stokes adds the non-gradient viscous impulse (6.2). | EXACT; smallness CONDITIONAL |
| C113 | The projected quadratic nonlinearity and its tangent vanish inside one homochiral shell; active departure produces the mode-leakage problem. | Dark tangent EXACT; stage closure OPEN |
