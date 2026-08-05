# Exact Navier--Stokes endpoint collars around steady Euler bubbles

**Date:** 2026-08-03

**Status:** uniform short-time Sobolev perturbation lemma and exact scaling
ledger; all-order Gevrey endpoint matching remains open.

**Scope:** the bare-endpoint obstruction in C38.  This note constructs the
canonical local Navier--Stokes collar launched from a steady Euler bubble.
It does not connect that collar to the amplified three-phase packet or
construct an infinite cascade.

## 1. Outcome

A nonzero steady Euler bubble cannot be held fixed on a Navier--Stokes
endpoint interval.  It need not be.  Use it only as the initial value of
the exact unforced Navier--Stokes evolution.

After rescaling one bubble to unit size and one turnover to unit time, the
only new parameter is
\[
                         \mu={\nu\over a\ell}
                            =\operatorname{Re}^{-1}.            \tag{1.1}
\]
For every fixed smooth steady Euler seed, the exact Navier--Stokes solution
exists on a parameter-uniform interval and remains \(O(\mu)\)-close to the
seed there.  Its first time derivative is exactly \(\mu\Delta V\), and all
higher viscous endpoint jets are generated automatically by the equation.

Thus the endpoint bath/wake required by C38 has a canonical realization of
relative size \(O(\operatorname{Re}^{-1})\).  The remaining issue is not
existence of the endpoint jet.  It is matching the active oscillatory
transition to this exact nonstationary collar with cascade-uniform
Gevrey bounds.

## 2. Dimensionless collar equation

Let \(V\) be a smooth divergence-free steady Euler field on
\(\mathbb R^3\) or a periodic box:
\[
                  \mathbb P(V\cdot\nabla V)=0.                 \tag{2.1}
\]
For a bubble of velocity size \(a\) and length scale \(\ell\), write
\[
 u(t,x)=aU(\tau,y),\qquad
 y={x-x_0\over\ell},\qquad
 \tau={a\over\ell}(t-t_0).                                   \tag{2.2}
\]
The unforced Navier--Stokes equation becomes
\[
 \partial_\tau U+\mathbb P(U\cdot\nabla U)=\mu\Delta U,
 \qquad U(0)=V,                                                \tag{2.3}
\]
with \(\mu\) given by (1.1).

Put \(U=V+W\).  Using (2.1),
\[
 \partial_\tau W
 +\mathbb P\big(V\cdot\nabla W+W\cdot\nabla V
                    +W\cdot\nabla W\big)
 =\mu\Delta W+\mu\Delta V,
 \qquad W(0)=0.                                                \tag{2.4}
\]
This is the exact endpoint-collar equation.  The source is proportional to
\(\mu\); it is not an order-one seam.

## 3. Uniform short-time perturbation lemma

Fix \(s>5/2\) and a smooth seed \(V\).  The standard \(H^s\) commutator
estimate applied to (2.4) gives, while \(\|W\|_{H^s}\le1\),
\[
 {d\over d\tau}\|W\|_{H^s}
 \le C_{s,V}\|W\|_{H^s}
     +C_s\|W\|_{H^s}^2
     +\mu\|\Delta V\|_{H^s}.                                 \tag{3.1}
\]
The term \(-\mu\Delta W\) is dissipative and does not worsen the
constant.  A bootstrap and Gronwall therefore give numbers
\(\tau_0,C>0\), depending on \(s,V\) but independent of
\(0\le\mu\le1\), such that
\[
 \sup_{0\le\tau\le\tau_0}
                 \|U(\tau)-V\|_{H^s}
 \le C\mu\tau.                                                 \tag{3.2}
\]
The same estimate supplies a viscosity-uniform existence interval, since
the Euler equation is the endpoint \(\mu=0\) and viscosity only contributes
dissipation.

For physical variables, use the scale-adapted norm
\[
 \|f\|_{H^s_\ell}^2
 =\sum_{m=0}^{\lfloor s\rfloor}
       \ell^{2m-3}\|\nabla^mf\|_2^2
   +\text{the corresponding fractional term}.                 \tag{3.3}
\]
Then a rescaled velocity \(aF((x-x_0)/\ell)\) has norm
\(a\|F\|_{H^s}\), and (3.2) becomes
\[
 \sup_{0\le t-t_0\le\tau_0\ell/a}
 \|u(t)-aV((\mathord\cdot-x_0)/\ell)\|_{H^s_\ell}
 \le Ca\operatorname{Re}^{-1}{a(t-t_0)\over\ell}.             \tag{3.4}
\]
In particular, the relative deformation over one fixed fraction of a
turnover is \(O(\operatorname{Re}^{-1})\).

## 4. The complete endpoint jet is canonical

Differentiating (2.3) at \(\tau=0\) and using (2.1) gives
\[
                         U_1:=\partial_\tau U(0)=\mu\Delta V.  \tag{4.1}
\]
The next derivative is
\[
 U_2
 =\mu^2\Delta^2V
  -\mu\mathbb P\big(V\cdot\nabla\Delta V
                    +\Delta V\cdot\nabla V\big).              \tag{4.2}
\]
More generally, if \(U_n=\partial_\tau^nU(0)\), then
\[
 U_{n+1}
 =\mu\Delta U_n
  -\mathbb P\sum_{a=0}^n{n\choose a}
                   U_a\cdot\nabla U_{n-a}.                    \tag{4.3}
\]
This exact recurrence is the full forward Navier--Stokes jet.  No Borel
guess or separately designed pressure tensor is needed to determine it.

Every nonzero \(U_n\), \(n\ge1\), contains at least one factor of \(\mu\).
For each fixed \(n\), the collar jet is therefore a perturbation of the
steady Euler endpoint.  Uniform control when \(n\) grows with the cascade
index is a separate Gevrey estimate and is not asserted here.

## 5. Energy, momentum, and reflection symmetry

The exact collar satisfies
\[
 {1\over2}\|U(\tau)\|_2^2
 +\mu\int_0^\tau\|\nabla U(s)\|_2^2\,ds
 ={1\over2}\|V\|_2^2.                                        \tag{5.1}
\]
For a fixed turnover interval, its relative energy loss is \(O(\mu)\).
Indeed, in physical units
\[
 {\nu(\ell/a)a^2\ell\over a^2\ell^3}
 ={\nu\over a\ell}=\mu.                                      \tag{5.2}
\]

Spatial mean is preserved because the nonlinear and viscous terms have
zero integral.  If \(V\) is equivariant under an orthogonal lattice
symmetry, uniqueness preserves that symmetry throughout the collar.  In
particular, an orientation-reversing reflection pair retains exactly zero
total helicity at every time: the solution is reflection equivariant while
helicity is a pseudoscalar.

The heat equation has infinite propagation speed, so a compact initial
bubble immediately develops a global tail.  On the torus this is exactly
the periodic wake of C86; it is part of the solution, not an external
error.

## 6. What this repairs and what remains

The bare endpoint in C38 failed because a constant-in-time Euler bubble has
residual \(-\nu\Delta V\).  Equations (2.3) and (4.1) cancel that residual
exactly.  The collar therefore supplies:

1. an exact unforced Navier--Stokes endpoint evolution;
2. the complete compatible pressure/velocity jet;
3. a canonical global wake;
4. \(O(\operatorname{Re}^{-1})\) deformation and dissipation over a
   turnover; and
5. preservation of the transition's mean and reflection constraints.

It does not yet supply a smooth splice.  The active transition has to land
on the full jet (4.3), not merely on the value \(V\).  The load-bearing
endpoint theorem is now:

> Construct the three-phase transition with terminal data equal to the
> exact collar jet through order \(M_j\asymp j^2/\log j\), with a
> Gevrey-2 bound uniform in the high-Reynolds scaling, while carrying the
> C86 periodic wake into the next stage.

This is sharper than asking for unspecified "viscous endpoint correctors":
their target, size, recurrence, energy loss, and symmetry are all fixed by
the exact local Navier--Stokes flow.
