# Gavrilov bubbles cannot serve as pressure-multipole actuators

Date: 2026-07-29

## Status and verdict

This note tests whether the kinematic curl packets in
`2026-07-29-pressure-multipole-control-audit.md` can be replaced by
rotated, translated, and scaled copies of Gavrilov's compact steady Euler
bubble.

The answer is a sharp **no** for multipole actuation, together with a
useful positive endpoint observation.

* Every compact steady Euler bubble with compact pressure is already
  pressure-dark to **all** harmonic orders.
* Its integrated covariance is exactly a positive scalar multiple of the
  identity.
* Rotations, translations, spatial scales, and amplitude scales preserve
  complete pressure darkness.
* Therefore the orbit of one Gavrilov bubble has rank zero in the
  exterior-plus-interior pressure chart.  A disjoint collection of such
  bubbles cannot cancel even one nonzero pressure moment of another field.
* A differentiable path that remains, at every time, a disjoint sum of
  exact steady bubbles can solve unforced Euler only if the physical
  velocity is constant in time.

Thus an exact static Gavrilov wake is an excellent **already-dark storage
endpoint**, but it has no control directions for darkening an anisotropic
active packet.  The active packet must be dynamically transmuted into the
bubble wake; adding static bubbles does not perform that transmutation.

The primary seed used here is Gavrilov's nonzero smooth compactly
supported steady Euler flow:
[A. V. Gavrilov, *A steady Euler flow with compact support*](https://arxiv.org/abs/1810.08020).

---

## 1. A universal complete-darkness identity

Let

\[
W\in C_c^\infty(\mathbb R^3;\mathbb R^3),
\qquad
\operatorname{div}W=0,
\]

and suppose

\[
(W\cdot\nabla)W+\nabla P=0
\tag{1.1}
\]

with \(P\) constant outside a compact set.  Subtract the exterior constant,
so \(P\) is compactly supported.  Since

\[
(W\cdot\nabla)W_i
=\partial_j(W_iW_j),
\]

divergence of (1.1) gives

\[
\partial_i\partial_j(W_iW_j)=-\Delta P.
\tag{1.2}
\]

Let \(H\) be **any** harmonic polynomial.  Two integrations by parts give

\[
\boxed{
\int_{\mathbb R^3}
(W\otimes W):D^2H\,dx
=-\int_{\mathbb R^3}H\Delta P\,dx
=-\int_{\mathbb R^3}P\Delta H\,dx
=0.}
\tag{1.3}
\]

Equation (1.3) contains every exterior pressure multipole of every
degree.  It also follows directly from the stronger physical fact that
the normalized pressure \(P\) vanishes identically outside its compact
support.

If a bubble is placed in an annulus and the active core lies in a
component of the complement of the joint velocity-pressure support, the
bubble pressure is constant there.  Hence every nonconstant interior
pressure jet is zero as well.  In the regular/Kelvin chart from the
previous audit, one bubble contributes

\[
(\text{all exterior moments},\text{all interior jets})=(0,0).
\tag{1.4}
\]

This is true to infinite order, not merely to the truncation depth \(M\).

---

## 2. The covariance orbit is a singleton

The degree-two part of (1.3) already forces an unexpectedly rigid
identity.  Define

\[
C=\int_{\mathbb R^3}W\otimes W\,dx.
\]

Hessians of harmonic quadratic polynomials span the trace-free symmetric
matrices, so

\[
C=cI.
\tag{2.1}
\]

The scalar is strictly positive for a nonzero bubble:

\[
c=\frac13\|W\|_2^2>0.
\tag{2.2}
\]

There is also a direct steady-stress virial proof.  Equation (1.1) is

\[
\partial_j(W_iW_j+P\delta_{ij})=0.
\]

Multiplication by \(x_k\) and integration give

\[
\int W_iW_k\,dx
=-\delta_{ik}\int P\,dx.
\tag{2.3}
\]

Thus \(c=-\int P\).

For amplitude \(A\), spatial scale \(\delta>0\), rotation \(R\in SO(3)\),
and translation \(z\), set

\[
\begin{aligned}
W_{A,\delta,R,z}(x)
&=A R
W\!\left(\frac{R^T(x-z)}{\delta}\right),\\
P_{A,\delta,R,z}(x)
&=A^2
P\!\left(\frac{R^T(x-z)}{\delta}\right).
\end{aligned}
\tag{2.4}
\]

Every pair in (2.4) is again steady Euler, and

\[
\int W_{A,\delta,R,z}\otimes W_{A,\delta,R,z}\,dx
=A^2\delta^3RCR^T
=A^2\delta^3cI.
\tag{2.5}
\]

Therefore:

* the rotation orbit of the covariance is the singleton \(\{cI\}\);
* allowing amplitude and spatial scale enlarges it only to the isotropic
  ray;
* its trace-free projection has rank zero and has no interior in the
  five-dimensional quadrupole-control space.

Higher translated or rotated raw stress tensors may look different, but
their contraction with every harmonic Hessian remains zero by (1.3).
Consequently the full pressure-chart orbit, not only its degree-two
projection, is the zero vector.

This answers the proposed orbit question: the bubble supplies a positive
isotropic **baseline**, but no pressure-control directions around that
baseline.

---

## 3. Disjoint bubble sums remain rank zero

Let \((W_a,P_a)\) be finitely many transformed Gavrilov pairs whose joint
velocity-pressure supports are pairwise disjoint.  Then

\[
U=\sum_aW_a,\qquad \Pi=\sum_aP_a
\tag{3.1}
\]

is an exact compact steady Euler solution.  Pointwise disjointness gives

\[
U\otimes U=\sum_aW_a\otimes W_a,
\]

and (1.3) gives

\[
\int(U\otimes U):D^2H=0
\tag{3.2}
\]

for every harmonic polynomial \(H\).

Now let \(V\) be a separate active field with support disjoint from all
bubble velocities.  At the instantaneous Poisson-source level,

\[
{\mathfrak m}_H\!\left(V+\sum_aW_a\right)
={\mathfrak m}_H(V)
+\sum_a{\mathfrak m}_H(W_a)
={\mathfrak m}_H(V).
\tag{3.3}
\]

No choice of bubble rotations, centres, scales, signs, or amplitudes
changes (3.3).  A disjoint static bubble bath can store arbitrary energy,
but it cannot cancel a nonzero multipole of the active packet.

There is an additional exact-Euler compatibility issue.  If \(V\) has an
algebraic pressure tail, then its **joint** velocity-pressure support is
not disjoint from the bubble region.  The usual disjoint-superposition
argument does not even make \(V+\sum W_a\) an Euler solution.  Complete
disjointness is available only after the active field has itself become
pressure-localized.

---

## 4. Gaussian overlap cannot rescue the orbit

An uncut Gaussian core never has literally disjoint velocity support, so
one might try to use its exponentially small overlap with a remote bubble.
Let \(V\) be the core and \(B_R\) a bubble supported near radius \(R\).
For a harmonic quadratic \(H\),

\[
\begin{aligned}
{\mathfrak m}_H(V+B_R)
={}&{\mathfrak m}_H(V)
+2\int
\operatorname{sym}(V\otimes B_R):D^2H,
\end{aligned}
\tag{4.1}
\]

because \({\mathfrak m}_H(B_R)=0\).  If

\[
|V(x)|\le Ce^{-a|x|^2},
\]

then

\[
\left|
\int\operatorname{sym}(V\otimes B_R):D^2H
\right|
\le
Ce^{-a'R^2}\|B_R\|_2.
\tag{4.2}
\]

Cancelling an order-one core quadrupole by (4.2) requires

\[
\|B_R\|_2\gtrsim e^{a'R^2},
\qquad
\|B_R\|_2^2\gtrsim e^{2a'R^2}.
\tag{4.3}
\]

At the proposed radius \(R\asymp j\), this is an
\(e^{cj^2}\) energy cost and destroys the cascade ledger.  If the
Gaussian is first cut so that supports are exactly disjoint, the cross
term in (4.1) is exactly zero and cancellation becomes impossible rather
than expensive.

Putting the bubble in the order-one core avoids the small overlap, but
then the sum of two overlapping steady/nonsteady fields has order-one
cross interactions and is not an exact Euler superposition.  Finding a
special interacting compact family would be a new construction, not a
consequence of Gavrilov's disjoint bubbles.

---

## 5. Rigidity of time-dependent paths through steady bubbles

There is a general dynamic obstruction stronger than checking amplitude,
translation, rotation, and scale separately.

### Proposition 5.1: a path of instantaneous steady states is stationary

Let \(u(t)\) be a differentiable, decaying divergence-free field.  Suppose
that for every \(t\) there is a decaying pressure \(p_{\rm st}(t)\) such
that

\[
(u(t)\cdot\nabla)u(t)+\nabla p_{\rm st}(t)=0.
\tag{5.1}
\]

If \(u\) also solves unforced Euler with pressure \(p(t)\),

\[
\partial_tu+(u\cdot\nabla)u+\nabla p=0,
\tag{5.2}
\]

then subtraction gives

\[
\partial_tu+\nabla q=0,
\qquad q=p-p_{\rm st}.
\tag{5.3}
\]

Since \(\operatorname{div}\partial_tu=0\),

\[
\Delta q=0.
\]

Decay at infinity forces \(\nabla q=0\).  Therefore

\[
\boxed{\partial_tu=0.}
\tag{5.4}
\]

The hypothesis applies to every differentiable path of pairwise disjoint
transformed Gavrilov bubbles, because each instantaneous sum is steady.
Consequently the only unforced Euler paths in this finite-dimensional
manifold are physically constant fields.

Parameter motion that leaves the physical field unchanged is harmless:
for example, rotating an axisymmetric bubble within its stabilizer merely
changes its coordinates.  Apart from such redundancies, (5.4) forces:

\[
A'(t)=0,\qquad
\delta'(t)=0,\qquad
z'(t)=0,\qquad
R'(t)=0.
\tag{5.5}
\]

The pressure-moment path is therefore the constant zero path, and every
raw covariance and energy is constant as well.

This also shows why time-dependent amplitude interpolation in the
Gavrilov transition ledger requires a Reynolds stress.  For example,
\(u=A(t)W\) leaves the uncancelled residual

\[
A'(t)W,
\tag{5.6}
\]

which is divergence-free and cannot be absorbed into pressure unless it
vanishes.

### Galilean motion is not a localized exception

For a constant vector \(V\),

\[
u(t,x)=V+W(x-Vt),
\qquad
p(t,x)=P(x-Vt)
\tag{5.7}
\]

is an exact Galilean translate.  But the uniform background \(V\) is not
in \(L^2(\mathbb R^3)\) and is not a compact wake.  Removing \(V\) outside
a neighbourhood of the bubble introduces a transition layer and destroys
the exact identity.  Thus (5.7) does not provide a localized moving-bubble
control.

For Navier--Stokes, even a static bubble has the additional residual
\(-\nu\Delta W\).  The rigidity result is therefore an Euler best case,
not a viscous construction.

---

## 6. The corrected dynamic target

Rotated Gavrilov bubbles do not replace the curl-atom multipole controls.
The distinction is structural:

* the arbitrary curl atoms have a full finite pressure-moment chart but
  are not steady Euler;
* exact compact steady Euler bubbles have perfect pressure localization
  but zero pressure-moment chart.

There are only three plausible ways forward.

1. **Full transmutation.**  Drain the anisotropic active packet completely
   and land on a disjoint Gavrilov-bubble state.  At the endpoint, pressure
   localization is then exact to all orders without a separate multipole
   correction.
2. **Active Euler--Reynolds bridge.**  Use the full-rank curl/Gaussian
   moment controls during the transition and arrange for their Reynolds
   stress to vanish precisely when the state reaches the dark bubble
   manifold.
3. **New interacting compact family.**  Construct a non-steady or
   multi-bubble exact Euler family whose cross interactions supply
   nonzero pressure moments while retaining compact joint support.

The first route is the cleanest revised endpoint theorem.  It also states
the hard part honestly: a path within the steady-bubble orbit cannot do
the work.  The transition must leave that orbit, carry order-one
nonstationary stress, and return to it only at the endpoint.
