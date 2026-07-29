# Pressure-multipole control for an analytic Gaussian return stage

Date: 2026-07-29

## Status and verdict

This note attacks the pressure-tail gate in the near-identity analytic
Gaussian proposal.

There is a genuine finite-dimensional positive result:

> For every fixed \(M\), every collection of pressure multipoles of
> degrees \(2,\ldots,M\) can be cancelled by a smooth, divergence-free
> velocity supported in a disjoint annulus.  The construction uses
> \(O(M^2)\) localized curl packets, has a strictly positive isotropic
> baseline, and can preserve the prescribed velocity exactly on the
> inner core.

At the stress level there is an explicit spherical right inverse with
quantitative bounds.  For a Gaussian core, putting the control shell at
normalized radius \(R\asymp j\) controls all moments through
\(M_{\rm p}=\mu j^2\) with order-one stress, for sufficiently small fixed
\(\mu>0\).  Cancelling to this deeper order, rather than merely
\(j^2/\log j\), makes the pressure tail \(e^{-cj^2}\) at a fixed larger
radius.  This depth is still compatible with \(K_j=j^A,\ A>4\).

There are also three sharp limitations.

1. Any fixed finite list of scalar sideband controls fails as
   \(M\to\infty\): the moment space has dimension \(M^2+2M-3\).
2. Sidebands confined to the preferred single-fast-direction class
   cannot even cancel the degree-two pressure quadrupole when \(K\gg M\).
   An annular wake with a longitudinal component, or at least one
   additional fast direction, is necessary.
3. A disjoint quadrupole corrector is generally nonperturbative in energy,
   and a nonzero analytic corrector cannot be compactly supported or leave
   an open inner core *exactly* unchanged.

Thus pressure multipoles are not a kinematic obstruction.  Their dynamic
generation as a recurrent, analytic or Gevrey-tame outgoing wake remains
an unsolved transition theorem.

---

## 1. The exact pressure moments

Let \(u\) be a sufficiently decaying real divergence-free field on
\(\mathbb R^3\), set

\[
T=u\otimes u,
\qquad
f=\partial_i\partial_jT_{ij},
\]

and solve

\[
-\Delta p=f.
\tag{1.1}
\]

For every polynomial \(H\), two integrations by parts give

\[
\int_{\mathbb R^3}Hf
=\int_{\mathbb R^3}T:D^2H.
\tag{1.2}
\]

Only harmonic polynomials enter the Newtonian far-field expansion.
Let \({\cal H}_m\) be the real homogeneous harmonic polynomials of degree
\(m\).  The pressure multipole associated with \(H\in{\cal H}_m\) is

\[
{\mathfrak m}_H(u)
:=\int_{\mathbb R^3}(u\otimes u):D^2H.
\tag{1.3}
\]

The degrees zero and one vanish automatically.  The number of independent
conditions through degree \(M\) is

\[
\begin{aligned}
D_M
&=\sum_{m=2}^{M}\dim{\cal H}_m
=\sum_{m=2}^{M}(2m+1)\\
&=(M+1)^2-4
=M^2+2M-3.
\end{aligned}
\tag{1.4}
\]

If the stress is supported in \(B_R\) and (1.3) vanishes for
\(2\le m\le M\), Taylor expansion of the Newton kernel gives, for
\(|x|\ge qR\) with fixed sufficiently large \(q>1\),

\[
|\partial^\alpha p(x)|
\le C_{\alpha,q}\|u\|_2^2
|x|^{-3-|\alpha|}
\left(\frac Cq\right)^{M-1}.
\tag{1.5}
\]

The precise constant is unimportant here; the important feature is the
geometric factor \(q^{-M}\).

For \(m=2\), \(D^2H\) ranges over all constant trace-free symmetric
matrices.  Hence the degree-two multipole vanishes exactly when

\[
C(u):=\int u\otimes u=cI
\tag{1.6}
\]

for some \(c\ge0\).

---

## 2. Three obstructions that any control proposal must respect

### 2.1 A fixed finite scalar control list cannot work to all orders

A family with \(N\) independently adjustable real scalar weights has
linearized moment rank at most \(N\).  Generic cancellation through degree
\(M\) therefore requires

\[
N\ge D_M=M^2+2M-3.
\tag{2.1}
\]

A fixed number of carrier *families* can still work if their envelopes
contain \(O(M^2)\) independent coefficients.  What is impossible is a
stage-independent fixed list of scalar sideband amplitudes.

### 2.2 The preferred one-carrier class fails already at degree two

Take a unit fast direction \(w\) and a mean-zero one-phase profile \(W\)
whose slow Fourier band is at most \(M\).  Exact incompressibility gives

\[
K\,w\cdot W
=-\partial_\vartheta^{-1}\operatorname{div}_xW.
\tag{2.2}
\]

Consequently

\[
\|w\cdot W\|_2
\le C\frac{M}{K}\|W\|_2.
\tag{2.3}
\]

Let

\[
C_W=\int W\otimes W.
\]

If \(C_W=cI\), then

\[
\frac13\operatorname{tr}C_W
=w^TC_Ww
\le C^2\frac{M^2}{K^2}\operatorname{tr}C_W.
\tag{2.4}
\]

For \(K>\sqrt3\,CM\), (2.4) forces \(W=0\).  Thus no number of
zero-target sidebands that remain inside the same transverse
single-carrier class can make a nonzero bath quadrupole-dark.  The
pressure corrector must contain an order-one component in the missing
direction.  It must be a separate low/annular wake or use another fast
direction.

This does not damage the all-chain \(K\)-null identity in the *active
core* if the extra directions are confined to the outgoing annulus.  It
does rule out treating pressure cancellation as a harmless extension of
the existing three-sideband chart.

### 2.3 A disjoint quadrupole corrector costs order-one energy

Suppose \(u\) and an annular correction \(v\) have disjoint supports.
Then

\[
C(u+v)=C(u)+C(v),
\qquad C(v)\ge0.
\]

Let the eigenvalues of \(C(u)\) be
\(\lambda_1\ge\lambda_2\ge\lambda_3\).  Quadrupole cancellation requires

\[
C(v)=cI-C(u),\qquad c\ge\lambda_1.
\]

The smallest possible correction energy is therefore

\[
\boxed{
\min\operatorname{tr}C(v)
=3\lambda_1-\operatorname{tr}C(u).}
\tag{2.5}
\]

For a nearly rank-one packet this is approximately twice the original
energy.  Terminal-only cancellation is possible, but it is not a small
Gaussian seam correction.  The added energy has to be included in the
wake and global energy ledger.

---

## 3. An explicit positive shell-stress right inverse

The following lemma is the central positive result.

### Proposition 3.1: positive shell stress prescribes every finite jet

Fix \(M\ge2\), \(R>0\), and desired numbers

\[
\beta_{m,k},
\qquad
2\le m\le M,\quad 1\le k\le2m+1.
\]

For each \(m\), choose a basis \(H_{m,k}\) of \({\cal H}_m\) orthonormal
for the rotation-invariant Hessian inner product

\[
\langle H,G\rangle_m
=\int_{\mathbb S^2}D^2H(\omega):D^2G(\omega)\,d\omega.
\tag{3.1}
\]

This is positive definite for \(m\ge2\).  Put

\[
E_{m,k}(\omega)=D^2H_{m,k}(\omega).
\]

Every component of \(E_{m,k}\) is a spherical harmonic of degree \(m-2\).
Different \(m\)'s are therefore orthogonal.  Define the trace-free tensor
field

\[
Q(\omega)
=\sum_{m=2}^{M}\sum_{k=1}^{2m+1}
R^{-(m-2)}\beta_{m,k}E_{m,k}(\omega).
\tag{3.2}
\]

Choose

\[
\lambda\ge2\|Q\|_{L^\infty(\mathbb S^2)}
\]

and set

\[
T(\omega)=\lambda I+Q(\omega).
\tag{3.3}
\]

Then \(T(\omega)\) is positive definite and

\[
\boxed{
\int_{\mathbb S^2}
T(\omega):D^2H_{m,k}(R\omega)\,d\omega
=\beta_{m,k}.}
\tag{3.4}
\]

Indeed,

\[
D^2H_{m,k}(R\omega)=R^{m-2}E_{m,k}(\omega),
\]

the \(Q\) term gives (3.4) by orthonormality, and the isotropic term gives
zero because

\[
I:D^2H_{m,k}=\Delta H_{m,k}=0.
\]

The addition theorem gives the quantitative bound

\[
\|Q\|_\infty
\le
C\sum_{m=2}^{M}
R^{-(m-2)}\sqrt m\,
|\beta_m|_{\ell^2}.
\tag{3.5}
\]

Thus the isotropic positive margin costs only the right-hand side of
(3.5).  It changes energy, but none of the pressure moments.

### Finite positive discretization

The integrands in (3.4) are restrictions of polynomials of degree at most
\(2M-4\).  A completely elementary positive cubature is obtained by:

* a trapezoidal rule in the azimuthal angle with more than \(2M-4\)
  nodes; and
* Gauss--Legendre quadrature in \(z=\cos\theta\) with enough nodes to
  integrate degree \(2M-4\).

It uses \(N\le CM^2\) positive nodes and is exact for (3.4).
At each node the positive matrix \(T(\omega_a)\) can be decomposed into
positive rank-one tensors.

For the divergence-free realization below, it is convenient instead to
use matrices

\[
I-b\otimes b.
\tag{3.6}
\]

Nine fixed unit directions

\[
\begin{gathered}
e_1,e_2,e_3,\\
\frac{e_1\pm e_2}{\sqrt2},\qquad
\frac{e_1\pm e_3}{\sqrt2},\qquad
\frac{e_2\pm e_3}{\sqrt2}
\end{gathered}
\tag{3.7}
\]

span all symmetric matrices through (3.6), and their equal positive sum
is isotropic.  Consequently every matrix sufficiently close to
\(\lambda I\) is a positive combination of (3.6).  Increasing \(\lambda\)
in (3.3) provides the required uniform positivity margin.

This gives a finite positive stress measure with \(O(M^2)\) atoms, exact
for every pressure moment through degree \(M\).

### Exterior moments are not the whole pressure-matching problem

An annular correction that cancels the exterior multipoles also creates a
harmonic pressure field inside the annulus.  Let

\[
H_{m,k}(x)=|x|^mY_{m,k}(x/|x|)
\]

and let its Kelvin transform be

\[
K_{m,k}(x)
=|x|^{-2m-1}H_{m,k}(x)
=|x|^{-m-1}Y_{m,k}(x/|x|).
\tag{3.8}
\]

The degree-\(m\) Taylor coefficient of the pressure produced in the inner
core is, up to a fixed normalization,

\[
\int T(y):D^2K_{m,k}(y)\,dy.
\tag{3.9}
\]

The \(m=0\) coefficient is a pressure gauge.  The \(m=1\) coefficients
control the pressure gradient and should also be set to zero; they add
only three more linear conditions.  The paired exterior/interior channels
below concern \(m\ge2\).

Therefore exterior cancellation alone does not preserve the active core's
pressure jet.  At radius \(R\),

\[
D^2H_{m,k}\sim R^{m-2},
\qquad
D^2K_{m,k}\sim R^{-m-3}.
\tag{3.10}
\]

The positive-atom argument extends to the combined exterior/interior
system.  Indeed, suppose a linear combination \(\Psi\) of the regular
harmonics \(H_{m,k}\) and singular harmonics \(K_{m,k}\) satisfies

\[
b^TD^2\Psi(y)b=0
\]

for every direction \(b\) and every \(y\) in an open annulus.  Then
\(D^2\Psi=0\) there, so \(\Psi\) is affine.  The distinct positive and
negative homogeneities in the chosen list force every coefficient to
vanish.  Hence localized stress atoms span the dual of the combined
finite jet space.  Isotropic orientation sums again give a strictly
positive null vector.

For quantitative separation, use two concentric control shells with
radii \(R\) and \(\sigma R\), where \(\sigma>1\) is fixed.  In each
angular channel the two radial homogeneities give the nonzero normalized
determinant

\[
\sigma^{-(m+3)}-\sigma^{m-2}\ne0.
\tag{3.11}
\]

After the angular Hessian maps are normalized, this is a two-shell radial
Vandermonde.  It permits one to prescribe the exterior multipole while
setting the corresponding interior pressure jet to zero.  The number of
coordinates and packets only doubles and remains \(O(M^2)\).

This two-sided cancellation is required if “preserving the core endpoint”
includes its pressure gradient, rather than only its velocity.

---

## 4. Exact smooth divergence-free annular realization

The stress atoms above are not merely formal.  Let
\(\psi\in C_c^\infty(B_1)\) be radial and define

\[
V_b(z)=\nabla\psi(z)\times b.
\tag{4.1}
\]

Then \(\operatorname{div}V_b=0\), and radial symmetry gives

\[
\int_{\mathbb R^3}V_b\otimes V_b\,dz
=c_\psi\left(|b|^2I-b\otimes b\right).
\tag{4.2}
\]

At a shell point \(y\), use the energy-normalized packet

\[
v_{y,b,\varepsilon}(x)
=\varepsilon^{-3/2}
V_b\!\left(\frac{x-y}{\varepsilon}\right).
\tag{4.3}
\]

For \(H\in{\cal H}_m\),

\[
\begin{aligned}
{\mathfrak m}_H(v_{y,b,\varepsilon})
&=
\int V_b(z)\otimes V_b(z):
D^2H(y+\varepsilon z)\,dz\\
&\longrightarrow
-c_\psi\,b^TD^2H(y)b
\qquad(\varepsilon\downarrow0),
\end{aligned}
\tag{4.4}
\]

because \(\Delta H=0\).

Choose the cubature nodes from Section 3, split the nine orientation
packets near each node by arbitrarily small distinct displacements, and
make all supports disjoint.  In the zero-radius limit their moment matrix
has:

1. full row rank \(D_M\); and
2. a strictly positive null vector coming from the isotropic orientation
   sums.

Both properties persist for sufficiently small positive
\(\varepsilon\).  The weights can then be perturbed to solve the finite
linear moment system **exactly**, while remaining positive.  Since the
supports are disjoint, square-rooting the weights and summing the packets
introduces no cross stresses.

This proves:

### Proposition 4.1: fixed-order annular velocity corrector

Let \(u\in C_c^\infty(\mathbb R^3;\mathbb R^3)\) be divergence-free and
supported strictly inside the inner boundary of a prescribed annulus.
For every finite \(M\) there is a real

\[
v\in C_c^\infty(\mathbb R^3;\mathbb R^3),
\qquad \operatorname{div}v=0,
\]

supported in that annulus, such that

\[
{\mathfrak m}_H(u+v)=0
\qquad
\text{for all }H\in{\cal H}_m,\quad2\le m\le M.
\tag{4.5}
\]

Using the two-shell version, \(v\) may simultaneously be required to have
zero interior harmonic pressure jet through any prescribed finite order.
The field \(v\) leaves \(u\) exactly unchanged on the inner core.  It does
not leave the *global* endpoint unchanged: the annular field is a
nontrivial outgoing wake.

For a compact Gevrey-\(s\) bump, the scaled packets obey

\[
\|\partial^qv_{y,b,\varepsilon}\|_\infty
\le
C^{q+1}(q!)^s\varepsilon^{-q-3/2}.
\tag{4.6}
\]

Using the explicit product cubature, all atoms can be separated with
\(\varepsilon\ge R M^{-C_0}\) for some fixed, nonoptimized \(C_0\).
Thus fixed-order derivative costs are polynomial in \(M/R\).  A
quasi-uniform positive spherical design improves the geometric exponent,
but is not needed for the existence statement.

### Gaussian overlap

An uncut Gaussian core is not literally disjoint from the annulus.  At
radius \(R\), its cross stress with (4.3) is \(O(e^{-cR^2})\).  The moment
map near the strictly positive solution above has a surjective Jacobian.
The finite-dimensional implicit-function theorem therefore absorbs these
cross terms and still gives exact moment cancellation for sufficiently
large \(R\).

Replacing \(\psi\) by Gaussian vector potentials gives real-analytic,
divergence-free packets.  Their mutual overlaps and their tails on the
inner core are exponentially small, and the same finite-dimensional
argument gives exact cancellation of the finite list of moments.  But an
analytic correction cannot vanish on an open core unless it vanishes
identically.  Therefore:

* compact Gevrey corrections preserve inner endpoint data exactly;
* analytic Gaussian corrections preserve it only to exponentially small
  accuracy.

This distinction cannot be removed by a better moment calculation.

### Time-dependent moment paths

The shell right inverse is linear in the desired moment vector.  If
\(\beta(t)\) is analytic in time and remains in a bounded set, choose one
fixed isotropic margin \(\lambda\) larger than the supremum in (3.5).
The nine-direction coefficients at the cubature nodes can then be chosen
analytic in \(t\) and uniformly positive.  After the finite-size
perturbation, the same conclusion follows from the parameter-dependent
implicit-function theorem.  Thus the **kinematic** moment corrector can be
made analytic in time with the same quantitative bounds.

This time-dependent version is necessary, not cosmetic.  A Gaussian
velocity with an uncancelled quadrupole has an algebraic pressure gradient.
The Euler acceleration therefore acquires an algebraic far field
immediately, even if the velocity at one instant is Gaussian.  A
terminal-only correction does not justify a Gaussian packet class during
the preceding stage.  One must either:

1. cancel the required pressure moments throughout the active interval;
   or
2. evolve the pressure-generated algebraic tail as part of the outgoing
   wake.

The linear right inverse solves the moment bookkeeping in option 1, but it
does not prove that its time-dependent amplitudes satisfy the fluid
equations without an order-one force.

There is one further compatibility condition.  If the same constant
trace-free affine strain \(S\) acts on the core and the corrector, making
the total covariance isotropic gives

\[
S:C(u+v)=0
\]

and removes the total affine work.  Simultaneous pressure darkness and
amplification therefore require the annular control to lie outside the
region of nearly constant active strain (or require a genuinely
non-affine parent).  Spatial separation is doing essential work here, not
merely simplifying the moment algebra.

---

## 5. Quantitative Gaussian scaling

Let the normalized active field satisfy

\[
|u(y)|\le C e^{-a|y|^2}.
\tag{5.1}
\]

For the Hessian-orthonormal basis in Section 3, the addition theorem and a
radial Gaussian moment estimate give

\[
|\beta_m|_{\ell^2}
\le
C^m\sqrt m\,(m!)^{1/2},
\qquad
\beta_{m,k}=-{\mathfrak m}_{H_{m,k}}(u).
\tag{5.2}
\]

Place the stress control at

\[
R_j=\kappa j.
\tag{5.3}
\]

Combining (3.5) and (5.2), the degree-\(m\) contribution is bounded
schematically by

\[
\operatorname{poly}(m)
\left(\frac{C\sqrt m}{\kappa j}\right)^m.
\tag{5.4}
\]

Choose

\[
M_{{\rm p},j}=\lfloor\mu j^2\rfloor
\tag{5.5}
\]

with \(\mu>0\) sufficiently small relative to \(a\kappa^2\).  Then the
sum in (3.5) is uniformly bounded.  Hence:

* the annular stress and its energy remain order one in normalized
  variables;
* the number of scalar moment coordinates is
  \(D_{M_{{\rm p},j}}=O(j^4)\);
* the high-degree control coefficients decay geometrically, so the
  stress-level spherical field \(Q\) has a uniform analytic norm.

Now cut the Gaussian velocity at a larger radius \(Lj\), and repair the
cutoff divergence by a compact annular right inverse.  Its norm is
polynomially controlled and its source is already Gaussian-small.  Include
this seam in the exact finite moment system.  If \(L\) is sufficiently
large (and then \(\mu\) is chosen sufficiently small), the seam is
\(e^{-cL^2j^2}\), including its normalized moments relative to the
control shell.  After cancellation through
\(M_{{\rm p},j}\), cut the pressure only at a further fixed multiple of
\(Lj\).  Equation (1.5) gives

\[
|\partial^\alpha p_{\rm tail}|
\le
\operatorname{poly}(j)e^{-cM_{{\rm p},j}}
\le e^{-c'j^2}.
\tag{5.6}
\]

The depth \(j^2/\log j\) used as a minimal WKB truncation elsewhere is not
enough at a fixed radial ratio: it gives only
\(e^{-cj^2/\log j}\).  The pressure moments should be cancelled to the
deeper order \(M_{{\rm p},j}\asymp j^2\).

This deeper order is compatible with the polynomial carrier.  If the
coupled coefficient hierarchy is Gevrey two,

\[
C^M(M!)^2K_j^{-M},
\qquad K_j=j^A,
\]

then at \(M=\mu j^2\)

\[
\log\!\left(C^M(M!)^2K_j^{-M}\right)
=-(A-4+o(1))\mu j^2\log j.
\tag{5.7}
\]

Thus \(A>4\) still has more than enough margin.  Polynomial deterioration
from packet separation can be absorbed by increasing \(A\).

---

## 6. What is and is not closed

The pressure-tail problem is kinematically finite-dimensional and
surjective.  It does not require a miraculous pressure-localized Gaussian
subclass.  A concrete viable architecture is:

1. perform the short analytic single-carrier amplification in the core;
2. retain its anisotropic stress during the energy-transfer interval;
3. evolve an order-one annular pressure-control wake containing the
   missing carrier direction, or explicitly retain the pressure-generated
   algebraic tail during the interval;
4. cancel harmonic pressure moments through
   \(M_{{\rm p},j}\asymp j^2\), while using a second control radius to
   cancel the pressure jet fed back into the core;
5. cut velocity and then pressure at two successively larger fixed
   multiples of \(j\ell_j\).

An isotropic common-helicity, common-shell Beltrami bath is a natural
pressure-dark positive baseline for step 3: before localization it is a
stationary Euler field, and small modulations can encode the trace-free
stress \(Q\).  This is preferable to trying to isotropize the active
single-carrier bath itself.

What remains unproved is decisive:

* the annular corrector must be generated by the Euler/Navier--Stokes
  dynamics, not inserted by an order-one terminal force;
* its extra carrier directions must not feed an \(O(K)\) interaction back
  into the active core;
* the order-one correction energy and subsequent wake interactions must
  fit the global recurrence;
* an analytic realization cannot preserve the inner core exactly, while a
  compact realization is not analytic;
* the endpoint map must retain a uniformly tame right inverse for
  \(O(j^4)\) moment coordinates together with the five active strain
  coordinates.

The result is therefore a real reduction, not a solution: pressure
multipoles admit a controlled finite-dimensional cure, but that cure
forces a nonperturbative, multi-directional outgoing wake.
