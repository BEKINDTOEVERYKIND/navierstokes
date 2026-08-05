# Multi-colour/core-buffer audit for the BAS transition

Status: internal leading-cell construction and obstruction, corrected after
independent audit.  This is not a Navier--Stokes singularity proof.  The
frozen-symbol high sum is pressure only before a transverse envelope is
introduced; the correction is recorded in Section 4.  The full formulas and
viscous-capture theorem are in
`affine_core_viscous_capture.md`.

## 1. A special flank value gives a two-colour strain dictionary

For radial tilt `l=0`, the Gavrilov-column BAS formula gives

```
H^2 = -(1/2) [3h/2+3+1/(2h)].
```

The descending-flank root with `H=-1` is

```
h_* = (-5-sqrt(22))/3 = -3.230138...
```

and lies strictly in the Leibovich--Stewartson unstable region `h<-2`.
For a small nonzero tilt, continuity permits a nearby flank point with
`H=-1` exactly.

Let `S` be any symmetric trace-free `3 x 3` matrix with eigenbasis
`(m_1,m_2,n)` and eigenvalues `(lambda_1,lambda_2,lambda_3)`, where
`lambda_3=-lambda_1-lambda_2`.  For colour `j=1,2`, choose

```
t_j x n = m_j,
e_rj = (n-m_j)/sqrt(2),
e_hj = -(n+m_j)/sqrt(2).
```

At `H=-1`, the mean growing polarization is common to both colours:

```
e_rj+H e_hj = e_rj-e_hj = sqrt(2)n.
```

For `q_j=Q_j e_rj`, `delta_j=Q_j/(2 Lambda_j)`, define

```
K_(j,s) = Lambda_j t_j+s q_j/2,
a_(j,s) = e_rj-s delta_j t_j-e_hj
        = sqrt(2)n-s delta_j t_j,             s in {+1,-1}.
```

Then `K_(j,s).a_(j,s)=0` exactly.  Within each colour, the positive-positive
sum is pressure, while the difference beat projects to

```
P_q [(a_+ . (-K_-))a_-+(a_- . K_+)a_+] = -2 Q_j e_hj.
```

The symmetric gradient atom of the written child is

```
Sym(e_hj tensor e_rj) = (1/2)(m_j tensor m_j-n tensor n).
```

Consequently

```
S = 2 lambda_1 Sym(e_h1 tensor e_r1)
  + 2 lambda_2 Sym(e_h2 tensor e_r2).
```

Relative carrier phase supplies either sign.  A single atom has eigenvalues
`(mu,0,-mu)` and determinant zero, so two colours are minimal for a generic
strain.

At order `Lambda`, all carriers have the same velocity direction `n`, and
all central wavevectors are perpendicular to `n`.  Their leading sum is an
exact shear `n f(x perpendicular n)`, for which `(u.grad)u=0`.  Thus every
order-`Lambda` self and cross interaction vanishes or is pressure.

## 2. Exact mixed-colour obstruction at order Q

The finite tilts needed to write a child produce mixed bands near
`Lambda(t_1+t_2)` and `Lambda(t_1-t_2)`.  In the eigenframe
`(m_1,m_2,n)=(e_1,e_2,e_3)`, take `t_1=e_2,t_2=-e_1` and equal `Lambda,Q`.
For positive-positive cross interactions, the `n` component is

```
(Q/sqrt(2))(s+t),       s,t in {+1,-1},
```

up to the common amplitude product.  Hence the same-sign pairs are nonzero
at order `Q`.  For positive-negative interactions the corresponding leading
component is

```
(Q/sqrt(2))(s-t),
```

so the opposite-sign pairs are also nonzero at order `Q`.  The four output
wavevectors are distinct.  Phases cannot cancel unique Fourier outputs.

If the carrier amplitudes are `x_(j,+),x_(j,-)`, then

```
|x_(1,+)x_(2,+)| |x_(1,-)x_(2,-)|
 = |x_(1,+)x_(1,-)| |x_(2,+)x_(2,-)|.
```

Both right-hand products must be nonzero to write both strain atoms, so at
least one same-sign parasitic is comparable to the geometric mean of the two
desired products.  The difference channels give the complementary
obstruction.  Therefore, within this four-primary-mode ansatz:

* two colours cancel every principal `O(Lambda)` interaction;
* no simultaneous finite-frequency two-colour cell supported only on these
  four primary modes cancels all `O(Q)` mixed bands while retaining both
  desired low beats.

The mixed bands have frequency about `sqrt(2)Lambda`.  If one nevertheless
uses simultaneous pulses, spectral damping requires roughly

```
nu Lambda^2 = vartheta g,       1/2 < vartheta < 1,
```

so the carrier can grow but a mixed band with comparable strain is damped.
Its forced amplitude is still comparable to the child before damping, so
sequential pulses or spatially disjoint active supports are much safer.

## 3. Nonlinear relative phase writes an exact affine core profile

A plane-wave child gives a cosine strain and recreates the sideband problem.
Instead, use a nonlinear slow relative phase.  For one colour, set

```
s = e_r . (x-x_0),
phi_+ = Lambda t.x+psi(s)/2,
phi_- = Lambda t.x-psi(s)/2,
delta(s)=psi'(s)/(2 Lambda),
a_+(s)=w-delta(s)t,
a_-(s)=w+delta(s)t,
w=e_r+H e_h.
```

The complex carrier waves `a_pm(s) exp(i phi_pm)` are exactly divergence
free: `a_pm.grad(phi_pm)=0` and `div(a_pm)=0`.

Their cross high sum has phase `2 Lambda t.x`.  The polarization-derivative
terms cancel, and the raw phase term is parallel to `t`, but its coefficient
depends on `s`.  It is therefore pressure only at frozen-symbol level.  After
subtracting the gradient, the exact surviving term can be written

```
P[ (psi' psi''/Lambda^2) e_r exp(i 2 Lambda t.x) ],
```

and has size `O(A^2 Q^3/Lambda^2)` for a `Q`-scale phase.

The low difference has phase `psi(s)` and, after pressure removes its
`e_r` component, equals

```
-i 2 H psi'(s) C exp(i psi(s)) e_h,
```

where `C=A_+ conjugate(A_-)` is the carrier product.  Choosing the relative
amplitude phase so `C=i c` gives the real low force

```
F_low(s) = 4 H c psi'(s) cos(psi(s)) e_h
         = 4 H c d/ds[sin(psi(s))] e_h.
```

Therefore any desired one-dimensional polynomial profile `f(s)e_h` can be
written on a core by setting

```
sin(psi(s)) = (1/(4 H c)) integral_0^s f(r) dr.
```

For an affine shear `f(s)=gamma s`, take

```
sin(psi(s)) = gamma s^2/(8 H c).
```

The primitive must stay strictly inside `(-1,1)` on the core.  A Gevrey
extension can flatten `psi` in a surrounding buffer.

More generally, finitely many such ridge shears span every divergence-free
polynomial jet.  Indeed a homogeneous divergence-free polynomial velocity
`U_m` has a homogeneous polynomial vector potential of degree `m+1`.
Scalar homogeneous polynomials are spanned by powers `(q.x)^(m+1)`; taking
curl turns each such term into

```
a (q.x)^m,       a.q=0,
```

which is exactly a ridge-shear profile handled by the phase law above.  Thus
constant or quadratic strain jets (affine or cubic velocity jets) admit a
finite sequential colour dictionary.  The two-colour formula in Section 1
is the sharp special case for a constant symmetric strain.

The individual self harmonics at phases `2phi_pm` do not vanish when
`psi'' != 0`.  Their leading Leray-projected radial term is

```
(A_pm^2/(4 Lambda^2))
  [psi' psi'' -/+ i psi'''] e_r exp(i 2 phi_pm).
```

Thus a phase varying uniformly on scale `Q` has a relative
`(Q/Lambda)^2` defect, two orders below the desired beat.  The `psi'''` term
must not be omitted.  This is the high-frequency corrector/wake to estimate;
it is not an order-`Lambda` principal sum.

Each modulated wave also has a zero-phase self term proportional to
`+/- |A_pm|^2 psi'' t/Lambda`.  The signs are opposite, so these terms cancel
when the two carrier magnitudes are equal.  Reciprocal child feedback does
*not* preserve that equality.  In the captured orbit its integrated size is
`O(Q/Lambda)` relative to the desired beat.  An exactly quadratic chirp makes
`psi''` constant, turning this term into a removable Galilean acceleration
with zero strain.

## 4. The envelope-dark equation is not full pressure darkness

For scalar envelopes `A(x),B(x)` multiplying the two carrier waves with
constant local polarizations

```
a=e_r-delta t+H e_h,
b=e_r+delta t+H e_h,
```

one explicit envelope-derivative contribution to their common high sum
vanishes when

```
A (a.grad B)+B (b.grad A)=0.                     (envelope-dark equation)
```

For `A=B=f`, this is

```
(e_r+H e_h).grad f=0.
```

At `H=-1` the direction in this equation is the common direction `n`.
However, this scalar equation is **not** sufficient for the entire high sum
to be a gradient.  The frozen phase contribution is parallel to the central
high frequency only when its scalar coefficient is constant in transverse
directions.

An exact counterexample is obtained by choosing a unit vector `m` orthogonal
to both `n` and `t` and setting `A=B=f(m dot x)`.  Then
`a dot grad(B)=b dot grad(A)=0`, so the displayed envelope-dark equation
holds.  Nevertheless the remaining high interaction is

```
N_high = -i (Q^2/Lambda) f^2 t exp(2 i Lambda t dot x),
```

and

```
curl(N_high)
 = -i (Q^2/Lambda) (f^2)' (m cross t)
   exp(2 i Lambda t dot x),
```

which is nonzero for every nonconstant lateral cutoff.  For a lateral scale
`R`, the Leray residual has leading size
`A_0^2 Q^2/(Lambda^2 R)`.  Thus lateral localization creates a small
high-frequency correction/wake; it cannot be declared exactly pressure.
Exact darkness in this symmetric example would force `f^2` to vary only
along `t`, which precludes compact lateral support.

A longitudinal cutoff also violates the scalar envelope equation.  Placing
it in remote endcaps of a long buffer tube still gives the previously
identified `O((Q L_b)^(-1))` contribution relative to the low beat, before
viscous gain.  The corrected core-buffer lemma must retain both that endcap
wake and the lateral Leray residual.
The full necessary-and-sufficient pressure criterion, classification of
real scalar envelopes, and characteristic-crossing replacement are in
`2026-08-02-characteristic-envelope-pressure-ledger.md`.

Simultaneous colours have the unavoidable mixed-band obstruction above.
A conservative protocol is sequential:

1. activate colour 1 and write its affine shear;
2. transport/damp its `Lambda` carrier and endcap wake;
3. activate colour 2 in the now constant-affine core;
4. include the first shear as a constant perturbation of the second frozen
   BAS matrix.

Constant affine strain changes carrier covectors but does not create a
discrete sideband ladder.  Four sequential pulses can be used if one wants a
pure symmetric velocity gradient rather than only its symmetric part: add
the transpose shear `e_rj tensor e_hj` for each `j`.  Then

```
(e_hj tensor e_rj)+(e_rj tensor e_hj)
 = m_j tensor m_j-n tensor n,
```

and the endpoint core is exactly `U(x)=S(x-x_0)`.  Since `S` is symmetric
and trace-free,

```
(U.grad)U=S^2(x-x_0)=grad[ (1/2)(x-x_0).S^2(x-x_0) ],
```

so it is an exact steady Euler flow on the core.

## 5. Localized core and a Gavrilov/packed-bubble seed

Let `chi` be a Gevrey cutoff equal to one on `B_R` and zero outside
`B_(2R)`.  The raw field `chi S(x-x_0)` fails to be divergence-free only in
the annulus.  Because `trace(S)=0`, its divergence has zero integral.  A
Bogovskii correction supported in the annulus gives a compact divergence-free
field that equals `S(x-x_0)` exactly on the core.  The same construction
works for any divergence-free quadratic Taylor polynomial.  All localization
and Euler defects are in the buffer.

Place a much smaller rescaled Gavrilov bubble inside the affine core.  After
a Galilean shift, the core strain supplies the load-bearing first jet while
the bubble supplies localized vorticity and the desired next-stage topology.
Alternatively cover the active region of a target bubble by finitely many
disjoint cores, match its affine or quadratic Taylor jet in each, and put the
Bogovskii/taper wakes in the gaps.  The Taylor error is `O(eta^2)` in velocity
on cells of radius `eta` times the bubble scale; the number of cells is
`O(eta^(-3))`.

There is an important orientation caveat.  A single frozen Gavrilov column
has one resonant central direction; it does not automatically provide the two
independently oriented BAS colours used for an arbitrary `S`.  Those require
an oriented/packed parent cluster or separate spatial-temporal blocks.
Fortunately the load-bearing strain of a next Gavrilov column is itself only
one atom.  In its resonant frame its velocity-gradient matrix has the form

```
B = [[0,-g,D],
     [g, 0,0],
     [-sigma H,0,0]],
```

so `Sym(B)` has only an `e_r--e_h` entry and eigenvalues `(mu,0,-mu)`.
Thus one nonlinear-phase colour can match the strain jet needed to amplify a
preseeded next Gavrilov block.  Two colours are needed only when the endpoint
must realize a generic symmetric trace-free matrix, and four transpose-paired
colours only when its full gradient must equal a generic symmetric matrix.

## 6. What remains to prove

The exact algebra supports a minimal two-colour strain dictionary and a
nonlinear-phase affine writer at frozen-envelope level.  It also excludes
simultaneous cancellation at the first nontrivial `O(Q)` order for the
specific four-primary-mode support in Section 2.  It does not exclude a
richer finite support with shared output modes and correctors.  Sequential
activation is therefore the safest current architecture, not a forced
theorem.  A viable PDE lemma must quantify:

* perturbation of the second BAS exponent by the first affine shear;
* viscous removal of each carrier before the next pulse;
* the `O((Q/Lambda)^2)` self harmonics;
* the lateral `O(Q^2/(Lambda^2 R))` Leray residual of the common high sum;
* remote endcap wake and Bogovskii-correction norms;
* preservation of a preseeded localized Gavrilov bubble inside the core.

These are analytic transition estimates, not a GPU search problem.
