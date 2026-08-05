# Exact audit of the symmetric `l=0` BAS doublet

Status: internal derivation, not a Navier--Stokes regularity or blow-up proof.

## 1. Geometry and frozen BAS block

Use an orthonormal frame `(e_r,e_t,e_h)` and put

```
q   = Q e_r,
K_s = Lambda e_t + s q/2 = Lambda(e_t+s delta e_r),
delta = Q/(2 Lambda),                         s in {+1,-1},
p_s = e_r-s delta e_t,
S = |p_s|^2 = 1+delta^2.
```

The divergence-free carrier plane at `K_s` is spanned by `p_s,e_h`.
For the resonant frozen-column symbol, evenness in the radial tilt gives the
same `sigma,H` at `+delta` and `-delta`, and

```
L_s p_s = sigma H e_h,
L_s e_h = (sigma/H) p_s.
```

Thus `p_s+H e_h` is the growing polarization and `p_s-H e_h` the decaying
polarization.  On the descending Gavrilov flank `H<0`.

Write

```
u_hat(K_s) = x_s p_s + y_s e_h,
u_hat(q)   = z e_h,
```

with the coefficients at negative wavevectors fixed by reality.

## 2. Exact Euler convolutions involving the central modes

With the Fourier convention

```
N_k = -i P_k sum_{p+r=k} (u_hat(p).r) u_hat(r),
```

the central interactions are

```
N_q = -i Q (x_+ conjugate(y_-) + conjugate(x_-) y_+) e_h,
N_K+ = -i Q z x_- e_h,
N_K- = +i Q conjugate(z) x_+ e_h.
```

At the high sum `K_++K_-=2 Lambda e_t`, the unprojected vector equals

```
-2 Q delta x_+ x_- e_t
+Q (x_- y_+ - x_+ y_-) e_h.
```

The first term is pressure.  Hence the high sum is killed exactly iff
`x_- y_+ - x_+ y_-=0`.

There is an exact exchange-symmetric phase family

```
x_+ = exp(i theta) x,       y_+ = exp(i theta) y,
x_- = exp(i(theta+phi)) x,  y_- = exp(i(theta+phi)) y,
z   = i exp(-i phi) Z,
```

with real `x,y,Z`.  If one truncates to the central modes and writes the
frozen growth rate as `g=dP`, it gives

```
x' = (g/H)y,
y' = g H x + Q Z x,
Z' = -2 Q x y.                                    (central truncation)
```

The central triad energy

```
E_c = 2(S x^2+y^2)+Z^2
```

is conserved by the `QZ` terms, while the frozen column changes it by
`4g(S+H^2)xy/H`.

This already disproves the advertised fixed-polarization scalar gate:

```
(y-Hx)' = Q Z x
```

on `y=Hx`.  Once a child is present, the decaying BAS polarization is forced.
The frozen BAS matrix and the child-feedback nilpotent matrix have no common
nonzero invariant line.

## 3. The decisive omitted interactions: outward sidebands

Reality supplies both `+q` and `-q`.  Therefore the same child also forces

```
K_+ + q = Lambda e_t + 3q/2,
K_- - q = Lambda e_t - 3q/2,
```

with projected amplitudes

```
-i Q z x_+ e_h,
+i Q conjugate(z) x_- e_h,
```

respectively.  These are nonzero whenever the desired carrier radial
component and child are nonzero.  Each is a unique convolution at its output
wavevector, so a phase choice cannot cancel it.

Consequently the three-wavevector set `{q,K_+,K_-}` is not invariant.  The
four-scalar system above is only a central-cell truncation.

## 4. General real/helical child does not repair the endpoints

Let the `+q` child polarization be `c=A e_t+B e_h`, with the `-q` coefficient
its conjugate.  For the BAS polarizations `p_++H e_h,p_-+H e_h`, killing the
positive outward sideband requires

```
A(1-3 delta^2)=0,      A Lambda H+Q B=0.
```

Killing the negative outward sideband requires

```
A(1-3 delta^2)=0,      A Lambda H-Q B=0
```

after conjugation.  If `delta != 1/sqrt(3)`, even one endpoint forces
`A=B=0`.  At the exceptional value the two endpoints still force `AH=0`;
`H=0` removes the desired low beat and the BAS mechanism.  Complex/helical
phase does not change this algebra.  Dropping the `-q` coefficient violates
reality and still produces an infinite one-sided ladder.

## 5. Exact carrier ladder and finite-support no-go

Put

```
K_n = Lambda e_t + (n+1/2)q,
p_n = e_r-(n+1/2)(Q/Lambda)e_t,
u_hat(K_n)=a_n p_n+b_n e_h,                    n in Z.
```

The exact child--carrier shift part is

```
b_n'|child = -i Q z a_(n-1)+i Q conjugate(z) a_(n+1).
```

For `z=iZ` and the reflection-real subspace this becomes

```
a_n' = (g_n/H_n)b_n-kappa_n a_n,
b_n' = g_n H_n a_n+QZ(a_(n-1)+a_(n+1))-kappa_n b_n.   (ladder)
```

The reflection `n -> -n-1` is preserved when `g_n,H_n,kappa_n` are even in
the radial tilt; it cancels the induced `e_t` component of the `q` mode, but
does not cancel outward sidebands.

If a nonzero sequence has finite support and `N` is an extreme occupied
index, then `q+K_N` (or `-q+K_N`) forces the next sideband by `Q z a_N e_h`.
If `a_N=0` but `b_N!=0`, the BAS block immediately produces
`a_N'=(g_N/H_N)b_N!=0`, after which the sideband is forced.  Thus no nonzero
finite-support carrier ladder is invariant when `Z!=0`.

Attempting exact cancellation with the next neighbor requires

```
a_(n+2)=(z/conjugate(z)) a_n,
```

whose modulus is constant on each parity class.  It has no nonzero `ell^2`
solution.  There is no finite-energy exactly dark ladder.

The full Euler nonlinearity also generates low harmonics `(n-m)q` and a
high band near `2 Lambda`; the ladder is therefore a necessary subsystem,
not a complete Fourier closure.

## 6. Bloch form and the plausible repair

At leading frozen order take `g_n=g`, `H_n=H`, and
`kappa_n=kappa`.  The Bloch transform in the sideband index has matrix

```
M(theta) = [[-kappa, g/H],
            [gH+2QZ cos(theta), -kappa]],
```

and eigenvalues

```
lambda_pm(theta) = -kappa
  +/- sqrt(g^2+(2gQZ/H)cos(theta)).
```

For `g>0,Z>0,H<0`, a packet near `theta=0` has the intended quenching/write
sign, while `theta=pi` is enhanced and has the wrong sign.  When `P` and
therefore `g=dP` changes sign, these favorable and unfavorable arcs swap.
This rules out a stationary, stable sign-flip landing.

The remaining plausible object is not a finite gate but an infinite
Gevrey sideband packet whose Bloch support is localized near the favorable
arc.  Frozen multiplication by `cos(theta)` preserves such support.  The
`n^2 Q^2` part of viscosity becomes heat flow in `theta`; over gain time
`G/g`, leakage across an order-one angular gap is of scale

```
exp(-c (Lambda/Q)^2/G).
```

The leading `2 Lambda` output is a Wronskian-type term
`Q(a b_theta-a_theta b)e_h`; it vanishes initially when `b/a` is constant,
but child-induced theta dependence regenerates it.  A narrow Bloch packet,
Gevrey tail estimates, and deliberate viscous damping at `2 Lambda` are the
remaining possible controls.

There is a useful rigidity statement behind this warning.  First freeze the
polarization ratio across the ladder and let `A(phi)=sum_n a_n exp(i n phi)`.
At high output index `s=n+m`, the nonpressure coefficient is proportional to

```
S_s = sum_(n+m=s) (n-m)^2 a_n a_m.
```

Only `s=-1` has zero radial output and can be pressure.  Therefore a
pressure-pure envelope must satisfy

```
2[(A')^2-A A''] = C exp(-i phi).
```

Put `y=exp(i phi/2)A`.  This becomes `y'^2-y y''=C`.  Differentiation gives
`(y''/y)'=0`; antiperiodicity of `y` then gives

```
y=c_+ exp(i mu phi)+c_- exp(-i mu phi),   mu in Z+1/2.
```

Thus `A` contains at most two mirror modes whose indices sum to `-1`.  To
write the first child harmonic their separation must be one, so the unique
nontrivial finite-energy equal-polarization writer is the original central
doublet `{n=-1,n=0}`.  An infinite localized envelope does not evade the
high-sum problem without varying polarization and adding explicit
`2 Lambda` correctors.

A complementary physical-space reduction reaches the same boundary.  If a
general child shear is allowed and one insists on a single pressure-pure
envelope, the leading system has the form

```
a_t=(g/H)b,
b_t=(gH-s(phi,t))a,
s_t=2 partial_phi^2 Re(a conjugate(b)).
```

Pressure-purity forces `b=rho(t)a`.  Invariance then forces `s` to be
spatially constant wherever `a` is nonzero, whereas the last equation forces
`partial_phi^2 |a|^2` to be constant there.  A nonzero quadratic intensity
cannot be smooth and compactly supported with a flat boundary.  Hence a
single exact compact envelope is also excluded at this leading level.  The
remaining escape needs several colours and/or a core--buffer construction
whose taper defect is carried by a Gevrey wake and damped high-band
correctors.

One useful spectral-filter window is to set `kappa=nu Lambda^2=theta_0 g`
with `1/4<theta_0<1`: the `Lambda` carrier can still grow, while a `2Lambda`
error with comparable strain is damped because `4kappa>g`.  The child loses
only `nu Q^2/kappa=(Q/Lambda)^2` on that time scale.  This is a proposed
estimate to prove, not an established PDE construction.

For the Palasek scale ledger, write

```
g ~ N^beta,       Q ~ N^b,
Lambda=(theta_0 g/nu)^(1/2),
epsilon=Q/Lambda ~ N^(b-beta/2).
```

The already required inequality `2b<beta` gives `epsilon -> 0`.  During a
net logarithmic gain `G`, sideband viscosity spreads a Bloch envelope by
`O(epsilon sqrt(G))`.  If its initial angular width is `eta`, the formal
narrow-packet hierarchy

```
epsilon sqrt(G) << eta << 1/G
```

is nonempty whenever `epsilon G^(3/2) << 1`.  This holds for
`G=(log Lambda)^2`, because `epsilon` is a negative power of `N`.  Leakage
across the order-one gap to the enhanced Bloch arc is then formally

```
exp[-c/(epsilon^2 G)].
```

These inequalities only clear the *linear* spectral-filter ledger.  They do
not bound the forced `2 Lambda` corrector, curvature, localization pressure,
or the nonlinear endpoint map.  In particular, damping at `2 Lambda` does
not by itself make a large forced response small.

## 7. Conditional central sign-flip calculation (not an exact closure)

If, only for diagnosis, one appends a single parent `P` with energy `P^2`
and imposes exact reciprocal energy exchange, put `C=S+H^2` and `g=dP`.
Then necessarily

```
P' = -(2dC/H)xy.
```

The inviscid central truncation has

```
x^2+P^2/(2C) = P0^2/(2C),
Z-(QH/(dC))P = -(QH/(dC))P0.
```

Writing `r=Q/(d sqrt(C))`, its first turning point is

```
P*/P0 = (r^2-1)/(r^2+1),
Z*/P0 = 2r(-H/sqrt(C))/(1+r^2),
x*/P0 = sqrt(2)r/(sqrt(C)(1+r^2)),
y*=0.
```

This is a transient sign flip, not a terminal state: the orbit reverses and
unwrites the child.  Invariants forbid a carrier-free endpoint with nonzero
child.  Moreover, after the parent changes sign the complementary BAS
polarization becomes growing.  A repair must stop parent--carrier overlap
before/near the sign change and use viscosity or spatial exit to damp the
remaining high-frequency packet.

Nor is the first omitted order automatically perturbative in the intended
cascade.  If the target child/parent ratio is

```
rho=N^(-(b-1)(alpha-beta)),
```

then the formal stable-polarization/sideband leakage is `O(rho^2)`, but a
nonpressure high-frequency residual carries a derivative leverage
`Lambda/Q`.  At the coarser dimensional level its ratio to the parent
amplification contains

```
N^((b-1)[1-2(alpha-beta)]).
```

The feasible window `beta>2b`, `alpha<=5/2`, `b>1` implies
`alpha-beta<1/2`, so this power grows.  Any sign-flip normal form would need
explicit cancellation beyond first omitted order, through the `P=0`
spectral-gap crossing.

Finally, all wavevectors in the exact `l=0` ladder lie in the rank-two
lattice generated by `q` and `Lambda e_t`.  Its nonlinear completion is a
2D3C system, which is globally regular.  Such a block can at most be one
local transition colour; a recurrent three-dimensional construction must
combine differently oriented, spatially controlled blocks without losing
the pressure and endpoint estimates.  The corresponding unlocalized
fixed-helical-phase ladder also sits in a globally regular invariant helical
subspace; compact localization breaks that global symmetry, so this is a
boundary on the module, not a no-go for every localized Gavrilov transition.

The strict conclusions are therefore:

1. the central high sum really is pressure under exact exchange symmetry;
2. the fixed-polarization three-scalar gate is not exact;
3. no real phase/helical choice or finite sideband ladder closes the system;
4. equal-polarization pressure-purity uniquely selects the central doublet,
   which is precisely the set that reality forces out of;
5. an approximate multi-colour Bloch/Gevrey packet with explicit damped
   high-band correctors is the only surviving version of this angle.

## Primary-source context

- S. Palasek, *Finite-time blow-up in an elementary model of the 3D
  Navier--Stokes equations*: https://arxiv.org/abs/2605.13827
- A. V. Gavrilov, *A steady Euler flow with compact support*:
  https://arxiv.org/abs/1810.08020
- N. Kishimoto and T. Yoneda, *Characterization of three-dimensional Euler
  flows supported on finitely many Fourier modes*:
  https://arxiv.org/abs/2110.08039
- A. Mahalov, E. S. Titi, and S. Leibovich, *Invariant helical subspaces for
  the Navier--Stokes equations*: https://doi.org/10.1007/BF00381234
