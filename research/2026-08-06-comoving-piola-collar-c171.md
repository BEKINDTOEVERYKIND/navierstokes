# C171: co-moving Piola collars remove the stationary drift loss, but not the wake block

**Date:** 2026-08-06
**Status:** exact curl/Piola/material-derivative and Helmholtz identities;
conditional quantitative bound under an explicit material-chart hypothesis;
exact same-charge reality-paired counterexample to an automatic extra
interaction factor
**Checker:**
[checks/comoving_piola_collar_c171.py](../checks/comoving_piola_collar_c171.py)

## 0. Claim boundary

This note stays inside the C107/C143 one-cell geometry.  It tests the
specific suggestion that the compact affine converter should move with the
parent rather than remain stationary in Eulerian coordinates.

There are two distinct conclusions.

1. The leading stationary-cutoff term really can be removed.  A curl
   potential transported as a one-form gives an exactly divergence-free
   material tube.  Its velocity is the contravariant Piola transform, and
   its exact parent cross-residual is

   \[
       2(U\mathbin\cdot\nabla)V.                         \tag{0.1}
   \]

   Thus the large parent velocity no longer multiplies the collar gradient.
   Under a scale-uniform Piola-chart bound, (0.1) gains one factor
   \(r/\ell=q^{-1}\) relative to the stationary estimate.  In backward
   \(L^2\) scale the result is

   \[
       hq^{-5/2}\log h=q^{-1}\log h
       =12n^{-8}\log n=o(n^{-6}).                         \tag{0.2}
   \]

   The gain in (0.2) comes from \(|U|\asymp\lambda r\), whereas
   the stationary transport used \(|\nabla U|\asymp\lambda\).  It is not
   an unproved second-derivative or pressure cancellation.

2. This does **not** close the wake-to-active block.  On the unchanged
   \(A_2\) support there is an explicit real, zero-normal-charge triad from
   the C140 wake frequency to the child frequency.  Reality partners do not
   cancel it.  A fixed linear collar block is homogeneous in the wake
   amplitude, so a nonzero block maps a \(b^2\) wake to \(b^2\) times its
   integrated operator, not automatically to \(b^3\).  Charge and ordinary
   Fourier reality may enforce exact darkness in a specified block, but by
   themselves they do not manufacture a continuous factor \(b\) in every
   allowed block.

The full material construction also has an explicit unresolved hypothesis:
the Piola chart must not acquire excessive distortion over the gain window.
A pointwise affine-jet bound alone gives only a Gronwall estimate and does
not prove this.  Moreover, the Helmholtz projector is an \(L^2\) contraction,
but the backward active evolution used by BAFL need not be.  These two
requirements are combined below into one named
**material-collar kernel closure (MCKC)** obligation.  Consequently C171 is
not an unforced one-cell stage, LCE, BAFL, a cascade, or a blow-up claim.

## 1. Exact curl and Piola localization

Let \(V\) be a smooth incompressible parent flow and let

\[
 \partial_tX(t,a)=V(t,X(t,a)),\qquad
 F(t,a)=D_aX(t,a),\qquad \det F(t,a)=1.             \tag{1.1}
\]

Take a fixed cutoff \(\chi\) equal to one on \(B_1\) and supported in
\(B_2\).  For a trace-free matrix \(S\), \(\|S\|\lesssim\lambda\), put

\[
 {\cal A}_0(a)=\chi(a/r){(Sa)\times a\over3},
 \qquad w_0=\nabla_a\times{\cal A}_0.              \tag{1.2}
\]

The elementary identity

\[
 \nabla\times\left({(Sa)\times a\over3}\right)
 =Sa-{\operatorname {tr}S\over3}a=Sa              \tag{1.3}
\]

shows that \(w_0=Sa\) on \(B_r\).  Globally, \(w_0\) is smooth,
divergence free, and supported in \(B_{2r}\).

Transport the potential as a one-form:

\[
 {\cal A}(t,X(t,a))=F(t,a)^{-T}{\cal A}_0(a),
 \qquad U(t,x)=\nabla_x\times{\cal A}(t,x).        \tag{1.4}
\]

The curl-Piola identity for a volume-preserving map gives exactly

\[
 \boxed{U(t,X(t,a))=F(t,a)w_0(a).}                 \tag{1.5}
\]

It follows at once that

\[
 \nabla_x\cdot U=0,
 \qquad
 \operatorname {supp}U(t)=X(t,\operatorname {supp}w_0).
                                                               \tag{1.6}
\]

Thus (1.4) is a compact material tube in \(\mathbb R^3\), or a compact
local tube in a torus chart before periodic copies meet.  No Leray correction
is needed to make the velocity divergence free.

The potential and velocity material derivatives are also exact.  Writing
\(A=\nabla V\),

\[
 D_t^V{\cal A}=-A^T{\cal A},
 \qquad
 \boxed{D_t^VU=AU=(U\mathbin\cdot\nabla)V.}        \tag{1.7}
\]

The second equality in (1.7) is the Cauchy law for a transported
divergence-free vector.  In particular, Piola transport removes
\(V\cdot\nabla\chi\) from the residual, but it does not make a transported
velocity solve the linearized Euler equation: the strain term in (1.7)
has the same sign as the other linearized cross term.

## 2. Exact Navier--Stokes and pressure residual

Suppose \((V,p_V)\) solves unforced Navier--Stokes on the time interval in
question.  Substitution of \(V+U\), followed by (1.7), gives the exact raw
residual

\[
 \boxed{
 G=2(U\mathbin\cdot\nabla)V
      +(U\mathbin\cdot\nabla)U-\nu\Delta U.}        \tag{2.1}
\]

On \(\mathbb R^3\) with decay, or on the torus after subtracting the
irrelevant constant, define

\[
 \pi=-\Delta^{-1}\nabla\cdot G,
 \qquad R=G+\nabla\pi=\mathbb P G.                 \tag{2.2}
\]

Then \((V+U,p_V+\pi)\) solves Navier--Stokes with residual \(R\).
The raw residual has zero spatial mean: each advective term is a divergence
because \(U\) is divergence free, and the Laplacian also integrates to
zero.  The Helmholtz decomposition is orthogonal in \(L^2\), hence

\[
 \boxed{
   \|R\|_2^2+\|\nabla\pi\|_2^2=\|G\|_2^2,
   \qquad \|R\|_2,\ \|\nabla\pi\|_2\leq\|G\|_2.} \tag{2.3}
\]

Therefore Leray projection and its pressure complement do not enlarge this
\(L^2\) residual.  Leray projection nevertheless destroys compact support:
although \(G\) is supported on the material tube, \(\nabla\pi\) is generally
global.  Spatial export therefore still owes a separated pressure-tail and
periodic-image estimate.  Equation (2.3) also does **not** say that a
localized active projection followed by a backward gain propagator is
contractive.  That non-normal response is precisely part of LCE/BAFL.

## 3. The quantitative \(q^{-1}\) gain and its hypothesis

The reference cutoff obeys, with constants depending only on finitely many
seminorms of \(\chi\),

\[
 \|\nabla^mw_0\|_2\lesssim\lambda r^{5/2-m}
 \quad(0\leq m\leq2),
 \qquad \|w_0\|_\infty\lesssim\lambda r.           \tag{3.1}
\]

For the full material tube, introduce the dimensionless chart quantity

\[
 {\cal D}_{2,r}(X)=\sup_{[0,T]\times B_{2r}}
 \left(
  |F|+|F^{-1}|+r|D_aF|+r^2|D_a^2F|
 \right).                                         \tag{3.2}
\]

The chain and product rules applied to (1.5) imply

\[
 \begin{aligned}
 \|U\|_2&\leq M_F\lambda r^{5/2},&
 \|U\|_\infty&\leq M_F\lambda r,\\
 \|\nabla U\|_2&\leq M_F\lambda r^{3/2},&
 \|\Delta U\|_2&\leq M_F\lambda r^{1/2},
\end{aligned}                                    \tag{3.3}
\]

To make the hidden chart dependence explicit, write

\[
 f=\|F\|_\infty,\quad g=\|F^{-1}\|_\infty,
 \quad d_1=r\|D_aF\|_\infty,
 \quad d_2=r^2\|D_a^2F\|_\infty.                 \tag{3.3a}
\]

After enlarging one fixed cutoff constant if necessary, an admissible common
factor in (3.3) is

\[
 M_F=C_\chi\left[
 f+g(f+d_1)+g^2(f+2d_1+d_2)+g^3d_1(f+d_1)
 \right]
 \leq C_\chi(1+{\cal D}_{2,r}(X))^5.             \tag{3.3b}
\]

The four terms respectively control \(U\), the first chain rule, the
second derivative of \(Fw_0\), and the derivative of \(F^{-1}\) in the
second spatial chain rule.  Thus \(M_F\) is a polynomial in the full
Piola/F-jet quantity, not merely in the condition number of \(F\).  Volume
preservation prevents a separate volume factor.

Assume on the tube that

\[
             \|\nabla V\|_\infty\leq C_V{a\over\ell}.
                                                               \tag{3.4}
\]

The three pieces in (2.1) then satisfy

\[
 \begin{aligned}
 \|2(U\cdot\nabla)V\|_2
    &\leq 2C_VM_F{a\lambda\over\ell}r^{5/2},\\
 \|(U\cdot\nabla)U\|_2
    &\leq M_F^2\lambda^2r^{5/2},\\
 \|\nu\Delta U\|_2
    &\leq \nu M_F\lambda r^{1/2}.                 \tag{3.5}
\end{aligned}
\]

The last two lines are the transported versions of C143's self-Euler and
viscous residual estimates.  C171 does not close their backward
active-propagator bounds; it only shows that transport introduces no new
raw \(L^2\) scale in those pieces.  Their stage-level treatment remains an
inherited C143/BAFL obligation.

Normalize by the parent \(L^2\) scale \(a\ell^{3/2}\), take
\(r=\ell/q\), and let \(\lambda T=\log h\).  The time-integrated parent
cross term is

\[
 {T\|2(U\cdot\nabla)V\|_2\over a\ell^{3/2}}
 \leq C M_Fq^{-5/2}\log h.                        \tag{3.6}
\]

After the full backward focus \(h=q^{3/2}\),

\[
 \boxed{
 h{T\|2(U\cdot\nabla)V\|_2\over a\ell^{3/2}}
 \leq C M_Fq^{-1}\log h.}                         \tag{3.7}
\]

This is one factor \(q^{-1}\) below C143's stationary bound.  The
stationary calculation uses

\[
 \|V\cdot\nabla U\|_2\sim a\lambda r^{3/2},
\]

whereas (3.5) uses

\[
 \|(U\cdot\nabla)V\|_2
 \lesssim (\lambda r)(a/\ell)r^{3/2}.
\]

Thus the factor is exactly the amplitude ratio \(r/\ell=q^{-1}\).  It is
not credited to a mysterious curvature cancellation.

On the factorial stage

\[
 q=n^8,\qquad h=n^{12},
\]

and (3.7) is

\[
 12CM_Fn^{-8}\log n.
\]

It is \(o(n^{-6})\) provided

\[
               M_F\log n=o(n^2)
 \quad\left(\text{equivalently }M_F\log q=o(q^{1/4})\right). \tag{3.8}
\]

In particular, a uniform material chart suffices.

The terminal active chart must also be charged exactly once.  Under the
stage hypothesis

\[
 \kappa_{\rm ch}\leq Cn^2=Cq^{1/4},
\]

(3.7) becomes

\[
 \kappa_{\rm ch}\,C M_Fq^{-1}\log h
 \leq C M_Fq^{-3/4}\log h
 =12CM_Fn^{-6}\log n=o(n^{-4}).                  \tag{3.8a}
\]

The same strict condition (3.8) justifies the last comparison.  The
logarithm means that the charted bound is not \(O(n^{-6})\); the required
conclusion is the weaker and sufficient \(o(n^{-4})\).  If instead one
pays C142's full two-polarization condition number \(h^2=q^3\), the crude
charted scale is \(M_Fq^2\log h\) and fails.  Thus (3.8a) is only a
prescribed-ray estimate under the stated \(q^{1/4}\) stage-chart
hypothesis.

### 3.1 What an affine jet does and does not improve

Let \(c'=V(t,c)\), \(A(t)=\nabla V(t,c)\), and let
\(Y(t,a)=c(t)+F_A(t)a\), \(F_A'=AF_A\), be the volume-preserving affine
jet flow.  Write

\[
 V(t,x)=c'(t)+A(t)(x-c(t))+E(t,x).                 \tag{3.9}
\]

The affine-Piola field \(U_A(t,Y)=F_Aw_0\) obeys the exact cross identity

\[
 D_t^VU_A+(U_A\cdot\nabla)V
 =2AU_A+E\cdot\nabla U_A+(\nabla E)U_A.            \tag{3.10}
\]

If \(|E|\lesssim(a/\ell^2)|x-c|^2\) and
\(|\nabla E|\lesssim(a/\ell^2)|x-c|\), the last two terms in (3.10) have
one further factor \(q^{-1}\) relative to \(2AU_A\).  The dominant
\(q^{-1}\) term is nevertheless the **full affine gradient**
\(2AU_A\).  Curvature supplies only a smaller remainder unless a separate
algebraic/pressure cancellation removes that affine term.

For comparison, the rigid center-following curl field

\[
 U_c(t,x)=w_0(x-c(t))                              \tag{3.11}
\]

is exactly divergence free and has

\[
 D_t^VU_c+(U_c\cdot\nabla)V
 =[V(t,x)-V(t,c)]\cdot\nabla U_c+(\nabla V)U_c.   \tag{3.12}
\]

A local Lipschitz bound gives (3.6)--(3.7) directly, without a Piola
distortion factor.  This proves that following the central Lagrangian
trajectory is already enough for the parent-cross **residual estimate**.
It is not a material collar: parent particles may cross its boundary, so
it does not by itself prove wake slaving or spatial export.

### 3.2 The material-chart gate is real

From (1.1), the basic Gronwall estimate is only

\[
 |F(t)|+|F(t)^{-1}|
 \lesssim\exp\left(\int_0^t\|\nabla V(s)\|_\infty ds\right), \tag{3.13}
\]

with analogous differentiated estimates involving \(D^2V,D^3V\).
A controlled instantaneous affine jet therefore does not imply (3.8) over
a long gain window.  The exact remaining hypothesis for the fully material
version is:

This is the first clause of MCKC:

> **MCKC(i), subcritical chart/boundary crossing.**  Prove (3.8), including
> the first two label derivatives in (3.2), for the actual A2 converter
> tube through its entire gain window; or replace the full material tube by
> (3.11) while controlling boundary crossing and its retained wake.

MCKC(i) is separate from the parent-cross scale arithmetic.  Equation (3.7)
settles that arithmetic conditional on MCKC(i); it does not prove that
clause.  In particular, a borderline chart polynomial
\(M_F\asymp q^{1/4}\) is insufficient: (3.7) then equals
\(q^{-3/4}\log h=n^{-6}(12\log n)\), which misses \(o(n^{-6})\).

## 4. Charge and reality do not supply the missing \(b\)

Retain the exact C114/C140 frequencies

\[
 \begin{aligned}
 N&=(1,1,1),& k_1&=(1,-1,0),\\
 k_c&=(1,0,-1),& e_1&=(2,-1,-1)=k_1+k_c.
 \end{aligned}                                    \tag{4.1}
\]

All three frequencies have zero normal charge: \(N\cdot k_1=N\cdot
k_c=N\cdot e_1=0\).  Choose the real divergence-free fields

\[
 V_*(x)=A_*\cos(k_1\cdot x),\qquad
 W_*(x)=N\cos(e_1\cdot x),
 \qquad A_*=(1,1,-2).                              \tag{4.2}
\]

Indeed \(A_*\cdot k_1=0\), \(N\cdot e_1=0\), and
\(A_*\cdot e_1=3\).  Direct differentiation gives

\[
 \begin{aligned}
 (W_*\cdot\nabla)V_*&=0,\\
 (V_*\cdot\nabla)W_*
  &=-{3\over2}N\left[
     \sin((e_1+k_1)\cdot x)+\sin(k_c\cdot x)
     \right].                                     \tag{4.3}
 \end{aligned}
\]

Since \(N\cdot k_c=0\), Leray projection leaves the displayed child
coefficient unchanged.  In Fourier language the only ordered paths to
\(k_c\) are

\[
 (-k_1,e_1),\qquad(e_1,-k_1).                     \tag{4.4}
\]

The first contributes \(3iN/4\); the second is zero because
\(N\cdot k_1=0\).  The path at \(-k_c\) is its complex conjugate, not a
cancellation.  Thus zero normal charge, real Fourier pairing, and the
unchanged \(A_2\) support allow a direct wake-to-child edge.

The bright ordering in (4.3) is \((V_*\cdot\nabla)W_*\), whereas the
Piola residual (0.1) with the particular assignment \(U=W_*\), \(V=V_*\)
contains \(2(W_*\cdot\nabla)V_*=0\).  Thus (4.3) is not a lower bound on
that particular Piola residual; it is a counterexample to a charge/reality
kernel claim for the full Navier--Stokes wake block.  Even the Piola
ordering has no charge/reality darkness on the whole allowed polarization
class.  Indeed, with

\[
 \widetilde W_*(x)=B_*\cos(e_1\cdot x),
 \qquad B_*=(0,-1,1),
\]

one has \(B_*\cdot e_1=0\), \(B_*\cdot k_1=1\), and

\[
 2(\widetilde W_*\cdot\nabla)V_*
 =-A_*\sin((e_1+k_1)\cdot x)+A_*\sin(k_c\cdot x), \tag{4.4a}
\]

where

\[
 P_{k_c}A_*=(-1/2,1,-1/2)\ne0.                  \tag{4.4b}
\]

This second example is an allowed-class test, not an identification of
\(\widetilde W_*\) with the actual C140 wake polarization.

For that actual paired-gate family, C141 supplies the sharper historical
boundary.  Its first fixed-projector wake-fed return is cubic and has

\[
 D_{3,\mathrm{wake}}(k_c)
 =-{27\over\sqrt2}(|d_1|^2+|d_2|^2)\,T N\ne0      \tag{4.4c}
\]

whenever the terminal coefficient \(T\ne0\).  The two reality paths add.
Thus the homogeneous fixed-projector gate obtains one factor \(b\) only
because one additional \(b\)-sized interaction is required; reality does
not cancel the first support-allowed return or give such a factor to an
arbitrary localized linear block.

There is also a general homogeneity obstruction.  Once the geometry,
active projection, and background trajectory are fixed, the linearized
collar block \({\cal L}_{WA}\) is linear in the wake.  Therefore

\[
       {\cal L}_{WA}(b^2w)=b^2{\cal L}_{WA}w.       \tag{4.5}
\]

A support or symmetry rule may put \(w\) exactly in the kernel of
\({\cal L}_{WA}\).  If it does not, (4.5) contains no additional factor
\(b\).  A stage-dependent family could still satisfy
\(\|{\cal L}_{WA}^{(b)}w\|=O(b)\); that is precisely a quantitative kernel
estimate to prove, not a consequence of linear homogeneity.  With a
generator of size \(O(\lambda)\) acting for
\(T=\log h/\lambda\), the direct Duhamel product bound remains

\[
                   O(b^2\log h),                  \tag{4.6}
\]

whose ratio to the required \(b^3\) is \(\log h/b=12n^2\log n\).
Equation (4.3) proves that charge and reality do not force the kernel on
the allowed class.  It does not prove a lower bound for every engineered
time-dependent converter: oscillatory cancellation, a slaving graph, or
spatial export may still kill its integrated response.

Exactly disjoint velocity supports kill the local bilinear cross term,
which is the C107 spatial-export branch.  Overlap alone, however, supplies
no universal factor \(b\); a quantitative gain would require an additional
small-overlap, oscillatory, polarization, or kernel estimate.  Pressure and
heat tails must still be controlled if spatial separation is used.

## 5. Exact surviving obstruction

C171 closes the narrow C143 parent-cross arithmetic:

* center following replaces the backward stationary \(O(\log h)\) bound by
  \(O(q^{-1}\log h)=o(n^{-6})\) in \(L^2\);
* full material Piola transport has the same scale conditional on MCKC(i);
* Leray projection and its pressure complement do not enlarge the raw
  \(L^2\) residual.

It does not close LCE.  The exact named remaining obstruction is the single
two-clause closure:

> **Material-collar kernel closure (MCKC).**  (i) Establish the subcritical
> material-chart bound (3.8), or the rigid-collar boundary-crossing
> substitute; and (ii) prove for the actual retained C140 wake \(w\) that
> \(\|{\cal L}_{WA}w\|\leq Cb\|w\|\) (for example via controlled distance
> to the integrated wake-to-active kernel), or export it with stage-uniform
> pressure/heat tails.

The explicit triad (4.3) shows why this must be a dynamic slaving/export
theorem rather than a bare charge or reality-parity assertion.  The needed
localized active projection and backward evolution can amplify (2.1) even
though \(\mathbb P\) itself is an \(L^2\) contraction.  The full unforced
stage map, the inherited self/viscous response bounds, and BAFL remain
conditional.
