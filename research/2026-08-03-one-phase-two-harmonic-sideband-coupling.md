# One-phase two-harmonic bath versus the rank-five sideband chart

**Date:** 2026-08-03
**Status:** exact quadratic charge ledger and finite-support no-go; the
forward charged hierarchy remains open
**Scope:** joint audit of C89 and the preferred single-parent construction
in [the forward multiphase parametrix](2026-07-29-forward-multiphase-parametrix.md).
This note does not construct a Navier--Stokes singularity.

## 1. Verdict

The two pieces are compatible at the level of the principal endpoint
Jacobian, but they do not form a finite Fourier cell.

* Rotate the old common parent onto the compressive C89 phase. The
  fundamental bath harmonic and the three partner families retain exactly
  the old rank-five zero-charge strain map.
* The second covariance harmonic does not alter that instantaneous
  zero-charge symbol. It creates three additional charged outputs from each
  partner sign. Every coefficient is \(O(1)\), not \(O(K)\), by the exact
  one-phase incompressibility cancellation.
* Those extra outputs cannot be removed by choosing a partner polarization.
  For the three checked low directions, the second-harmonic/partner map is
  two-dimensionally invertible for all sufficiently large \(K\).
* More strongly, no nonzero finite collection of corrector modes can make
  the linearized two-harmonic source purely zero charge. An extreme charge
  always leaks through one of the second harmonics. Harmonic separation
  labels the wake but does not close it.
* Sequential partner pulses preserve the rank-five Jacobian and remove the
  direct simultaneous partner--partner products. They do not erase the
  bath--partner ladder or interactions with previously retained charged wake.
* Choosing every wave in sine phase gives the C89 covariance, preserves the
  rank-five sine child, and imposes central-odd symmetry. Total helicity is
  then exactly zero for the full Navier--Stokes evolution, including
  viscosity.

The resulting viable architecture is

\[
 \boxed{\text{two-harmonic one-phase bath}
        +\text{ sequential rank-five partners}
        +\text{ a retained infinite charged wake}.}
\tag{1.1}
\]

It is a concrete principal control design, not a finite-dimensional closed
transition. All rational identities below are checked in
[the exact checker](../checks/one_phase_two_harmonic_sideband.py).

## 2. Put both constructions in the same coordinates

Use the compressive phase direction and the two C89 covariance columns

\[
 w=e_1,\qquad a=e_2,\qquad d=\frac12e_3 .                 \tag{2.1}
\]

The centrally odd two-harmonic bath is

\[
 W_{\rm bath}(\theta)
 =\sqrt2\,a\sin\theta+\sqrt2\,d\sin(2\theta),
 \qquad \theta=K w\cdot x.                               \tag{2.2}
\]

It has covariance

\[
 \langle W_{\rm bath}\otimes W_{\rm bath}\rangle_\theta
 =a\otimes a+d\otimes d
 =\operatorname{diag}\left(0,1,\frac14\right).           \tag{2.3}
\]

Apply the cyclic orthogonal rotation

\[
 (x_1,x_2,x_3)\longmapsto(x_3,x_1,x_2)                   \tag{2.4}
\]

to the old chart. Its low directions become

\[
 \begin{aligned}
 q_1&=(20,-45,-36),\\
 q_2&=(9,-4,-5),\\
 q_3&=(1,1,1).
 \end{aligned}                                           \tag{2.5}
\]

For each \(j\), introduce a partner

\[
 r_j=q_j-Kw,\qquad b_j\cdot r_j=0.                       \tag{2.6}
\]

A bounded parameterization of its two polarization coordinates is

\[
 b_K(\beta)
 =\beta+\frac{\beta\cdot q}{K-q\cdot w}w,
 \qquad \beta\in w^\perp.                               \tag{2.7}
\]

Thus the physical amplitudes stay \(O(1)\) as \(K\to\infty\).

## 3. Complete first quadratic charge ledger

For divergence-free plane waves, suppressing the harmless Fourier factor
\(i\), write

\[
 {\cal B}_{p+r}(A,B)
 =\mathbb P_{p+r}\big[(A\cdot r)B+(B\cdot p)A\big].       \tag{3.1}
\]

Let a bath harmonic have charge \(m\in\{-2,-1,1,2\}\), wave vector
\(mKw\), and transverse polarization \(d_m\), where

\[
 d_{\pm1}\parallel a,\qquad d_{\pm2}\parallel d.          \tag{3.2}
\]

The sign and imaginary scalar imposed by the sine convention multiply the
whole coefficient and do not affect its charge, rank, or size. Since

\[
 Kb\cdot w=b\cdot q,                                     \tag{3.3}
\]

the interaction with the charge \(-1\) partner is exactly

\[
 \boxed{
 C_{m,q,K}b
 =\mathbb P_{q+(m-1)Kw}
   \left[(d_m\cdot q)b+m(b\cdot q)d_m\right].}
                                                               \tag{3.4}
\]

There is no \(K\) in the unprojected bracket. The projection has norm one,
and (2.7) is uniformly bounded. Hence every first bath--partner coefficient
is \(O(1)\). One partner sign gives:

| Bath charge \(m\) | Output charge \(m-1\) | Output wave | Role |
|---:|---:|---|---|
| \(-2\) | \(-3\) | \(q-3Kw\) | extra charged wake |
| \(-1\) | \(-2\) | \(q-2Kw\) | reality outward sideband |
| \(1\) | \(0\) | \(q\) | desired low child |
| \(2\) | \(1\) | \(q+Kw\) | extra covariance-harmonic sideband |

The reality-conjugate partner has slow frequency \(-q\) and gives the
mirrored charges \(-1,0,2,3\). Thus the complete first extra-charge set is

\[
 \{-3,-2,-1,1,2,3\}.                                    \tag{3.5}
\]

The bath has no self-interaction at frozen slow coefficients. Any two bath
waves are collinear in wave number and both polarizations are perpendicular
to \(w\), so both transport contractions in (3.1) vanish exactly. This is
true even though \(a\) and \(d\) are not parallel.

## 4. The zero-charge rank survives

At \(m=1\), (3.4) is precisely

\[
 L_{q,K}b
 =\mathbb P_q\left[(a\cdot q)b+(b\cdot q)a\right],        \tag{4.1}
\]

the orthogonal rotation of the checked single-parent child map. Therefore,
for the checked large carrier and throughout its asymptotic large-\(K\)
regime, each of the three blocks has rank two, and
the six symmetric-gradient columns

\[
 \operatorname{sym}(L_{q_j,K}b\otimes q_j)               \tag{4.2}
\]

have combined rank five. The second bath harmonic cannot contribute directly
to charge zero from an initial charge \(-1\) partner, because \(2-1=1\).
Thus the instantaneous principal zero-charge map is exactly unchanged.

This also gives a finite-time local statement. Linearize the forward
Galerkin system about the two-harmonic bath, start with partner control \(b\)
and zero low child, and denote the low response by \(z_0(t)\). Then

\[
 z_0(T)=T L_{q,K}b+O(T^2).                                \tag{4.3}
\]

Consequently the checked five-by-five endpoint minor is nonzero for every
sufficiently small \(T>0\). If the controlled endpoint is the deformation
rather than the child velocity, one more integration gives

\[
 \delta F(T)
 =\frac{T^2}{2}\nabla L_{q,K}b+O(T^3),                   \tag{4.4}
\]

and the same rank conclusion holds. The \(O(T^2)\) terms include feedback
paths through the new charged sidebands; they alter but cannot instantly
destroy a nonzero minor. This is a fixed-band, local-in-control statement,
not a uniform inverse for the full infinite hierarchy.

The new \(m=2\) output cannot be suppressed by a clever partner
polarization. Section 6 gives its exact determinant. For every \(q_j\) and
all sufficiently large \(K\),

\[
 b\longmapsto
 \mathbb P_{q+Kw}\left[(d\cdot q)b+2(b\cdot q)d\right]    \tag{4.5}
\]

is an isomorphism between the two transverse planes. Setting the charge
\(+1\) output to zero therefore forces \(b=0\) and destroys the desired
child.

## 5. Partner--partner products and sequential pulses

A single exact partner wave has no self-interaction, with either itself or
its conjugate, because \(b_j\cdot r_j=0\). Two different simultaneous
partners do interact. Exact transversality gives

\[
 \begin{aligned}
 b_j\cdot r_\ell
 &=b_j\cdot(q_\ell-q_j),\\
 b_j\cdot(-r_\ell)
 &=b_j\cdot(q_j-q_\ell).
 \end{aligned}                                           \tag{5.1}
\]

These coefficients are again \(O(1)\), with no \(K\). The first line creates
charge \(-2\) at slow frequency \(q_j+q_\ell\); the second creates charge
zero at the unwanted difference frequency \(q_j-q_\ell\). Conjugates give
the opposite charges and slow frequencies.

Disjoint sequential partner pulses remove these direct simultaneous
partner--partner products. At the zero-control point, interactions of an
earlier pulse's generated wake with a later partner are at least quadratic
in the controls, so the rank-five Jacobian remains the direct triangular
one for a sufficiently short pulse schedule. This justifies sequential
pulsing as a principal finite-dimensional design.

It does **not** give an exact nonlinear decoupling. The earlier charged wake
is still physical state. Unless it is exported or included in the next
forward solve, it interacts with the later bath and partner. Sequentiality
removes one source family; it does not reset the hierarchy.

One may instead turn off the second harmonic during the three endpoint
pulses. That deletes the new charges \(-3\) and \(1\) from (3.4), but it also
degenerates the rank-two covariance during the pulse, and the temporal ramp
has its own residual. The negative fundamental still creates the ordinary
outward charge \(-2\). Temporal segregation is a possible splice design,
not a finite-cell cure.

## 6. Harmonic separation does not close the charge ladder

Put the controlled parent at an arbitrary harmonic \(N>0\), and take its
matched partner at

\[
 r=q-NKw,\qquad NKb\cdot w=b\cdot q.                      \tag{6.1}
\]

Interaction with bath harmonic \(m\) is exactly

\[
 \mathbb P_{q+(m-N)Kw}
 \left[(d_m\cdot q)b+\frac mN(b\cdot q)d_m\right].         \tag{6.2}
\]

Increasing \(N\) suppresses only the second term. The first term
\((d_m\cdot q)b\) stays order one. Since the three directions (2.5) span
\(\mathbb R^3\), no nonzero bath polarization is perpendicular to all three.
Mere harmonic separation therefore cannot decouple the covariance bath from
a full rank-five partner family.

There is also an exact finite-support obstruction. Fix one of (2.5), write

\[
 q=(x,y,\zeta),\qquad p_h=q+hKw,                           \tag{6.3}
\]

and consider a perturbation \(Z_h\perp p_h\), including a possible
charge-zero component. The positive second
harmonic contributes to output charge \(h+2\) through

\[
 T_hZ_h
 =\mathbb P_{p_{h+2}}
 \left[(d\cdot q)Z_h+2K(Z_h\cdot w)d\right].               \tag{6.4}
\]

Put \(X=x+hK\) and use the basis

\[
 u_1=(-y,X,0),\qquad u_2=(-\zeta,0,X)                     \tag{6.5}
\]

of \(p_h^\perp\). Before final projection, let \(R_h\) denote the bracketed
map in (6.4). Direct expansion gives

\[
 \boxed{
 \det(R_hu_1,R_hu_2,p_{h+2})
 =\frac{\zeta^2}{4}\,X
   \left(|p_h|^2-4K^2\right).}                            \tag{6.6}
\]

Thus \(T_h:p_h^\perp\to p_{h+2}^\perp\) is invertible unless

\[
 X=0\quad\hbox{or}\quad |p_h|=2K,                         \tag{6.7}
\]

provided \(\zeta\ne0\). All three \(q_j\)'s have \(x>0\) and
\(\zeta\ne0\). If

\[
 K>|q|,
 \qquad
 K>\frac{|q|^2}{4x},                                     \tag{6.8}
\]

then (6.7) fails for every integer \(h\):

* for \(|h|\ge3\), \(|p_h|>2K\);
* for \(|h|\le1\), \(|p_h|<2K\);
* at \(h=2\), the shell gap is \(4Kx+|q|^2>0\);
* at \(h=-2\), it is \(-4Kx+|q|^2<0\).

For the three vectors (2.5), the worst threshold is \(q_1\), with
\(|q_1|=61\) and \(|q_1|^2/(4q_{1,1})=3721/80\). Hence every integer
\(K\ge62\) satisfies (6.8) simultaneously for all three blocks.

This proves a finite-support no-go.

> **Proposition 6.1.** Under (6.8), let
> \(\{Z_h:h\in S\subset\mathbb Z\}\) be a nonempty finite perturbation.
> Its linearized interaction with the real
> two-harmonic bath cannot have only charge-zero output.

Let \(H=\max S\). If \(H\ne-2\), output charge \(H+2\ne0\) lies above the
support. Only the positive second harmonic acting on \(Z_H\) can reach it,
so (6.6) forces \(Z_H=0\), a contradiction. If \(H=-2\), use
\(L=\min S\): output \(L-2\ne0\) lies below the support, and the identical
argument with the negative second harmonic forces \(Z_L=0\).

The bath has zero slow wave number, so its linearization preserves each
slow-\(q\) sector. On any nonaliased finite band the argument therefore
applies sector by sector; modes with different \(q\)'s cannot cancel the
extreme leakage.

The proposition concerns a static finite Fourier corrector. It does not
obstruct the forward Galerkin hierarchy, an infinite Gevrey profile, or a
finite band whose terminal charged state is explicitly exported. It proves
that one of those non-finite-cell options is mandatory.

## 7. Zero helicity is fully compatible

Equation (2.2) has zero averaged helicity. Curl converts every sine harmonic
to a cosine harmonic, and

\[
 \langle\sin(m\theta)\cos(n\theta)\rangle_\theta=0         \tag{7.1}
\]

for all \(m,n\). Distinct harmonics also have zero cross covariance, giving
(2.3).

The partners may all be placed in sine phase. The product of a parent sine
with a partner derivative contains a sine at their sum frequency, so the
desired low child remains in sine phase and the real two-polarization
rank-five chart loses no control coordinate.

More importantly, every sine wave is odd under central inversion. With an
even localization envelope and an odd low affine core, the complete initial
field can satisfy

\[
 u(0,-x)=-u(0,x).                                        \tag{7.2}
\]

Navier--Stokes preserves this symmetry: velocity, convection, pressure
gradient, and viscosity all remain odd, while pressure itself is even. Curl
of an odd velocity is even. Therefore the helicity density is odd and

\[
 \int u(t,x)\cdot\operatorname{curl}u(t,x)\,dx=0           \tag{7.3}
\]

for every time for which the symmetric solution exists. This is stronger
than the initial Fourier average calculation and survives unequal viscous
damping of the first and second harmonics. The wake router and collars must
be chosen reflection-equivariantly to retain this conclusion globally.

## 8. What is and is not gained

The audit resolves the immediate joint-coupling question.

1. The C89 bath does **not** reintroduce an \(O(K)\) cross-colour chain.
   Every new coefficient belongs to the existing one-phase \(K\)-null
   hierarchy.
2. The checked zero-charge rank-five derivative survives exactly at pulse
   onset and for short finite pulses.
3. Zero total helicity is compatible without adding another carrier or a
   compensating wake.
4. The two-harmonic bath makes a finite charged closure impossible. The
   positive result must be proved in the forward/export framework, not as a
   finite Fourier return cell.

The next load-bearing estimate is a uniform endpoint inverse for the
sequential pulses **inside** the complete one-phase material hierarchy,
with all charges in (3.5) and their descendants retained through order
\(M_j\), followed by the carried-wake contraction required by C90. The
present exact ledger says that this route has no naked carrier-frequency
loss, but it does not prove the needed Gevrey-2 semigroup or terminal wake
bound.
