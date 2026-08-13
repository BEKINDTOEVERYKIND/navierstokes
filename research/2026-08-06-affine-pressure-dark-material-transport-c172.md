# C172: affine pressure-dark fibers cannot form the localized A2 focus packet

**Date:** 2026-08-06

**Status:** exact constant-matrix fiber classification, localized-kernel
no-go for rank at least two, fixed rank-one shear exception, universal
zero-order affine-multiplier uniqueness, and translation-invariant
finite-order symbol no-go; full one-cell dynamics remain conditional

**Checker:**
[checks/affine_pressure_dark_material_transport_c172.py](../checks/affine_pressure_dark_material_transport_c172.py)

## 0. Claim boundary

C171 leaves the exact affine part of the co-moving parent cross residual

\[
                         2AU,                       \tag{0.1}
\]

where \(A=\nabla V\) is trace free.  This note asks whether (0.1) can be
made pressure-dark, so that Leray projection removes it before the
backward active propagator sees it.

For one Fourier fiber the answer has a complete algebraic classification.
For \(k\ne0\) and \(k\cdot u=0\),

\[
             {\mathbb P}_kAu=0
             \quad\Longleftrightarrow\quad
             Au\in{\mathbb R}k.                    \tag{0.2}
\]

The dark fiber is the kernel of a two-by-two compression.  Its determinant
is

\[
 {k^T\operatorname {adj}(A)k\over |k|^2}.          \tag{0.3}
\]

Thus an invertible \(A\) has a one-dimensional dark polarization exactly
on the quadratic cone \(k^TA^{-1}k=0\).  A rank-two \(A\) has dark
frequencies on the union of two planes determined by its left and right
kernels.  These sets have three-dimensional measure zero.  Consequently:

> If \(\operatorname {rank}A\ge2\), the only \(L^2(\mathbb R^3)\)
> divergence-free field satisfying \({\mathbb P}(AU)=0\) is \(U=0\).
> In particular there is no nonzero compact/curl packet with exact affine
> pressure darkness.

Real global Fourier pairs are compatible with the dark cone, and a
periodic torus may contain individual dark lattice modes.  The obstruction
is localization: a spatially localized packet has a continuum of
frequencies and cannot live on a measure-zero characteristic set.

There is one genuine local exception.  If
\(A=a\otimes b\) has rank one and trace zero, then \(a\cdot b=0\), and

\[
                       U=\nabla\psi\times b          \tag{0.4}
\]

is compact, divergence free, and lies in \(\ker A\) for every compact
\(\psi\).  But \(A^2=0\), and its Piola image has no amplitude strain on
this kernel.  This shear exception cannot provide the scheduled
exponential focus.

For the C142 selector in its displayed initial affine frame, \(A\) is
symmetric rank two with eigenvalues \(0,\lambda,-\lambda\).  The C140
child and both named wakes start with \(N\)-polarization and are not dark.
Making those three \(N\)-polarized fibers dark for one constant selector
would force \(AN=0\), which makes their exact Kelvin velocity amplitude
\(N\) constant and hence removes the C142 velocity-amplitude gain.

Finally, among datum-independent zero-order affine matrix multipliers,
divergence preservation forces the Piola map up to one scalar factor.  The
scalar can only replace \(2A\) by \(2A+\gamma I\); for the selector this
matrix still has rank at least two for every \(\gamma\).  A universal
translation-invariant finite-order differential rule acting directly on
the velocity cannot reproduce the degree-zero rational Leray symbol on an
open frequency cone.  The exact Fourier-multiplier alternative is a
nonlocal Hodge/Biot--Savart correction (or a specially engineered
support-dependent antidivergence), whose global tail remains in BAFL.

C172 therefore rules out pressure darkness as a local cure of C171's
affine term on the existing focus geometry.  It does not rule out the
small \(q^{-1}\) estimate already obtained in C171, a dynamic integrated
cancellation, spatial export, MCKC, LCE, BAFL, or the full unforced stage.

## 1. Exact dark-fiber classification

Let

\[
 {\cal K}_A(k)=
 \{u\in k^\perp:{\mathbb P}_kAu=0\},
 \qquad
 {\mathbb P}_k=I-{k\otimes k\over|k|^2}.           \tag{1.1}
\]

This is the kernel of

\[
 T_A(k)={\mathbb P}_kA|_{k^\perp}:k^\perp\to k^\perp.
                                                               \tag{1.2}
\]

Choose an oriented orthonormal basis \(e,f,\widehat k\), where
\(\widehat k=k/|k|\).  The determinant of (1.2) is the oriented area
multiplier

\[
 \begin{aligned}
 \det T_A(k)
 &=\widehat k\cdot(Ae\times Af)\\
 &=\widehat k^T\operatorname {cof}(A)\widehat k\\
 &=\boxed{{k^T\operatorname {adj}(A)k\over|k|^2}.} \tag{1.3}
 \end{aligned}
\]

Here \(\operatorname {cof}(A)_{ij}=(-1)^{i+j}
\det A_{\widehat i\widehat j}\) and
\(\operatorname {adj}(A)=\operatorname {cof}(A)^T\).  The area identity is
\((Ae)\times(Af)=\operatorname {cof}(A)(e\times f)\); only the final
quadratic is unchanged by transposition.  Hence

\[
 \dim{\cal K}_A(k)=
 \begin{cases}
 0,&k^T\operatorname {adj}(A)k\ne0,\\
 1,&k^T\operatorname {adj}(A)k=0
       \text{ and }{\mathbb P}_kA{\mathbb P}_k\ne0,\\
 2,&{\mathbb P}_kA{\mathbb P}_k=0.
 \end{cases}                                      \tag{1.4}
\]

This is a complete classification for every real \(3\times3\) matrix; the
trace-free condition is not needed until the dynamical consequences below.

### 1.1 Invertible matrices

If \(A\) is invertible, (0.2) gives

\[
 u=\alpha A^{-1}k,\qquad
 0=k\cdot u=\alpha k^TA^{-1}k.                    \tag{1.5}
\]

Thus

\[
 {\cal K}_A(k)=
 \begin{cases}
 \operatorname {span}\{A^{-1}k\},&
          k^TA^{-1}k=0,\\
 \{0\},&k^TA^{-1}k\ne0.
 \end{cases}                                      \tag{1.6}
\]

The two-dimensional case in (1.4) cannot occur because the restriction of
an invertible map to a two-dimensional subspace is injective.

For an invertible symmetric trace-free \(A\), the eigenvalues have both
signs, so the real cone in (1.6) is nonempty.  It still has Lebesgue
measure zero.

### 1.2 Rank two

Let \(r\) span \(\ker A\), let \(\ell\) span \(\ker A^T\), and normalize
neither.  With the convention fixed above, the adjugate has rank one and
the cofactor has the transposed orientation:

\[
 \operatorname {adj}(A)=c\,r\otimes\ell,
 \qquad \operatorname {cof}(A)=c\,\ell\otimes r.              \tag{1.7}
\]

Indeed, the columns of the adjugate lie in \(\ker A\), while its rows lie
in \((\ker A^T)^T\); hence \(c\ne0\).  Consequently

\[
 k^T\operatorname {adj}(A)k
       =c(k\cdot r)(k\cdot\ell).                  \tag{1.8}
\]

The characteristic set of frequencies with a nontrivial dark fiber is
exactly

\[
                  r^\perp\cup\ell^\perp.          \tag{1.9}
\]

The generic dark fiber on these planes is one dimensional; the
two-dimensional exceptional directions are exactly those satisfying the
last line of (1.4).

For the symmetric selector, \(r=\ell\), so the characteristic set reduces
to the single plane \(r^\perp\).  In an eigenbasis with

\[
                         A=\operatorname {diag}(0,\lambda,-\lambda),
                                                               \tag{1.10}
\]

every \(k\) with a nonzero kernel-direction component has no dark
polarization.  On \(r^\perp\), the kernel vector \(r\) is always one dark
polarization.  A second appears only on the two null directions
\(k_2=\pm k_3\).

### 1.3 Rank one and zero

Every real rank-one matrix is \(A=a\otimes b\).  For the incompressible
affine gradients under study, trace freedom additionally says

\[
                         a\cdot b=0.               \tag{1.11}
\]

For \(k\not\parallel a\), (0.2) forces \(b\cdot u=0\); hence

\[
 {\cal K}_A(k)=k^\perp\cap b^\perp.                \tag{1.12}
\]

It is one dimensional except when \(k\parallel b\), when it is two
dimensional.  If \(k\parallel a\), every \(u\in k^\perp\) is dark, so the
dimension is again two.  For \(A=0\), every transverse polarization is
dark.

## 2. Reality is compatible; localization is not

If \(A\) and \(k\) are real and \(Au=\alpha k\), then

\[
 A\overline u=\overline\alpha k
              =-\overline\alpha(-k).              \tag{2.1}
\]

Thus the coefficient at \(-k\) required by Fourier reality is dark whenever
the coefficient at \(k\) is dark.  There is no reality-pair obstruction to
one global mode.

Now let \(U\in L^2(\mathbb R^3)\) be divergence free and suppose

\[
                         {\mathbb P}(AU)=0.        \tag{2.2}
\]

For almost every \(\xi\ne0\),

\[
 \widehat U(\xi)\in{\cal K}_A(\xi).                \tag{2.3}
\]

If \(\operatorname {rank}A=3\), then \(\operatorname {adj}(A)\) is
invertible.  Its quadratic form cannot vanish identically: that would make
its symmetric part zero, hence make it a real skew-symmetric matrix in odd
dimension, which is singular.  If \(\operatorname {rank}A=2\), (1.8) is
a nonzero polynomial.  In either case the determinant polynomial is
nonzero and its zero set has three-dimensional measure zero.  Equation
(2.3) therefore forces \(\widehat U=0\) almost everywhere:

\[
 \boxed{\operatorname {rank}A\ge2,\quad
        U\in L^2,\quad\nabla\cdot U=0,\quad
        {\mathbb P}(AU)=0
        \ \Longrightarrow\ U=0.}                 \tag{2.4}
\]

This is stronger than a compact-support no-go.  Distributional plane waves
or surface-supported Fourier data can live on the characteristic set, but
an \(L^2\) packet cannot.

For rank one, choose any \(\psi\in C_c^\infty\) and put

\[
                  U=\nabla\psi\times b
                   =\nabla\times(\psi b).          \tag{2.5}
\]

Then

\[
 \nabla\cdot U=0,\qquad b\cdot U=0,\qquad AU=0.    \tag{2.6}
\]

This is an exact compact/curl pressure-dark packet for every rank-one
matrix.  For the fixed trace-free affine matrix relevant here, (1.11)
also gives \(A^2=0\).  For the affine flow
\(F(t)=e^{tA}=I+tA\), its Piola image satisfies

\[
        U(t,F(t)a)=F(t)U_0(a)=U_0(a).             \tag{2.7}
\]

It is sheared in position but receives no material, \(L^p\), or pointwise
amplitude gain (\(\det F=1\)); the fixed nilpotent flow and its inverse grow
only linearly in time.  Thus the only nonzero rank admitting dark fibers
over an open frequency set loses the exponential-in-time strain mechanism
required by the displayed C142 stage.  This conclusion is only for one
fixed rank-one trace-free \(A\); it does not classify time-dependent shears
whose kernel direction changes.

On a torus, integer points of (1.6) or (1.9) may give real smooth global
dark modes.  A rank-two characteristic plane carries a rank-two lattice
only when its normal is rational relative to the torus lattice; an
irrational plane can have a smaller intersection.  The displayed C142
plane is rational, so existence of periodic dark modes is not the issue.
Multiplying one by a genuinely three-dimensional spatial cutoff introduces
normal-frequency content off the characteristic set and reintroduces the
projected affine residual.  A cutoff constant in the normal direction can
stay on the plane, but is not localized in three dimensions.  No stronger
assertion about an exotic characteristic Fourier series is needed here:
the C143 curl-localized packet under test is the three-dimensional cutoff
construction just described.

## 3. The C140 start-frame child and wakes are bright for the C142 selector

Use the unchanged \(A_2\) vectors

\[
 \begin{aligned}
 N&=(1,1,1),&r&=k_c=(1,0,-1),\\
 d&=(-1,2,-1),&
 e_1&=(2,-1,-1),\qquad e_2=(1,1,-2).
 \end{aligned}                                    \tag{3.1}
\]

They are mutually arranged as in C142:

\[
 r\perp d,\quad r\perp N,\quad d\perp N,\qquad
 e_{1,2}={3\over2}r\mp{1\over2}d.                 \tag{3.2}
\]

At the initial affine frame used to launch the C142 Kelvin modes, the
symmetric selector gradient is a nonzero scalar multiple of

\[
                         A_*=d\otimes N+N\otimes d. \tag{3.3}
\]

More precisely, in C142's notation
\(B=v\otimes w+w\otimes v=A_*/(3\sqrt2)\) and
\(A=-\rho'B=\sigma A_*\), where
\(\sigma=-\rho'/(3\sqrt2)\).  The normalization \(A_*\) is used below.

It is trace free and symmetric, with eigenvalues
\(0,3\sqrt2,-3\sqrt2\), and

\[
 A_*r=0,\qquad A_*N=3d.                           \tag{3.4}
\]

The C140 child coefficient at \(k_c=r\) and the two named wake coefficients
at \(e_1,e_2\) are all \(N\)-directed at this launch frame.  The projected
one-\(A_*\) factors have exact squared sizes

\[
 \begin{aligned}
 |{\mathbb P}_rA_*N|^2&=54,\\
 |{\mathbb P}_{e_1}A_*N|^2
 =|{\mathbb P}_{e_2}A_*N|^2&={81\over2}.           \tag{3.5}
 \end{aligned}
\]

Thus, for \(A=\sigma A_*\), the full C171 cross residual \(2AU\) has
squared projected factors

\[
 |{\mathbb P}_r(2A)N|^2=216\sigma^2,
 \qquad
 |{\mathbb P}_{e_i}(2A)N|^2=162\sigma^2.          \tag{3.6}
\]

Equivalently these are \(12(\rho')^2\) and \(9(\rho')^2\), respectively.
These are the symbol sizes per displayed \(N\)-polarized coefficient; a
mode coefficient \(\alpha N\) contributes the additional factor
\(|\alpha|^2\).
None is pressure-dark when \(\sigma\ne0\).  This agrees with (1.10): the
child frequency is the kernel direction itself, which lies outside the
characteristic plane, and each wake has a nonzero component in that
direction.

There is also a selector-independent incompatibility.  If one constant
matrix \(A\) made the \(N\)-polarized child and both \(N\)-polarized wakes
dark, then

\[
 AN\in\operatorname {span}r
       \cap\operatorname {span}e_1
       \cap\operatorname {span}e_2=\{0\}.          \tag{3.7}
\]

Thus \(AN=0\).  The exact affine Kelvin system is

\[
 k'=-A^Tk,\qquad
 a'=-Aa+2{k\over|k|^2}(k^TAa),\qquad k\cdot a=0.  \tag{3.8}
\]

Starting with \(a(0)=N\) and \(k(0)\cdot N=0\), the condition \(AN=0\)
gives

\[
               a(t)=N,\qquad k(t)\cdot N=0         \tag{3.9}
\]

for the whole affine interval.  This last conclusion requires the same
constant \(A\) (or, more generally, \(A(t)N=0\) throughout the interval):
darkness imposed only at one time does not propagate by itself.  Hence
simultaneous pressure darkness of the three start-frame passive
\(N\)-polarized components removes their affine velocity gain.  A constant
selector can darken one specially aligned fiber, but not the nonparallel
C140 start-frame child/wake family while retaining the C142 velocity gain.

## 4. Piola is the unique universal zero-order affine multiplier

Let \(F(t)\) be a volume-preserving affine flow map.  Consider a zero-order
material rule with a spatially constant matrix multiplier

\[
                    U(t,F(t)a)=M(t)u_0(a),         \tag{4.1}
\]

where \(M\) is independent of the material label and of the particular
datum.  Apply (4.1) to every divergence-free plane wave
\(u_0=v e^{i\xi\cdot a}\),
\(\xi\cdot v=0\).  Its Eulerian divergence is zero exactly when

\[
        \xi\cdot(F^{-1}M)v=0
        \quad\hbox{for every }v\perp\xi.           \tag{4.2}
\]

Thus \((F^{-1}M)^T\xi\) is parallel to \(\xi\) for every \(\xi\).
Every vector is an eigenvector of this matrix, so it is scalar:

\[
                         \boxed{M(t)=c(t)F(t).}     \tag{4.3}
\]

Piola transport is therefore unique among universal zero-order affine
matrix multipliers, up to a scalar amplitude.  This does not classify a
label-dependent multiplier, whose derivative contributes additional
terms, or a rule tailored to one prescribed datum.

On any interval on which the transported field is nontrivial and
\(c\ne0\), write \(\gamma=c'/c\).  Since \(F'=AF\),

\[
 D_t^VU=(A+\gamma I)U.                            \tag{4.4}
\]

The linearized Euler cross residual becomes

\[
                         (2A+\gamma I)U.           \tag{4.5}
\]

For the C142 selector with eigenvalues \(0,\lambda,-\lambda\), the three
eigenvalues in (4.5) are

\[
                  \gamma,\quad\gamma+2\lambda,\quad
                  \gamma-2\lambda.                \tag{4.6}
\]

When \(\lambda\ne0\), at most one can vanish for any \(\gamma\).  Hence
\(2A+\gamma I\) has rank at least two.  The localized no-go (2.4) applies
for every scalar reweighting \(c(t)\).  A scalar-Piola alternative cannot
produce a nonzero compact exact linearized-Euler packet for this selector.

For a different affine matrix with a repeated eigenvalue, a special
\(\gamma\) can make (4.5) rank one, and the compact kernel construction
(2.5) for that rank-one residual matrix may apply.  (That residual need
not itself be trace free or nilpotent, so the no-gain conclusion after
(2.7) does not automatically transfer.)  C172 does not rule out such a
degenerate repeated-eigenvalue case abstractly.  It is not the C142
\(0,\lambda,-\lambda\) selector, and it does not rescue the C140 family
analyzed in Section 3.

## 5. Why a translation-invariant finite-order rule still meets Hodge

At one instant, the exact material-coordinate amplitude generator in
(3.8) is

\[
 {\cal L}_A(k)u=-Au+
       2{k(k^TAu)\over|k|^2},\qquad u\perp k.       \tag{5.1}
\]

The second term is the order-zero rational Leray/Riesz symbol.  Suppose a
universal finite-order, constant-coefficient differential operator acting
directly on velocity agreed with (5.1) on every transverse polarization
in a nonempty open frequency set \(O\).  Write its polynomial symbol as
\(P(k)=\sum_{m=0}^M P_m(k)\), with \(P_m\) homogeneous of degree \(m\).

For each \(k\in O\), all \(sk\) with \(s\) in a small interval about one
remain in \(O\), while both \(u\perp sk\) and
\({\cal L}_A(sk)u={\cal L}_A(k)u\).  The polynomial identity in \(s\)
therefore gives

\[
 P_m(k)u=0\quad(m>0),\qquad P_0u={\cal L}_A(k)u
 \quad (u\perp k).                                \tag{5.2}
\]

The positive-degree symbols need not be zero as full matrices: they may
annihilate transverse inputs and act only on longitudinal inputs.  Modulo
precisely such invisible symbols, however, the induced operator on
divergence-free data is the constant matrix \(C=P_0\).

If such a \(C\) agreed with (5.1), then

\[
 (C+A)u\in\operatorname {span}k
 \quad\hbox{for }k\in O,\ u\perp k.                \tag{5.3}
\]

Put \(D=C+A\).  For each coordinate vector \(e_j\), the polynomial
\(k\times D(k\times e_j)\) vanishes on \(O\), hence everywhere.  Thus
\(D(k^\perp)\subseteq\mathbb Rk\) for every \(k\); testing the three
coordinate axes gives \(D=0\), so \(C=-A\).  Equation (5.1) would then
require \(k^TAu=0\) for every \(u\perp k\) on \(O\).  The same polynomial
continuation makes this true for every \(k\), equivalently
\(A^Tk\parallel k\) for every \(k\).  Thus \(A\) is scalar.  If it is
trace free, \(A=0\).

Therefore a nonzero trace-free affine gradient has no universal
translation-invariant finite-order differential operator on velocity equal
to the exact instantaneous Kelvin generator on a full open packet class.
Writing the operator as a curl does not change this conclusion if the
vector potential is itself obtained by finitely many constant-coefficient
derivatives of the velocity: the combined symbol is still polynomial.
The exact generator uses the Hodge multiplier \(|k|^{-2}\).  Implementing
it as a Fourier multiplier requires a pseudodifferential
pressure/Biot--Savart operation, or a support-dependent antidivergence
construction whose collar and exterior tail must be estimated separately.

This is not a no-go for variable-coefficient local systems with auxiliary
fields, for a potential already obtained nonlocally, or for a rule tailored
to finitely many modes.  It is about a universal, translation-invariant
finite-order rule on a full open packet class.  It does not contradict
individual dark modes on (1.6)/(1.9), the rank-one compact kernel (2.5), or
an engineered finite-dimensional cancellation.

## 6. Exact surviving obstruction

The pressure-dark attack ends at one precise incompatibility:

> **Pressure-dark localization/focus incompatibility (PDFI).**  The
> rank-two C142 selector has only measure-zero pressure-dark Fourier fibers,
> so it admits no nonzero \(L^2\) localized dark packet; its C140
> child/wake polarizations are bright in the displayed launch frame; and
> making all three start-frame fibers dark for one constant selector kills
> their Kelvin velocity gain.  Every universal zero-order affine matrix
> multiplier reduces to scalar Piola and retains rank at least two, while a
> universal translation-invariant finite-order velocity rule cannot replace
> the nonlocal Hodge symbol on an open frequency cone.

Thus pressure cannot erase the affine factor in C171 without either losing
localization/focus or paying the global pressure/collar tail.  The viable
route remains C171's small \(q^{-1}\) residual plus MCKC, not an exact local
pressure cancellation.  The self/viscous residuals, integrated wake kernel,
full unforced stage, LCE, and BAFL remain conditional.
