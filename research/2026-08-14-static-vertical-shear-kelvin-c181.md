# C181: the exact static vertical-shear Kelvin map is bounded on a fixed cone, but can have long-sweep gain

**Date:** 2026-08-14

**Status:** exact affine/principal Kelvin propagator, compact-cone bound, and
explicit long-sweep lower bound; no finite-frequency, localized, viscous, or
nonlinear terminal-converter theorem

**Checker:**
[checks/static_vertical_shear_kelvin_c181.py](../checks/static_vertical_shear_kelvin_c181.py)

## 0. Claim boundary

C179 gives an exact unforced 2D3C reservoir and an invertible instantaneous
source--reservoir edge.  An instantaneous edge determinant is not a
propagator.  This note computes the simpler geometric-optics propagator for
the static vertical shear

\[
             U(x)=\Theta(x_H)n,\qquad n\cdot\nabla\Theta=0,       \tag{0.1}
\]

along one trajectory.  The horizontal base point is fixed, so
\(g=\nabla_H\Theta\) is constant on the trajectory.  Equivalently, this is
the exact Kelvin system for the affine shear \(U=(g\cdot x_H)n\), and the
principal/WKB system for a smooth static shear.

There are two complementary conclusions.

1.  If the initial covector lies in one fixed cone
    \(|k_H(0)|\le C|m|\), \(m=n\cdot k\ne0\), then the Kelvin velocity map has
    an operator norm bounded by a constant depending only on \(C\).  It
    cannot supply C161's missing \(\sqrt q\) terminal gain on such a fixed
    projective input tube, even if the later covector leaves it.
2.  If the initial cone degenerates, long-sweep gain is possible. A
    covector starting with horizontal coordinate \(-X\) and fixed normal
    charge, with a nonzero transverse horizontal coordinate, can be swept
    to \(+X\) and has an exact polarization whose gain is asymptotic to
    \((\pi/2)X\).  The initial covector itself then lies outside every fixed
    normal cone as \(X\to\infty\), and the construction must also pay its
    frequency, residence, localization, and wake ledgers.

The formula sometimes written as \(P_{F^{-T}k}F^{-T}a\) is a
Piola--Leray transport, not the incompressible Euler Kelvin velocity map for
this nonsymmetric shear.  It must not be substituted for the calculation
below.

Nothing here proves a full finite-frequency \(L^2\) estimate for a
non-affine shear, an A2-advected passive scalar, a physical \(q\)-star,
C125/RIGM/BAFL, or a one-cell stage.

## 1. Exact covector and Kelvin equations

Rotate horizontal coordinates so that

\[
       n=e_3,\qquad g=G e_1,\qquad G=|\nabla_H\Theta|.       \tag{1.1}
\]

Write the covector and its conserved normal charge as

\[
                 k(t)=(x(t),y,m),\qquad m\ne0.              \tag{1.2}
\]

With the convention \((\nabla U)_{ij}=\partial_jU_i\), the velocity
gradient is \(A=n\otimes g\).  The Kelvin system is

\[
 \dot k=-A^Tk,\qquad
 \dot a=-Aa+2k\,{k\cdot Aa\over |k|^2},\qquad k\cdot a=0. \tag{1.3}
\]

Thus

\[
             \dot x=-mG,\qquad \dot y=\dot m=0,             \tag{1.4}
\]

and, writing \(a=(u,v,w)\),

\[
 \begin{aligned}
 R^2&=x^2+y^2+m^2,\
 \dot u&={2mGxu\over R^2},\qquad
 \dot v={2mGyu\over R^2},\qquad
 w=-{xu+yv\over m}.
 \end{aligned}                                             \tag{1.5}
\]

The last relation is the divergence-free constraint.  It also recovers the
third Kelvin equation after differentiating and using (1.4)--(1.5).

## 2. Closed-form propagator

Put

\[
 d^2=y^2+m^2,\qquad R_0^2=x_0^2+d^2,
 \qquad
 I_d(x)={x\over2d^2(x^2+d^2)}
       +{1\over2d^3}\arctan{x\over d}.                    \tag{2.1}
\]

Since \(x'=-mG\), the first equation in (1.5) gives the exact invariant

\[
                         u(t)R(t)^2=u_0R_0^2.               \tag{2.2}
\]

Moreover \(I_d'(x)=(x^2+d^2)^{-2}\).  Therefore

\[
 \boxed{
 \begin{aligned}
 x(t)&=x_0-mGt,\qquad y(t)=y,\qquad m(t)=m,\\
 u(t)&=u_0{R_0^2\over R(t)^2},\\
 v(t)&=v_0-2y u_0R_0^2\{I_d(x(t))-I_d(x_0)\},\\
 w(t)&=-{x(t)u(t)+yv(t)\over m}.
 \end{aligned}}                                           \tag{2.3}
\]

This is the complete two-polarization Kelvin velocity map.  The value of
\(G\) enters only through the swept horizontal coordinate \(x(t)\).

### 2.1 A fixed-cone bound

Assume at the initial endpoint

\[
                   |x_0|,|y|\le C|m|.                     \tag{2.4}
\]

One has \(d\ge|m|\), \(R_0^2\le(1+2C^2)m^2\), and

\[
 |I_d(x)-I_d(x_0)|
 =\left|\int_{x_0}^{x}{ds\over(s^2+d^2)^2}\right|
 \le {\pi\over2|m|^3}.                                   \tag{2.5}
\]

The apparent terminal \(x(t)\) in \(w(t)\) is also harmless because
\(|x|/(x^2+d^2)\le1/(2d)\).  Consequently every coefficient of the
two-by-two map \((u_0,v_0)\mapsto a(t)\) is bounded by a finite explicit
polynomial in \(C\), uniformly for every \(x(t)\in\mathbb R\).  Since
\(|a_0|\ge (u_0^2+v_0^2)^{1/2}\), one may for example enlarge the elementary
component estimates to get

\[
                 \boxed{|a(t)|\le K_C|a_0|}                \tag{2.6}
\]

holds with

\[
 K_C=4(1+C)(1+2C^2)(1+\pi C).                              \tag{2.7}
\]

The constant is deliberately crude.  Its independence of \(q,G,t\) is
the useful statement.  Any fixed projective neighborhood of the C159 ray
has \(|m|/|k(0)|\ge\kappa>0\), hence is contained in (2.4) for some fixed
\(C\).  A terminal static-shear interval starting in that neighborhood has
only an \(O_C(1)\) principal multiplier, not the missing \(\sqrt q\), even
if the shear later moves the covector out of the cone.

By the exact Kelvin composition law, this conclusion may be appended to
the A2 amplifier if the covector entering the terminal shear satisfies
(2.4).  It says that this terminal factor cannot be
counted as an independent \(\sqrt q\) gain.  It does not bound an
interleaved A2/shear cocycle whose projective path leaves the cone.

## 3. Degenerating input cones can grow

The uniform input-cone hypothesis is essential.  Take

\[
 m=y=1,\qquad x_0=-X,\qquad x(T)=X,qquad v_0=0,
 \qquad u_0=(1+X^2)^{-1/2}.                                \tag{3.1}
\]

Then \(|a_0|=1\).  Since \(R(T)=R_0\), (2.3) gives \(u(T)=u_0\), and

\[
 v(T)=-2(X^2+2)u_0\{I_{\sqrt2}(X)-I_{\sqrt2}(-X)\}.       \tag{3.2}
\]

Using

\[
 I_{\sqrt2}(X)-I_{\sqrt2}(-X)
        \longrightarrow {\pi\over4\sqrt2},                \tag{3.3}
\]

and \(w(T)=-Xu(T)-v(T)\), one obtains

\[
                       \boxed{|a(T)|\sim {\pi\over2}X.}    \tag{3.4}
\]

Thus static vertical shear has no input-cone-independent amplification
bound.  This family starts with \(|m|/|k_0|\asymp X^{-1}\).  Obtaining a
factor \(\sqrt q\) by this channel therefore requires an input cone at
least as degenerate as \(|m|/|k_0|\lesssim q^{-1/2}\), outside C176's fixed
projective tube.  The exact formula identifies that cost; it does not prove
that a separately redesigned correlated packet, frequency, coherence, or
wake ledger is impossible.

## 4. Stationary inviscid scalar profiles on A2

For the A2 planar flow \(v=N\times\nabla f\), every

\[
                              \Theta=H(f)                   \tag{4.1}
\]

is an exact inviscid passive scalar because
\(v\cdot\nabla\Theta=H'(f)(N\times\nabla f)\cdot\nabla f=0\).
This supplies a static vertical profile in the common amplitude clock of
the inviscid A2 orbit.  Under viscosity,

\[
 \partial_tH(f)-\nu\Delta H(f)
 =-\nu\{H''(f)|\nabla f|^2+H'(f)\Delta f\}                 \tag{4.2}
\]

unless an independently chosen time dependence cancels the right side.
On C127, a fixed carrier-frequency profile has normalized heat drift
\(O(\mu_j T)\), factorially small on every polynomial or logarithmic
window.  A terminal profile with planar frequency \(qK_j\) instead pays
\(O(\mu_jq^2T)\); this is still factorially small for the polynomial
\(q=n^8\) schedule, but it does not preserve the exact form \(H(f)\) and
does not solve terminal preparation.

## 5. What remains open

C181 is an affine/principal calculation only.  A viable PPRG terminal
converter must still provide at least one of:

1. an interleaved noncommuting A2/passive-shear cocycle with the required
   \(\sqrt q\) gain while controlling its projective excursion; or
2. a finite-frequency focusing effect absent from the affine principal
   endpoint factor.

Either route must include the actual passive advection--diffusion orbit,
reality-complete Fourier spreading, reverse blocks, depletion, physical
point coherence, localization, C125/RIGM/BAFL, and the unforced endpoint
map.  The long-sweep example prevents upgrading (2.6) to a global
no-amplification theorem.

## 6. Verification boundary

The dependency-free checker verifies with exact rational arithmetic the
Kelvin reduction (1.4)--(1.5), invariant (2.2), and rational part of the
propagator; checks the antiderivative derivative numerically with a strict
margin; validates (2.3) against direct RK4 integration; checks the crude
fixed-cone bound on an exhaustive rational grid; and verifies convergence
of the long-sweep ratio to \(\pi/2\).  It does not certify a non-affine
finite-frequency propagator or any one-cell obligation listed above.
