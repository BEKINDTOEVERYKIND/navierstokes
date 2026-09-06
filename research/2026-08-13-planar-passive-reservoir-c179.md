# C179: an exact planar passive reservoir has a full-rank one-edge symbol, but no focusing theorem

**Date:** 2026-08-13

**Status:** exact 2D3C unforced embedding, Fourier-edge determinant,
compact-cone inverse, static-shear Piola--Leray identity, and translation-coset
arithmetic; terminal preparation, the physical star, C125, RIGM, BAFL, and
the stage remain open

**Checker:**
[checks/planar_passive_reservoir_c179.py](../checks/planar_passive_reservoir_c179.py)

## 0. Claim boundary

C177 found an exact unforced homogeneous gate by putting the pump and
reservoir in one curl eigenspace. Exact darkness then confines the gate to
the pump shell, while C161's terminal chart asks for a broad family of
translations. This note tests the other exact unforced reservoir already
present in the C121 geometry.

Let \(N=(1,1,1)\) and \(H=N^\perp\). The C121 field is independent of the
\(N\)-coordinate and splits into an in-plane velocity plus an \(N\)-directed
passive scalar. The latter may contain arbitrary planar frequencies, not
only one Beltrami shell. It has no self-interaction and no feedback into
the planar pump. This gives a genuine positive replacement for prescribed
gate pulses at the level of the **unforced background**.

For one planar reservoir mode \((g,h)\), \(g\in H\), \(h\parallel N\), and
one off-plane source mode \((p,a)\), the exact Leray-projected tangent edge
to \(k=p+g\) is

\[
 {\cal A}_{p,g}a=-iP_k\{(h\cdot p)a+(a\cdot g)h\}.       \tag{0.1}
\]

As a map \(p^\perp\to k^\perp\), its determinant, in oriented orthonormal
area conventions, is

\[
 \boxed{\det A_{p,g,h}
 ={(h\cdot p)^2(|p|^2-|g|^2)\over |p|\,|p+g|}.}          \tag{0.2}
\]

Here \(A\) is the real symbol with the common factor \(-i\) omitted. That
factor, and a change of orientation, change the
displayed sign, not its zeros or singular values. Thus a separated compact
cone with nonzero normal charge, nonzero output, and unequal source/gate
radii has a uniform two-polarization inverse. The one-helicity rank defect
which motivated C177 is not a universal one-edge obstruction here.

This is **not** a focusing theorem.  The static-shear Piola--Leray proxy is
exactly reducible after the matching covector is transported, while C181
computes the distinct Euler Kelvin propagator and proves that it is uniformly
bounded on the fixed C159/C176 input cone.  The passive-scalar maximum
principle also prevents the total reservoir itself from being amplified.  A
complete terminal theorem must analyze the full finite-frequency
source--daughter--reverse evolution, not multiply one-edge determinants.

The surviving **planar passive reservoir gate (PPRG)** is:

> Prepare, by one unforced advection--diffusion orbit, a reality-complete
> planar scalar whose full Leray propagator takes the C176 \(q^2\) packet to
> a coherent normalized \(q^3\) endpoint, including reverse edges,
> collisions, depletion, endpoint profile, C125/RIGM, and BAFL.

No such theorem, one-cell stage, or singularity is claimed here.

## 1. Exact unforced triangular reservoir

Write a field independent of the \(N\)-coordinate as

\[
                      B(x,t)=v(x,t)+\Theta(x,t)N,
 \qquad v\cdot N=0.                                  \tag{1.1}
\]

Both \(v\) and \(\Theta\) depend only on the two planar variables. Since
\(N\cdot\nabla=0\), substitution into unforced Navier--Stokes gives exactly

\[
 \begin{aligned}
  \partial_tv+{\mathbb P}_{H}(v\cdot\nabla v)&=\nu\Delta v,\\
  \partial_t\Theta+v\cdot\nabla\Theta&=\nu\Delta\Theta.
 \end{aligned}                                           \tag{1.2}
\]

There is no term from \(\Theta\) in the first line. In particular, take
the planar part \(v_{A_2}(t)\) of the exact C121/C140 pump and let
\(\Theta\) be the unique smooth solution of the second line for arbitrary
smooth periodic planar data. Then

\[
                  B(t)=v_{A_2}(t)+\Theta(t)N               \tag{1.3}
\]

is an exact unforced three-dimensional Navier--Stokes solution. The
reservoir may have \(qK\)-scale planar modes even though the pump has
\(K\)-scale modes. Gate--gate interaction and gate-to-pump feedback vanish
identically, not perturbatively.

For inviscid transport, if \(X(t,a)\) is the Lagrangian flow of
\(v_{A_2}\), then

\[
 \Theta(t,X(t,a))=\Theta_0(a),\qquad
 \Theta_0=\Theta_T\circ X(T,\cdot)                         \tag{1.4}
\]

prepares any smooth terminal profile \(\Theta_T\) exactly. With
\(\nu>0\), forward evolution from chosen initial data is still exact, but
arbitrary exact backward terminal preparation is a parabolic range
problem. C178's finite polynomial heat preparation suggests an
approximation mechanism; it does not prove the required
advection--diffusion preparation in the analytic buffer or the C125 norm.

The total scalar obeys

\[
 \|\Theta(t)\|_\infty\le\|\Theta(0)\|_\infty,
 \qquad
 {1\over2}{d\over dt}\|\Theta\|_2^2
       =-\nu\|\nabla\Theta\|_2^2.                         \tag{1.5}
\]

These identities do not forbid transient growth of the **linearized
velocity perturbation** about (1.3); they do forbid counting rearrangement
of a large hidden scalar as newly created reservoir amplitude.

## 2. Exact one-edge Leray symbol

Use arbitrary Euclidean coordinates for the calculation. Let

\[
 p\ne0,\quad g\ne0,\quad g\cdot h=0,\quad
 a\cdot p=0,\quad k=p+g\ne0,                              \tag{2.1}
\]

and specialize below to \(h\parallel N\), \(g\in N^\perp\). The
linearized Euler interaction between the source
\(a e^{ip\cdot x}\) and reservoir \(h e^{ig\cdot x}\) is (0.1). Define
the real symbol without the common factor \(-i\) by

\[
       A_{p,g,h}a=P_k\{\alpha a+(a\cdot g)h\},
       \qquad \alpha=h\cdot p.                            \tag{2.2}
\]

For oriented area forms on \(p^\perp\) and \(k^\perp\), its determinant is

\[
 \det A={k\over|k|}\cdot
      \left(Ae_1\times Ae_2\right),                       \tag{2.3}
\]

where \((e_1,e_2,p/|p|)\) is an oriented orthonormal basis. Projection
does not change the \(k\)-normal triple product, so

\[
 \begin{aligned}
 k\cdot\bigl[\{\alpha e_1+(e_1\cdot g)h\}
              \times
              \{\alpha e_2+(e_2\cdot g)h\}\bigr]
 &=\alpha^2{k\cdot(p-g)\over|p|}\\
 &=\alpha^2{|p|^2-|g|^2\over|p|}.                         \tag{2.4}
 \end{aligned}
\]

Here one uses the exact oriented identity

\[
 \{\alpha e_1+(e_1\cdot g)h\}\times
 \{\alpha e_2+(e_2\cdot g)h\}
 ={\alpha^2\over|p|}(p-g),                                \tag{2.4a}
\]

which follows from
\((e_1\cdot g)e_2\times h+(e_2\cdot g)h\times e_1
=-\{(p\cdot g)/|p|\}\,h\), \(g\cdot h=0\), and
\(h\cdot p=\alpha\). Dividing (2.4) by \(|k|\) proves (0.2).

The zero set has a clean interpretation:

1. \(h\cdot p=0\): the source has zero normal charge, and one
   polarization is not transported by the vertical reservoir;
2. \(|p|=|g|\): the two-dimensional symbol has rank at most one; or
3. \(p+g=0\): the output fiber in (2.3) is undefined.

The equal-radius rank loss is exact. It is unrelated to an unlucky choice
of source basis or Leray coordinates.

### 2.1 A uniform compact-cone inverse

Fix constants \(0<c_0<C_0<\infty\), \(c_1,c_2>0\). On the normalized
compact set

\[
 \begin{aligned}
 {\cal K}=\{(p,g):\;&c_0\le|p|,|g|\le C_0,\quad
 |h\cdot p|\ge c_1,\\
 &\big||p|^2-|g|^2\big|\ge c_2,\quad |p+g|\ge c_1\},       \tag{2.5}
 \end{aligned}
\]

with fixed unit \(h\perp g\), (0.2) is bounded away from zero. The matrix
entries are continuous, so its largest singular value is bounded above;
therefore its smallest singular value is bounded below uniformly on
\({\cal K}\). This is a genuine two-polarization chart for each edge.

For the C176 application, source and gate frequencies are both of order
\(qK_j\), so the derivative in (2.2) is order \(qK_j\). To make each of
\(q\) orthogonal daughter edges have normalized strength
\(b/\sqrt q\), the scalar Fourier coefficient of each reservoir mode has
the scale

\[
                         \boxed{|\theta_a|\asymp
                         {b\over q\sqrt q}.}               \tag{2.6}
\]

The complete \(q\)-mode scalar has coefficient-\(\ell^2\) scale \(b/q\).
This is only the forward row normalization. It neither proves a skew
reverse row nor makes the nonlinear endpoint energy-preserving.

## 3. The planar shifts fit the C176 short direction

Let \(u\) be the C176 long-tube direction, so its reciprocal slab has
half-widths \(q,q,q/J\) in a frame whose short covector direction is
\(u/|u|\). The vector \(u\) is not planar: C152 gives

\[
                u=-T'(0)U_0-\sqrt2T(0)N.                  \tag{3.1}
\]

Choose instead a primitive planar integer vector
\(d\in\mathbb Z^3\cap N^\perp\) with \(u\cdot d\ne0\); for example
\(d=r_3=(-1,0,1)\), since \(U_0\parallel r_3\) and \(T'(0)\ne0\).
Translation by \(ad\) therefore changes the short C176 coordinate at a
nonzero fixed rate. It is false, and unnecessary, that \(d\parallel u\).

There is enough quotient arithmetic for \(q^2\) sources. Let
\({\cal B}_{q,J}\) be the C176 slab, with \(q/J\to\infty\). Its lattice
count is \(\Omega(q^3/J)\). A line \(k+\mathbb Zd\) has nonzero projection
on the short axis; hence its intersection with the slab has
\(O_d(q/J+1)=O_d(q/J)\) points. Therefore the slab meets

\[
                         \boxed{\Omega_d(q^2)}             \tag{3.2}
\]

distinct cosets of \(\mathbb Zd\). On the cofinal even-\(n\) schedule
already used in C161, \(q=n^8\) is even. Select \(q^2\) representatives
\(S\), one per coset, and a reality-compatible set

\[
 {\cal A}=\{a\in\mathbb Z:0<|a|\le q/2\},\qquad
 G=\{aK_jd:a\in{\cal A}\}.                               \tag{3.3}
\]

Impose the explicit separation hypothesis

\[
 (S+G)\cap\{-S,-G,-(S+G)\}=\varnothing,                  \tag{3.3a}
\]

after placing representatives and the shift band away from zero, and use
physical representatives \(K_jS\). The translation map on the selected
half-lattice representatives

\[
                       S\times G\longrightarrow S+G       \tag{3.4}
\]

is injective: equality first puts the two representatives in the same
\(\mathbb Zd\) coset, hence makes them equal, and then makes the shift
labels equal. Thus \(|S+G|=q^3\), and (3.3a) lets the reality partners be
adjoined without identifying an active representative with its conjugate.
The existence of one selection satisfying both the actual C176 slab and
(3.3a) is a hypothesis here, not a consequence of the quotient count
alone. This is the correct replacement for
identifying the planar direction with \(u\) or assuming the C161 product
chart inside the rotated C176 slab.

Equations (3.2)--(3.4) are support arithmetic only. They do not say that
every selected pair lies in one compact inverse cone (2.5), that the
passively evolved reservoir retains its Fourier support, or that repeated
sums avoid all active/wake collisions.

## 4. Why the determinant is not yet a focus

The geometric-optics boundary can be seen exactly on a static vertical
shear. In unit-normal coordinates write

\[
              W(x)=\Theta(x_H)n,\qquad n\cdot\nabla\Theta=0,
 \qquad A=\nabla W=n\otimes\nabla\Theta.                  \tag{4.1}
\]

Since \(A^2=0\), the physical flow derivative and covector map are

\[
        F(t)=I+t\,n\otimes\nabla\Theta,\qquad
        G(t)=F(t)^{-T}=I-t\,\nabla\Theta\otimes n.         \tag{4.2}
\]

Let a covector/amplitude pair satisfy \(p\cdot a=0\), and let the normal
charge \(m=n\cdot p\ne0\). The transported covector and the associated
Piola--Leray amplitude are

\[
 p(t)=G(t)p,\qquad
 a_{\rm PL}(t)=P_{p(t)}G(t)a.                             \tag{4.3}
\]

Indeed,

\[
 G(t)a-{n\cdot a\over m}G(t)p
        =a-{n\cdot a\over m}p                             \tag{4.4}
\]

is constant. Since subtracting a multiple of \(p(t)\) does not change
the projection,

\[
 \boxed{P_{G(t)p}G(t)a
       =P_{G(t)p}\left(a-{n\cdot a\over m}p\right).}       \tag{4.5}
\]

Thus the same phase/shear admits an exactly reducible transported
Piola--Leray response rather than an independent gain on every Fourier
edge. On compact time/covector sets away from \(G(t)p=0\), (4.5) is
bounded; over longer shear times the moving covector and projection must
be tracked explicitly. The actual incompressible Euler Kelvin amplitude
has an additional pressure/stretching law; (4.3) is not asserted to solve
that law. C181 computes that distinct propagator exactly. This
Piola--Leray identity is nevertheless the exact microlocal warning
needed here: instantaneous edge rank alone does not control the transported
propagator. It parallels C164/C169's phase-rearrangement and
maximum-principle warnings.

The A2-advected scalar is not static, and the exact full finite-frequency
propagator can mix these shear episodes noncommutatively. Equation (4.5)
does not rule that out. It says precisely what must be beaten: the
terminal mechanism has to produce coherent physical focus after the
covector/phase transport, not merely an invertible instantaneous edge.

## 5. Exact surviving obligations

C179 supplies an autonomous high-frequency reservoir class and removes a
local two-polarization rank objection. PPRG still requires all of:

1. **terminal-profile preparation:** an exact inviscid pullback or a
   quantitatively controlled viscous advection--diffusion range theorem on
   the C176 tube, with localization and C125-relative error;
2. **uniform packet symbol:** one reality-complete selection of \(S,G\)
   lying in compact separated cones, with coefficients/phases giving the
   \(b/\sqrt q\) forward row simultaneously for all \(q^2\) sources;
3. **full evolution:** reverse blocks, depletion of the source/reservoir,
   repeated translations, gate--daughter collisions, Leray pressure, and
   the transported-shear cancellation (4.5);
4. **physical focus:** a lower bound for one real velocity component at the
   endpoint, not a coefficient-\(\ell^1\) upper scale or a product of local
   determinants; and
5. **closure:** the actual wake in the C175 admissible graph/export class,
   followed by C125, RIGM, BAFL, LCE, and the pump-to-pump endpoint map.

The stage is therefore narrower but still open.

## 6. Verification boundary

The dependency-free checker verifies:

* exact 2D3C triangular advection algebra on Fourier modes;
* the one-edge Leray determinant by exact rational arithmetic on many
  nondegenerate examples, plus both exact rank-loss surfaces;
* a concrete compact separated family with a positive rational determinant
  margin;
* the \(b/(q\sqrt q)\) forward-row and reservoir-norm ledgers;
* quotient-coset injectivity for the planar \(r_3\) shifts and the
  nonparallel \(u\cdot r_3\ne0\) identity; and
* the static vertical-shear Piola--Leray identity (4.4)--(4.5).

It cannot verify viscous terminal preparation, the physical \(q\)-star,
uniformity over the actual C176 packet, reverse/collision/depletion
control, localization, C125, RIGM, BAFL, or a one-cell stage.
