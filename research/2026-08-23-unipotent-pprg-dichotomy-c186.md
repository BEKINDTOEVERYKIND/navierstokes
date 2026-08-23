# C186: abstract unipotent episode algebra does not force a common flag

**Date:** 2026-08-23

**Status:** conditional PPRG algebra: exact two-unipotent hyperbolic witness,
explicit Lyapunov and robustness constants, exact bounded passive 2D3C
rotating-gradient example, and a conditional common-conjugacy alignment
theorem; C183 does not prove that Kelvin polarization blocks are unipotent,
and realization of the two episode maps by one PPRG orbit remains open

**Checker:**
[checks/unipotent_pprg_dichotomy_c186.py](../checks/unipotent_pprg_dichotomy_c186.py)

## 0. Dichotomy verdict

C183 makes the **planar flow return** on a regular closed streamline
unipotent: \(L=I+N\), \(N^2=0\).  This is not the Kelvin polarization
monodromy \(\Phi\).  For a returning covector C183 proves only
\(\Phi\in SL(2,\mathbb R)\), not \((\Phi-I)^2=0\).  The classification below
is therefore a conditional algebraic lemma for candidate polarization
blocks whose unipotence has been established separately.

Unipotence, determinant one, and uniform block bounds alone do not force a
common flag.  The exact bounded pair

\[
 U=\begin{pmatrix}1&1\\0&1\end{pmatrix},\qquad
 V=\begin{pmatrix}1&0\\1&1\end{pmatrix}                 \tag{0.1}
\]

has no common invariant flag and its two-episode product is hyperbolic with
an explicit per-episode Lyapunov exponent greater than \(12/25\).  The
conclusion survives entrywise perturbations of size \(1/100\), with
per-episode exponent greater than \(9/20\), provided the perturbed blocks
remain in \(SL(2,\mathbb R)\).

This gives a positive **conditional finite-dimensional algebraic** PPRG
witness.  It is not the requested PDE dichotomy or a PPRG realization
theorem.  C183 does not prove unipotence of the pressure/polarization
connection, nor that one smooth passive 2D3C orbit can realize neighborhoods
of both matrices in (0.1) with a common return-fiber identification.  The
robustness result below does apply to general \(SL(2)\) blocks near this
specific pair; it is not a universal neighborhood theorem for arbitrary
unipotents.

There is also an exact PDE-level structural sublemma: a smooth bounded
passive scalar on an unforced steady 2D Euler background can rotate its
gradient from \(e_1\) to \(e_2\) along one trajectory.  Thus a blanket
claim that scalar passivity forces gradient alignment is false.  This still
does not identify the full Kelvin episode matrices with (0.1).

## 1. Exact two-episode witness

If unipotence of the two candidate polarization blocks is separately
proved, the complete two-by-two classification reduces their algebra to
one scalar.  Let \(N_1,N_2\ne0\) be real two-by-two matrices with
\(N_1^2=N_2^2=0\), and put

\[
                         \tau=\operatorname {tr}(N_1N_2). \tag{1.0}
\]

Each nilpotent has the unique flag
\(\ker N_i=\operatorname {im}N_i\).  Writing
\(N_i=u_i\otimes\alpha_i\), with
\(\alpha_i(u_i)=0\), gives

\[
 \tau=\alpha_1(u_2)\alpha_2(u_1).                        \tag{1.0a}
\]

In two dimensions, \(\ker\alpha_i=\langle u_i\rangle\).
Consequently

\[
 \boxed{\tau=0\quad\Longleftrightarrow\quad
   N_1,N_2\text{ have a common flag}
   \quad\Longleftrightarrow\quad N_1,N_2\text{ are proportional}.} \tag{1.0b}
\]

Moreover

\[
 \det\{(I+N_1)(I+N_2)\}=1,qquad
 \operatorname {tr}\{(I+N_1)(I+N_2)\}=2+\tau,            \tag{1.0c}
\]

so the discriminant is \(\tau(\tau+4)\).  The pair is positive
hyperbolic for \(\tau>0\), elliptic for \(-4<\tau<0\), and negative
hyperbolic for \(\tau<-4\).  Thus, after unipotence is established, a
common-frame interval calculation can certify hyperbolicity through
\(\tau\ge\tau_0>0\) (or \(\tau\le-4-\tau_0\)).  Without that premise the
correct test is directly
\(|\operatorname{tr}(\Phi_2\Phi_1)|>2\).  A common-flag identity plus an
explicit subcritical growth theorem for the complete admissible Kelvin
class would trigger the pre-registered PPRG no-go; merely excluding this
two-block hyperbolic witness would not.

Put

\[
 N_U=U-I=\begin{pmatrix}0&1\\0&0\end{pmatrix},\qquad
 N_V=V-I=\begin{pmatrix}0&0\\1&0\end{pmatrix}.           \tag{1.1}
\]

Then

\[
 N_U^2=N_V^2=0,qquad \det U=\det V=1.                  \tag{1.2}
\]

The unique invariant flags of the two nonidentity unipotents are

\[
                    \ker N_U=\langle e_1\rangle,qquad
                    \ker N_V=\langle e_2\rangle.         \tag{1.3}
\]

They are distinct.  Nevertheless each block and its inverse has a fixed
explicit Euclidean bound:

\[
 \|U^{\pm1}\|_2=\|V^{\pm1}\|_2
       ={1+\sqrt5\over2}<{13\over8}.                    \tag{1.4}
\]

The product

\[
 P=UV=\begin{pmatrix}2&1\\1&1\end{pmatrix},qquad
 \det P=1,qquad \operatorname {tr}P=3                 \tag{1.5}
\]

has eigenvalues

\[
                    \lambda_\pm={3\pm\sqrt5\over2}.     \tag{1.6}
\]

For a completely rational load-bearing lower bound, take

\[
                         w=(13/8,1)^T.                   \tag{1.7}
\]

Direct multiplication gives

\[
 Pw=(17/4,21/8)^T\ge {34\over13}w                    \tag{1.8}
\]

componentwise.  Therefore

\[
              \rho(P)\ge {34\over13}.                   \tag{1.9}
\]

For \(x>1\), with \(z=(x-1)/(x+1)\),

\[
 \log x=2\left(z+{z^3\over3}+{z^5\over5}+\cdots\right).
                                                               \tag{1.10}
\]

At \(x=34/13\), \(z=21/47\), and exact rational arithmetic gives

\[
 2\left({21\over47}+{1\over3}\left({21\over47}\right)^3
          +{1\over5}\left({21\over47}\right)^5\right)
       ={24\over25}+{1185042\over5733625175}>{24\over25}. \tag{1.11}
\]

Thus the exponent per individual episode satisfies

\[
 \boxed{\gamma={1\over2}\log\rho(P)>{12\over25}.}       \tag{1.12}
\]

After \(r\) two-episode cycles, the Perron vector grows by at least
\((34/13)^r\).  In particular, using natural logarithms,

\[
 r\ge\left\lceil{25\over64}\log q\right\rceil
       \Longrightarrow \|P^rw\|\ge q^{3/8}\|w\|,       \tag{1.13}
\]

and

\[
 r\ge\left\lceil{25\over48}\log q\right\rceil
       \Longrightarrow \|P^rw\|\ge q^{1/2}\|w\|.       \tag{1.14}
\]

Consequently the C182 power deficit is compatible with a bounded
noncommuting two-polarization cocycle over \(O(\log q)\) episodes.  The
algebra does not prove that the passive PDE supplies those episodes.

## 2. Explicit robustness box

Let \(\widetilde U,\widetilde V\in SL(2,\mathbb R)\) and assume that every
entry differs from the corresponding entry in (0.1) by at most

\[
                              \varepsilon={1\over100}.    \tag{2.1}
\]

The trace of a two-by-two product is a sum of four scalar products.  Three
base products in \(\operatorname {tr}(UV)\) have both factors equal to one;
the fourth has both factors zero.  Hence

\[
 |\operatorname {tr}(\widetilde U\widetilde V)-3|
   \le6\varepsilon+4\varepsilon^2={151\over2500}.        \tag{2.2}
\]

It follows that

\[
 \operatorname {tr}(\widetilde U\widetilde V)
       \ge {7349\over2500}>{29\over10}.                  \tag{2.3}
\]

An \(SL(2,\mathbb R)\) matrix with trace greater than \(29/10\) has
expanding eigenvalue greater than \(5/2\), because

\[
                  {5\over2}+{2\over5}={29\over10}.       \tag{2.4}
\]

Using (1.10) at \(x=5/2\), now with \(z=3/7\),

\[
 \log{5\over2}>
 2\left({3\over7}+{1\over3}\left({3\over7}\right)^3\right)
 ={312\over343}>{9\over10}.                             \tag{2.5}
\]

Therefore the robust per-episode exponent obeys

\[
                         \boxed{\widetilde\gamma>{9\over20}.} \tag{2.6}
\]

This margin is deliberately much wider than a floating eigenvalue check.
It is suitable for a future interval enclosure of two actual episode maps.

## 3. Exact passive 2D3C rotation of the gradient

On \(\mathbb T^2\), define

\[
                 v(x,y)=(-\sin y,\sin x).                 \tag{3.1}
\]

It is divergence free.  Its scalar vorticity and pressure are

\[
 \omega=\cos x+\cos y,qquad p=-\cos x\cos y.             \tag{3.2}
\]

Direct calculation gives

\[
              v\cdot\nabla\omega=0,qquad
              (v\cdot\nabla)v+\nabla p=0,                \tag{3.3}
\]

so \(v\) is a smooth steady unforced two-dimensional Euler solution.  Let
\(\Phi_t\) be its flow and prescribe

\[
       \Theta_0(x,y)=\sin x,qquad
       \Theta(t,x,y)=\Theta_0(\Phi_{-t}(x,y)).             \tag{3.4}
\]

Then

\[
        \partial_t\Theta+v\cdot\nabla\Theta=0,qquad
        \|\Theta(t)\|_\infty=1.                          \tag{3.5}
\]

Consequently

\[
                         {\cal U}=(v_1,v_2,\Theta)        \tag{3.6}
\]

is an exact smooth bounded unforced inviscid 2D3C Euler solution.

The origin is a fixed trajectory.  There

\[
 \nabla v(0)=\begin{pmatrix}0&-1\\1&0\end{pmatrix}=:J,
 \qquad D(t)=e^{tJ}=:R_t.                                \tag{3.7}
\]

The exact passive-gradient law gives

\[
        g(t)=\nabla\Theta(t,0)=D(t)^{-T}e_1=R_te_1.       \tag{3.8}
\]

In particular,

\[
                         g(0)=e_1,qquad
                         g(\pi/2)=e_2.                   \tag{3.9}
\]

Thus scalar \(L^\infty\) passivity does not force its gradient to preserve
one fixed line.  This says nothing by itself about invariant lines of the
pressure-coupled Kelvin polarization connection.  Equations (3.1)--(3.9)
are an actual PDE example, but the full Kelvin maps over the two
quarter-period episodes have not been shown to lie in the robust boxes
around \(U,V\).

## 4. The conditional common-conjugacy boundary

There is a precise additional hypothesis under which a common flag *is*
forced.  Let

\[
 C=I+sE_{12},\quad s\ne0,\qquad
 U_j=C^{-j}U_0C^j,\qquad U_0=I+M,                         \tag{4.1}
\]

where

\[
 M=\begin{pmatrix}a&b\\c&-a\end{pmatrix},\qquad M^2=0. \tag{4.2}
\]

Exact multiplication gives

\[
 C^{-j}MC^j=
 \begin{pmatrix}
 a-jsc&b+2jsa-j^2s^2c\\
 c&-a+jsc
 \end{pmatrix}.                                         \tag{4.3}
\]

If these entries are uniformly bounded for every integer \(j\ge0\), the
quadratic coefficient forces \(c=0\), and the linear coefficient then
forces \(a=0\).  Hence \(M=bE_{12}\), and all \(U_j\) have the common flag
\(\langle e_1\rangle\).

There is also an explicit finite-window version.  Suppose \(J\) is even
and the entries in (4.3) have absolute value at most \(B\) at
\(j=0,J/2,J\).  Taking the second finite difference of the upper-right
entry and then its first difference gives

\[
                   |c|\le {8B\over s^2J^2},qquad
                   |a|\le {5B\over |s|J}.                \tag{4.4}
\]

Thus a genuinely common passive conjugacy produces quantitative flag
alignment over long bounded windows.

C183 proves polynomial return formulas for \(g_j\) and \(p_j\).  It does
not prove (4.1) for the complete pressure/polarization monodromy.  Applying
(4.3)--(4.4) to PPRG without first deriving that conjugacy would therefore
insert the desired no-go as an extra premise.

## 5. Architecture verdict and next admissibility test

The pre-registered negative trigger does not fire at this conditional
algebraic level: bounded candidate unipotents need not share a flag.  The
actual PPRG test must use the Kelvin maps themselves:

> On one smooth bounded unforced passive 2D3C orbit, enclose two consecutive
> return-fiber polarization maps \(\Phi_1,\Phi_2\in SL(2,\mathbb R)\) in the
> respective \(1/100\) boxes of Section 2, or directly prove
> \(|\operatorname{tr}(\Phi_2\Phi_1)|\ge2+\delta\) for explicit
> \(\delta>0\), while satisfying C183's affine covector resonance and a
> witness-specific finite-frequency PDE error bound through the required
> logarithmic window.

A rigorous common-flag/subcritical-growth theorem or exhaustive Lyapunov
upper bound for that full admissible class invokes `FRONTIER.md`'s
architecture trigger; it must not be replaced by a sixth named gate.
Failure of one candidate search, or exclusion of only the canonical boxes,
is not such a theorem.  Box success gives the explicit exponent (2.6);
afterward the result must be inserted into the UVSR residual rather than
credited as a standalone endpoint.

## 6. Verification boundary

The checker verifies with exact integer and rational arithmetic:

1. nilpotence, determinants, distinct flags, product, and rational cone
   growth for (0.1);
2. the logarithmic-series bounds (1.11) and (2.5);
3. the \(q^{3/8}\) and \(q^{1/2}\) logarithmic schedules;
4. the entrywise robustness estimate (2.2)--(2.6);
5. the steady Euler, passive transport, and rotating-gradient identities in
   Section 3 at the level of exact symbolic coefficients; and
6. the conjugacy and finite-difference identities in Section 4.

It does not produce the two canonical matrices from the PDE example, bound
a finite-frequency propagator, prove that C183's Kelvin blocks are
unipotent, or verify viscosity, localization, reverse edges, depletion, the
full UVSR residual, or a Navier--Stokes singularity.
