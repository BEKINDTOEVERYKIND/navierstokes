# An explicit AO profile with a certified global `b` maximum

**Date:** 2026-08-02
**Status:** exact self-contained calculus certificate, with dependency-free
rational interval checker; not independently cross-audited
**Scope:** the straight Euler vortex column and its principal-symbol/BAS
growth edge.  This note does not prove a uniform semigroup bound, a curved
vortex-ring quasimode, or a Navier--Stokes transition.

## 1. Why this closes one previously listed gate

Albritton--Ożański's Assumption A does not itself require that the selected
radius maximize `b` globally.  Their Appendix A.1 proves the stronger fact
for the Batchelor vortex family.  The calculation below selects one member
of that family by exact formulas and repeats the global-maximality argument
in a form suitable for the BAS comparison.

The primary source is D. Albritton and W. Ożański, *Linear and nonlinear
instability of vortex columns*, arXiv:2310.20674v3, Assumption A, equations
(1.10)--(1.16), (1.28), (1.31), and Appendix A.1:

<https://arxiv.org/abs/2310.20674v3>

## 2. Exact profile and parameters

Put

\[
 X:=\log {9\over5},\qquad r_0:=\sqrt X,
 \qquad \beta:={1\over\sqrt{4-X}},                 \tag{2.1}
\]

and define the swirl parameter

\[
 Q:={\beta X^2\over 4/5-X}.                         \tag{2.2}
\]

Consider the Batchelor column

\[
 V(r)={Q\over r}(1-e^{-r^2}),\qquad W(r)=e^{-r^2}.
                                                               \tag{2.3}
\]

Writing `x=r^2`, `y=e^{-x}`, and `Gamma=rV`, its basic quantities are

\[
 \Gamma=Q(1-y),\qquad
 \Omega={Q(1-y)\over x},\qquad
 q(r)=-{\Gamma'(r)\over W'(r)}=Q,                  \tag{2.4}
\]

where the quotient extends continuously at the origin.  The functions
`Omega,W` are smooth even functions of `r`; all their odd derivatives vanish
at zero and all derivatives tend to zero at infinity.  Thus the regularity
part of Assumption A holds and `0<inf q=sup q=Q<infinity`.

The elementary bounds

\[
 {4\over7}<X<{3\over5}                              \tag{2.5}
\]

will be useful.  The lower bound follows from
`log(1+t)>2t/(2+t)` at `t=4/5`; the upper bound follows from
`e^(3/5)>1+3/5+(3/5)^2/2+(3/5)^3/6=227/125>9/5`.
In particular, `4/5-X>0`.  Moreover

\[
 \beta Q={X^2\over(4-X)(4/5-X)}<1,                 \tag{2.6}
\]

because

\[
 (4-X)(4/5-X)-X^2={8\over5}(2-3X)>0.               \tag{2.7}
\]

## 3. The AO stationarity conditions

For this column,

\[
 \Lambda(x)=\beta e^{-x}-Q{1-e^{-x}\over x},       \tag{3.1}
\]

and direct differentiation shows that `Lambda'(r)=0` away from the axis
exactly when

\[
 e^x=1+x+{\beta\over Q}x^2.                         \tag{3.2}
\]

By (2.2), `(beta/Q)X^2=4/5-X`; hence (3.2) at `x=X`
is the identity `9/5=1+X+4/5-X`.

To verify the minimum and the uniqueness of its level, set

\[
 c={\beta\over Q}={4/5-X\over X^2},\qquad
 F(x)=e^x-1-x-cx^2.                                 \tag{3.3}
\]

Equation (2.5) gives `c>1/2`.  Since `F(0)=F'(0)=0` and
`F''(x)=e^x-2c` is strictly increasing from a negative value, `F'`
first decreases and then increases through zero exactly once.  Therefore
`F` decreases and then increases, and it has exactly one positive zero.
That zero is `X`.  Also

\[
 {d\Lambda\over dx}={Qe^{-x}\over x^2}F(x),        \tag{3.4}
\]

so `Lambda` has a unique strict global minimum at `r=r_0`.  In particular,
there is no other radius with `Lambda(r)=Lambda(r_0)`.  At the minimum,

\[
 \Lambda''(r_0)
 ={8\beta\over9}{7X-4\over4/5-X}>0,                \tag{3.5}
\]

where positivity follows from `X>4/7`.

Next,

\[
 b(x)=4\beta Q(1-\beta Q)
 {e^{-x}(1-e^{-x})\over1+\beta^2x}.                 \tag{3.6}
\]

The prefactor is positive by (2.6).  For an arbitrary `a=beta^2>0`, the
sign of the derivative of the last fraction is the sign of

\[
 H_a(x)=e^{-x}[2(1+ax)+a]-(1+ax+a).                 \tag{3.7}
\]

Its derivative is

\[
 H_a'(x)=e^{-x}(a-2-2ax)-a<0.                       \tag{3.8}
\]

Indeed, if the parenthesis is negative this is immediate; if it is
non-negative, multiplication by `e^{-x}<=1` still makes the result at most
`a-2-2ax-a<0`.  Thus `H_a` decreases strictly from `H_a(0)=1` to a negative
limit and has one zero.  With `a=1/(4-X)` and `e^{-X}=5/9`,

\[
 (4-X)H_a(X)={5\over9}(8+1)-(4+1)=0.               \tag{3.9}
\]

Consequently `r_0` is the unique global maximizer of `b`.  This supplies
all of the geometric conditions in AO Assumption A, including the stronger
global statements used in their Appendix A.1.

## 4. Certified constants and the fixed-sector resonant BAS ratio

At the selected radius,

\[
 b_0={32X^2(2-3X)\over
       81(4-X)(4/5-X)^2}.                            \tag{4.1}
\]

The rational-interval checker accompanying this note certifies

\[
\begin{aligned}
 0.58778665&<X<0.58778667,\\
 0.76667&<r_0<0.76668,\\
 0.54135&<\beta<0.54136,\\
 0.88134&<Q<0.88136,\\
 0.21018&<b_0<0.21020,\\
 0.45845&<\sqrt{b_0}<0.45848.
\end{aligned}                                       \tag{4.2}
\]

For a resonant helical bicharacteristic with radial covector `ell`, the
column BAS reduction is

\[
 \sigma^2(\ell;r)=b(r)
 {n^2/r^2+\alpha^2\over
  \ell^2+n^2/r^2+\alpha^2}.                         \tag{4.3}
\]

Equations (3.6)--(3.9) imply `sigma^2<=b(r)<=b_0`, while equality is
attained at `(r,ell)=(r_0,0)`.  Hence the global resonant BAS exponent at
this fixed helical ratio is exactly

\[
 \Lambda_{\rm BAS}=\sqrt{b_0}.                      \tag{4.4}
\]

AO's fixed-`m` ring eigenvalues therefore satisfy

\[
 {\operatorname{Im}\omega_{m,n}\over\Lambda_{\rm BAS}}
 =1-D_m n^{-1/2}+O(n^{-1+\delta}),                  \tag{4.5}
\]

where

\[
 D_m=(2m-1)
 \left({\Lambda''(r_0)\over8p_0\sqrt{b_0}}\right)^{1/2},
 \qquad p_0={4\over X(4-X)}.                        \tag{4.6}
\]

For the ground state, the checker gives

\[
             0.1883<D_1<0.1885.                    \tag{4.7}
\]

Thus the asymptotic spectral/resonant-BAS ratio is exactly one, with an
explicit first correction.  AO's remainder constant is not effective, so (4.5)
does not yield a certified numerical wavenumber at which a prescribed
finite-`n` ratio is crossed.

## 5. What remains open

1. The profile (2.3) is smooth and decaying but not compactly supported:
   `V(r)~Q/r`.  AO Assumption A does not require compact support.  An exact
   compactification preserving the simultaneous critical point and global
   inequalities is a separate localization problem.
2. The resonant BAS maximum is a fixed-helical-ratio certificate.  The full
   reduction in `2026-08-02-ao-batchelor-full-bas-cocycle.md` proves that
   non-resonant orbits have zero Lyapunov exponent, but refutes unrestricted
   edge matching by finding a faster resonance at another ratio.  In any
   event, the principal-symbol calculation does not supply the Fourier-sector
   semigroup estimate with a frequency-uniform polynomial prefactor needed by
   the long-gain ledger.  AO explicitly leave mode dependence in the constants
   following their equation (3.32).
3. No estimate here controls the residual created by bending the straight
   column into a finite-curvature ring or coupling it to the viscous
   transition.

The first previously open item in C54--global maximization of `b` for an
admissible explicit AO profile--is therefore closed.  The uniform
propagator and curvature obligations remain load-bearing.
