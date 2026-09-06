# C190: the C186 rotating-gradient orbit cannot realize the two-block PPRG witness

**Date:** 2026-08-25

**Status:** exact orbit-specific principal-Kelvin obstruction; this is
outcome (b) of the pre-registered PPRG witness test, not a no-go for the
complete admissible passive 2D3C class

**Checker:**
[checks/rotating_gradient_return_obstruction_c190.py](../checks/rotating_gradient_return_obstruction_c190.py)

## 0. Result and claim boundary

C186 left one explicit smooth bounded unforced passive 2D3C orbit as the
chosen realization candidate:

\[
 v(x,y)=(-\sin y,\sin x),\qquad
 \Theta_0(x,y)=\sin x,
 \tag{0.1}
\]

with \(\Theta\) passively transported by \(v\).  Along the fixed trajectory
at the origin, its passive gradient rotates from \(e_1\) to \(e_2\) in one
quarter period.  The pre-registered test asked for two consecutive Kelvin
polarization maps in one C183 return frame and fixed success as either

\[
 \Phi_1\in B_{1/100}(U),\quad
 \Phi_2\in B_{1/100}(V),
 \qquad
 U=I+E_{12},\quad V=I+E_{21},
 \tag{0.2}
\]

entry by entry, or

\[
             |\operatorname {tr}(\Phi_2\Phi_1)|\ge 2+\delta
             \quad\hbox{for an explicit }\delta>0.
\tag{0.3}
\]

Neither criterion is available on this orbit.  The obstruction is exact:

1. a genuine coefficient return forces the conserved vertical Kelvin
   charge to be \(m=0\); and
2. on that only return-coherent branch every episode map belongs to the
   same one-parameter unipotent subgroup.

For the explicit returning covector selected below, the two consecutive
quarter-episode maps in the canonical co-rotating orthonormal fiber frame
are

\[
 \boxed{\Phi_1=\Phi_2=I+{\pi\over2}E_{12}.}
\tag{0.4}
\]

The genuine full-return block is \(S=I+2\pi E_{12}\), and consecutive
full-return blocks are also identical.

Consequently

\[
 (\Phi_i-I)^2=0,\qquad
 \operatorname {tr}(\Phi_2\Phi_1)=2,
 \qquad
 \|S^N\|_2\le 1+{44\over7}N.
\tag{0.5}
\]

Thus the square-zero gate is kept honestly, but it yields a common flag and
zero Lyapunov exponent rather than the C186 hyperbolic pair.

This is the theorem-grade negative outcome for the **chosen C186 orbit**.
It does not prove that all passive 2D3C Kelvin cocycles have a common flag,
does not fire the architecture-wide PPRG trigger in `FRONTIER.md`, and does
not assert a finite-frequency, viscous, nonlinear, UVSR, or singularity
theorem.  No successor gate is introduced.

## 1. Exact coefficient return and the vertical-charge obstruction

Put

\[
 J=\begin{pmatrix}0&-1\\1&0\end{pmatrix},\qquad
 R_t=e^{tJ},\qquad n=e_3.
\tag{1.1}
\]

At the origin, which is a fixed planar trajectory, C186 gives

\[
 M=\nabla v(0)=J,\qquad D(t)=R_t,\qquad
 g(t)=\nabla\Theta(t,0)=R_te_1.
\tag{1.2}
\]

The full velocity gradient along the orbit is therefore

\[
 A(t)=\begin{pmatrix}J&0\\g(t)^T&0\end{pmatrix}.
\tag{1.3}
\]

For a Kelvin covector \(k=(p,m)\), C183's exact covector formula specializes
to

\[
              \boxed{p(t)=R_t(p_0-mt e_1),\qquad m(t)=m.}
\tag{1.4}
\]

The first positive return of the coefficient pair \((M,g)\) is
\(T=2\pi\).  At that time

\[
                    p(T)-p_0=-mT e_1.
\tag{1.5}
\]

Hence \(k(T)=k(0)\) if and only if \(m=0\).  Allowing only a projective
return does not enlarge the class.  Indeed, if
\(k(T)=\lambda k(0)\) and \(m\ne0\), the unchanged vertical component
forces \(\lambda=1\), after which (1.5) again forces \(m=0\), a
contradiction.  Thus for \(m\ne0\) the endpoint planes
\(k(T)^\perp\) and \(k(0)^\perp\) are distinct, and a matrix obtained by an
arbitrary identification of them is not the common-return-fiber object
pre-registered in (0.2)--(0.3).

There is a second formulation that directly addresses the two quarter
episodes used in C186.  Put \(h=\pi/2\) and write \(p_0=(a,b)\).  From
(1.4),

\[
 |k(t)|^2=(a-mt)^2+b^2+m^2.
\tag{1.6}
\]

C183's exact area formula gives

\[
 \det\Phi_1={|k(0)|\over|k(h)|},\qquad
 \det\Phi_2={|k(h)|\over|k(2h)|}.
\tag{1.7}
\]

Thus the pre-registered requirement \(\Phi_1,\Phi_2\in SL(2,\mathbb R)\)
forces \(|k(0)|=|k(h)|=|k(2h)|\).  If \(m\ne0\), the first equality gives

\[
                         a={mh\over2},
\tag{1.8}
\]

whereas the second gives

\[
                         a={3mh\over2},
\tag{1.9}
\]

which is impossible.  Hence even the two consecutive quarter maps can
both pass the honest determinant-one gate only when \(m=0\).

The same conclusion is independent of the equal-quarter choice.  For any
three distinct consecutive endpoint times \(t_0<t_1<t_2\), determinant one
for both episode maps forces
\(|k(t_0)|=|k(t_1)|=|k(t_2)|\).  If \(m\ne0\), (1.6) minus this common
value is a nonzero quadratic in \(t\) with leading coefficient \(m^2>0\),
but it would have the three distinct roots \(t_0,t_1,t_2\).  Therefore any
two consecutive determinant-one episodes on this orbit force \(m=0\), not
only the two quarters highlighted by C186.

This is the first half of the obstruction: the vertical charge that could
make the passive shear enter the pressure coupling prevents both a genuine
C183 return and the requested pair of determinant-one quarter blocks on
this particular rotating-gradient orbit.

## 2. Exact polarization map on the only coherent branch

Now set \(m=0\) and choose the unit returning horizontal covector

\[
                     p_0=-e_2,qquad p(t)=R_t(-e_2).
\tag{2.1}
\]

Use the oriented orthonormal return frame

\[
 E_1(t)=n,qquad
 E_2(t)={k(t)\times n\over |p(t)|}=-R_te_1.
\tag{2.2}
\]

It is \(2\pi\)-periodic and satisfies
\(E_1\times E_2=k/|k|\).  Write the Kelvin velocity amplitude as
\(a=E(t)z\).  The exact Kelvin generator is

\[
 K(t)=-A(t)+2{k(t)k(t)^TA(t)\over |k(t)|^2},
 \qquad z'=B(t)z,
 \quad B=E^TKE-E^TE'.
\tag{2.3}
\]

This connection is constant.  To see it without integrating an ODE, put
\(Q_t=\operatorname {diag}(R_t,1)\).  Then

\[
 A(t)=Q_tA_0Q_t^T,\qquad k(t)=Q_tk_0,qquad E(t)=Q_tE_0,
\tag{2.4}
\]

where

\[
 A_0=\begin{pmatrix}0&-1&0\\1&0&0\\1&0&0\end{pmatrix},
 \quad
 k_0=\begin{pmatrix}0\\-1\\0\end{pmatrix},
 \quad
 E_0=\begin{pmatrix}0&-1\\0&0\\1&0\end{pmatrix},
 \quad
 \Omega=Q_t^TQ_t'=
 \begin{pmatrix}0&-1&0\\1&0&0\\0&0&0\end{pmatrix}.
\tag{2.5}
\]

Direct exact multiplication gives

\[
 K_0=-A_0+2k_0k_0^TA_0
 =\begin{pmatrix}0&1&0\\1&0&0\\-1&0&0\end{pmatrix},
\tag{2.6}
\]

and therefore

\[
 \boxed{B=E_0^TK_0E_0-E_0^T\Omega E_0
       =\begin{pmatrix}0&1\\0&0\end{pmatrix}=E_{12}.}
\tag{2.7}
\]

Since \(B^2=0\), the exact propagator over every interval is

\[
                    \Phi(t_1,t_0)=I+(t_1-t_0)E_{12}.
\tag{2.8}
\]

For \(h=\pi/2\), (2.8) gives the two consecutive maps (0.4).  The
coefficient, covector, and frame all repeat after \(T=2\pi\), so every
genuine full-return block is \(S=I+2\pi E_{12}\).

For completeness, an arbitrary unit horizontal returning direction
\(h_0=p_0/|p_0|\) gives the same conclusion with

\[
 B=\beta E_{12},\qquad
 \beta=e_1\cdot Jh_0,\qquad |\beta|\le1.
\tag{2.9}
\]

Thus changing the returning horizontal direction changes only the shear
coefficient, not its invariant flag.

## 3. Evaluation of the pre-registered criteria

### 3.1 The two entrywise boxes

In the displayed common frame, every return block has lower-left entry
zero.  Hence

\[
            |(\Phi_2)_{21}-V_{21}|=|0-1|=1>{1\over100}.
\tag{3.1}
\]

The failure persists under a single fixed change of basis applied at all
three endpoint coordinate spaces.  Such a common change preserves
\(\Phi_1=\Phi_2\).  But the two closed entrywise boxes in (0.2) are
disjoint, since one matrix in both would imply

\[
       1=|U_{12}-V_{12}|\le {1\over100}+{1\over100}
         ={1\over50},
\tag{3.2}
\]

which is false.  Therefore the requested ordered \(U,V\) box pair cannot
be obtained in the canonical C183 co-rotating frame or under a single
fixed common conjugacy.  Independent endpoint gauges are excluded for the
reason stated in Section 3.3.

### 3.2 The direct trace branch and the square-zero gate

For the quarter blocks \(Q=I+(\pi/2)E_{12}\),

\[
 (Q-I)^2=0,
 \qquad
 \Phi_2\Phi_1=Q^2=I+\pi E_{12},
 \qquad
 |\operatorname {tr}(\Phi_2\Phi_1)|=2.
\tag{3.3}
\]

Thus the square-zero premise is proved rather than inferred from the
planar return, but the trace criterion fails for every \(\delta>0\).
For an arbitrary return-coherent horizontal direction, (2.9) instead gives

\[
 Q_\beta=I+{\pi\over2}\beta E_{12},\qquad
 \operatorname {tr}(Q_\beta^2)=2
 \quad (|\beta|\le1),
\tag{3.3a}
\]

so the conclusion is exhaustive over the \(m=0\) covectors, not only the
displayed choice \(p_0=-e_2\).

Moreover, for every integer \(N\ge0\),

\[
 S^N=I+2\pi N E_{12},
 \qquad
 \boxed{\|S^N\|_2\le1+2\pi N\le1+{44\over7}N.}
\tag{3.4}
\]

Here \(\|E_{12}\|_2=1\) and the last inequality uses the elementary
Archimedean bound \(\pi<22/7\).  In particular this return cocycle has
Lyapunov exponent exactly zero.

### 3.3 Quarter-period ambiguity and finite frequency

The gradient rotation highlighted in C186 occurs over quarter periods.
Those quarters are not genuine C183 coefficient returns because \(g\) has
not returned, but (1.6)--(1.9) show that requiring both canonical
orthonormal-fiber maps to have determinant one already forces \(m=0\).
In the natural common co-rotating trivialization (2.2), each quarter
propagator is then the same matrix

\[
                       I+{\pi\over2}E_{12},
\tag{3.5}
\]

so their product still has trace two and a common flag.  Independent
endpoint gauges could manufacture different coordinate matrices, but then
the trace and the boxes would be gauge artifacts rather than the
pre-registered common-fiber test.

Because the exact principal criterion has no positive margin, outcome (a)
is not reached and no finite-frequency error is charged as a success.
This checkpoint does not claim that finite-frequency effects on another
orbit are impossible.  It records the required outcome (b) for this orbit
and leaves the unchanged full admissible-class PPRG dichotomy open.

## 4. Verification boundary

The dependency-free checker verifies with exact rational and formal
polynomial arithmetic:

1. the matrices (2.5)--(2.7), including the moving-frame term;
2. the full-return covector defect \((-mT,0)\);
3. the consecutive-quarter determinant-one contradiction for \(m\ne0\),
   determinant one on \(m=0\), exact square-zero, equality of consecutive
   blocks, product trace two, and the formula for every tested integer
   power;
4. disjointness of the two \(1/100\) boxes; and
5. the explicit \(1+(44/7)N\) return-growth majorant.

It does not enclose a finite-frequency PDE propagator, prove a viscous
estimate, exclude any other passive 2D3C orbit, close the architecture-wide
PPRG branch, construct UVSR, or prove a Navier--Stokes singularity.
