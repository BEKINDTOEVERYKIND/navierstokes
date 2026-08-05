# Positive-stress charts: explicit three-phase formula and two-phase repair

**Date:** 2026-08-01
**Status:** corrected after independent audit; exact principal algebra
**Scope:** leading averaged covariance only.  No multiphase evolution or
Navier--Stokes transition is constructed.

**Correction (2026-08-02).**  The original note claimed that three phases
were minimal near isotropy.  That conclusion fixed the two kernel directions
while taking the tangent space.  Pointwise directions are algebraic variables.
Once their variations are included, the two-phase sum map is a submersion
at a positive decomposition of \(qI\).  Two pointwise directions therefore
admit a smooth local positive covariance chart; they are not thereby realized
as curl-free gradients of transported material phases.

## 1. Corrected result

A single incompressible fast phase with gradient \(k\) has amplitudes in
\(k^\perp\), so its averaged covariance belongs to

\[
 \mathcal S(k)=\{R\in\operatorname{Sym}_3:Rk=0\},
 \qquad \dim\mathcal S(k)=3.                             \tag{1.1}
\]

This note sharpens the one-carrier obstruction, but not to three phases.

> **Positive covariance chart near isotropy.**  One phase is impossible,
> while two variable pointwise kernel directions are algebraically sufficient
> and hence minimal for covariance decomposition.  Three fixed orthogonal
> phase gradients give a particularly
> simple explicit formula.  In fact, every symmetric matrix \(Q\) with
> \(\|Q-qI\|_{\rm op}<q/4\) has an explicit decomposition
> \[
>                         Q=R_1+R_2+R_3,                  \tag{1.2}
> \]
> where \(R_i\) is positive definite on \(e_i^\perp\) and
> \(R_ie_i=0\).

The two-phase local section below is smooth but implicit.  The three-phase
formula is explicit and may still be convenient.  Either choice has
cross-colour interactions at leading fast frequency, so neither replaces
the forward material-phase theorem.

## 2. What is true for two fixed phase directions

Let \(k,l\) be nonparallel.  Since

\[
 \dim(\mathcal S(k)\cap\mathcal S(l))=1,
\]

the linear sum has dimension

\[
              \dim(\mathcal S(k)+\mathcal S(l))=3+3-1=5. \tag{2.1}
\]

The trace-free stress space also has dimension five, so dimension counting
alone does not rule out two generic nonorthogonal phases.  Positivity around
isotropy supplies the missing obstruction.

Suppose an isotropic baseline decomposes as

\[
                  cI=A+B,\qquad
                  A\in\mathcal S(k),\quad B\in\mathcal S(l),
                  \quad c\ne0.                            \tag{2.2}
\]

Using symmetry, \(Ak=0\), and \(Bl=0\),

\[
 c\,k\cdot l
 =k\cdot(A+B)l
 =(Ak)\cdot l+k\cdot(Bl)=0.                              \tag{2.3}
\]

Therefore \(k\cdot l=0\) is necessary.  When \(k\perp l\), the isotropic
matrix does lie in the five-dimensional sum (2.1).  Projection of that sum
to \(\operatorname{Sym}_0^3\) then has the isotropic line as kernel, and
hence rank only four.  Thus, **with the directions held fixed**:

* nonorthogonal phases cannot supply a nonzero isotropic positive baseline;
* orthogonal phases have such a baseline but only four trace-free tangent
  directions after the phase directions are artificially frozen.

The final conclusion in the original version was false: it omitted tangent
directions obtained by rotating the kernels.

## 2A. Two variable pointwise directions give a local smooth chart

Fix \(0<a<q\) and take

\[
 k_0=e_1,\qquad \ell_0=e_2,
\]

\[
 A_0=\operatorname{diag}(0,q,a),\qquad
 B_0=\operatorname{diag}(q,0,q-a).                    \tag{2.4}
\]

Then \(A_0+B_0=qI\); both matrices are positive definite on the plane
orthogonal to their respective kernels.  Parameterize \(A\) by its positive
\(2\times2\) block and a nearby unit kernel \(k\), and parameterize \(B\)
analogously using \(\ell\).  Let

\[
                    F(A,k,B,\ell)=A+B.                \tag{2.5}
\]

At (2.4), fixed-block variations already supply the five symmetric
coordinates

\[
 E_{11},\ E_{22},\ E_{33},\ E_{13}+E_{31},\
 E_{23}+E_{32}.                                      \tag{2.6}
\]

If \(\delta k=xe_2\), differentiating \(Ak=0\) gives

\[
             \delta A_{12}=\delta A_{21}=-qx.         \tag{2.7}
\]

This supplies the missing \(E_{12}+E_{21}\) coordinate.  Hence \(DF\) has
rank six, the full dimension of \(\operatorname{Sym}_3\).  Selecting these
six domain coordinates gives an invertible minor, so the implicit-function
theorem provides a smooth local right inverse of \(F\) near \(qI\).
Positivity on the two transverse planes persists by openness.

Consequently every \(Q\) in some neighborhood of \(qI\) has

\[
 Q=A(Q)+B(Q),\qquad A(Q)k(Q)=0,\quad
 B(Q)\ell(Q)=0,                                      \tag{2.8}
\]

with both covariances uniformly positive on their transverse planes.  Since
one phase has rank at most two, two variable pointwise directions are minimal
for this algebraic covariance problem.  For a spatial field `Q(x)`, the
selected `k(Q(x))` and `ell(Q(x))` need not be curl-free and hence need not be
gradients of scalar phases.

For pointwise existence without a smooth choice, diagonalize a positive
definite \(Q=\sum_i\lambda_i e_i\otimes e_i\) and choose \(0<t<1\):

\[
 A=\lambda_2e_2\otimes e_2+t\lambda_3e_3\otimes e_3,
 \qquad
 B=\lambda_1e_1\otimes e_1+(1-t)\lambda_3e_3\otimes e_3. \tag{2.9}
\]

Then \(\ker A=\operatorname{span}(e_1)\) and
\(\ker B=\operatorname{span}(e_2)\).  Eigenvectors cannot be chosen
smoothly through isotropy, which is why the submersion argument, rather
than (2.9), is needed for the local chart claim.

## 3. Explicit fixed-direction three-phase decomposition

Write \(Q=(Q_{ij})=Q^T\).  Define

\[
R_1=
\begin{pmatrix}
0&0&0\\
0&Q_{22}/2&Q_{23}\\
0&Q_{23}&Q_{33}/2
\end{pmatrix},                                           \tag{3.1}
\]

\[
R_2=
\begin{pmatrix}
Q_{11}/2&0&Q_{13}\\
0&0&0\\
Q_{13}&0&Q_{33}/2
\end{pmatrix},
\qquad
R_3=
\begin{pmatrix}
Q_{11}/2&Q_{12}&0\\
Q_{12}&Q_{22}/2&0\\
0&0&0
\end{pmatrix}.                                           \tag{3.2}
\]

Then

\[
 R_ie_i=0,\qquad R_1+R_2+R_3=Q.                         \tag{3.3}
\]

If \(\|Q-qI\|_{\rm op}<q/4\), every diagonal entry of \(Q\) is greater
than \(3q/4\), while every off-diagonal entry has magnitude less than
\(q/4\).  Each nonzero \(2\times2\) block in (3.1)--(3.2) therefore has
diagonal entries greater than \(3q/8\) and off-diagonal magnitude less
than \(q/4\).  Gershgorin's bound gives its smallest eigenvalue greater
than \(q/8\).  Hence each \(R_i\) is uniformly positive on
\(e_i^\perp\).

At the exact isotropic baseline,

\[
                 R_i={q\over2}(I-e_i\otimes e_i),         \tag{3.4}
\]

and the three matrices sum to \(qI\).

## 4. Realization by transverse oscillations

Let the two positive eigenpairs of the restriction of \(R_i\) to
\(e_i^\perp\) be \((\lambda_{i,1},v_{i,1})\) and
\((\lambda_{i,2},v_{i,2})\).  The mean-zero one-phase amplitude

\[
 a_i(\theta)
 =\sqrt{2\lambda_{i,1}}v_{i,1}\cos\theta
  +\sqrt{2\lambda_{i,2}}v_{i,2}\sin\theta               \tag{4.1}
\]

satisfies

\[
 e_i\cdot a_i=0,qquad
 {1\over2\pi}\int_0^{2\pi}a_i(\theta)\otimes a_i(\theta)\,d\theta
 =R_i.                                                    \tag{4.2}
\]

Thus six real polarization coordinates realize this explicit chart.  The
averaged covariance is exactly \(Q\), and incompressibility holds at
principal order for every colour.  A two-phase realization uses four real
polarization coordinates plus the two variable kernels.  Calling those
kernels phase gradients requires a separate integrability and transport
construction.

## 5. Relation to the active-transition program

The result closes only the positive principal-symbol algebra left open by
the invalid one-carrier proposal.  The minimal local algebra uses two
variable pointwise directions; the displayed three-phase formula remains a
convenient fixed-direction alternative.

It does not address four dynamic issues.

1. Coexisting phases produce cross-colour sum and difference frequencies.
2. The IFT-selected direction fields need not be curl-free gradients.  Even
   if their integrability is arranged, material transport and lower-order
   divergence/pressure correctors remain.
3. Spatially separating the colours removes cross interactions but also
   prevents their pointwise covariances from adding without a further
   intermittent partition/gluing argument.
4. The endpoint must retain the nonlocal viscous pressure wake and meet the
   \(C^M(M!)^2\) majorant through \(M\asymp j^2/\log j\).

The next useful theorem should decide dynamically between a genuine
integrable/transported realization of the two-direction submersion
(2.4)--(2.8), the explicit three-phase formula (3.1)--(4.2), and
a sequential one-colour protocol.  Static covariance rank alone no longer
selects the architecture.
