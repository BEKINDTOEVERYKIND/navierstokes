# C118--C120: the finite A2 hexagon, depletion orbits, and an enclosed ladder eigenvalue

**Date:** 2026-08-05  
**Status:** finite-dimensional statements exact and arithmetically checked;
identification with an invariant Navier--Stokes Fourier subsystem is open  
**Scope:** the existing Beltrami/A2 stage geometry only.  This is not a
regularity or blow-up theorem.

## 1. Verdict

The six first A2 sidebands do support an exact energy-conserving pump
normal form.  Every nonzero singular value of its hexagon coefficient matrix
has a nonlinear depletion orbit, not merely a linearly growing mode.  The
three equal leaves admit an explicit heteroclinic orbit which empties the
pump at a finite section.  A one-edge-weighted hexagon has characteristic
polynomial

\[
 q(\lambda)=\lambda^6-9\lambda^4+18\lambda^2-9
\tag{1.1}
\]

and a simple positive eigenvalue \(\sigma_*\) with the rigorous rational
enclosure

\[
 \boxed{\frac{633}{250}<\sigma_*<\frac{2533}{1000}}.
\tag{1.2}
\]

Those are exact ODE facts.  They do **not** yet prove that the corresponding
six Fourier modes are invariant under the full Euler or Navier--Stokes
equation.  Support addition alone does not close: the pump additions which
make a desired hexagon edge also allow a leaf to reach the centre and an
outer second sideband.  This support calculation does not prove those
coefficients are nonzero.  The physical theorem therefore still needs
either an exact polarization cancellation of those outputs or a
quantitative off-ladder slaving estimate.

The exact arithmetic is checked by
[`checks/hexagon_depletion_eigenmode_c118_c120.py`](../checks/hexagon_depletion_eigenmode_c118_c120.py).

## 2. Exact A2 support geometry

Put

\[
 r_1=(1,-1,0),\qquad r_2=(0,1,-1),\qquad
 r_3=(-1,0,1).
\tag{2.1}
\]

Then

\[
 r_1+r_2+r_3=0,\qquad |r_i|^2=2,\qquad
 r_i\mathbin\cdot r_j=-1\quad(i\ne j).
\tag{2.2}
\]

Thus \({\cal R}=\{\pm r_1,\pm r_2,\pm r_3\}\) is the A2 root hexagon.
For integers \(m\) and \(K>0\), put \(q_m=m(1,1,1)\).  Every first leaf

\[
 k_{i,\pm}=q_m\pm K r_i
\tag{2.3}
\]

lies on the same shell,

\[
 |k_{i,\pm}|^2=3m^2+2K^2.
\tag{2.4}
\]

A positive vertex \(q_m+Kr_i\) and a negative vertex
\(q_m-Kr_j\) differ by a pump root precisely when \(i\ne j\).  Hence the
support graph is the bipartite six-cycle with incidence matrix

\[
 B_0=J-I=
 \begin{pmatrix}
 0&1&1\\
 1&0&1\\
 1&1&0
 \end{pmatrix}.
\tag{2.5}
\]

This is a support statement, not a coefficient statement.  Leray
projection and the two transverse polarizations put matrix weights on the
six edges.  Replacing those physical blocks by the scalar entries of
\(B_0\), or by the weighted example below, is the finite-dimensional
normal-form hypothesis whose PDE realization remains to be proved.

## 3. C118: every singular direction has an exact depletion orbit

Let \(B\) be any real \(d\)-by-\(d\) matrix.  The pump--leaf Hamiltonian
normal form is

\[
\begin{aligned}
 \dot p&=-2\sum_{a=1}^N x_a^T B y_a,\\
 \dot x_a&=pB y_a,\\
 \dot y_a&=pB^T x_a.
\end{aligned}
\tag{3.1}
\]

It conserves

\[
 \mathcal E=p^2+\sum_{a=1}^N\bigl(|x_a|^2+|y_a|^2\bigr),
\tag{3.2}
\]

because the pump work cancels the leaf work exactly.

Suppose

\[
 B\eta=\sigma\xi,\qquad B^T\xi=\sigma\eta,
 \qquad |\xi|=|\eta|=1,\qquad \sigma>0.
\tag{3.3}
\]

The submanifold

\[
 x_a=A_a\xi,\qquad y_a=A_a\eta
\tag{3.4}
\]

is invariant and reduces (3.1) to

\[
 \dot p=-2\sigma\sum_{a=1}^N A_a^2,\qquad
 \dot A_a=\sigma p A_a.
\tag{3.5}
\]

All nonzero ratios \(A_a/A_b\) are constant.  With
\(R^2=\sum_aA_a^2\), the entire \(N\)-leaf problem is therefore the
two-dimensional system

\[
 \dot p=-2\sigma R^2,\qquad \dot R=\sigma pR,\qquad
 p^2+2R^2=A^2.
\tag{3.6}
\]

Its nontrivial heteroclinic orbit is

\[
 p(t)=-A\tanh z,\qquad
 R(t)=\frac{A}{\sqrt2}\operatorname{sech}z,\qquad
 z=\sigma A(t-t_0).
\tag{3.7}
\]

Thus the growing eigenmode of the frozen pump is accompanied by a complete
nonlinear pump-depletion trajectory.  Notice the boundary: after the
depletion section the orbit reverses the pump sign and returns the leaf
energy.  It is not by itself a one-way cascade stage.

## 4. C119: the equal three-leaf orbit

For \(B_0\), let \({\bf 1}=(1,1,1)^T\).  Since
\(B_0{\bf1}=2{\bf1}\), the symmetric hexagon state

\[
 x=y=a(t){\bf1}
\tag{4.1}
\]

obeys

\[
 \dot p=-12a^2,\qquad \dot a=2pa,\qquad
 p^2+6a^2=A^2.
\tag{4.2}
\]

The exact orbit is

\[
 \boxed{
 p(t)=-A\tanh(2A(t-t_0)),\qquad
 a(t)=\frac{A}{\sqrt6}\operatorname{sech}(2A(t-t_0)).}
\tag{4.3}
\]

At \(t=t_0\), \(p=0\) and all the normal-form energy is in the three
equal leaf pairs.  This proves depletion inside the finite hexagon model.
It does not prove that the three physical pump Fourier coefficients remain
locked to one scalar \(p(t)\) after leaf feedback.

## 5. C120: one weighted hexagon eigenmode with a rational enclosure

To ensure that the instability is not an artefact of the double singular
value of \(B_0\), strengthen one hexagon edge and set

\[
 B_*=
 \begin{pmatrix}
 0&1&1\\
 1&0&2\\
 1&1&0
 \end{pmatrix},\qquad
 A_*=B_*B_*^T=
 \begin{pmatrix}
 2&2&1\\
 2&5&1\\
 1&1&2
 \end{pmatrix}.
\tag{5.1}
\]

The inviscid frozen-pump leaf matrix is

\[
 L_*=
 \begin{pmatrix}0&B_*\\B_*^T&0\end{pmatrix}.
\tag{5.2}
\]

Direct expansion gives

\[
 \det(\mu I-A_*)=p(\mu)=\mu^3-9\mu^2+18\mu-9,
 \qquad
 \det(\lambda I-L_*)=p(\lambda^2),
\tag{5.3}
\]

which is (1.1).  Since \(p(5)<0\), \(p'(\mu)>0\) for
\(\mu\ge5\), and

\[
 p\!\left((633/250)^2\right)<0,\qquad
 p\!\left((2533/1000)^2\right)>0,
\tag{5.4}
\]

there is exactly one root \(\mu_*>5\), and (1.2) holds for
\(\sigma_*=\sqrt{\mu_*}\).  An exact eigenvector of \(A_*\), evaluated
at any root \(\mu\) of \(p\), is

\[
 \xi(\mu)=
 \begin{pmatrix}
 \mu-3\\ \mu\\ \mu^2-7\mu+6
 \end{pmatrix}.
\tag{5.5}
\]

Indeed the first two rows of \((A_* -\mu I)\xi\) vanish identically and
the last row equals \(-p(\mu)\).  After normalizing \(\xi\), taking
\(\eta=B_*^T\xi/\sigma_*\) gives \(|\eta|=1\) and produces an
eigenvector \((\xi,\eta)\) of \(L_*\) with positive eigenvalue
\(\sigma_*\).  The corresponding C118 orbit is exact.

For comparison, the unweighted matrix has the complete exact spectrum

\[
 \operatorname{spec}
 \begin{pmatrix}0&B_0\\B_0&0\end{pmatrix}
 =\{2,-2,1,1,-1,-1\},
\tag{5.6}
\]

and characteristic polynomial
\((\lambda^2-4)(\lambda^2-1)^2\).

## 6. The exact obstruction exposed by the same geometry

The six-cycle does not close under pump addition.  For every root \(r_i\),

\[
 (q_m+Kr_i)-Kr_i=q_m,\qquad
 (q_m+Kr_i)+Kr_i=q_m+2Kr_i.
\tag{6.1}
\]

Both outputs lie outside the six first leaves.  They are allowed after one
pump interaction, exactly the same support-interaction order as the desired
edge

\[
 (q_m+Kr_i)+Kr_j=q_m-Kr_k,\qquad \{i,j,k\}=\{1,2,3\}.
\tag{6.2}
\]

Support separation therefore supplies no small parameter.  A full one-cell
theorem must prove one of the following on the actual helical polarizations:

1. the centre and second-sideband coefficients cancel exactly;
2. they are slaved to the unstable leaf with a quantitatively smaller
   Duhamel norm; or
3. they are included in a larger invariant/slaved state whose return map
   retains the required gain.

Until one of these is established, C118--C120 are exact finite-dimensional
modules and not a Navier--Stokes stage map.
