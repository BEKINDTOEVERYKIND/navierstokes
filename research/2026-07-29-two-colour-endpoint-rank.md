# Two-colour Kelvin beats: an exact child-rank calculation and the minimum three-beat repair

## Status

This note tests the smallest genuinely non-collinear Fourier/Kelvin
transition suggested by the one-phase cell audit.  It does **not** prove a
Navier--Stokes singularity.

There are four exact conclusions.

1. Two divergence-free carrier modes at \(k\) and \(l\), with
   \(q=k+l\) low, have a projected high--high-to--low symbol of rank two
   on \(q^\perp\) if and only if
   \[
   |k|\ne |l|.
   \]
   The equal-radius doublet, although geometrically natural, has rank only
   one.  This is an exact algebraic obstruction, not a numerical failure.
2. There is an explicit integer unequal-radius family for which the two
   child singular values stay bounded away from zero as the carrier
   frequency tends to infinity.  In one child direction the leading
   \(k-l\) wake derivative vanishes exactly; in the other, its size relative
   to the child is exactly
   \[
   \frac{|q|}{|k-l|}=O(|k|^{-1}).
   \]
3. The full generated Fourier lattice has an exact charge decomposition
   into one-dimensional sideband chains.  All nonzero charges may be
   declared outgoing wake **at the endpoint**, but they are not dynamically
   decoupled from the child: opposite charges feed charge zero.
4. One doublet can create only a one-phase shear child
   \(Z(q\cdot x)\).  Its velocity gradient has rank at most one, so it
   cannot be the nonsingular affine Kelvin pump
   \(\operatorname{diag}(-\alpha,-\beta,\alpha+\beta)\).
   At least three independent beat directions are necessary.  They are
   also sufficient at the principal algebraic level: an explicit rational
   three-beat family gives both the exact
   \(\gamma=5/4\) affine strain and a rank-five chart of
   \(\mathrm{Sym}_0^3\).

Thus the smallest doublet is a positive, uniformly ranked **daughter
charge**, but not a return cell.  The next prize-relevant target is a
three-beat, forward material-phase construction in which the nonzero
charge chains are retained through the active interval and routed to the
annular wake only afterward.

The identities are regression-checked in
[`checks/two_colour_endpoint_rank.py`](../checks/two_colour_endpoint_rank.py).

---

## 1. The exact two-wave Euler symbol

Use the velocity Fourier convention

\[
 u(x)=\sum_{n\in\mathbb Z^3}\widehat u_n e^{in\cdot x},
 \qquad
 n\cdot\widehat u_n=0,
 \qquad
 \widehat u_{-n}=\overline{\widehat u_n}.
\]

The Navier--Stokes equation (with Euler obtained by setting \(\nu=0\)) is

\[
 \partial_t\widehat u_n
 + i\mathbb P_n
 \sum_{p+r=n}
 (\widehat u_p\cdot r)\widehat u_r
 + \nu |n|^2\widehat u_n
=0.
\tag{1.1}
\]

Take two positive-frequency coefficients

\[
 \widehat u_k=a,\qquad
 \widehat u_l=b,\qquad
 a\cdot k=b\cdot l=0.
\]

At the sum \(q=k+l\), the symmetrized quadratic coefficient, with the
common factor \(-i\) suppressed, is

\[
 {\cal B}_q(a,b)
 :=
 \mathbb P_q\left((a\cdot l)b+(b\cdot k)a\right).
\tag{1.2}
\]

Since \(a\cdot k=b\cdot l=0\),

\[
 a\cdot l=a\cdot q,\qquad b\cdot k=b\cdot q.
\tag{1.3}
\]

The same real waves also generate the difference \(d=k-l\).  Its
coefficient is

\[
 {\cal B}_d(a,\overline b)
 :=
 \mathbb P_d\left(
 (a\cdot(-l))\overline b
+(\overline b\cdot k)a
\right).
\tag{1.4}
\]

Self-interactions at \(2k\) and \(2l\) vanish because a single transverse
plane wave is nonlinearly dark.  Thus \(q\) and \(d\) are the complete
second-order output of the real two-wave doublet.

---

## 2. Rank theorem and the equal-radius obstruction

Let

\[
 Q=|q|,\qquad e_3=\frac qQ.
\]

Assume \(k\) and \(l\) are not collinear.  After rotating about \(q\), write

\[
\begin{aligned}
 k&=K e_1+\kappa e_3,\\
 l&=-K e_1+(Q-\kappa)e_3,
\end{aligned}
\qquad K>0,
\tag{2.1}
\]

and set

\[
 e_2=e_3\times e_1,\qquad
 \delta=2\kappa-Q.
\tag{2.2}
\]

The difference vector is

\[
 d=2K e_1+\delta e_3.
\tag{2.3}
\]

Use the normalized in-plane transverse polarizations

\[
\begin{aligned}
 p&=\frac{\kappa e_1-Ke_3}{n_k},
 &n_k&=\sqrt{K^2+\kappa^2}=|k|,\\
 r&=\frac{(Q-\kappa)e_1+Ke_3}{n_l},
 &n_l&=\sqrt{K^2+(Q-\kappa)^2}=|l|.
\end{aligned}
\tag{2.4}
\]

Every complex transverse coefficient has the unique form

\[
 a=x p+y e_2,\qquad b=z r+w e_2.
\tag{2.5}
\]

Direct substitution into (1.2) gives

\[
 {\cal B}_q(a,b)
 =
 \frac{KQ\delta}{n_kn_l}xz\,e_1
+KQ\left(\frac{zy}{n_l}-\frac{xw}{n_k}\right)e_2.
\tag{2.6}
\]

The missing \(e_3\) component is removed by \(\mathbb P_q\).

Now

\[
 |k|^2-|l|^2
 =
Q(2\kappa-Q)
=Q\delta.
\tag{2.7}
\]

Equation (2.6) proves the exact dichotomy

\[
 \operatorname{span}_{a\perp k,\ b\perp l}
 {\cal B}_q(a,b)
 =
 \begin{cases}
 q^\perp,& |k|\ne |l|,\\[2mm]
 \operatorname{span}\{k\times l\},& |k|=|l|.
 \end{cases}
\tag{2.8}
\]

### Consequence

An equal-radius near-opposite doublet can create only the polarization
normal to the \(k,l\) plane.  Varying phases, amplitudes, or the two
transverse polarizations does not recover the missing in-plane child
direction.  The loss follows from \(|k|=|l|\) itself.

This is relevant to symmetric two-chain searches: pairing two carriers on
the same Fourier shell builds the rank defect into the ansatz.

---

## 3. A uniformly ranked integer family

Take

\[
 q=(0,0,1),\qquad
 k_N=(N,0,N),\qquad
 l_N=(-N,0,1-N).
\tag{3.1}
\]

Then \(q=k_N+l_N\), while

\[
 d_N=k_N-l_N=(2N,0,2N-1).
\tag{3.2}
\]

The two carriers lie on asymptotically equal shells but are not exactly
equal in length.  For the polarizations (2.4),

\[
 {\cal B}_q(p,r)
 =
 c_Ne_1,
\qquad
 c_N
 =
\frac{N(2N-1)}
{\sqrt{2}\,N\sqrt{2N^2-2N+1}}.
\tag{3.3}
\]

In particular,

\[
 c_N\longrightarrow1.
\tag{3.4}
\]

The in-plane daughter coefficient therefore has a uniform inverse.  The
relative shell mismatch tends to zero, so rank two does not require
placing the carriers at parametrically different frequencies.

### 3.1 Exact leading difference leakage

For the base pair \(p,r\),

\[
 {\cal B}_{d_N}(p,r)
 =
\frac{KQ^2}{n_kn_l}\mathbb P_{d_N}e_1.
\tag{3.5}
\]

Since

\[
\left|\mathbb P_de_1\right|
=\frac{|\delta|}{|d|},
\tag{3.6}
\]

comparison with (2.6) gives the exact identity

\[
 \frac{|{\cal B}_d(p,r)|}{|{\cal B}_q(p,r)|}
 =
\frac{Q}{|d|}.
\tag{3.7}
\]

For (3.1), this is \(O(N^{-1})\).  Hence the amplitude-product control
changes the in-plane child by order one while changing the leading
difference wake by only order \(N^{-1}\).

### 3.2 A child control with zero difference derivative

At the base point \(a=p,b=r\), choose

\[
 \delta a=\frac{n_l}{2KQ}e_2,\qquad
 \delta b=-\frac{n_k}{2KQ}e_2.
\tag{3.8}
\]

Differentiating (1.2) and (1.4) gives

\[
\begin{aligned}
 D{\cal B}_q[\delta a,\delta b]&=e_2,\\
 D{\cal B}_d[\delta a,\delta b]&=0.
\end{aligned}
\tag{3.9}
\]

Together with the amplitude-product variation, the projected two-by-two
child Jacobian is

\[
 \begin{pmatrix}
 c_N&0\\
 0&1
 \end{pmatrix}.
\tag{3.10}
\]

Its determinant tends to one.  This is the positive endpoint-rank result
which a same-shell doublet misses.

For complex coefficients, (2.6) has complex rank two, hence real rank four
onto the complete complex \(q\)-coefficient.  On the odd, sine-phase
subspace relevant to a child strain, (3.10) gives the two real transverse
strain coordinates.

---

## 4. What this says about the endpoint map

Start the \(q\)-mode at zero and let the two carriers be nonzero.  For each
fixed carrier \(K\), local smooth Euler evolution gives

\[
 \widehat u_q(T)
 =
-iT\,{\cal B}_q(a,b)+O_K(T^2).
\tag{4.1}
\]

Therefore the projected endpoint derivative at the base pair has the form

\[
 -iT
 \begin{pmatrix}
 c_N&0\\
 0&1
 \end{pmatrix}
 +O_N(T^2).
\tag{4.2}
\]

For every fixed \(N\), it is onto for all sufficiently small \(T>0\).
With small carrier heat, the Navier--Stokes Duhamel weights perturb these
nonzero entries and do not change their rank.

This statement must not be overstated.  The error constant in (4.1) is
not uniform in \(N\) in fixed Eulerian Fourier coordinates.  Once the
low child appears, its transport of a carrier differentiates the carrier
phase and produces factors of \(N\).  Over an order-one interval this
creates the full sideband chain described below.

Thus (3.10) is a scale-uniform **principal material-phase symbol**, not yet
a scale-uniform Eulerian endpoint theorem.  A Kelvin/WKB coordinate must
resum the transport chain instead of estimating it term by term.

---

## 5. Exact charge decomposition of every generated chain

Because \(k\) and \(l\) are linearly independent, every generated Fourier
mode has a unique label

\[
 \lambda_{r,s}=rk+sl,\qquad (r,s)\in\mathbb Z^2.
\tag{5.1}
\]

Define its charge

\[
 h=r-s.
\tag{5.2}
\]

Since \(l=q-k\),

\[
 \lambda_{r,s}=h k+s q.
\tag{5.3}
\]

For each fixed \(h\), the modes form a one-dimensional \(q\)-sideband
chain

\[
 \Lambda_h
 =
\{hk+nq:n\in\mathbb Z\}.
\tag{5.4}
\]

The quadratic Euler convolution respects the exact grading

\[
 \Lambda_h+\Lambda_g\subseteq\Lambda_{h+g}.
\tag{5.5}
\]

The principal pieces are:

\[
\begin{array}{c|c|c}
\text{object}&(r,s)&h\\ \hline
k&(1,0)&1\\
l&(0,1)&-1\\
q=k+l&(1,1)&0\\
d=k-l&(1,-1)&2.
\end{array}
\tag{5.6}
\]

Reality supplies the opposite charges.  At cubic order the new labels
include

\[
(\pm2,\pm1),\qquad(\pm1,\pm2),
\tag{5.7}
\]

as well as feedback into the parent modes.  Iteration generates all
allowed labels as an algebraic closure; special polarization cancellations
can leave individual coefficients zero, while a generic doublet populates
the chains.

### Endpoint wake split

The charge-zero chain is precisely

\[
 \Lambda_0=\{nq:n\in\mathbb Z\},
\tag{5.8}
\]

the one-phase child/centre profile.  It is therefore natural to define

\[
\begin{aligned}
 X_{\rm child}&=\bigoplus_{|n|\le M}\widehat u_{nq},\\
 X_{\rm wake}&=\bigoplus_{h\ne0}\Lambda_h
 \ \oplus\
 \bigoplus_{|n|>M}\widehat u_{nq}.
\end{aligned}
\tag{5.9}
\]

All nonzero charges, including the leading \(h=\pm2\) difference chain,
can be assigned to the outgoing wake at the endpoint.

They cannot be deleted during the active interval.  Equation (5.5) also
gives

\[
 \Lambda_h+\Lambda_{-h}\subseteq\Lambda_0,
\tag{5.10}
\]

so opposite wake charges feed the child.  The correct transition problem
is a forward evolution of the complete truncated charged lattice followed
by the projection (5.9).  It is not a finite triad and not a block-triangular
system.

Spatial export is a further issue.  A Gavrilov wake is disjoint in
physical space, whereas (5.9) is only a spectral endpoint split.  Moving
the nonzero charges into disjoint bubbles without changing the prescribed
charge-zero jet requires a localized router and pressure-tail estimate.
The algebra here identifies exactly what the router must accept; it does
not construct it.

---

## 6. Decisive limitation of one doublet

Even retaining every charge-zero harmonic does not make a three-dimensional
affine child.  A field supported on \(\Lambda_0\) has the form

\[
 U_0(x)=Z(q\cdot x),\qquad q\cdot Z=0.
\tag{6.1}
\]

Its gradient is

\[
 \nabla U_0(x)=Z'(q\cdot x)\otimes q,
\tag{6.2}
\]

so

\[
 \operatorname{rank}\nabla U_0(x)\le1,\qquad
 \det\nabla U_0(x)=0.
\tag{6.3}
\]

The Kelvin pump required by the cascade is

\[
 S=\operatorname{diag}(-\alpha,-\beta,\alpha+\beta),
\qquad \alpha,\beta>0,
\tag{6.4}
\]

which is nonsingular.  No one-phase profile (6.1) can agree with \(Sx\)
to first order in a core.

More generally, the sum of \(m\) distinct one-phase low outputs has
gradient rank at most \(m\) at a point.  Hence

\[
 \boxed{\text{at least three independent low beat directions are
 necessary for a full-rank affine child.}}
\tag{6.5}
\]

The rank-two result (3.10) is therefore a daughter-charge submersion, not
the complete affine-return submersion.

This lower bound applies to the bare beat output.  Arbitrary
space-dependent envelopes can introduce additional low directions through
their gradients, but then the endpoint variables are the envelope and
localization wake, not the two finite Fourier coefficients tested here.

---

## 7. Three beats are algebraically sufficient

The lower bound (6.5) is sharp at the level of a prescribed affine jet.
For the admissible cascade exponent

\[
 \gamma=\frac54,
\tag{7.1}
\]

set

\[
 S_*=\operatorname{diag}\left(-1,-\frac54,\frac94\right).
\tag{7.2}
\]

Take the three integer low directions

\[
\begin{aligned}
q_1&=(-45,-36,20),\\
q_2&=(-4,-5,9),\\
q_3&=(1,1,1),
\end{aligned}
\tag{7.3}
\]

and the rational transverse coefficients

\[
\begin{aligned}
c_1&=\left(\frac1{13},-\frac5{56},\frac9{728}\right),\\
c_2&=\left(-\frac4{13},\frac{25}{56},\frac{81}{728}\right),\\
c_3&=\left(\frac{16}{13},-\frac{125}{56},\frac{729}{728}\right).
\end{aligned}
\tag{7.4}
\]

Exact rational arithmetic gives

\[
 c_j\cdot q_j=0
\quad(j=1,2,3),
\tag{7.5}
\]

and

\[
 \sum_{j=1}^3c_j\otimes q_j
 =
S_*.
\tag{7.6}
\]

Consequently, the divergence-free low field

\[
 U_{\rm low}(x)
 =
\sum_{j=1}^3c_j\sin(q_j\cdot x)
\tag{7.7}
\]

has

\[
 U_{\rm low}(0)=0,\qquad
 \nabla U_{\rm low}(0)=S_*.
\tag{7.8}
\]

This is an exact local affine-jet witness, not an approximation.

There is also room for a robust strain submersion.  For each \(q_j\), let
\(c\) range over the two-dimensional plane \(q_j^\perp\) and consider

\[
 \operatorname{sym}(c\otimes q_j).
\tag{7.9}
\]

The six resulting columns span the five-dimensional space
\(\mathrm{Sym}_0^3\).  With the transverse bases

\[
(-q_{j,2},q_{j,1},0),\qquad
(-q_{j,3},0,q_{j,1}),
\tag{7.10}
\]

the first five columns have determinant

\[
-1\,214\,003\,700
\tag{7.11}
\]

in the coordinates
\((M_{11},M_{22},M_{12},M_{13},M_{23})\).
Thus the rank is exactly five over \(\mathbb Q\).

### 7.1 How to realize the three low symbols

For each \(q_j\), choose a high integer vector \(k_j\) with

\[
 k_j\not\parallel q_j,\qquad
2k_j\cdot q_j\ne |q_j|^2,
\tag{7.12}
\]

and put

\[
 l_j=q_j-k_j.
\tag{7.13}
\]

The second condition is exactly \(|k_j|\ne|l_j|\).  By (2.8), the matched
doublet \((k_j,l_j)\) controls both transverse coefficients at \(q_j\).
Taking

\[
 k_j=Nw_j
\tag{7.14}
\]

with \(w_j\cdot q_j\ne0\) and \(w_j\not\parallel q_j\) keeps the matched
symbol nondegenerate as \(N\to\infty\).

The six high carriers can be organized into two spectral colour classes,
\(\{k_j\}\) and \(\{l_j\}\).  If the \(w_j\) are distinct, every unmatched
sum

\[
 k_i+l_j=N(w_i-w_j)+q_j,\qquad i\ne j,
\tag{7.15}
\]

stays high.  Therefore the **low projected** quadratic derivative is the
direct sum of the three matched rank-two symbols and is onto
\(\mathrm{Sym}_0^3\).

This is a formal principal-level positive submersion.  The unmatched high
products in (7.15) are new wake charges and need not be small.  They are
the next analytic gate.  Three sequential doublet pulses avoid simultaneous
unmatched products but then the earlier low modes deform the later carrier
phases.  Either implementation requires forward Kelvin/material
propagation and a wake router.

---

## 8. Why no exact finite Fourier closure follows

The support generated by two non-collinear real carriers is the lattice
(5.1), not the four modes \(\{\pm k,\pm l\}\).  Except for special
stationary two-dimensional-like or Beltrami configurations, a
three-dimensional finite Fourier set cannot close under Euler evolution.
The rank-producing unequal-radius doublet is not one of those stationary
exceptions.

The multi-frequency Craik--Criminale literature does not supply a global
shortcut.  It does provide the correct material-phase lesson:
incommensurate Kelvin waves can be added sequentially along selected
Lagrangian trajectories.  But that construction is explicitly evaluated
along a single flowline for a general base, and the authors emphasize that
simultaneously adding the waves is not in general an exact global
Navier--Stokes solution.  It therefore supports a local convective WKB
description, not the compact global endpoint cell required here.

The finite Fourier calculation should be read in exactly that way:

* (2.6) is the principal child symbol;
* (5.4) is the material sideband state that must be resummed;
* (5.9) is the outgoing child/wake projection;
* no finite-mode truncation is an exact stage.

---

## 9. Revised breakthrough theorem

The algebra now specifies a narrower and testable theorem.

> **Three-beat forward Kelvin--Reynolds submersion.**  Construct a compact
> Gevrey Euler--Reynolds transition, using three unequal-radius matched
> doublets and a finite-band charged material-phase state, such that:
>
> 1. the matched low symbols remain uniformly close to the rank-five chart
>    in Section 7;
> 2. the three low outputs create a localized full-rank strain core with
>    \(\gamma\in(1,3/2)\);
> 3. every charge chain through the stage-dependent truncation is solved
>    forward, not discarded or inverted by a stationary operator;
> 4. unmatched high products and nonzero endpoint charges are transported
>    into the disjoint annular wake;
> 5. opposite-charge feedback into the child is included in the endpoint
>    derivative;
> 6. the projected child derivative has a scale-uniform tame right inverse;
>    and
> 7. the discarded charge/harmonic tail and the physical seams are
>    \(e^{-cj^2}\), hence terminally \(C^\infty\)-flat after physical
>    rescaling.

Sections 2, 3, and 7 prove the principal rank and give an explicit target
point.  Sections 5 and 8 identify why that is not yet the theorem:
charge feedback and physical wake routing remain.

### Computational decision

A generic GPU optimization is not the next decisive step.  The equal-shell
ansatz is already ruled out exactly, and optimizing its gain cannot restore
the missing rank.  The useful next computation, after a concrete localized
three-beat parametrix is written, would be a **charged endpoint Jacobian**
test:

* retain the full finite charge band \(|h|,|n|\le M\);
* work in material/Kelvin coordinates;
* measure the five singular values of the child strain projection;
* measure unmatched-charge growth in a weighted wake norm; and
* verify convergence as \(M\), carrier \(N\), and spatial localization are
  increased.

Until that parametrix exists, the remaining gate is analytic construction,
not raw compute.

---

## Primary references

* R. T. Craik and W. O. Criminale,
  [Evolution of wavelike disturbances in shear flows: a class of exact
  solutions of the Navier--Stokes
  equations](https://doi.org/10.1098/rspa.1986.0061).
* B. R. Fabijonas and D. D. Holm,
  [Multi-frequency Craik--Criminale solutions of the Navier--Stokes
  equations](https://arxiv.org/abs/nlin/0304049).
* H. R. Dullin and J. Worthington,
  [Stability Theory of the 3-Dimensional Euler
  Equations](https://arxiv.org/abs/1903.09970).
* N. Kishimoto and T. Yoneda,
  [Characterization of three-dimensional Euler flows supported on finitely
  many Fourier modes](https://arxiv.org/abs/2110.08039).
