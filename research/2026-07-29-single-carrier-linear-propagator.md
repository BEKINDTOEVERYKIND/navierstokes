# The single-carrier linear propagator: one derivative, one phase gauge, and no \(K\)-chain

Date: 2026-07-29

## Status and verdict

This note isolates the **linearized finite-charge question** left between
the exact single-carrier algebra in
`2026-07-29-single-carrier-shear-sideband-audit.md` and the nonlinear
Cheverry gate in `2026-07-29-cheverry-material-phase-gate.md`.

It does not construct a Navier--Stokes singularity.

The conclusion is positive, but deliberately narrower than a nonlinear
one-phase theorem.

* For
  \[
  k=Ke_3,\qquad a=e_1,\qquad r_j=q_j-Ke_3,
  \]
  the three published low directions retain an exact rank-five strain
  chart.  No helicity or second carrier direction is needed.
* On the full one-direction charge lattice, transversality removes \(K\)
  from every nonzero-charge/nonzero-charge Euler interaction.
* The unrestricted **quadratic** map on a band of radius \(M\) really can
  have size \(M^2\).  In contrast, its **linearization about any fixed
  finite-label star** has norm \(O(M)\) in weighted Wiener and
  coefficient \(\ell^2\) norms.  The distinction is essential.
* The only elementary \(O(K)\) block is the exact length-two shear chain
  \[
  e_3\longmapsto \partial_3W_K\longmapsto0.
  \]
  It is a phase-translation mode.  A material eikonal moves it into the
  phase variable; even without doing so its cost is only polynomial in
  \(K\), not \(K^M\).
* In affine material coordinates the inviscid propagator therefore costs
  at most
  \[
  K^D M^D\exp(CMT)
  \tag{0.1}
  \]
  for a fixed \(D\).  A terminal heat inverse adds
  \[
  \exp(C\theta M^2),\qquad \theta=\varepsilon K^2.
  \tag{0.2}
  \]
  Consequently, if
  \[
  K_j=j^A,\qquad
  M_j\asymp\frac{j^2}{\log j},\qquad
  \theta_jM_j^2\longrightarrow0,
  \tag{0.3}
  \]
  then the logarithm of the linear propagator cost is \(o(j^2)\).

This passes the requested **linear propagator** ledger.  It does not prove
that the order-one endpoint map stays onto, nor that the nonlinear
material-profile trajectory has a uniform analytic radius.  Cheverry's
hidden instability is precisely a one-slow-derivative effect; on the
retained band it costs \(e^{CM}\), which is acceptable here, while in an
untruncated Gevrey-two Cauchy problem it remains a genuine obstruction.

The finite-label and nilpotent algebra below is checked independently in
[`checks/single_carrier_linear_propagator.py`](../checks/single_carrier_linear_propagator.py).

---

## 1. Exact child chart, in the form needed here

Let

\[
q=(\alpha,\beta,\delta),\qquad
k=Ke_3,\qquad r=q-k,\qquad a=e_1.
\tag{1.1}
\]

For \(b\cdot r=0\), the symmetrized Euler coefficient at \(q=k+r\) is

\[
\begin{aligned}
{\cal C}_{q,K}b
&=
P_q\big[(a\cdot r)b+(b\cdot k)a\big]\\
&=
P_q\big[\alpha b+(b\cdot q)e_1\big].
\end{aligned}
\tag{1.2}
\]

Both equalities are exact.  When \(K\ne\delta\), parametrize
\(r^\perp\) by

\[
\begin{aligned}
b_1&=\left(1,0,\frac{\alpha}{K-\delta}\right),\\
b_2&=\left(0,1,\frac{\beta}{K-\delta}\right).
\end{aligned}
\tag{1.3}
\]

Let \(V_s\) denote the unprojected vector in (1.2) associated with
\(b_s\).  Since adding a multiple of \(q\) to either column does not
change its triple product with \(q\),

\[
\det[P_qV_1,P_qV_2,q]
=\det[V_1,V_2,q].
\tag{1.4}
\]

A direct calculation gives

\[
\boxed{
\det[V_1,V_2,q]
=
\frac{\alpha^2(2K\delta-|q|^2)}{K-\delta}.}
\tag{1.5}
\]

The pole at \(K=\delta\) belongs only to this normalized basis.  A basis
valid at every \(K\) is

\[
\begin{aligned}
d_1&=(-\beta,\alpha,0),&
d_2&=(K-\delta,0,\alpha)
\end{aligned}
\tag{1.6}
\]

for \(r^\perp\), together with

\[
\begin{aligned}
c_1&=(-\beta,\alpha,0),&
c_2&=(-\delta,0,\alpha)
\end{aligned}
\tag{1.7}
\]

for \(q^\perp\).  In these bases the exact matrix is

\[
[{\cal C}_{q,K}]
=
\begin{pmatrix}
\alpha&-2\alpha\beta K/|q|^2\\
0&\alpha(|q|^2-2K\delta)/|q|^2
\end{pmatrix}.
\tag{1.8}
\]

Thus the exact rank-two condition is

\[
\boxed{\alpha\ne0,\qquad 2K\delta\ne|q|^2.}
\tag{1.9}
\]

For fixed \(q\), the limiting oriented area is \(2\alpha^2\delta\).
Hence

\[
\alpha\delta\ne0
\tag{1.10}
\]

is the uniform large-\(K\) condition.

The three directions

\[
\begin{aligned}
q_1&=(-45,-36,20),\\
q_2&=(-4,-5,9),\\
q_3&=(1,1,1)
\end{aligned}
\tag{1.11}
\]

satisfy (1.7).  Since each child map is onto \(q_j^\perp\), the six
strain columns are the same six transverse columns whose rank is five in
the earlier exact chart.  Therefore the single-carrier construction keeps
the rank-five principal endpoint symbol.

This is only an instantaneous principal chart.  It is not yet an
order-one endpoint submersion.

---

## 2. Exact \(K\)-null identity on the full lifted lattice

Let \(Q\) have columns \(q_1,q_2,q_3\), and set

\[
\xi_n=Qn,\qquad n\in\mathbb Z^3.
\tag{2.1}
\]

The one-direction lifted modes are

\[
p_\alpha=Kh\,w+\xi_n,\qquad
\alpha=(h,n),\qquad w=e_3.
\tag{2.2}
\]

Take

\[
\begin{aligned}
p&=Kh\,w+\xi,\\
r&=Kg\,w+\eta,
\end{aligned}
\qquad h,g\ne0,
\tag{2.3}
\]

and amplitudes

\[
A\cdot p=0,\qquad B\cdot r=0.
\tag{2.4}
\]

Transversality gives

\[
Kh\,A\cdot w=-A\cdot\xi,\qquad
Kg\,B\cdot w=-B\cdot\eta.
\tag{2.5}
\]

Consequently,

\[
\boxed{
\begin{aligned}
A\cdot r
&=A\cdot\eta-\frac ghA\cdot\xi,\\
B\cdot p
&=B\cdot\xi-\frac hgB\cdot\eta.
\end{aligned}}
\tag{2.6}
\]

The projected bilinear coefficient is therefore

\[
\boxed{
\begin{aligned}
{\cal B}_{p+r}(A,B)
=P_{p+r}\bigg[
&\left(A\cdot\eta-\frac ghA\cdot\xi\right)B\\
&+
\left(B\cdot\xi-\frac hgB\cdot\eta\right)A
\bigg].}
\tag{2.7}
\]

There is no \(K\) in (2.7).  No helicity assumption was used.

The same identity survives an affine material deformation.  If

\[
\kappa_\alpha(t)
=F(t)^{-T}(Kh\,w+\xi_n),
\tag{2.8}
\]

write

\[
\omega(t)=F^{-T}w,\qquad
\zeta_n(t)=F^{-T}\xi_n.
\tag{2.9}
\]

Then (2.6) holds verbatim with
\((w,\xi,\eta)\) replaced by
\((\omega,\zeta_n,\zeta_m)\).  It is an algebraic consequence of
\(A\cdot\kappa_\alpha=0\), not a Euclidean-shell coincidence.

---

## 3. The quadratic map is \(M^2\), but the fixed-star derivative is \(M\)

This is the main distinction of the note.

### 3.1 The unrestricted bilinear map

On

\[
|h|+|g|+|n|_1+|m|_1\le M,
\tag{3.1}
\]

one has

\[
|\xi_n|+|\xi_m|\le C_QM.
\tag{3.2}
\]

Formula (2.7) only gives

\[
|{\cal B}_{p+r}(A,B)|
\le CM^2|A||B|.
\tag{3.3}
\]

This is sharp.  For example, taking

\[
h=1,\quad \xi=Mq_1,\qquad
g=M,\quad\eta=0
\tag{3.4}
\]

produces a projected coefficient asymptotic to \(45M^2\) for suitable
unit transverse amplitudes.  Thus an unrestricted quadratic
one-direction profile is not a one-derivative analytic vector field in
the naive isotropic charge norm.

### 3.2 Linearization about a fixed star

Let \({\cal S}_*\) be a fixed finite set of background labels

\[
\beta=(g,m),\qquad g\ne0,
\tag{3.5}
\]

independent of \(K\) and \(M\).  Assume

\[
\max_{\beta\in{\cal S}_*}
\left(
|g|+|g|^{-1}+|m|_1+\|B_\beta\|_{L^\infty_t}
\right)
\le C_*.
\tag{3.6}
\]

This includes the parent and any fixed finite collection of initial
sidebands.

Let the perturbation label be

\[
\alpha=(h,n),\qquad
h\ne0,\qquad |h|+|n|_1\le M.
\tag{3.7}
\]

In (2.7), make one input \(\beta\in{\cal S}_*\) and the other input
\(\alpha\).  The slow frequency and charge of the background input are
bounded.  Therefore

\[
\begin{aligned}
\left|
A_\alpha\cdot\eta_\beta
-\frac{g}{h}A_\alpha\cdot\xi_n
\right|
&\le C_*M|A_\alpha|,\\
\left|
B_\beta\cdot\xi_n
-\frac{h}{g}B_\beta\cdot\eta_\beta
\right|
&\le C_*M|B_\beta|.
\end{aligned}
\tag{3.8}
\]

Since \(\|P_\kappa\|_{\mathbb C^3\to\mathbb C^3}=1\), every nonzero
matrix edge of the linearized charged operator is \(O(M)\), uniformly in
\(K\).

Each label is connected to only \(2|{\cal S}_*|\) background shifts,
including reality conjugates.  Schur's test consequently gives:

> **Theorem 3.1 (fixed-star charge bound).**  
> Let \({\cal L}_{K,M,*}(t)\) be the Galerkin linearized Euler operator
> about a fixed finite-label one-direction star satisfying (3.6), acting
> on the nonzero-charge band (3.7).  Suppose the material deformation
> \(F,F^{-1}\) is uniformly bounded.  Then
> \[
> \boxed{
> \|{\cal L}_{K,M,*}(t)\|_{\ell^2\to\ell^2}
> \le C(1+M),}
> \tag{3.9}
> \]
> where \(C\) is independent of \(K\) and \(M\).

The same estimate holds in a weighted Wiener norm

\[
\|Z\|_{\rho,1}
=\sum_\alpha e^{\rho(|h|+|n|_1)}|Z_\alpha|.
\tag{3.10}
\]

Indeed, every background shift has fixed degree, so the ratio of the
weight at the target to the weight at the source is bounded by a constant
depending only on \(\rho\) and \({\cal S}_*\).

This is the precise sense in which the off-shell **linearized** operator
costs one charge derivative, even though the unrestricted quadratic map
costs two.

### 3.3 Affine strain and moving Leray projection

For an affine low field \(V=S(t)x\), the Kelvin amplitude term is

\[
-S A_\alpha
+2\kappa_\alpha
\frac{\kappa_\alpha\cdot SA_\alpha}{|\kappa_\alpha|^2}.
\tag{3.11}
\]

Its norm is at most \(3\|S\|\), independently of
\(|\kappa_\alpha|\).  The material wave vector solves

\[
\kappa_\alpha'=-S^T\kappa_\alpha.
\tag{3.12}
\]

Since the Leray matrix

\[
P_\kappa=I-\frac{\kappa\otimes\kappa}{|\kappa|^2}
\tag{3.13}
\]

is homogeneous of degree zero, its time derivative under (3.12) is
\(O(\|S\|)\), again independently of \(K\).  Thus affine transport and
pressure add no carrier derivative to (3.9).

---

## 4. Charge zero and the exact phase-translation Jordan block

The hypothesis \(h\ne0\) in Theorem 3.1 is necessary.  A charge-zero
velocity can have an order-one component in the fast normal direction.
Its advection of the carrier is then \(O(K)\).

This \(O(K)\) term has an exact geometric meaning.  Let

\[
W_K(x)=U(Kx_3)e_1
\tag{4.1}
\]

be any smooth periodic shear, and define

\[
L_{W_K}z
=P\big[(W_K\cdot\nabla)z+(z\cdot\nabla)W_K\big].
\tag{4.2}
\]

Then

\[
\boxed{
L_{W_K}e_3=\partial_3W_K,\qquad
L_{W_K}(\partial_3W_K)=0.}
\tag{4.3}
\]

The first identity follows because \(e_3\cdot\nabla W_K=\partial_3W_K\),
which is already divergence-free.  The second follows because both
\(W_K\) and \(\partial_3W_K\) point in the \(e_1\) direction and depend
only on \(x_3\).

Thus the unmodulated linearized Euler evolution contains the exact
length-two chain

\[
e_3\longmapsto\partial_3W_K\longmapsto0,
\qquad
\|\partial_3W_K\|\asymp K.
\tag{4.4}
\]

It gives a polynomial \(K\) loss, not an exponential \(e^{cK}\) and not
a nilpotent chain of length \(M\).

It is exactly a phase translation.  For every constant \(c\),

\[
u_c(t,x)
=ce_3+U(K(x_3-ct))e_1
\tag{4.5}
\]

is an exact Euler solution.  Differentiating (4.5) at \(c=0\) gives

\[
e_3-t\partial_3W_K,
\tag{4.6}
\]

the solution generated by (4.3).

For a general charge-zero velocity \(V\), introduce a phase

\[
\partial_t\phi+V\cdot\nabla\phi=0.
\tag{4.7}
\]

Writing the oscillatory field as \(W(t,x,K\phi(t,x))\) cancels

\[
K(\partial_t\phi+V\cdot\nabla\phi)\partial_\vartheta W
\tag{4.8}
\]

exactly.  Physical incompressibility gives, on the mean-zero profile,

\[
K\nabla\phi\cdot W
=-\partial_\vartheta^{-1}\operatorname{div}_xW.
\tag{4.9}
\]

Therefore the high--high fast term becomes

\[
K(W\cdot\nabla\phi)\partial_\vartheta W
=-
\left(
\partial_\vartheta^{-1}\operatorname{div}_xW
\right)\partial_\vartheta W.
\tag{4.10}
\]

Equations (4.7)--(4.10) are the material-phase normal form.  They remove
the only naked \(K\) from the inviscid profile equation.

The change of variables is not uniformly conditioned in the ordinary
physical velocity norm: differentiating a phase shift produces
\(K\partial_\vartheta W\).  On a finite band its reconstruction cost is
at worst a fixed polynomial \(K^DM^D\).  That polynomial is harmless for
the logarithmic ledger below.

---

## 5. Propagator and terminal inverse ledger

Consider the material-phase finite-band system

\[
\partial_tZ={\cal L}_{K,M,*}(t)Z+G.
\tag{5.1}
\]

Assume:

1. the background star has fixed finite lifted support as in (3.6);
2. the charge-zero velocity is included in the material phase (4.7);
3. the material deformation and the phase remain nondegenerate with
   uniformly bounded normalized derivatives;
4. pressure is either the exact affine Fourier projection or an
   order-zero material-phase projector with uniform finite-band bounds.

Then Theorem 3.1, the affine estimate (3.11), and the polynomial phase
coordinate change imply

\[
\boxed{
\|{\cal U}_{K,M}(t,s)\|
+\|{\cal U}_{K,M}(s,t)\|
\le
C K^DM^D e^{CM|t-s|}.}
\tag{5.2}
\]

This is an ordinary finite-dimensional evolution estimate; it does not
use a sign-definite energy--helicity sector.  The one-carrier shear is not
a nonzero curl eigenfield, so the Beltrami second-variation argument is
not the relevant symmetrizer.

Add normalized viscosity \(\varepsilon\).  Forward heat is contractive.
If a complete retained terminal state is prescribed, the backward heat
factor obeys

\[
\exp\left(
\varepsilon\int_s^t
\max_{\alpha\in{\cal I}_M}|\kappa_\alpha(\tau)|^2\,d\tau
\right)
\le
\exp(C_T\theta M^2),
\qquad
\theta=\varepsilon K^2.
\tag{5.3}
\]

Combining (5.2) and (5.3) gives

\[
\log C_{\mathrm{prop}}
\le
CM+C\theta M^2+D\log K+D\log M+O(1).
\tag{5.4}
\]

For (0.3),

\[
CM_j
=O\left(\frac{j^2}{\log j}\right)=o(j^2),
\tag{5.5}
\]

and

\[
\log K_j+\log M_j=O(\log j)=o(j^2).
\tag{5.6}
\]

Hence

\[
\boxed{\log C_{\mathrm{prop},j}=o(j^2).}
\tag{5.7}
\]

This is the desired finite-band linear-propagator threshold.

### Endpoint qualification

The exact rank-five calculation proves that the instantaneous
sideband-to-child symbol is onto.  Estimate (5.2) shows that transporting
an already-surjective control map through the stage costs at most
\(\exp(o(j^2))\).

It does **not** prove that the projected order-one endpoint Jacobian
remains onto.  A singular value can cross zero while the propagator stays
bounded.  One still needs either:

* a direct endpoint determinant computation along the selected path;
* a short-pulse argument with accumulated order-one child output; or
* an invariant graph/shooting theorem which keeps the rank-five block
  transverse.

Thus (5.7) closes the stiffness estimate, not the endpoint rank theorem.

---

## 6. Relation to Cheverry's phase cascade

Cheverry's primary paper contains all three warnings relevant here.

1. Its large-amplitude discussion says that a fixed ordinary one-phase
   BKW ansatz on an order-one interval is generically inadequate and that
   an infinite hierarchy of phase corrections appears.
2. In the profile divergence constraint, the apparently singular normal
   component is rewritten with
   \(\partial_\vartheta^{-1}\operatorname{div}_x\), exactly the mechanism
   in (4.9)--(4.10).
3. After this rewrite, the inviscid operators are first order but not all
   skew-symmetric.  Cheverry identifies the resulting derivative loss as
   the source of hidden instabilities and proves stability only after
   adding a sufficiently strong compatible profile viscosity.

These facts neither invalidate nor prove (5.7).

On a retained slow band \(|\eta|\le M\), one first-order hidden
instability can grow like

\[
e^{c|\eta|t}\le e^{cMt}.
\tag{6.1}
\]

This can be sharp for unstable hydrostatic profiles, but

\[
\log e^{cM_j}=O(j^2/\log j)=o(j^2).
\tag{6.2}
\]

Therefore Cheverry's one-derivative loss is compatible with the present
finite-band propagator ledger.  It is incompatible with a generic
untruncated Gevrey-two Cauchy theorem, because a Gevrey-two weight cannot
absorb \(e^{c|\eta|t}\) on a fixed interval.

There is also no applicable parabolic shortcut.  Cheverry's stability
theorem uses a compatible profile viscosity bounded below by a sufficiently
large constant.  Here

\[
\theta_j=\varepsilon_jK_j^2\longrightarrow0
\tag{6.3}
\]

on the fast direction, while the slow damping is even smaller.  His
stability theorem cannot be imported into this regime.

The material phase transported by the complete charge-zero flow does
formally resum the geometrical phase hierarchy, and adjusting phases can
be represented as angular profile translations.  What remains unproved
is uniform control of that exact \(K\)-dependent phase/profile system.

---

## 7. Pressure and projection losses

There are three different pressure statements.

### 7.1 Affine material phases

For affine \(V=Sx\), every mode has the exact Kelvin wave vector (2.8).
The Leray projector is the block-diagonal matrix \(P_{\kappa_\alpha}\).
It has norm one in every coefficient norm, and its time variation costs
only \(O(\|S\|)\).  There is no \(K\)- or \(M\)-loss from pressure.

### 7.2 Variable nondegenerate phases

For a variable phase, the pressure operator is

\[
\Delta_K
=K^2|\nabla\phi|^2\partial_\vartheta^2
+K\left(
2\nabla\phi\cdot\nabla_x+\Delta\phi
\right)\partial_\vartheta
+\Delta_x.
\tag{7.1}
\]

On nonzero angle charges and under
\(|\nabla\phi|\ge c>0\), its semiclassical inverse has an expansion in
\(K^{-1}\).  The coefficient of \(K^{-n}\) costs at most \(n\) slow
derivatives: the first-order term advances the expansion by one, while
\(\Delta_x\) advances it by two.  Thus pressure does not double the
Gevrey order.

This is a slow-envelope parametrix, not an inverse on every
angle-mean-zero profile.  Exact gauge modes

\[
H(\vartheta-K\phi(x))
\tag{7.2}
\]

lie in its kernel.  A microlocal separation from that gauge set, followed
by the ordinary physical pressure solve for the exponentially small
residual, is required.

### 7.3 Charge-zero pressure

The angle mean is the ordinary low pressure and belongs in the material
flow.  On the torus, nonzero physical integer frequencies have no
small-divisor loss.  Exact lifted aliases occur only at charge order
\(O(K)\) for the published \(Q\); the retained regime \(M\ll K\) stays
separated from them.

These pressure facts support the \(O(M)\) fixed-star linear estimate.
They do not turn the \(O(M^2)\) unrestricted quadratic map into an
\(O(M)\) map.

---

## 8. Sharp conclusion

No sign-definite Beltrami sector is needed for the single-carrier repair,
and no fatal \(O(K)\) unstable or length-\(M\) nilpotent chain appears in
the material-phase fixed-star linearization.

What is proved is:

\[
\boxed{
\text{fixed-star linearized cost }
\le K^DM^D e^{CM+C\theta M^2}
=\exp(o(j^2)).}
\tag{8.1}
\]

What is not proved is:

1. an \(O(M)\) bound for the unrestricted quadratic charge map;
2. a uniform analytic lifespan for the nonlinear hydrostatic profile;
3. preservation of the rank-five endpoint determinant for order-one
   time;
4. convergence or Gevrey summability of the full phase/profile hierarchy;
5. localization, alias routing, or wake closure.

The next decisive theorem is therefore not another shell calculation.
It is a selected-path result: construct the actual three-sideband
material-profile trajectory and prove either

\[
\sup_{0\le t\le T}
\|D{\cal N}(W_*(t))\|_M\le CM
\tag{8.2}
\]

or an equivalent slow-analytic estimate with a radius that remains
positive, while tracking the rank-five endpoint block.  Failure of
(8.2) through a background slope of order \(M\) would restore the
\(M^2\) obstruction.

---

## Primary references

* Christophe Cheverry,
  [*Cascade of phases in turbulent flows*](https://arxiv.org/abs/math/0402408),
  especially the large-amplitude discussion and Section 6's first-order
  reduction and hidden-instability estimate.
* A. D. D. Craik and W. O. Criminale,
  [*Evolution of wavelike disturbances in shear flows: a class of exact
  solutions of the Navier--Stokes equations*](https://doi.org/10.1098/rspa.1986.0061),
  for exact material/Kelvin phases on affine backgrounds.
