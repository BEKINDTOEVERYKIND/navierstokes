# Exact audit of the common-carrier shear sideband repair

## Executive verdict

The proposed replacement of the common-helicity Beltrami star by one
common **shear** carrier is a real algebraic improvement.

Take

\[
 w=e_3,\qquad a=e_1,\qquad k=Kw,
\]

and, for each prescribed low direction \(q\), take

\[
 r=q-Kw,\qquad b\in r^\perp.
\]

Then the matched high--high interaction at \(q=k+r\) is exactly

\[
 L_{q,K}b
 =
 P_q\big[(a\cdot q)b+(b\cdot q)a\big].
\tag{0.1}
\]

This is not merely its \(K\to\infty\) limit.  For the three published
directions

\[
\begin{aligned}
q_1&=(-45,-36,20),\\
q_2&=(-4,-5,9),\\
q_3&=(1,1,1),
\end{aligned}
\tag{0.2}
\]

the map (0.1) has rank two for **every positive integer \(K\)**.  Hence
the six partner-polarization variables retain the previously verified
rank-five chart onto \(\mathrm{Sym}_0^3\) at every integer carrier
frequency.  The exact inherited five-by-five minor is

\[
 -\frac{644\,815\,080\,000}{3721}
 (3721-40K)(61-9K),
\tag{0.3}
\]

which never vanishes for integer \(K\).

The one common fast direction also removes the former multi-colour
off-shell obstruction much more strongly than an initial-shell
calculation suggests.  On the complete charged lattice, transversality
cancels the nominal factor \(K\) in **every charged--charged
interaction**, at every generation.  This remains exactly true after
pullback by the full low-flow material map.

There are, however, four qualifications that a proof cannot skip.

1. In an unrestricted charge ball of radius \(M\), the surviving
   multiplier is \(O(M^2)\), not \(O(M)\).  A direct family below has
   projected symbol asymptotic to \(45M^2\).  Ordinary charge-analytic
   norms therefore do not turn the system into a one-derivative
   nonlinearity.  The natural raw majorant is Gevrey two.  An \(O(M)\)
   bound requires an additional invariant slope cone or a new weighted
   unknown.
2. Full material phases absorb the only remaining \(O(K)\) inviscid
   term, namely low-mode transport of charged modes.  They do not absorb
   viscosity, whose principal size remains
   \(\nu K^2|F^{-T}w|^2\).
3. Under a material deformation \(F\), matched-child invertibility
   requires a metric transversality condition.  It is safe along the
   desired diagonal strain flow, but is not automatic under an arbitrary
   spatially varying low flow.
4. On the torus the lifted charge grading is not globally injective:
   \(Q(1,9,81)=182e_3\).  The first alias is at order \(O(K)\), so it does
   not affect any fixed-order audit, but an exact infinite construction
   must route it rather than pretend it is absent.

Thus the common-carrier shear is, at present, the cleanest principal
transition ansatz in this project.  The remaining gate is analytic and
geometric, not the finite-dimensional child rank or an immediate
\(K\)-stiff high--high chain.

The exact arithmetic and the quadratic multiplier example are checked
in
[`checks/single_carrier_shear_audit.py`](../checks/single_carrier_shear_audit.py).

---

## 1. Exact finite-\(K\) child map

For divergence-free Fourier amplitudes \(u_p\perp p\) and
\(u_r\perp r\), write the symmetrized projected Euler coefficient,
without its harmless factor \(i\), as

\[
 \mathcal B_{p+r}(u_p,u_r)
 =
 P_{p+r}\big[
   (u_p\cdot r)u_r+(u_r\cdot p)u_p
 \big].
\tag{1.1}
\]

Put

\[
 k=Kw,\qquad r=q-Kw,\qquad
 w=e_3,\qquad a=e_1.
\tag{1.2}
\]

Since \(a\cdot w=0\),

\[
 a\cdot r=a\cdot q.
\tag{1.3}
\]

Since \(b\cdot r=0\),

\[
 b\cdot k=K b\cdot w=b\cdot q.
\tag{1.4}
\]

Substitution in (1.1) proves the exact formula (0.1).

Write

\[
 q=(\alpha,\beta,\delta),\qquad
 Q_q=|q|^2.
\tag{1.5}
\]

The published directions all have \(\alpha\ne0\).  Convenient bases of
\(r^\perp\) and \(q^\perp\) are respectively

\[
\begin{aligned}
d_1&=(-\beta,\alpha,0),&
d_2&=(K-\delta,0,\alpha),\\
c_1&=(-\beta,\alpha,0),&
c_2&=(-\delta,0,\alpha).
\end{aligned}
\tag{1.6}
\]

Direct projection gives

\[
 L_{q,K}d_1=\alpha c_1
\tag{1.7}
\]

and

\[
 L_{q,K}d_2
 =
 -\frac{2\alpha\beta K}{Q_q}c_1
 \alpha\left(1-\frac{2K\delta}{Q_q}\right)c_2.
\tag{1.8}
\]

Thus, in the bases (1.6),

\[
 [L_{q,K}]
 =
 \begin{pmatrix}
 \alpha&-2\alpha\beta K/Q_q\\
 0&\alpha(Q_q-2K\delta)/Q_q
 \end{pmatrix},
\tag{1.9}
\]

and

\[
 \det [L_{q,K}]
 =
 \frac{\alpha^2(Q_q-2K\delta)}{Q_q}.
\tag{1.10}
\]

In oriented orthonormal bases of the two planes, the same determinant is

\[
 \det_\perp L_{q,K}
 =
 \frac{\alpha^2(Q_q-2K\delta)}
 {|q|\,|q-Kw|}.
\tag{1.11}
\]

Consequently the exact condition is

\[
 \boxed{\alpha\ne0,\qquad Q_q\ne2K\delta.}
\tag{1.12}
\]

The simpler condition
\((a\cdot q)(w\cdot q)\ne0\) is the correct large-\(K\) condition, but
at a fixed finite \(K\) one must also exclude the single exceptional
value \(K=Q_q/(2\delta)\).

For (0.2), the exceptional values are

\[
 \frac{3721}{40},\qquad
 \frac{61}{9},\qquad
 \frac32.
\tag{1.13}
\]

None is an integer.  Moreover (1.11) tends to the nonzero limit

\[
 -\frac{2\alpha^2\delta}{|q|}
\tag{1.14}
\]

as \(K\to\infty\).  The two-dimensional inverse is therefore uniformly
conditioned for all sufficiently large \(K\); the apparent growth in
the second column of (1.9) is only the growth of the unnormalized domain
basis vector \(d_2\).

---

## 2. Exact inherited strain rank

For each \(q_j\), let

\[
\begin{aligned}
c_{j,1}&=(-q_{j,2},q_{j,1},0),\\
c_{j,2}&=(-q_{j,3},0,q_{j,1}).
\end{aligned}
\tag{2.1}
\]

The six columns

\[
 \operatorname{sym}(c_{j,s}\otimes q_j),
 \qquad j=1,2,3,\quad s=1,2,
\tag{2.2}
\]

were previously checked to have rank five in
\(\mathrm{Sym}_0^3\).  In the coordinates

\[
 (M_{11},M_{22},M_{12},M_{13},M_{23}),
\tag{2.3}
\]

their first five columns have determinant

\[
 D_0=-1\,214\,003\,700.
\tag{2.4}
\]

The partner bases \(d_{j,s}\) in (1.6) map block-diagonally to the
child bases (2.1).  Hence the first-five-partner minor is

\[
 D_K
 =
 D_0\,
 \det[L_{q_1,K}]\,
 \det[L_{q_2,K}]\,
 (L_{q_3,K})_{11}.
\tag{2.5}
\]

Using

\[
\begin{aligned}
|q_1|^2&=3721,&
\det[L_{q_1,K}]
&=\frac{2025(3721-40K)}{3721},\\
|q_2|^2&=122,&
\det[L_{q_2,K}]
&=\frac{16(122-18K)}{122},\\
(L_{q_3,K})_{11}&=1,
\end{aligned}
\tag{2.6}
\]

gives exactly (0.3).  Therefore

\[
 \boxed{\operatorname{rank}
 \left[
 (b_{j,s})\mapsto
 \sum_j\operatorname{sym}
 \big(L_{q_j,K}b_{j,s}\otimes q_j\big)
 \right]=5}
\tag{2.7}
\]

for every positive integer \(K\).

In particular, the exact diagonal target

\[
 S_*=\operatorname{diag}\left(-1,-\frac54,\frac94\right)
\tag{2.8}
\]

can be produced by solving three independent two-dimensional child
systems.  No second carrier direction or helicity constraint is needed
for this endpoint.

---

## 3. Every charged--charged interaction is \(K\)-null

Let \(Q\) be the matrix with columns \(q_1,q_2,q_3\), and label a general
lifted Fourier mode by

\[
 (h,n)\in\mathbb Z\times\mathbb Z^3,\qquad
 p_K(h,n)=Kh\,w+\xi,\qquad \xi=Qn.
\tag{3.1}
\]

Take two nonzero fast charges,

\[
 p=Kh\,w+\xi,\qquad
 r=Kg\,w+\eta,\qquad h,g\ne0,
\tag{3.2}
\]

with amplitudes

\[
 A\cdot p=0,\qquad B\cdot r=0.
\tag{3.3}
\]

Dividing the two divergence identities by \(h\) and \(g\) gives

\[
\boxed{
\begin{aligned}
A\cdot r
 &=A\cdot\left(\eta-\frac gh\xi\right),\\
B\cdot p
 &=B\cdot\left(\xi-\frac hg\eta\right).
\end{aligned}}
\tag{3.4}
\]

There is no \(K\) on the right.  Therefore

\[
\begin{aligned}
\mathcal B_{p+r}(A,B)
=P_{p+r}\bigg[
&A\cdot\left(\eta-\frac gh\xi\right)B\\
&+
B\cdot\left(\xi-\frac hg\eta\right)A
\bigg].
\end{aligned}
\tag{3.5}
\]

This is an all-generation identity.  It includes:

* parent--partner matched beats;
* parent--conjugate-partner remote high modes;
* partner--partner high sums;
* partner--conjugate-partner low differences; and
* every later charged--charged product.

Thus the first- and second-generation \(O(K)\) chain found for the
two-carrier Beltrami star is structurally impossible here.  It relied on
different leading carrier directions.  With one common fast charge, the
divergence constraint cancels the fast derivative before any special
polarization or helicity identity is used.

---

## 4. The effective charge multiplier is generally quadratic

Projection has norm at most one, so (3.5) implies

\[
\begin{aligned}
|\mathcal B_{p+r}(A,B)|
\le |A||B|\bigg(
&|\eta|+\frac{|g|}{|h|}|\xi|\\
&+|\xi|+\frac{|h|}{|g|}|\eta|
\bigg).
\end{aligned}
\tag{4.1}
\]

If

\[
 |h|+|g|+|n|_1+|m|_1\le M,
\tag{4.2}
\]

then (4.1) is \(O(M^2)|A||B|\), and this is sharp in an unrestricted
charge ball.

Indeed, let \(q=q_1=(-45,-36,20)\), choose \(K=M^3\), and set

\[
\begin{aligned}
p&=Kw+Mq
&& (h=1,\ \xi=Mq),\\
r&=MKw
&& (g=M,\ \eta=0).
\end{aligned}
\tag{4.3}
\]

Take

\[
 A_M
 =
 e_1-\frac{M q_1^{(1)}}{K+M q_1^{(3)}}e_3
 =
 e_1+\frac{45M}{K+20M}e_3,
\qquad
B=e_1.
\tag{4.4}
\]

Then \(A_M\cdot p=0\), \(B\cdot r=0\), and
\(|A_M|,|B|\to1\).  The leading term in the raw interaction is

\[
 (A_M\cdot r)B
 =
 \frac{45M^2K}{K+20M}e_1.
\tag{4.5}
\]

The other term is only \(O(M)\), and the output direction is
asymptotic to \(e_3\), so Leray projection does not remove (4.5).
Consequently

\[
 \frac{
 |\mathcal B_{p+r}(A_M,B)|
 }{M^2}
 \longrightarrow45.
\tag{4.6}
\]

The exact checker gives, for \(M=20,40,80,160\),

\[
 40.605275,\quad
 43.319413,\quad
 44.297313,\quad
 44.683621.
\tag{4.7}
\]

In one-phase physical variables the same issue is visible from

\[
 w\cdot V_h
 =
 -\frac1K\partial_\theta^{-1}
 \operatorname{div}_{\!y}V_h,\qquad h\ne0.
\tag{4.8}
\]

The nominal fast term becomes schematically

\[
 (V\cdot w)K\partial_\theta V
 =
 -
 \big(\partial_\theta^{-1}
 \operatorname{div}_{\!y}V\big)
 \partial_\theta V.
\tag{4.9}
\]

The factor \(K\) is gone, but one slow derivative and one
\(\theta\)-derivative remain, distributed across the two inputs.
Ordinary isotropic analytic Cauchy estimates therefore lose two radii.

An \(O(M)\) bound would follow on a slope cone

\[
 |\xi|\le C|h|,\qquad |\eta|\le C|g|,
\tag{4.10}
\]

because the ratios in (4.1) would then be harmless.  But that cone is
not invariant under the real charge algebra.  With initial parent charge
\(P=(1,0)\) and partner charge \(S_1=(-1,e_1)\),

\[
 (M+1)P+M S_1=(1,Me_1),
\tag{4.11}
\]

whose slope grows like \(M\).  A polarization theorem could conceivably
show that all such algebraically allowed coefficients vanish, but no
such theorem is presently available, and generic partner polarizations
give no reason to expect it.

The safe current conclusion is therefore:

\[
\boxed{
\text{exactly \(K\)-null, but order two in unrestricted charge.}}
\tag{4.12}
\]

A Gevrey-two majorant is consistent with this conclusion.  A genuinely
analytic or one-derivative closure needs an additional normal form,
weighted unknown, or invariant support condition.

---

## 5. Full low-flow material phases

Let \(U(t,x)\) be the complete charge-zero velocity and let

\[
\partial_tX(t,y)=U(t,X(t,y)),\qquad F=D_yX.
\tag{5.1}
\]

For incompressible \(U\), \(\det F=1\).  Write \(y=X^{-1}(t,x)\).
Every lifted charged phase

\[
 \exp\big(i(Kh\,w+\xi)\cdot y\big)
\tag{5.2}
\]

satisfies

\[
 (\partial_t+U\cdot\nabla_x)
 \exp\big(i(Kh\,w+\xi)\cdot y\big)=0.
\tag{5.3}
\]

Thus the entire \(O(K)\) term

\[
 U\cdot\nabla V_{\ne0}
\tag{5.4}
\]

is absorbed, not just transport by one selected low Fourier mode.

The physical covector corresponding to (3.1) is

\[
 \pi_{h,n}=F^{-T}(Kh\,w+\xi).
\tag{5.5}
\]

Using the Piola amplitude \(FA\), one has

\[
 (FA)\cdot\pi_{h,n}
 =
 A\cdot(Kh\,w+\xi).
\tag{5.6}
\]

Hence the all-chain cancellation is invariant under arbitrary material
deformation:

\[
\begin{aligned}
(FA)\cdot F^{-T}(Kg\,w+\eta)
&=
A\cdot\left(\eta-\frac gh\xi\right),\\
(FB)\cdot F^{-T}(Kh\,w+\xi)
&=
B\cdot\left(\xi-\frac hg\eta\right).
\end{aligned}
\tag{5.7}
\]

The other low--high ordering,

\[
 V_{\ne0}\cdot\nabla U,
\tag{5.8}
\]

differentiates only the slow field and is \(O(\|\nabla U\|)\), not
\(O(K)\).  Variable-coefficient pressure and Piola factors introduce
condition numbers and slow derivatives of \(F\), but no new principal
fast derivative in the inviscid equation.

This identifies the correct second-chain verdict:

> after one common material phase is attached to all charged modes,
> neither high--high nor low--high inviscid coupling carries an
> uncancelled factor \(K\).

Viscosity is different.  Its material-coordinate principal symbol is

\[
 \nu\left|F^{-T}(Kh\,w+\xi)\right|^2,
\tag{5.9}
\]

so for fixed nonzero \(h\) and bounded slow charge it remains

\[
 \nu K^2h^2|F^{-T}w|^2+O(\nu K|\xi|).
\tag{5.10}
\]

The material transform does not solve the heat-clock problem.

---

## 6. Child rank after material deformation

There is one further condition hidden by evaluating the child map only
at \(F=I\).  Transport the parent and partner amplitudes by \(F\).
The physical child map is

\[
 b\longmapsto
 P_{F^{-T}q}\,
 F\big[\alpha b+(b\cdot q)a\big],
\qquad \alpha=a\cdot q.
\tag{6.1}
\]

Put

\[
 H=F^{-1}F^{-T}.
\tag{6.2}
\]

Using the cofactor identity for

\[
 T=\alpha I+a\otimes q,
\qquad
 \operatorname{cof}T
 =
 2\alpha^2I-\alpha q\otimes a,
\tag{6.3}
\]

the oriented area determinant of (6.1), up to strictly nonzero
normalizing factors, is

\[
 \alpha^2
 \left[
 (Hq)\cdot q-2K(Hq)\cdot w
 \right].
\tag{6.4}
\]

Therefore large-\(K\) uniform invertibility requires

\[
 |a\cdot q|\ge c,\qquad
 |(Hq)\cdot w|\ge c.
\tag{6.5}
\]

The original Euclidean condition \(q\cdot w\ne0\) does not imply the
second inequality for every positive-definite \(H\).  A sufficiently
rotated or sheared material metric can make \((Hq)\cdot w=0\).

For the desired diagonal strain flow, however, \(F\) and \(H\) remain
diagonal.  Then

\[
 (Hq_j)\cdot e_3=H_{33}q_{j,3},
\tag{6.6}
\]

which never vanishes at finite time because all three
\(q_{j,3}\ne0\).  The proposed target is therefore safe at the affine
centre.  A localized construction must keep (6.5) uniformly true
throughout the active support.

There is a related strain-chart condition.  Once the child maps are
onto, the physical child directions are

\[
 \widetilde q_j=F^{-T}q_j.
\tag{6.7}
\]

The symmetric chart has rank five if the only trace-free symmetric
matrix for which every \(\widetilde q_j\) is an eigenvector is zero.
A sufficient criterion is that the graph connecting nonorthogonal
pairs of the three \(\widetilde q_j\) be connected.

For diagonal positive \(H\),

\[
 q_1\cdot Hq_2
 =
 180(H_{11}+H_{22}+H_{33})>0.
\tag{6.8}
\]

Moreover \(q_3\) cannot be simultaneously \(H\)-orthogonal to both
\(q_1\) and \(q_2\): the two required equations would be

\[
\begin{aligned}
20H_{33}&=45H_{11}+36H_{22},\\
9H_{33}&=4H_{11}+5H_{22},
\end{aligned}
\tag{6.9}
\]

which are incompatible for positive \(H_{11},H_{22}\).  Hence the
rank-five strain chart also survives every positive diagonal material
deformation.  It can fail for a general \(F\), since some positive
metric can make any prescribed basis orthogonal.

---

## 7. The torus alias at order \(O(K)\)

The three slow directions satisfy the exact integer identity

\[
 Q(1,9,81)=182e_3.
\tag{7.1}
\]

Consequently the lifted map

\[
 (h,n)\longmapsto Kh e_3+Qn
\tag{7.2}
\]

cannot be globally injective.  If \(d=\gcd(K,182)\), its primitive
integer kernel is generated by

\[
 h_*=\frac{182}{d},\qquad
 n_*=-\frac Kd(1,9,81).
\tag{7.3}
\]

The \(\ell^1\) charge size of this alias is

\[
 |h_*|+|n_*|_1
 =
 \frac{182+91K}{d}.
\tag{7.4}
\]

Thus a fixed perturbative order \(M\ll K\) sees an injective grading
(with the precise safe range depending on \(d\)), but the complete
infinite series does not.  This matters because a lifted mode with
nonzero \(h\) can become a zero physical Fourier mode, and the exact
viscous symbol in (5.9) then cancels rather than behaving like
\(\nu K^2h^2\).

There are three honest ways to proceed:

1. keep the exact combined symbol \(Kh e_3+Qn\) and route aliases in the
   infinite charged system;
2. work first on the lifted four-dimensional phase space and prove that
   evaluation on \(\theta=Kx_3\) remains controlled despite the kernel;
3. move the localized construction to \(\mathbb R^3\), where continuous
   slow frequencies need not create an exact integer alias.

What is not legitimate is to use finite-order injectivity as if it were
an exact all-order grading on \(\mathbb T^3\).

---

## 8. What this does and does not establish

The audit establishes the following exact finite-dimensional and
symbolic module.

1. One real shear parent and six sideband-polarization variables create
   all five trace-free strain coordinates.
2. This remains true at every positive integer carrier frequency.
3. All charged--charged inviscid interactions are exactly free of the
   carrier factor \(K\), at every generation.
4. A single full low-flow material map removes the remaining \(O(K)\)
   low--high transport.
5. The same cancellations survive material deformation, provided the
   child transversality and strain-rank conditions remain uniform.

This is substantially stronger than the earlier Beltrami-star
parametrix.  It does **not** yet establish:

* a one-derivative analytic inverse;
* convergence or Borel summability of the Gevrey-two charged expansion;
* control of the torus aliases;
* a localized material flow preserving the metric transversality
  everywhere;
* the viscosity/heat-clock inequalities;
* routing of the physical wake and flat forcing; or
* a complete scale-recursive Navier--Stokes solution.

The most decisive next analytic experiment is therefore not another
finite shell optimization.  It is a small exact/Galerkin implementation
of the **lifted one-phase charged equations in material coordinates**,
with:

* the full \(h=0\) velocity used to update the material map;
* all \(h\ne0\) interactions evaluated through (3.5);
* exact physical symbols \(Kh e_3+Qn\), so aliases are visible;
* a comparison of analytic, slope-weighted, and Gevrey-two norms; and
* the viscous symbol (5.9) retained without asymptotic replacement.

That computation would test the remaining plausible breakthrough
claim: whether a weighted normal form improves the sharp \(O(M^2)\)
charge multiplier to a summable one-stage transition while preserving
the exact rank-five child chart.
