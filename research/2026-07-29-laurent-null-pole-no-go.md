# No irreducible Laurent pole divisor is null for the complexified Laplacian

Date: 2026-07-29

## Decision

The meromorphic-pole route in
`2026-07-29-flat-force-borel-attack.md` required an irreducible Laurent
polynomial \(q\) satisfying

\[
 q\ \big|\ \sum_{j=1}^{3}(D_jq)^2,
 \qquad D_j=z_j\partial_{z_j}.                    \tag{0.1}
\]

The earlier audit excluded binomials and irreducible trinomials.  The
result below excludes **every** nonunit irreducible Laurent polynomial,
with any finite number of monomials.  Thus there is no four-monomial
candidate to find: the finite-Laurent simple-pole search is closed before
the Navier--Stokes residue equations are reached.

This is not a Navier--Stokes regularity theorem.  It is an exact no-go for
one proposed singular ansatz.

## 1. The algebraic theorem

Let

\[
 R=\mathbb C[z_1^{\pm1},z_2^{\pm1},z_3^{\pm1}]
\]

and equip \(\mathbb C^3\) with the complex-bilinear extension of the
Euclidean dot product,

\[
 v\mathbin{\cdot}w=\sum_{j=1}^{3}v_jw_j.
\]

> **Theorem.**  If \(q\in R\) is a nonunit irreducible Laurent
> polynomial, then
> \[
> q\nmid \sum_{j=1}^{3}(D_jq)^2.                 \tag{1.1}
> \]

The vector

\[
 [D_1q(z):D_2q(z):D_3q(z)]
\]

at a smooth point of \(V(q)\subset(\mathbb C^*)^3\) is the logarithmic
Gauss map.  This standard interpretation is recalled, for example, in
Madani--Nisse, *Generalized logarithmic Gauss map and its relation to
(co)amoebas*, arXiv:1205.2917.  The proof below is self-contained.

### Step 1: divisibility makes the lifted hypersurface null

Assume for contradiction that

\[
 \sum_j(D_jq)^2=qh                              \tag{1.2}
\]

for some \(h\in R\).  Pull this identity to the logarithmic cover by
putting

\[
 Q(y)=q(e^{y_1},e^{y_2},e^{y_3}),\qquad
 H(y)=h(e^{y_1},e^{y_2},e^{y_3}).
\]

Then

\[
 |\nabla Q|^2=QH.                               \tag{1.3}
\]

At every smooth point \(y_0\) of \(Q=0\), the nonzero normal

\[
 n_0=\nabla Q(y_0)
\]

is null:

\[
 n_0\mathbin{\cdot}n_0=0.                       \tag{1.4}
\]

### Step 2: every null-normal trajectory is a straight complex line

Let \(y(s)\) solve

\[
 y'(s)=\nabla Q(y(s)),\qquad y(0)=y_0.
\]

Equation (1.3) gives

\[
 \frac{d}{ds}Q(y(s))=Q(y(s))H(y(s)),
\]

so a trajectory starting on \(Q=0\) stays on \(Q=0\).  Differentiating
(1.3) and then restricting to \(Q=0\) gives

\[
 2\,\operatorname{Hess}(Q)\nabla Q=H\nabla Q.    \tag{1.5}
\]

Consequently

\[
 \frac{d}{ds}\nabla Q(y(s))
 =\frac{H(y(s))}{2}\nabla Q(y(s)).               \tag{1.6}
\]

The normal direction is therefore constant along the trajectory.
After a local reparametrization, the trajectory is

\[
 y_0+t n_0.
\]

It follows first for small \(t\), and then for every \(t\in\mathbb C\) by
the identity theorem, that

\[
 Q(y_0+t n_0)=0.                                \tag{1.7}
\]

Thus the logarithmic lift of the hypersurface is ruled by affine null
lines.

### Step 3: finite Fourier support forces only finitely many normal directions

Write

\[
 Q(y)=\sum_{m\in A}c_m e^{m\cdot y},
\]

where \(A\subset\mathbb Z^3\) is finite and has at least two elements.
Restricting this exponential polynomial to the line (1.7) gives

\[
 \sum_{m\in A}
 c_m e^{m\cdot y_0}e^{(m\cdot n_0)t}=0
 \quad\hbox{for all }t.                          \tag{1.8}
\]

Exponentials with distinct exponents are linearly independent.  Since
every coefficient \(c_m e^{m\cdot y_0}\) is nonzero, (1.8) implies that
there are distinct \(m,m'\in A\) such that

\[
(m-m')\mathbin{\cdot}n_0=0.                    \tag{1.9}
\]

More explicitly, partition \(A\) into groups on which \(m\cdot n_0\) is
constant.  Linear independence says that the coefficient sum in every
group is zero; no singleton can have zero sum, so at least one tied pair
must occur.

Projectively, every value of the logarithmic Gauss map therefore lies in

\[
 \{[v]\in\mathbb P^2:v\cdot v=0\}
 \ \cap\
 \bigcup_{\substack{m,m'\in A\\m\ne m'}}
 \{[v]:(m-m')\cdot v=0\}.                       \tag{1.10}
\]

The first set is a nondegenerate conic.  Each set in the finite union is
a projective line, so (1.10) is finite.  The smooth locus of the
irreducible hypersurface \(V(q)\) is irreducible and Zariski dense.
Hence its logarithmic Gauss map, whose image is irreducible, must be
constant:

\[
 [D_1q:D_2q:D_3q]=[v]                           \tag{1.11}
\]

for one fixed nonzero null vector \(v\).

### Step 4: a constant logarithmic normal contradicts the integer lattice

On a local logarithmic sheet of \(V(q)\), (1.11) says that every tangent
plane is

\[
 v^\perp=\{w:v\cdot w=0\}.
\]

The sheet is therefore an open subset of an affine plane

\[
 \{y:v\cdot y=c\}.                               \tag{1.12}
\]

The restriction of \(Q\) vanishes identically on this plane.  Linear
independence of exponential functions on the two-dimensional plane again
implies that two distinct exponents \(m,m'\in A\) have the same
restriction to it.  Thus

\[
 d:=m-m'\in (v^\perp)^\perp=\mathbb C v.         \tag{1.13}
\]

But \(d\) is a nonzero real integer vector.  If \(d=\lambda v\), then

\[
 0=\lambda^2(v\cdot v)=d\cdot d
   =d_1^2+d_2^2+d_3^2,
\]

which is impossible for \(d\in\mathbb Z^3\setminus\{0\}\).  This
contradiction proves the theorem.

## 2. Consequence for the pole-collision attack

For a regular simple pole

\[
 u\sim a/q
\]

the audited leading equations gave \(a\cdot\nabla q=0\).  The nominal
most singular self-advection then vanishes, while a pressure pole is
normal to \(q=0\).  Unless another pole component supplies a same-order
tangential cancellation, viscosity requires the complex normal to be
null.  For a Laurent pole divisor this is precisely (0.1).

The theorem therefore rules out the entire **single irreducible
finite-Laurent simple-pole divisor** mechanism, regardless of the number
of monomials or their coefficients.  An exhaustive four-monomial
computer-algebra search would add no information.

The same conclusion applies to a reduced simple pole written with a
reducible Laurent denominator: at a generic point of any irreducible
factor, the other factors are nonzero units, so that factor is the local
pole divisor and must satisfy the same null-normal condition.

The result does not rule out:

1. a genuinely coupled collection of pole components producing a
   same-order tangential cancellation before the null condition is
   imposed;
2. higher-order poles whose preceding residue hierarchy changes the
   leading balance;
3. non-Laurent analytic or essential-singularity divisors; or
4. the packet-plus-wake/Gevrey return construction studied elsewhere in
   this repository.

The first three possibilities no longer have a concrete closure equation
or energy-dissipation ledger.  They should not displace the one-carrier
Gevrey cell problem unless such a structure is found.

## Reference

* F. Madani and M. Nisse, “Generalized logarithmic Gauss map and its
  relation to (co)amoebas,” arXiv:1205.2917,
  <https://arxiv.org/abs/1205.2917>.
