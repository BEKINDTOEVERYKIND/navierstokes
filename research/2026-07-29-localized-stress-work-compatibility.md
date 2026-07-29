# Localized high--high reset: an explicit six-wave cell and an unavoidable work-carrying wake

Date: 2026-07-29

## Claim boundary

This note addresses only the algebraic Reynolds-stress part of a proposed
localized Kelvin cascade.  It does not construct an Euler return orbit or a
Navier--Stokes singularity.

There are two complementary conclusions.

1. At one point, or in an affine flow with a boundary energy flux, six
   linearly polarized Kelvin-wave families give an explicit positive
   covariance which extracts energy from
   \[
   S=\operatorname{diag}(-\alpha,-\beta,\alpha+\beta),
   \qquad \alpha,\beta>0,
   \]
   represents an arbitrary sufficiently small daughter stress modulo
   pressure, and has zero helicity.  Six is minimal for an exact fixed-ray
   chart of a full \(\mathrm{Sym}^3\) neighborhood, and for a positive
   fixed-ray chart through the isotropic class modulo pressure.
2. A strictly localized stress which is required to have a prescribed
   Leray divergence cannot independently choose its work against the parent.
   Its work is fixed by that divergence.  In particular, the pure
   symmetric-gradient daughter stress used in the
   Cheskidov--Dai--Palasek inverse cascade does zero total work against an
   exactly constant trace-free strain.  Cutting off the extracting part
   necessarily creates a boundary source.  That source is the first
   mathematically unavoidable component of the outgoing wake.

Thus positivity, helicity, polarization transport, and quadratic resonance
do not kill the averaged cell.  The no-wake localization requirement does.

## 1. The stress--work identity

Let \(U\) be a smooth divergence-free vector field and let \(R=R^T\) be a
smooth symmetric tensor, either periodic or sufficiently decaying that no
boundary term remains.  Write

\[
f=\mathbb P\operatorname{div}R .
\]

Since

\[
\operatorname{div}R=f+\nabla q
\]

for some scalar \(q\), and since \(R\) is symmetric,

\[
\boxed{
\int R:\operatorname{sym}\nabla U
=-\int f\cdot U .
}
\tag{1.1}
\]

Indeed,

\[
\int R:\operatorname{sym}\nabla U
=\int R:\nabla U
=-\int(\operatorname{div}R)\cdot U
=-\int f\cdot U .
\]

Consequences of (1.1) are exact.

* If two localized stresses have the same Leray divergence, then they do
  exactly the same total work on \(U\).
* Adding an arbitrary pressure gauge \(\rho(x)I\) can make a stress positive
  without altering its Leray divergence or its work, since
  \(\operatorname{tr}\nabla U=0\).
* Adding a localized source-dark stress
  \(K\), with \(\mathbb P\operatorname{div}K=0\), cannot change the work.

In the Euler--Reynolds convention

\[
\partial_tU+\mathbb P\operatorname{div}(U\otimes U+R)=0,
\]

the low field receives energy at the rate

\[
\int R:\operatorname{sym}\nabla U,
\]

and the unresolved high field receives the opposite amount.  Hence the
high field extracts energy precisely when the left side of (1.1) is
negative.  This sign is a compatibility condition on \(f\), not a freely
adjustable feature of a particular rank-one factorization.

## 2. Constant-strain first-moment obstruction

Take the affine parent

\[
U_0(x)=Sx,\qquad \operatorname{tr}S=0,
\]

and a compactly supported stress on \(\mathbb R^3\).  Formula (1.1) becomes

\[
\int R:S=-\int f(x)\cdot Sx\,dx .
\tag{2.1}
\]

Thus the work is determined by the first moment of the prescribed daughter
source.

Now let \(V\) be compactly supported or periodic and consider a CDP-type
daughter stress

\[
R_{\mathrm{dau}}
=c\,\operatorname{sym}\nabla V+pI ,
\tag{2.2}
\]

with any constant \(c\).  Its contraction with the constant strain is

\[
\int R_{\mathrm{dau}}:S
=c\int\operatorname{sym}\nabla V:S
=0.
\tag{2.3}
\]

Equivalently, its Leray divergence is a constant multiple of
\(\Delta\mathbb PV\), whose zeroth and first moments vanish, so the
right-hand side of (2.1) is zero.

The precise CDP identity has the form

\[
\sum_j a_j^2\theta_j\otimes\theta_j
=2\mathcal D V+pI,\qquad
\operatorname{div}\mathcal D=\Delta\mathbb P .
\]

Therefore that rank-one technology can prescribe the daughter strain, but
the resulting pure daughter stress cannot simultaneously draw nonzero net
energy from a homogeneous affine strain.  This is not a positivity
obstruction; it is the duality identity (1.1).

There is a useful quantitative version.  Suppose the stress is supported in
a ball of diameter \(L\), \(S_U(x_0)=S\), and the daughter part has zero
contraction with \(S\).  Then

\[
\left|\int R_{\mathrm{dau}}:S_U(x)\,dx\right|
\leq
C L\|\nabla S_U\|_{L^\infty}
\|R_{\mathrm{dau}}\|_{L^1}.
\tag{2.4}
\]

At daughter scale \(L=r^{-1}\), work obtained only from variation of the
parent strain is suppressed by \(r^{-1}\).  An order-one handoff must either
amplify the stress by that missing factor, retain a boundary/wake source, or
abandon the locally homogeneous-strain model.

## 3. Explicit six-ray rank-one dictionary

Define

\[
\begin{array}{lll}
\theta_{12}^{+}=(e_1+e_2)/\sqrt2,&
\theta_{12}^{-}=(e_1-e_2)/\sqrt2,\\[2mm]
\theta_{13}^{+}=(e_1+e_3)/\sqrt2,&
\theta_{13}^{-}=(e_1-e_3)/\sqrt2,\\[2mm]
\theta_{23}^{+}=(e_2+e_3)/\sqrt2,&
\theta_{23}^{-}=(e_2-e_3)/\sqrt2 .
\end{array}
\tag{3.1}
\]

For a symmetric matrix \(R=(r_{ij})\), set

\[
\begin{aligned}
A&=r_{11}+r_{22}-r_{33},\\
B&=r_{11}+r_{33}-r_{22},\\
C&=r_{22}+r_{33}-r_{11},
\end{aligned}
\]

and

\[
\begin{array}{ll}
c_{12}^{\pm}=A/2\pm r_{12},&
c_{13}^{\pm}=B/2\pm r_{13},\\[1mm]
c_{23}^{\pm}=C/2\pm r_{23}.&
\end{array}
\tag{3.2}
\]

A direct calculation gives

\[
\boxed{
R=\sum_{ab\in\{12,13,23\}}\sum_{\sigma=\pm}
c_{ab}^{\sigma}
\theta_{ab}^{\sigma}\otimes\theta_{ab}^{\sigma}.
}
\tag{3.3}
\]

All six coefficients are positive in a neighborhood of \(I\).  More
generally, take the biased positive matrix

\[
R_*=\rho I+\kappa e_1\otimes e_1,\qquad 0<\kappa<\rho.
\tag{3.4}
\]

Then

\[
R_*:S=-\alpha\kappa<0,
\tag{3.5}
\]

and its coefficients are

\[
c_{12}^{\pm}=c_{13}^{\pm}=(\rho+\kappa)/2,\qquad
c_{23}^{\pm}=(\rho-\kappa)/2.
\tag{3.6}
\]

Consequently every tensor \(R_*+G\), for \(G\) sufficiently small, has a
smooth positive six-ray factorization and still extracts energy from \(S\).
This is an explicit version of the local rank-one geometric lemma.

The extraction in (3.5) is a pointwise affine-cell statement.  A spatially
constant \(R_*\) does zero total work on a periodic parent because
\(\int_{\mathbb T^3}\nabla U=0\).  To obtain nonzero global work one needs
an affine boundary flux or a localized envelope, and the latter is exactly
what produces the source in Section 7.

Six fixed rays are minimal for an exact chart of a full neighborhood in
\(\mathrm{Sym}^3\), simply by dimension.  They are also minimal for a
positive fixed-ray chart *through the isotropic class modulo pressure*:
the quotient \(\mathrm{Sym}^3/\mathbb RI\) has dimension five, and with
only five fixed rays a full-rank coefficient map to this quotient has zero
kernel.  It therefore cannot also have a strictly positive coefficient
vector mapping to the isotropic class.  This does not rule out a five-ray
quotient chart around a non-isotropic biased base such as \(R_*\).
Three space-dependent eigenvectors suffice for one positive matrix, but do
not give a fixed modulated-plane-wave dictionary.

## 4. Six linearly polarized waves and zero helicity

Choose carrier directions

\[
\begin{array}{lll}
n_{12}^{+}=e_3,&n_{12}^{-}=3e_3,\\
n_{13}^{+}=e_2,&n_{13}^{-}=3e_2,\\
n_{23}^{+}=e_1,&n_{23}^{-}=3e_1.
\end{array}
\tag{4.1}
\]

Each carrier is orthogonal to its polarization.  The six real waves

\[
W_{ab}^{\sigma}(x)
=\sqrt{2c_{ab}^{\sigma}}\,
\theta_{ab}^{\sigma}
\cos(\lambda n_{ab}^{\sigma}\cdot x+\phi_{ab}^{\sigma})
\tag{4.2}
\]

have

\[
\left\langle
\sum_{ab,\sigma}W_{ab}^{\sigma}\otimes W_{ab}^{\sigma}
\right\rangle_{\mathrm{low}}
=R.
\tag{4.3}
\]

The factor-three colour separates the two waves with parallel carrier
axes.  Spatially disjoint CDP colours can be used instead if equal carrier
radii are important.

Every wave in (4.2) has zero helicity pointwise:

\[
W_{ab}^{\sigma}\cdot
\operatorname{curl}W_{ab}^{\sigma}=0.
\]

The carrier supports are distinct, so all cross-helicity integrals vanish
by Fourier orthogonality.  Hence the total helicity is exactly zero, not
merely small.

For a modulated block, one may use the velocity-potential form

\[
W_{\lambda}
=-\mathbb P\Delta\left(
\lambda^{-2}a(x)\theta\sin(\lambda n\cdot x)
\right).
\tag{4.4}
\]

For each individual block the Fourier coefficient at wave number \(k\) is
a scalar multiple of the real vector \(\mathbb P_k\theta\).  Therefore

\[
\widehat W(k)^*\cdot(ik\times\widehat W(k))=0
\]

at every \(k\), and the block has exactly zero integrated helicity.  If the
six modulated carrier bands are disjoint, their sum also has zero integrated
helicity.  Compact spatial cutoffs destroy exact spectral disjointness, but
an improper reflection symmetry, such as \(u(-x)=-u(x)\), cancels total
helicity exactly.

Thus helicity does not obstruct the averaged reset.

## 5. Exact affine polarization transport

For a Kelvin wave

\[
w=\operatorname{Re}\{a(t)e^{ik(t)\cdot x}\},
\qquad k\cdot a=0,
\]

in the affine Euler field \(Sx\),

\[
k'=-S^Tk,\qquad
a'=-Sa+2k\frac{k\cdot Sa}{|k|^2}.
\tag{5.1}
\]

For the carriers in (4.1), \(k\) remains on a coordinate axis and
\(a\) remains in its orthogonal coordinate plane.  Hence the pressure term
in (5.1) vanishes.  If the initial covariance is (3.4), the six transported
waves retain the diagonal covariance

\[
R_*(t)=
\operatorname{diag}\left(
(\rho+\kappa)e^{2\alpha t},
\rho e^{2\beta t},
\rho e^{-2(\alpha+\beta)t}
\right).
\tag{5.2}
\]

Its extraction rate is

\[
-R_*(t):S
=
\alpha(\rho+\kappa)e^{2\alpha t}
+\beta\rho e^{2\beta t}
-(\alpha+\beta)\rho e^{-2(\alpha+\beta)t}.
\tag{5.3}
\]

At \(t=0\) this equals \(\alpha\kappa>0\), and its derivative is strictly
positive.  Polarization transport therefore strengthens, rather than
reverses, the extraction for this explicit cell.

Two carrier pairs in (4.1) move to higher wave number under the compressive
axes; the pair carried by \(e_3\) moves downward.  The latter pair is an
auxiliary tensor-shaping colour, not a demonstrated forward-cascade child.
Replacing it while retaining a six-dimensional positive stress chart is a
dynamical design problem, not an algebraic no-go.

## 6. Resonance and finite-Fourier closure

The six carrier lines in (4.1) have no quadratic interaction which lands
back on a selected carrier line at the selected radius: the same-axis
difference is \(2\lambda e_i\), and mixed-axis sums have two nonzero
coordinates.  Thus quadratic resonance can be kept outside the child band.
Alternatively, in three dimensions the physical supports of different
colours can be separated as in the CDP construction.

This does not make the six waves an exact finite-dimensional Euler
solution.  Their quadratic products create new carrier sums and
differences, and further interactions generate an expanding lattice.
Kishimoto--Yoneda's classification of finite-Fourier-support 3D Euler flows
shows why exact closure is the wrong target: apart from stationary
two-dimensional-like and Beltrami cases, finite-mode Euler dynamics of the
required kind do not exist.

The admissible object is therefore a modulated, infinite-tail
packet-plus-wake cell.  The six waves specify only its principal stress.

## 7. Why localization forces a wake

The extracting bias \(R_*\) in (3.4) is source-dark only while it is
spatially constant.  Localizing it with a cutoff gives

\[
\mathbb P\operatorname{div}(\chi^2R_*)
=\mathbb P\big(R_*\nabla\chi^2\big).
\tag{7.1}
\]

This boundary source is not optional.  If one adds a localized correction
to cancel (7.1) while preserving the pure prescribed daughter source, the
stress--work identity forces that correction to cancel the work of the
localized extracting bias as well.

Accordingly, a realizable localized principal stress has the form

\[
R_{\mathrm{loc}}
=\chi^2 R_*+R_{\mathrm{dau}}+\rho(x)I,
\tag{7.2}
\]

and its Leray divergence splits into

\[
\mathbb P\operatorname{div}R_{\mathrm{loc}}
=
\underbrace{\mathbb P(R_*\nabla\chi^2)}_{\text{work-carrying boundary/wake}}
+
\underbrace{\mathbb P\operatorname{div}R_{\mathrm{dau}}}_{\text{prescribed child}}.
\tag{7.3}
\]

The first term supplies exactly the first moment required by (2.1).  It must
be retained in the outgoing state.  Calling it an error and forcing it to
vanish removes the energy extraction.

This gives a sharper theorem target:

> Construct a time-dependent localized Euler packet whose high-frequency
> covariance has the six-ray principal part above, whose projected interior
> source is the next daughter, and whose boundary source is propagated as a
> finite-energy wake carrying the exact work and moment balances.

The wake must then be included in the renormalized endpoint map.  It cannot
be hidden in a flat external force.

## 8. Separate restriction on the parent

The affine field is only a local model.  A recent rigidity theorem of
Peralta-Salas and Slobodeanu proves that analytic localizable steady Euler
flows in bounded domains are axisymmetric, with a rotationally symmetric
domain whose transverse section is a disk or annulus with convex boundary.
Consequently a compact Gavrilov-style steady parent cannot be treated as a
freely designable realization of an arbitrary localized biaxial affine
strain together with this reset.  A viable parent may have to be
time-dependent, nonlocalizable, or explicitly packet-plus-wake.

## Conclusion

An explicit six-ray positive stress dictionary is available, and the usual
suspects do not kill it:

* positivity is handled by the biased positive base \(R_*\);
* helicity is zero by linear polarization or one improper symmetry;
* Kelvin transport preserves the strict extraction;
* quadratic resonances can be coloured away.

The exact obstruction is instead

\[
\int R:\operatorname{sym}\nabla U
=-\int(\mathbb P\operatorname{div}R)\cdot U.
\]

A pure localized CDP daughter source has no independent energy-extraction
degree of freedom.  The boundary/wake source created by localizing the
anisotropic drain is therefore necessary at leading order.

## Primary sources

* A. Cheskidov, M. Dai, and S. Palasek, *Instantaneous Type I blow-up and
  non-uniqueness of smooth solutions of the Navier--Stokes equations*:
  https://arxiv.org/abs/2511.09556
* N. Kishimoto and T. Yoneda, *Characterization of three-dimensional Euler
  flows supported on finitely many Fourier modes*:
  https://arxiv.org/abs/2110.08039
* D. Peralta-Salas and R. Slobodeanu, *A symmetry theorem for localizable
  steady solutions of the 3D Euler equations*:
  https://arxiv.org/abs/2606.13462
