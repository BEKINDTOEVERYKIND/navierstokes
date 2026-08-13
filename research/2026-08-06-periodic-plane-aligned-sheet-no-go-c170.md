# C170: the C159 periodic plane cannot carry a \(q^2\) aligned lattice sheet

**Date:** 2026-08-06

**Status:** exact periodic-plane geometry, projected-lattice capacity, and
fiber-shear obstruction; qualitative principal-cocycle perturbation only

**Checker:**
[checks/periodic_plane_aligned_sheet_no_go_c170.py](../checks/periodic_plane_aligned_sheet_no_go_c170.py)

## 0. Claim boundary

This note tests one narrowly defined repair of the C161 terminal charge
star. One might try to place \(q^2\) amplified sources on the C159
periodic-covector plane while making their horizontal projections parallel,
or confining them to a sector of aperture \(O(q^{-1})\). One fixed
horizontal pure-normal gate polarization would then be source-adapted, up to
\(O(q^{-1})\), across the whole source set.

At frequency \(O(q)\), this repair is impossible. Precisely:

1. horizontal projection is a linear isomorphism on the C152
   periodic-covector plane, so an exactly parallel horizontal family in
   that plane is a single real line;
2. for \(0<\delta\le\delta_0<\pi/2\), any angular \(\delta\)-sector of
   projected integer covectors of horizontal size \(O(q)\) contains
   \(O_{R,\delta_0}(\delta q^2+q)\) distinct points. Thus a
   \(C/q\)-sector contains only \(O_{R,C}(q)\) points. The same is true
   for the exact periodic-plane subset;
3. the parallel real ray through the C159 covector has exactly the same
   Kelvin monodromy and expanding physical-amplitude line at every nonzero
   scale. It contains a nonzero integer frequency if and only if an
   unresolved period-derivative ratio is rational, and even then contains
   only \(O(q)\) frequencies of size \(O(q)\);
4. using integer normal-charge lifts off the periodic plane can naively
   restore another factor \(q\), but C154 shears those lifts in the
   transverse horizontal direction. To remain in a \(C/q\)-sector through
   \(\Theta(\log q)\) returns, at most one integer lift over each projected
   lattice point is possible for all sufficiently large \(q\). The total
   compatible capacity is again \(O(q)\).

The exact-plane \(C/q\)-sector inherits a uniform expanding principal
cocycle from C159 for sufficiently large \(q\), by smooth perturbation of a
simple hyperbolic multiplier. Its principal gain varies by only
\(1+O((\log q)/q)\) through \(O(\log q)\) returns. This is a statement about
the periodic Kelvin cocycle, not a finite-frequency localized packet.
Nothing here is a no-go for a full two-polarization, noncommuting, unfolded
terminal map.

## 1. Exact coordinates for the periodic plane

Retain the C152 vectors

\[
 N=(1,1,1),\qquad r=r_3=(-1,0,1),\qquad
 d=(-1,2,-1)=-r_1+r_2.                              \tag{1.1}
\]

They are mutually orthogonal, with

\[
 |N|^2=3,\qquad |r|^2=2,\qquad |d|^2=6.             \tag{1.2}
\]

At the algebraic base point of C152, put

\[
 \gamma={\sqrt{21}\over5},\qquad
 g_0=\gamma d,\qquad U_0=3\gamma r.                 \tag{1.3}
\]

The rank-one physical return is \(F=I+u\otimes g_0\), where

\[
 u=-T'(0)U_0-\sqrt2T(0)N=a_0r-b_0N,
 \quad a_0=-3\gamma T'(0)>0,
 \quad b_0=\sqrt2T(0)>0.                            \tag{1.4}
\]

Define

\[
 \boxed{\sigma={2a_0\over3b_0}
 =-{\sqrt{42}\over5}{T'(0)\over T(0)}>0.}           \tag{1.5}
\]

Writing a vector in the orthogonal basis as \(xr+yd+zN\), the C152
periodic-covector plane is

\[
 \begin{aligned}
 \mathcal P=u^\perp
 &=\{xr+yd+zN:2a_0x-3b_0z=0\}\\
 &=\boxed{\operatorname{span}_{\mathbb R}
          \{d,r+\sigma N\}}.                       \tag{1.6}
 \end{aligned}
\]

Let

\[
 \pi_Hk=k-{N\cdot k\over3}N                       \tag{1.7}
\]

be projection onto the \(A_2\) plane \(N^\perp\). The kernel of \(\pi_H\)
is \(\mathbb RN\), while \(u\cdot N=-3b_0\ne0\). Therefore

\[
 \boxed{\pi_H|_{\mathcal P}:\mathcal P\longrightarrow N^\perp
        \text{ is a linear isomorphism}.}           \tag{1.8}
\]

Its inverse is

\[
 \mathcal I_\sigma(xr+yd)=x(r+\sigma N)+yd.         \tag{1.9}
\]

Thus a fixed horizontal direction determines the normal charge inside
\(\mathcal P\); the compatible set is one real line, not a sheet.

### 1.1 Arithmetic caveat

For an integer vector \(k\), its coefficients in the orthogonal
\((r,d,N)\) basis are rational. Hence

\[
 \sigma\notin\mathbb Q
 \quad\Longrightarrow\quad
 \mathcal P\cap\mathbb Z^3=\mathbb Zd.             \tag{1.10}
\]

Indeed, the plane condition is \(z=\sigma x\); irrationality forces
\(x=z=0\), and the primitive integer direction left is \(d\). If
\(\sigma=p/s\in\mathbb Q\) in lowest terms, then
\(sr+pN\in\mathcal P\cap\mathbb Z^3\), so the intersection is a rank-two
lattice. No rationality decision for the ratio in (1.5) is supplied by
C152, C159, or this note.

This distinction matters below. The projected ambient lattice
\(\pi_H(\mathbb Z^3)\) has rank two, but not every one of its points has an
integer lift in the exact plane \(\mathcal P\). Counts for the ambient
projected lattice are unconditional upper bounds for the exact-plane
subset, not assertions of exact-plane occupancy.

## 2. Projected-lattice sector capacity

The horizontal projections of integer frequencies form the weight lattice

\[
 \Lambda_H=\pi_H(\mathbb Z^3)\subset N^\perp.       \tag{2.1}
\]

For \(k=(k_1,k_2,k_3)\in\mathbb Z^3\), define the two integer coordinates

\[
 A(k)=r\cdot k=k_3-k_1,
 \qquad D(k)=d\cdot k=-k_1+2k_2-k_3.                \tag{2.2}
\]

Orthogonality gives the exact formulas

\[
 \boxed{\pi_Hk={A(k)\over2}r+{D(k)\over6}d,}
 \qquad
 \boxed{|\pi_Hk|^2={A(k)^2\over2}+{D(k)^2\over6}.} \tag{2.3}
\]

The pair \((A(k),D(k))\) determines \(\pi_Hk\). Its image consists exactly
of the integer pairs with \(A(k)\equiv D(k)\pmod2\): necessity follows from
\(D-A=2(k_2-k_3)\), and \((0,(A+D)/2,A)\) realizes every such pair. This
parity constraint can only reduce the upper count below. Formula (2.3) then
shows that distinct points of \(\Lambda_H\) are separated and that the least
nonzero squared norm is \(2/3\).

Fix constants \(0<r_-<R\), and let \(0<\delta\le\delta_0<\pi/2\). Consider
the two-sided sector annulus

\[
 \mathcal S_{q,\delta}=
 \left\{h\in\Lambda_H:
 r_-q\le|h|\le Rq,
 \ \angle(h,\mathbb Rr)\le\delta\right\}.          \tag{2.4}
\]

For every \(h=\pi_Hk\) in this set,

\[
 |A(k)|\le\sqrt2Rq,
 \qquad
 |D(k)|\le\sqrt6Rq\sin\delta.                     \tag{2.5}
\]

Counting the possible integer pairs and using
\(\sin\delta\le\delta\) gives

\[
 \boxed{\#\mathcal S_{q,\delta}
 \le (2\lfloor\sqrt2Rq\rfloor+1)
      (2\lfloor\sqrt6Rq\delta\rfloor+1)
 =O_R(\delta q^2+q).}                               \tag{2.6}
\]

The lower annular bound in (2.4) is needed when speaking about a direction
at scale \(q\), but not for the upper count. The constant is uniform for
\(0<\delta\le\delta_0\). In particular,

\[
 \boxed{\#\mathcal S_{q,C/q}=O_{R,C}(q).}            \tag{2.7}
\]

Because \(\pi_H|_{\mathcal P}\) is injective, (2.6)--(2.7) also bound
\(\mathcal P\cap\mathbb Z^3\). In the irrational case of (1.10), that
exact-plane subset is substantially sparser: for a sufficiently narrow
sector about \(\mathbb Rr\), it contains no nonzero point at all.

### 2.1 The fixed-direction C159 ray

Inside \(\mathcal P\), horizontal projection parallel to \(r\) forces the
vector to lie on

\[
 \mathcal R_*=\mathbb R(r+\sigma N).                \tag{2.8}
\]

The C159 covector lies on this ray. Indeed, C152 gives

\[
 k_*=U_0+\beta N=3\gamma(r+\sigma N),
 \qquad \beta=3\gamma\sigma.                        \tag{2.9}
\]

The ray \(\mathcal R_*\) contains a nonzero integer vector if and only if
\(\sigma\in\mathbb Q\). In that case its integer points are multiples of
one primitive vector and only \(O_R(q)\) have size at most \(Rq\). Thus
even a favorable answer to the open arithmetic question cannot provide a
\(q^2\) aligned sheet.

## 3. Uniformity on the surviving ray and exact-plane sector

The covector equation is linear in its initial covector. If \(k(t)\) is a
solution, so is \(\lambda k(t)\). The Kelvin-amplitude generator

\[
 \mathcal L(t,k)v=-\mathsf A(t)v
 +2{k\over|k|^2}\bigl(k^T\mathsf A(t)v\bigr)        \tag{3.1}
\]

is homogeneous of degree zero in every nonzero real \(k\). Consequently
all nonzero points of \(\mathcal R_*\) have exactly the same physical
amplitude equation, monodromy on the common plane \(k_*^\perp\), and
expanding physical-amplitude line. C159 therefore applies at every scale on
the ray:

\[
 \rho(M_*)>1,\qquad \det M_*=1,
 \qquad \operatorname{tr}M_*>2.                    \tag{3.2}
\]

There is also a qualitative uniform statement for a shrinking exact-plane
sector. Every initial covector in \(\mathcal P\) is exactly periodic under
the C152 return. By (1.9), an angular \(C/q\) perturbation of its horizontal
direction about \(\mathbb Rr\) is an \(O_C(q^{-1})\) perturbation of its
projective full direction about \([k_*]\). Choose a smooth transverse frame
near \([k_*]\). The one-period monodromy depends smoothly on this
projective direction. Since (3.2) has a simple expanding multiplier
\(\rho_* >1\), there exist constants \(q_0,C_M>0\) such that, for
\(q\ge q_0\),

\[
 |\log\rho(k)-\log\rho_*|\le {C_M\over q},
 \qquad
 \operatorname{dist}(E_+(k),E_+(k_*))\le {C_M\over q}. \tag{3.3}
\]

The constants depend on the fixed aperture constant \(C\), the C152 orbit,
and the chosen local frame; this note does not compute them. Along the
individual expanding line \(E_+(k)\), for \(0\le\ell\le L\log q\),

\[
 {\rho(k)^\ell\over\rho_*^\ell}
 =\exp\!\left(O_{C,L}\left({\ell\over q}\right)\right)
 =1+O_{C,L}\left({\log q\over q}\right).           \tag{3.4}
\]

A source-adapted horizontal polarization chosen on the central ray has only
\(O_C(q^{-1})\) kinematic angular mismatch in this sector, provided it stays
uniformly away from its rank-defect line. Equations (3.3)--(3.4) concern the
principal periodic cocycle only. They do not promote the lattice candidates
to localized finite-frequency solutions.

## 4. C154 removes the off-plane normal lifts

Could one recover \(q^2\) modes by taking \(O(q)\) integer normal-charge
lifts over each of the \(O(q)\) projected points in a \(C/q\)-sector? Such
lifts generally do not lie in \(\mathcal P\).

For a projected point \(h\in\Lambda_H\), let
\(k_{\mathcal P}=\mathcal I_\sigma(h)\) be its unique **real** periodic-plane
lift. It need not be an integer vector. All integer lifts over \(h\), if one
is chosen as an origin, form \(k_m=k_{m_0}+mN\). Put
\(\Delta_m=k_m-k_{\mathcal P}\). C154 gives

\[
 K^\ell=I-\ell g_0\otimes u,
 \qquad
 K^\ell k_m=k_{\mathcal P}+\Delta_m
             -\ell g_0(u\cdot\Delta_m).             \tag{4.1}
\]

Here \(g_0=\gamma d\) is horizontal and perpendicular to \(r\). Suppose
that at return zero and at some return
\(c_1\log q\le\ell_q\le c_2\log q\), the horizontal projections remain in
a tube of width \(W_q\le C_0\delta q\) about \(\mathbb Rr\). Taking signed
components in the \(d\)-direction in (4.1) gives the necessary tolerance

\[
 \boxed{|u\cdot\Delta_m|
 \le {2W_q\over\ell_q|g_0|}
 =O\!\left({\delta q\over\log q}\right).}           \tag{4.2}
\]

For the aligned aperture \(\delta=C/q\), this becomes

\[
 \boxed{|u\cdot\Delta_m|=O((\log q)^{-1}).}          \tag{4.3}
\]

Successive integer lifts have fixed nonzero spacing

\[
 u\cdot\Delta_{m+1}-u\cdot\Delta_m
 =u\cdot N=-3b_0\ne0.                               \tag{4.4}
\]

Thus the number of compatible lifts over one projected point is at most

\[
 O\!\left(1+{\delta q\over\log q}\right).           \tag{4.5}
\]

When \(\delta=C/q\), the interval in (4.3) is eventually shorter than the
spacing in (4.4), so there is at most one compatible integer lift over each
projected point. Combining with (2.7),

\[
 \boxed{\#\{\text{aligned, C154-compatible integer lifts}\}=O(q).} \tag{4.6}
\]

This argument compares exact discrete returns of the principal covector
cocycle. It does not control the between-return path, a neighborhood of the
base orbit, or a finite-frequency packet. Those remain separate analytic
obligations.

## 5. Surviving obligation

C170 closes only the periodic-plane aligned-sheet repair. The unresolved
arithmetic condition for the surviving ray is

\[
 -{\sqrt{42}\over5}{T'(0)\over T(0)}\in\mathbb Q,   \tag{5.1}
\]

but deciding (5.1) cannot restore the missing factor \(q\) in cardinality.
The one-cell theorem must instead use a genuinely two-directional source
sheet and overcome its polarization/synchronization problem, or retain both
source polarizations in a noncommuting or unfolded transfer. Any such
proposal still needs finite-frequency localization, relative leakage,
retained wake, and the full Navier--Stokes stage estimates absent here.
