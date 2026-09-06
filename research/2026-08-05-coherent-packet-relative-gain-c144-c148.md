# C144--C148: coherent \(A_2\) focus, endpoint shell no-go, and the relative-gain gate

**Date:** 2026-08-05
**Status:** exact kinematic packet, exact shell/cardinality obstructions,
exact scalar reservoir arithmetic, and exact axial-layer linearization;
localized broad-band gain and endpoint conversion remain open
**Checker:**
[checks/coherent_packet_relative_gain_c144_c148.py](../checks/coherent_packet_relative_gain_c144_c148.py)

## 0. Claim boundary

This note stays on the C114--C143 \(A_2\) stage. It tests the only two
surviving replacements for the failed generic affine focus:

1. phase-coherent concentration using many existing-lattice modes; and
2. a prelocalized microseed amplified by the C121 decaying pump.

The results narrow the stage substantially.

- A packet with \(q^3\) conjugate pairs, hence \(2q^3\) nonzero Fourier
  wavevectors, and the required \(q^{3/2}\) point gain exists
  kinematically on the existing lattice.
- It repeats in every parent period cell and has order-one relative radial
  bandwidth. It is not a localized next pump.
- No exact or relatively narrow C121 shell can supply the required point
  gain at the scheduled child-volume energy.
- Collapsing a three-dimensional coherent packet into six next-pump modes
  by one quadratic layer requires a comparably large family of gate and
  companion frequencies before polarization cancellation.
- A microseed/reservoir ledger and an unfolded \(b^{-1}\)-gate ledger both
  meet all scalar exponents, but C125 requires errors relative to an
  \(n^{-28}\) growing seed, not merely absolutely summable errors.
- The exact C121 linearization on the relevant broad axial layer is
  non-Toeplitz. C122's six-mode eigenvalue does not imply a coherent
  \(q\)-wide growing packet.

The remaining target is named the **localized broad-band relative-gain
theorem (LBRG)**. It is the active-focus/gain component of the full BAFL
stage, not a separate cascade architecture.

## 1. C144: an exact coherent packet on the existing lattice

Use

\[
 N=(1,1,1),\qquad k_1=(1,-1,0),\qquad k_2=(0,1,-1).
 \tag{1.1}
\]

For integers \(q,K\ge1\) and \(0\le a,b,c<q\), define

\[
 v_{abc}=(4q+a)N+bk_1+ck_2,\qquad
 k_{abc}=Kv_{abc}.
 \tag{1.2}
\]

The determinant of the basis \((N,k_1,k_2)\) is \(3\), so the \(q^3\)
positive wavevectors are distinct. They all have positive \(N\)-component;
adding their reality partners gives \(2q^3\) distinct nonzero modes.

Orthogonality of \(N\) to the root plane gives

\[
 |v_{abc}|^2
 =3(4q+a)^2+2(b^2+c^2-bc),
 \tag{1.3}
\]

and therefore

\[
 48q^2\le |v_{abc}|^2
 \le77q^2-34q+5\le77q^2.
 \tag{1.4}
\]

Let \(E=k_1\) and

\[
 p_k=P_kE=E-{(E\cdot k)k\over|k|^2}.
 \tag{1.5}
\]

Then \(k\cdot p_k=0\) and

\[
 |p_{k_{abc}}|^2
 =2-{(2b-c)^2\over|v_{abc}|^2}
 \ge {23\over12}.
 \tag{1.6}
\]

For completeness, (1.6) is uniform in \(q\). After multiplying by
\(12|v_{abc}|^2\), it reduces to
\[
 46b^2-46bc+10c^2\le3(4q+a)^2.
\]
The left side is convex in \(c\), and its two endpoint restrictions are
convex in \(b\); its maximum on \(0\le b,c\le q-1\) is
\(46(q-1)^2<48q^2\), while the right side is at least \(48q^2\).

On \(\mathbb T^3=(\mathbb R/2\pi\mathbb Z)^3\) with normalized Haar
measure, put

\[
W(x)=\sum_{a,b,c}
 \left(p_{k_{abc}}e^{ik_{abc}\cdot x}
       +p_{k_{abc}}e^{-ik_{abc}\cdot x}\right)
 =2\sum_{a,b,c}p_{k_{abc}}\cos(k_{abc}\cdot x),
 \quad
 S=\sum_{a,b,c}|p_{k_{abc}}|^2,
 \quad
 U={W\over\sqrt{2S}}.
 \tag{1.7}
\]

Fourier orthogonality gives \(\|U\|_2=1\). At every coherence point
\(x_0=2\pi m/K\),

\[
 {E\over|E|}\cdot U(x_0)=\sqrt S
 \ge\sqrt{23\over12}\,q^{3/2}.
 \tag{1.8}
\]

Moreover, on

\[
 B\left(x_0,{1\over2\sqrt{77}Kq}\right),
 \tag{1.9}
\]

all phases have absolute value at most \(1/2\), so
\(\cos(k\cdot(x-x_0))\ge7/8\). Hence

\[
 \int_B|U|^2\,d\mu
 \ge {1127\over36864\pi^2\,77^{3/2}}\,K^{-3}.
 \tag{1.10}
\]

More precisely, let
\[
 Q_m=x_0+[-\pi/K,\pi/K)^3,\qquad x_0={2\pi m\over K}.
 \tag{1.11}
\]
Every centered coordinate parent cube \(Q_m\) has energy exactly
\(K^{-3}\), because \(U\) is \(2\pi/K\)-periodic in each coordinate, and
the ball (1.9) lies in \(Q_m\). Thus the fixed constant in (1.10) is a
lower bound for the fraction of each centered parent-cell energy lying in
a child-scale ball.

This is a kinematic per-cell concentration lemma, not a stage. The peak
repeats in every parent and primitive period cell. Isolating one physical
child requires a parent-scale envelope/cutoff, with the Leray, pressure,
and non-shell costs excluded from (1.7). The packet also spans the
order-one band (1.4), so it is not a Beltrami shell.

## 2. C145: exact/narrow-shell endpoint no-go

The schedule asks a single localized parent of scale \(K^{-1}\), amplitude
\(a\), and energy \(E\asymp a^2K^{-3}\) to produce a child with

\[
 b=n^{-2},\qquad q=n^8,\qquad h=q^{3/2}=n^{12},
 \tag{2.1}
\]

so the child has energy \(b^2E\) before the next depletion and point
amplitude \(bha\).

Let a periodic vector field with \(M\) Fourier wavevectors have precisely
that \(L^2\) scale:

\[
 \|u\|_2^2=b^2a^2K^{-3}.
 \tag{2.2}
\]

Cauchy--Schwarz at a point gives

\[
 |u(x_0)|^2\le M\|u\|_2^2.
 \tag{2.3}
\]

Requiring \(|u(x_0)|\ge bha\) forces

\[
 \boxed{M\ge h^2K^3=q^3K^3=(qK)^3.}
 \tag{2.4}
\]

By contrast, an exact integer shell of radius \(R\) contains at most

\[
 2(2\lfloor R\rfloor+1)^2=O(R^2)
 \tag{2.5}
\]

wavevectors: after fixing the first two coordinates, there are at most two
choices for the third. At \(R\asymp qK\), (2.5) is short of (2.4) by one
full radial degree of freedom.

The same conclusion holds for a relatively narrow shell. Associate a unit
cube to each lattice point in a shell of radius \(R\) and thickness
\(\Delta\). The cubes are disjoint and lie in a Euclidean shell enlarged
by \(\sqrt3/2\), so their number is

\[
 O\!\left(R^2(\Delta+1)+R(\Delta+1)^2+(\Delta+1)^3\right).
 \tag{2.6}
\]

If \(\Delta=o(R)\), (2.6) is \(o(R^3)\) and cannot meet (2.4). Thus the
required coherent focus has order-one relative radial bandwidth.

An exact six-root C121 pump is an even more immediate mismatch: fixed
point amplitude \(a_{j+1}\) gives \(L^2\) size comparable to
\(a_{j+1}\), not to
\(a_{j+1}(qK)^{-3/2}\). C126 independently proves that a nonzero compact
pump cannot be exactly monochromatic. Therefore the honest endpoint can
only be a localized C121-like core plus a non-shell collar/wake; it cannot
be the exact global C121 orbit.

Heat evolution cannot collapse the broad support. If the broad packet is
made exactly homochiral on one shell, its Leray-projected Euler
nonlinearity vanishes, so it also cannot nonlinearly redistribute itself
into the six target roots. A converter must use heterochiral, off-shell,
or localized gates.

## 3. C146: the exact support cost and the unfolded scalar opening

Let \(S\) be a source set of \(M\) wavevectors and let \(T\) be the six
target pump wavevectors. A one-step quadratic assignment sends
\(s\in S\) to \(t\in T\) using

\[
 g=t-s.
 \tag{3.1}
\]

One fixed oriented gate \(g\) serves at most six sources, one for each
\(t\). Consequently at least \(M/6\) distinct oriented gate frequencies
are required. If the real pairs \(\{g,-g\}\) are counted as single
unoriented gates, this statement alone gives only \(M/12\) pairs.
Reality also supplies \(-g\), and the same source then has the
support-allowed companion

\[
 s-g=2s-t.
 \tag{3.2}
\]

For a fixed companion wavevector and target \(t\), (3.2) determines \(s\);
therefore the collision multiplicity is at most six. Before polarization
cancellation, at least \(M/6\) distinct oriented companion wavevectors are
support-allowed as well (or \(M/12\) unoriented reality pairs).

This is a cardinality statement, not a coefficient lower bound. The C116
dual-helicity gate can cancel its first companion coefficient, although
C141 then gives a nonzero cubic return. A simultaneous \(M\)-channel
converter also creates cross-channel sums not counted by (3.2).

There is one exact scalar opening for temporal unfolding. Put

\[
 J=b^{-1}=n^2,\qquad\theta=b=n^{-2}.
 \tag{3.3}
\]

If the \(J\) sequential microgate target contributions are aligned and
each has size \(\asymp b\theta\), while the retained wake and active return
from each gate are bounded by \(O(b\theta^2)\) and \(O(b\theta^3)\),
respectively, then triangle summation gives the conditional power ledger

\[
 \begin{aligned}
 \text{target}&\asymp Jb\theta=b,\\
 \text{wake}&\lesssim Jb\theta^2=b^2=n^{-4},\\
 \text{active return}&\lesssim Jb\theta^3=b^3=n^{-6},\\
 J\theta&=1.
 \end{aligned}
 \tag{3.4}
\]

Thus temporal unfolding is compatible with the C140/BAFL powers under
the displayed per-gate hypotheses. It is not a construction. For the
\(M\gtrsim(qK)^3\) coherent degrees in (2.4), each
temporal batch must still handle many channels or the stage would require
far more than \(J\) batches. The missing estimate is uniform control of
all cross-channel interactions, localization collars, and pressure return,
independent of the number of coherent channels.

## 4. C147: exact microseed/reservoir arithmetic

The scalar alternative is to place an extremely weak seed in the child
volume before gain. Set

\[
 s=n^{-16},\qquad q^{-3}=n^{-24},\qquad H=n^{26}.
 \tag{4.1}
\]

Then the incoming pointwise and \(L^2\) fractions, final pointwise and
\(L^2\) fractions, and final energy fraction are exactly

\[
 \begin{aligned}
 \varepsilon&=s q^{-3/2}=n^{-28},\\
 Hs&=n^{10}=g,\\
 H\varepsilon&=n^{-2}=b,\\
 (H\varepsilon)^2&=n^{-4}.
 \end{aligned}
 \tag{4.2}
\]

The raw child-frequency writer arithmetic is also favorable absolutely:

\[
 sq\log H=26n^{-8}\log n=o(n^{-6}).
 \tag{4.3}
\]

If it has the same child support, its \(L^2\) size is
\(26n^{-20}\log n\). Relative to the growing seed, however,

\[
 {26n^{-20}\log n\over n^{-28}}
 =26n^8\log n.
 \tag{4.4}
\]

Thus absolute summability is irrelevant to C125 unless the writer term is
part of the intended growing coordinate, is exactly dark, or satisfies a
new relative cancellation.

The ideal retained-line gain itself is compatible with the factorial
Reynolds number. In inertial time

\[
 \vartheta=aKt,\qquad \mu={\nu K\over a},
 \tag{4.5}
\]

let \(\lambda>0\) be the dimensionless physical growing coefficient for a
leaf \(qN+r_i\). Then

\[
 {dz\over d\vartheta}
 =\left[\lambda e^{-2\mu\vartheta}
        -\mu(3q^2+2)\right]z
 \tag{4.6}
\]

and

\[
 G_q(\vartheta)
 ={\lambda\over2\mu}(1-e^{-2\mu\vartheta})
  -\mu(3q^2+2)\vartheta.
 \tag{4.7}
\]

For

\[
 R={\lambda\over\mu(3q^2+2)}>1,
 \tag{4.8}
\]

the exact maximum is

\[
 \vartheta_*={\log R\over2\mu},\qquad
 G_q(\vartheta_*)={3q^2+2\over2}(R-1-\log R).
 \tag{4.9}
\]

On C127,

\[
 \mu_j=\nu(j!)^{-2},\qquad q_j=n^8,\qquad
 \mu_jq_j^2\log n\longrightarrow0.
 \tag{4.10}
\]

Conditionally on a physical retained eigenline with
\(\lambda\ge\lambda_0>0\), the first time at which
\(G_q=26\log n\) is

\[
 \vartheta_H={26\over\lambda}\log n\,[1+o(1)].
 \tag{4.11}
\]

Pump decay and leaf heat are then summable. Equations (4.6)--(4.11) do
not supply a localized broad-band eigenline.

The final packet is also not perturbatively small in a critical norm.
Its \(L^2\) fraction is \(n^{-2}\), but the child frequency is \(q=n^8\),
so the scale-invariant \(H^{1/2}\) factor is

\[
 q^{1/2}n^{-2}=n^2.
 \tag{4.12}
\]

C124's finite normal-form comparison therefore cannot control the
localized packet nonlinearity by smallness alone.

## 5. C148: exact broad axial-layer operator and LBRG

The gap between (4.7) and a coherent packet can be written exactly. Set
the pump frequency \(K=1\); restoring dimensions multiplies the inviscid
operator by \(aK\) and gives \(\mu=\nu K/a\). Let

\[
 r_1=k_1,\qquad r_2=k_2,\qquad r_3=-r_1-r_2,
 \tag{5.1}
\]

\[
 d_1=(1,0),\qquad d_2=(0,1),\qquad d_3=(-1,-1),
 \tag{5.2}
\]

and

\[
 t_i=N\times r_i,\qquad h_i=t_i+i\sqrt2N,\qquad
 u_{+i}=c_ih_i,\qquad u_{-i}=\overline{c_ih_i}.
 \tag{5.3}
\]

For an axial layer

\[
 \kappa_{bc}=mN+br_1+cr_2,\qquad
 D_{bc}=|\kappa_{bc}|^2
 =3m^2+2(b^2+c^2-bc),
 \tag{5.4}
\]

write \(v_{bc}\perp\kappa_{bc}\). The exact inertial-time
linearization about the heat-decaying C121 pump is

\[
 \boxed{
 \begin{aligned}
 (A_m(\vartheta)v)_{bc}
 ={}&-\mu D_{bc}v_{bc}\\
 &-ie^{-2\mu\vartheta}
 \sum_{\substack{i=1,2,3\\ \epsilon=\pm1}}
 P_{\kappa_{bc}}\bigl[
 (u_{\epsilon i}\cdot\kappa_{bc})
       v_{(b,c)-\epsilon d_i}\\
 &\hspace{43mm}
 +\epsilon(v_{(b,c)-\epsilon d_i}\cdot r_i)u_{\epsilon i}
 \bigr].
 \end{aligned}}
 \tag{5.5}
\]

The scalar weights are

\[
 \begin{aligned}
 u_{+i}\cdot\kappa_{bc}
 &=3c_i(\zeta_i+i\sqrt2m),\\
 u_{-i}\cdot\kappa_{bc}
 &=3\overline{c_i}(\zeta_i-i\sqrt2m),
 \end{aligned}
 \tag{5.6}
\]

where

\[
 \zeta_1=c,\qquad\zeta_2=-b,\qquad\zeta_3=b-c.
 \tag{5.7}
\]

Thus the exact interior operator is not translation invariant: both the
weights (5.6) and the transverse fibers/projections vary with \((b,c)\).
A nonzero constant-polarization Bloch ansatz is impossible. Orthogonality
at \((b,c)=(0,0),(1,0),(0,1)\) would force one constant polarization to be
orthogonal to \(N,r_1,r_2\), hence to all of \(\mathbb R^3\).

Freezing one fiber at a fixed \(\kappa\) gives only a diagnostic symbol.
For \(\theta\in\mathbb T^2\), define

\[
 \begin{aligned}
 U_\theta&=\sum_i\left(
 u_{+i}e^{-i\theta\cdot d_i}
 +u_{-i}e^{i\theta\cdot d_i}\right),\\
 C_\theta&=\sum_i\left(
 u_{+i}r_i^Te^{-i\theta\cdot d_i}
 -u_{-i}r_i^Te^{i\theta\cdot d_i}\right).
 \end{aligned}
 \tag{5.8}
\]

Here \(U_\theta\) is real and \(C_\theta\) is purely imaginary. If \(R_\kappa\)
is a real orthonormal frame for \(\kappa^\perp\), put

\[
 {\cal G}_\theta=-iR_\kappa^TC_\theta R_\kappa.
 \tag{5.9}
\]

The two frozen eigenvalues are

\[
 \begin{aligned}
 \lambda_\pm={}&-\mu|\kappa|^2
 -ie^{-2\mu\vartheta}\kappa\cdot U_\theta\\
 &+{e^{-2\mu\vartheta}\over2}\left[
 \operatorname{tr}{\cal G}_\theta
 \pm\sqrt{(\operatorname{tr}{\cal G}_\theta)^2
          -4\det{\cal G}_\theta}\right].
 \end{aligned}
 \tag{5.10}
\]

The \(O(m)\) principal term is purely imaginary transport. Possible Euler
growth is the \(O(1)\) stretching matrix \({\cal G}_\theta\), plus the
fiber-connection effects omitted by freezing; those connection terms are
also \(O(1)\), the same order as prospective growth. C121 allows arbitrary
pump coefficients \(c_i\); it therefore supplies no pump-independent
positive band. For example, with only one nonzero root coefficient the
frozen stretching block has rank at most one and its possible nonzero real
eigenvalue is a trigonometric function which vanishes for some \(\theta\).

There is also an exact obstruction for the natural equal-real three-root
pump. If \(c_1=c_2=c_3\in\mathbb R\), then

\[
 \sum_i r_i=0,\qquad \sum_i t_i=0
 \quad\Longrightarrow\quad
 U_0=0,\qquad C_0=0.
 \tag{5.11}
\]

At \(\theta=0\), its frozen symbol is therefore only
\(-\mu|\kappa|^2I\). It has no uniformly positive frozen Bloch band. This
does not rule out a non-frozen positive cocycle for an optimized
three-root pump.

The C144 packet has \(m\in[4q,5q)\) and tangential width \(q\). Across that
box, \(\zeta_i/m\) varies by order one; the fixed C120/C122 finite matrix
does not approximate (5.5) uniformly. A successful packet may use an
eikonal/Lagrangian rephasing to absorb the \(O(m)\) transport, but it still
needs a uniform positive \(O(1)\) amplitude cocycle, phase coherence,
boundary control, and the relative C125 estimate.

We call that missing statement **LBRG**. To keep its relative estimate
non-vacuous, let \(A^0_{H,j}(s)\) denote the intended retained broad-band
linear generator on \(0\le s\le T_j\). The theorem must produce a right
trajectory \(w_{+,j}(s)\), a left functional \(\Phi_j(s)\), and a real
exponent \(G_j(s)\) satisfying
\[
 \partial_sw_{+,j}=A^0_{H,j}w_{+,j},\qquad
 \|w_{+,j}(0)\|=1,\qquad G_j(0)=0,
 \tag{5.12}
\]
\[
 \|\Phi_j(s)\|=1,\qquad
 \Phi_j(s)w_{+,j}(s)=e^{G_j(s)},\qquad
 \partial_s\Phi_j+\Phi_jA^0_{H,j}=G_j'\Phi_j.
 \tag{5.13}
\]
The last identity says that \(e^{-G_j}\Phi_j\) is an adjoint trajectory.
Together with the unit-dual normalization and the right trajectory, it
ties \(G_j\) to physical norm growth and prevents rescaling an arbitrary
test functional to make the defect small.

> Construct, for the actual localized C121-like pump and every sufficiently
> large stage, a coherent growing packet/bundle covering the required
> three-dimensional Fourier degrees, together with data satisfying
> (5.12)--(5.13) and \(G_j(T_j)\ge26\log n\), such that the full C125
> weighted defect satisfies
> \[
> {1\over\varepsilon_j}
> \int_0^{T_j} e^{-G_j(s)}
> \left|\Phi_j(s)\!\left(
> P_HA v_\perp+E_Hv_H+P_H{\cal N}(v,v)
> \right)\right|ds
> \le\theta<1,
> \qquad \varepsilon_j=n^{-28},
> \tag{LBRG}
> \]
> uniformly in \(j\), while the packet remains coherent enough to enter the
> unfolded localized next-pump conversion and obeys the BAFL wake bounds.

This packages the exact relative obstruction rather than hiding it in
absolute summability. LBRG is open. No full one-cell map, singularity, or
Millennium conclusion is claimed.
