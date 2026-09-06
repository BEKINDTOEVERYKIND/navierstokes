# C195: quantitative off-ray \(A_2\) finite-horizon dominated cones

**Date:** 2026-08-30

**Status:** outward computer-assisted principal phase-space dominated-cone
field for a finite \(R_{\rm filt}\)-block horizon; no invariant or
canonical stable/unstable bundle, closed periodic phase, cone-field
derivatives, finite-frequency packet, viscosity, nonlinear return, UVSR,
or singularity theorem

**Checker:**
[checks/a2_offray_hyperbolic_bundle.py](../checks/a2_offray_hyperbolic_bundle.py)

## 0. Verdict

The C159/C192 polarization cone gain is not confined to one exact ray over
the finite filter horizon.  An explicit fat neighborhood of the returning
\(m\ne0\) orbit has strictly mapped forward and reflected-inverse cones,
with coefficient gain larger than \(3000\) per block in both directions.
For any positive forward-cone line \(\ell_+\) and any sign-reflected
backward-cone line \(\ell_-\), their physical angle and the corresponding
two-line oblique projectors obey

\[
 \sin\angle(\ell_+,\ell_-)\ge\frac{520}{569}>\frac9{10},
 \qquad
 \|P_{\ell_+\parallel\ell_-}\|,
 \|P_{\ell_-\parallel\ell_+}\|\le\frac{569}{520}.     \tag{0.1}
\]

Choose a common Euclidean lift of the initial torus neighborhood and put

\[
 r=\max\{\|X_0-X_{0,*}\|_2,\ \|k_0-k_{0,*}\|_2\}.       \tag{0.2}
\]

For the C193 filter horizon, this initial base/covector displacement is
certified whenever

\[
 r(1+TR_{\rm filt})^3\le\frac1{87{,}000{,}000{,}000{,}000}. \tag{0.3}
\]

Thus on \(q=n^8\), the radii used in the present packet ledgers are
certified at the explicit thresholds

\[
\begin{array}{c|c}
r&q=n^8\text{ threshold}\\ \hline
q^{-1/4}&n\ge10^{10},\\
q^{-1/12}&n\ge10^{39},\\
q^{-1/3}&n\ge10^8.
\end{array}                                               \tag{0.4}
\]

The large \(n\) values are an honest consequence of the coarse explicit
constants; they are not suppressed into asymptotic notation.  These facts
give a finite-horizon dominated cone field, not an invariant splitting or
canonical lines \(E^u,E^s\). The proposed anisotropic two-width box is
governed by its larger radius \(q^{-1/12}\), so the combined box is certified
at \(n\ge10^{39}\).

## 1. Correct off-level Kelvin block

Write \(n=N/\sqrt3\), split the covector as \(k=p+mn\), put

\[
 D=|p|^2,\qquad Q=D+m^2,\qquad t=p\times n,
\]

and write the velocity-gradient block as

\[
 A=\begin{pmatrix}S&0\\ l^T&0\end{pmatrix},\qquad
 \operatorname {tr}S=0,
\]

in the horizontal--axial splitting. The covector equation and the physical
Kelvin equation are

\[
 \dot p=-S^Tp-ml,\quad \dot m=0,\quad
 \dot a=-Aa+2k\,{k\cdot Aa\over Q}.                    \tag{1.1}
\]

In the moving orthogonal frame \(E_1=P_kn=n-mk/Q\),
\(E_2=k\times n=t\), one has

\[
 |E_1|^2={D\over Q},\quad |E_2|^2=D,\quad
 \dot E_1=-{m\dot k\over Q}+{m\dot Q\over Q^2}k,
 \quad \dot E_2=\dot k\times n.                       \tag{1.2}
\]

Writing \(a=z_1E_1+z_2E_2\), orthogonal projection gives, for
\(i,j\in\{1,2\}\),

\[
 B_{ij}={E_i\cdot[-AE_j+2k(k\cdot AE_j)/Q-\dot E_j]
                 \over |E_i|^2}.                       \tag{1.3}
\]

Substitution of (1.1)--(1.2), using only \(\operatorname{tr}S=0\), gives
all four coefficient entries

\[
\boxed{\begin{aligned}
B_{11}&={m\,l\cdot p\over D},\\
B_{12}&={2m\,p\cdot St-(l\cdot t)(D-m^2)\over D},\\
B_{21}&={m(t\cdot Sp-p\cdot St)-m^2l\cdot t\over QD},\\
B_{22}&={2p\cdot Sp+m\,l\cdot p\over D}.
\end{aligned}}                                             \tag{1.4}
\]

Away from \(f=0\), \(S\) remains trace free but is no longer symmetric.
Thus \(B_{21}\) contains the skew term

\[
                    \frac{m(t\cdot Sp-p\cdot St)}{QD}.   \tag{1.5}
\]

The checker proves (1.4) coefficient by coefficient in the exact symbolic
variables \(p_1,p_2,m,S_{11},S_{12},S_{21},l_1,l_2\), after independently
differentiating \(E_1,E_2\) and clearing every denominator. Rational
nonsymmetric specializations are retained as additional sentinels. For
the \(A_2\) pump, \(l\cdot t=-\sqrt2c_h\), where
\(c_h=p\cdot(N\times\nabla f)\), and the root-frame algebra gives

\[
                  t\cdot Sp-p\cdot St=2\sqrt3\,fD.       \tag{1.6}
\]

Consequently the extra \(B_{21}\) contribution is

\[
                         \frac{2\sqrt3\,mf}{Q}.           \tag{1.7}
\]

The certificate retains the invariant level

\[
                              |f|\le10^{-11}              \tag{1.8}
\]

instead of silently reusing the zero-level symmetric formula.

The remaining coordinates are also allowed to vary independently. Around
the already certified C159 reference tube, the interval calculation
fattens every cell by

\[
\begin{array}{c|c}
\text{coordinate}&\text{additional radius}\\ \hline
\text{phase}&10^{-7},\\
\gamma&10^{-5},\\
\beta=m/\sqrt3&10^{-7},\\
c_h=p\cdot(N\times\nabla f)&10^{-5},\\
f&10^{-11}.
\end{array}                                               \tag{1.9}
\]

It also retains the reference interpolation radii \(2\cdot10^{-6}\) in
phase and \(8\cdot10^{-4}\) in \(\gamma\).

## 2. Outward cone certificate

The checker regenerates the complete C159 reference path and subdivides it
into \(2048\) cells. On the coefficient cone

\[
                    {\cal K}=\{(z_1,z_2):z_1>0,\ 
                              0.137\le z_2/z_1\le0.2\},  \tag{2.1}
\]

directed interval Taylor propagation gives

\[
\begin{aligned}
 z_{1,\rm fwd}&>3009.757266,\\
 z_{1,\rm inv}&>3005.135765.                             \tag{2.2}
\end{aligned}
\]

The complete slope images satisfy

\[
\begin{aligned}
 {\cal R}_{\rm fwd}([0.137,0.2])
   &\subset[0.1405737,0.1898012],\\
 {\cal R}_{\rm inv}([0.137,0.2])
   &\subset[0.1402257,0.1903237].                        \tag{2.3}
\end{aligned}
\]

Both lie strictly inside (2.1). Since \(3000>e^8\), every admissible
forward block and every reflected inverse block has coefficient-cone gain
above \(e^8\). Induction through the declared horizon gives
\(z_1(jT)>3000^jz_1(0)\) for every forward-cone vector and
\(0\le j\le R_{\rm filt}\), with the analogous statement for the reflected
inverse evolution. No eigenvalue, invariant line, or canonical splitting is
assigned to an individual nonreturning block.

At each comparison-section time \(3<|P_nk|<4\), so the coordinate
determinant of a Kelvin block is below \(16/9<2\). Equation (2.2)
therefore makes that block's projective slope-map contraction smaller than

\[
 \operatorname{Lip}({\cal R}_{\rm block})
 <\frac2{3000^2}=\frac1{4{,}500{,}000}.                \tag{2.4}
\]

Every positive forward-cone coefficient slope and every sign-reflected
backward-cone slope convert to physical orthonormal magnitudes in
\([13/20,6/5]\). The elementary two-line angle formula then gives (0.1)
for any such pair; it does not construct preferred lines.

Between comparison-section times, C193's \(\|DU\|_{\rm op}\le6\) makes the physical
Kelvin generator norm at most \(18\). Since \(T<76/25\) and \(e<3\), each
one-block propagator and inverse has the explicit fixed bound

\[
                 \|\Phi(t,s)\|_{\rm op}<3^{55}
                 =174449211009120179071170507.           \tag{2.5}
\]

This controls continuous within-block overshoot by a fixed factor; it
does not provide derivatives of the cone field.

## 3. Initial-radius transport and thresholds

C194 gives the action-angle bounds

\[
 J_1(t)\le216(1+t),\qquad J_2(t)\le286{,}992(1+t)^2.
\]

All norms below are physical Euclidean norms on the common lift used in
(0.2). First,

\[
 |\delta X(t)|\le J_1(t)r\le216Lr,\qquad L=1+t.         \tag{3.1}
\]

For the Eulerian transported covector field C194 gives
\(\|D_xk(t,\cdot)\|\le |k_0|J_2(t)\). Comparing rays launched from distinct
initial base points first pays (3.1), while changing the initial covector
pays \(J_1r\). Since \(|k_{0,*}|+r<7\),

\[
 |\delta k(t)|
 \le216Lr+216\cdot286{,}992\cdot7\,L^3r
 <434{,}000{,}000L^3r.                                 \tag{3.2}
\]

The remaining radius rows are consequences, not extra premises:

* the roots have length \(\sqrt2\), so every sine/cosine phase coordinate
  differs by less than \(306Lr\);
* \(m=k\cdot n\) is invariant and \(\beta=m/\sqrt3\), hence
  \(|\delta\beta|\le r\);
* \(f\) is invariant and \(\sup|\nabla f|<4\), hence
  \(|\delta f|\le4r\);
* \(c_h=p\cdot(N\times\nabla f)=k\cdot U+\sqrt6\,mf\) is invariant because
  each term on the right is invariant. At the initial section,
  \(|N\times\nabla f|<7\), its derivative is below \(10\), and
  \(|p_{0,*}|<4\), so \(|\delta c_h|<47r<200r\);
* \(\gamma=(p\cdot g)/|g|^2\). On the annulus,
  \(|g|<4,\ \|Dg\|<6,\ |g|^2\ge3/2\), while the certified reference has
  \(|p_*|<12,\ |\gamma_*|<1\). Direct subtraction of the numerator and
  denominator gives
  \[
  |\delta\gamma|\le2|\delta k|+100|\delta X|.           \tag{3.3}
  \]

Equations (3.1)--(3.3) show that (0.3) places phase, \(\gamma,\beta,c_h,f\)
inside every additional box in (1.9): the respective upper bounds are
\(10^{-7},10^{-5},10^{-7},10^{-5},10^{-11}\). On the C193 clock,

\[
 R_{\rm filt}=\left\lceil\frac38\log n\right\rceil+1,
 \qquad 1+TR_{\rm filt}<14n^{1/14}.                     \tag{3.4}
\]

Raising the three resulting rational inequalities to their common
denominators proves the thresholds (0.4) with integer arithmetic only.

## 4. Claim boundary

C195 closes only the requested **finite-horizon \(C^0\) off-ray
dominated-cone robustness** problem. It does not close:

* an invariant or canonical stable/unstable splitting;
* a globally exact periodic eikonal phase or integer carrier;
* \(C^1\) dependence of the cone field or selected cone lines;
* a real, finite-frequency, periodic packet;
* a uniform multi-beam Fourier-integral parametrix or endpoint band;
* viscosity, nonlinear depletion/wake, UVSR, or singular-center tracking.

In particular, the discarded draft estimate for cone/projector first
derivatives is not part of this claim. C196 separately shows that even
successful shrinking-tube cone control does not by itself supply enough
endpoint Fourier phase space for raw \(q^{3/2}\) focus.
