# Aspect-uniform Piola transport of a compact straight packet

**Date:** 2026-08-03

**Status:** exact geometric commutator and norm lemma; conditional
frequency-\(p\) specialization. The uniform coalescing-edge scalar
pseudomode is not constructed here.

**Scope:** this note closes the purely geometric part of the ambient-Piola
obligation left open in C79. It transports the exact q-free
velocity--pressure lift from a straight periodic column to a thin embedded
torus, without comparing two global Hodge projectors. It does not turn the
resulting quasimode into an eigenmode or a Navier--Stokes singularity.

## 1. Geometry and the exact norm identity

Let \(D\Subset\mathbb R^2\) be a fixed bounded cross-section and put

\[
 Q_\varepsilon=D\times
 \mathbb R/(2\pi/\varepsilon)\mathbb Z,\qquad
 y=(y_1,y_2,z).
                                                               \tag{1.1}
\]

Assume \(|y_1|\le R_0\) on \(D\) and
\(\delta=\varepsilon R_0<1\). The normalized thin-torus map is

\[
 F_\varepsilon(y)=
 \left(
  (\varepsilon^{-1}+y_1)\cos(\varepsilon z),
  (\varepsilon^{-1}+y_1)\sin(\varepsilon z),
  y_2
 \right).                                                     \tag{1.2}
\]

In the moving orthonormal frame

\[
 E_1=e_R(\varepsilon z),\qquad E_2=e_Z,\qquad
 E_3=e_\phi(\varepsilon z),
\]

write

\[
 h=1+\varepsilon y_1.
\]

Then

\[
 DF_\varepsilon=\operatorname{diag}(1,1,h),\qquad
 J_\varepsilon=h,\qquad
 F_\varepsilon^*(dx^2)=dy_1^2+dy_2^2+h^2dz^2.                \tag{1.3}
\]

The contravariant Piola transform and scalar pressure transport are

\[
 {\cal P}_\varepsilon u
 =J_\varepsilon^{-1}DF_\varepsilon u
 =\left(\frac{u_1}{h},\frac{u_2}{h},u_3\right),
 \qquad
 \Pi=\pi\circ F_\varepsilon^{-1}.                             \tag{1.4}
\]

The physical Euclidean \(L^2\) norm has the exact pullback formula

\[
\begin{aligned}
 \|{\cal P}_\varepsilon u\|_{L^2(F_\varepsilon Q_\varepsilon)}^2
  &=
  \int_{Q_\varepsilon}
  \left(\frac{|u_1|^2+|u_2|^2}{h}
                  +h|u_3|^2\right)\,dy,\\
 (1-\delta)\|u\|_2^2
  &\le
  \|{\cal P}_\varepsilon u\|_2^2
  \le (1-\delta)^{-1}\|u\|_2^2.                              \tag{1.5}
\end{aligned}
\]

No constant in (1.5) depends on the aspect length
\(2\pi/\varepsilon\).

For a single axial Fourier mode whose cross-sectional amplitude is
\(u(\cdot,\cdot,0)\),

\[
 \|u\|_{L^2(Q_\varepsilon)}^2
 =\frac{2\pi}{\varepsilon}
   \|u(\cdot,\cdot,0)\|_{L^2(D)}^2.                           \tag{1.6}
\]

Thus an ambient unit packet carries the harmless amplitude factor
\((\varepsilon/2\pi)^{1/2}\). The same factor multiplies its residual, so
there is no aspect loss in the normalized residual ratio.

## 2. Exact divergence preservation and support

For a vector field \(v=(v_1,v_2,v_3)\) written in the moving frame,

\[
 \operatorname{div}_xv
 =\partial_1v_1+\partial_2v_2+h^{-1}\partial_zv_3
       +\frac{\varepsilon}{h}v_1.                             \tag{2.1}
\]

Substitution of (1.4) gives the exact identity

\[
 \operatorname{div}_x({\cal P}_\varepsilon u)
 =h^{-1}\operatorname{div}_yu.                               \tag{2.2}
\]

Suppose \(u\) and \(\pi\) are smooth and supported in
\(D_0\times\mathbb T_{2\pi/\varepsilon}\), where
\(D_0\Subset D\). The Piola field and transported pressure then vanish in
a full collar of the boundary of the tubular chart. Extension by zero
produces a smooth compactly supported ambient field on \(\mathbb R^3\);
(2.2) remains true distributionally. In particular, a straight
divergence-free packet becomes an exactly ambient divergence-free packet.

This is why scalar localization must precede q-free reconstruction and
bending. No new cutoff is inserted after the divergence equation has been
solved.

## 3. Exact first-order commutator

Let \(a\) be the straight base field and let \(A_\varepsilon\) be the
moving-frame pullback of the exact curved Gavrilov base. Define

\[
 b_\varepsilon=A_\varepsilon-{\cal P}_\varepsilon a.          \tag{3.1}
\]

For a reference vector field \(c\), write

\[
 D_c=c_1\partial_1+c_2\partial_2+c_3\partial_z.
\]

The physical covariant derivative in the moving frame is

\[
\begin{aligned}
 (\nabla^x_Av)_1
  &=A_1\partial_1v_1+A_2\partial_2v_1
    +\frac{A_3}{h}\partial_zv_1
    -\frac{\varepsilon A_3}{h}v_3,\\
 (\nabla^x_Av)_2
  &=A_1\partial_1v_2+A_2\partial_2v_2
    +\frac{A_3}{h}\partial_zv_2,\\
 (\nabla^x_Av)_3
  &=A_1\partial_1v_3+A_2\partial_2v_3
    +\frac{A_3}{h}\partial_zv_3
    +\frac{\varepsilon A_3}{h}v_1.                           \tag{3.2}
\end{aligned}
\]

Set

\[
 B_\varepsilon(c,u)
 =\nabla^x_{{\cal P}_\varepsilon c}
              ({\cal P}_\varepsilon u)
  -{\cal P}_\varepsilon(D_cu).                               \tag{3.3}
\]

A direct calculation gives

\[
\begin{aligned}
 (B_\varepsilon(c,u))_1
  &=-\frac{h-1}{h^2}D_cu_1
    -\frac{\varepsilon c_1u_1}{h^3}
    -\frac{\varepsilon c_3u_3}{h},\\
 (B_\varepsilon(c,u))_2
  &=-\frac{h-1}{h^2}D_cu_2
    -\frac{\varepsilon c_1u_2}{h^3},\\
 (B_\varepsilon(c,u))_3
  &=-\frac{h-1}{h}D_cu_3
    +\frac{\varepsilon c_3u_1}{h^2}.                         \tag{3.4}
\end{aligned}
\]

Transporting pressure as a scalar does not exactly intertwine gradients.
The defect is nevertheless explicit:

\[
\begin{aligned}
 C_\varepsilon(\pi)
  &=
 \nabla_x(\pi\circ F_\varepsilon^{-1})
 -{\cal P}_\varepsilon\nabla_y\pi\\
  &=\frac{h-1}{h}
     \left(\partial_1\pi,\partial_2\pi,-\partial_z\pi\right).
                                                                    \tag{3.5}
\end{aligned}
\]

Define the unprojected linearized Euler momentum expressions

\[
\begin{aligned}
 {\cal M}_0(a;u,\pi)
  &=\partial_tu+D_au+D_ua+\nabla_y\pi,\\
 {\cal M}_\varepsilon(A;v,\Pi)
  &=\partial_tv+\nabla^x_Av+\nabla^x_vA+\nabla_x\Pi.
                                                                    \tag{3.6}
\end{aligned}
\]

Because the Piola map is time independent, (3.3)--(3.5) imply the exact
identity

\[
\begin{aligned}
 {\cal M}_\varepsilon
 (A_\varepsilon;{\cal P}_\varepsilon u,\Pi)
 &={\cal P}_\varepsilon{\cal M}_0(a;u,\pi)
   +B_\varepsilon(a,u)+B_\varepsilon(u,a)
   +C_\varepsilon(\pi)\\
 &\quad
   +\nabla^x_{b_\varepsilon}({\cal P}_\varepsilon u)
   +\nabla^x_{{\cal P}_\varepsilon u}b_\varepsilon.           \tag{3.7}
\end{aligned}
\]

In particular, there is no curvature factor multiplying
\(\partial_tu\), hence no separate \(\varepsilon|\omega|\) commutator.
The frequency loss in (3.7) comes only from spatial first derivatives and
the pressure gradient.

## 4. Aspect-uniform physical \(L^2\) estimate

Assume

\[
 \varepsilon R_0\le\delta_0<1,\qquad
 \|a\|_{W^{1,\infty}(D)}\le A_0,\qquad
 \|A_\varepsilon-a\|_{W^{1,\infty}(D)}
       \le A_1\varepsilon.                                   \tag{4.1}
\]

The hollow Gavrilov comparison supplies the final bound in (4.1) on a
fixed active annulus. Since
\({\cal P}_\varepsilon a-a=O_{W^{1,\infty}}(\varepsilon)\),

\[
 \|b_\varepsilon\|_{W^{1,\infty}}
 \le C(R_0,\delta_0,A_0,A_1)\varepsilon.                     \tag{4.2}
\]

Equations (3.2)--(3.5), the exact weighted norm (1.5), and (4.2) yield

\[
\begin{aligned}
 &\|{\cal M}_\varepsilon
 (A_\varepsilon;{\cal P}_\varepsilon u,\Pi)\|_{L^2_x}\\
 &\qquad\le
 C_0\|{\cal M}_0(a;u,\pi)\|_{L^2_y}\\
 &\qquad\quad
 +C_1\varepsilon
 \left(\|\nabla_yu\|_{L^2_y}
       +\|u\|_{L^2_y}
       +\|\nabla_y\pi\|_{L^2_y}\right),                      \tag{4.3}
\end{aligned}
\]

where \(C_0,C_1\) depend on the fixed cross-section and the constants in
(4.1), but **not** on \(2\pi/\varepsilon\).

The derivative estimate used here is

\[
 \|\nabla_x{\cal P}_\varepsilon u\|_{L^2_x}
 \le C(\delta_0)
 \left(\|\nabla_yu\|_{L^2_y}
                +\varepsilon\|u\|_{L^2_y}\right).             \tag{4.4}
\]

Let \(\mathbb P_{\rm amb}\) be the one fixed Leray projector on
\(\mathbb R^3\). Since it is an \(L^2\) contraction, (4.3) immediately
gives the same upper bound for the projected ambient residual. No inverse
Laplacian on the long tube occurs.

## 5. Frequency-\(p\) specialization and winding

Let the straight q-free packet satisfy

\[
\begin{aligned}
 f_p&={\cal M}_0(a;u_p,\pi_p),\\
 \|\nabla u_p\|_2+\|\nabla\pi_p\|_2
   &\le C_{\rm env}p\|u_p\|_2+C_{\rm env}\|f_p\|_2.
                                                                    \tag{5.1}
\end{aligned}
\]

For the reconstruction in the companion q-free note, the first line is
the exact force

\[
 f_p=e^{i(\alpha z-n\theta-\omega t)}
       \frac{\gamma}{in^2}G_p\,e_r,                           \tag{5.2}
\]

so, when \(|\gamma|\le\gamma_1\) and \(n=p\),

\[
 \frac{\|f_p\|_2}{\|u_p\|_2}
 \le Cp^{-2}\frac{\|G_p\|_2}{\|H_p\|_2}.                     \tag{5.3}
\]

For the q-free reconstruction this is the natural first-order frequency
envelope: the spatial Fourier covectors and the radial WKB covector are
\(O(p)\), while the exact straight momentum equation also bounds
\(\nabla\pi_p\) from \(|\omega_p|=O(p)\).

More concretely, on the fixed annulus assume

\[
 0<\gamma_0\le|\gamma|\le\gamma_1,\qquad
 \|\partial_r^kH_p\|_2\le C_Hp^k\|H_p\|_2
 \quad(0\le k\le2),                                         \tag{5.4}
\]

with fixed \(C^2\) profile coefficients and
\(|\alpha|+|n|\lesssim p\). The exact q-free formulas give
\(\|u_p\|_2\gtrsim\|H_p\|_2\) from \(u_r=H_p/r\),
\(\|\nabla u_p\|_2\lesssim p\|H_p\|_2\), and
\(\|\nabla\pi_p\|_2\lesssim p\|H_p\|_2\). Hence (5.1) follows
from (5.4). The unproved part is the uniform construction of \(H_p\)
with its exponentially small \(G_p\), not a further velocity-pressure
or geometric inversion.

Combining (1.5), (4.3), and (5.1) gives

\[
\begin{aligned}
 \frac{\|\mathbb P_{\rm amb}
 {\cal M}_\varepsilon
 (A_\varepsilon;{\cal P}_\varepsilon u_p,\Pi_p)\|_2}
 {\|{\cal P}_\varepsilon u_p\|_2}
 &\le
 C\left(
 \frac{\|f_p\|_2}{\|u_p\|_2}
          +\varepsilon p+\varepsilon\right)\\
 &\le C\left(
 p^{-2}\frac{\|G_p\|_2}{\|H_p\|_2}
          +\varepsilon p+\varepsilon\right).
\end{aligned}                                                \tag{5.5}
\]

Thus the desired finite-curvature power is genuinely
\(O(\varepsilon p)\) in the physical ambient \(L^2\) residual ratio,
provided (5.1) holds uniformly for the scalar packet being used.

The axial period introduces no derivative larger than \(p\). For
cross-sectional charge \(n=p\), choose

\[
 m=\operatorname{round}\left(\frac{\beta_*p}{\varepsilon}\right),
 \qquad
 \alpha_\varepsilon=\varepsilon m,\qquad
 \beta_\varepsilon=\frac{\alpha_\varepsilon}{p}.             \tag{5.6}
\]

Then \(e^{i\alpha_\varepsilon z}\) is exactly periodic and

\[
 |\beta_\varepsilon-\beta_*|
 \le\frac{\varepsilon}{2p},\qquad
 |\alpha_\varepsilon|\le |\beta_*|p+\frac{\varepsilon}{2}.   \tag{5.7}
\]

Using \(\beta_\varepsilon\) directly in the q-free formulas creates no
lattice residual. The perturbation from the target pitch is lower than
the curvature term in (5.5).

## 6. Viscosity and the precise remaining obstruction

The same fixed-chart calculation for the Euclidean vector Laplacian gives

\[
 \|\Delta_x{\cal P}_\varepsilon u
       -{\cal P}_\varepsilon\Delta_yu\|_2
 \le C\varepsilon
 \left(\|u\|_{H^2}+\|u\|_{H^1}+\|u\|_2\right),               \tag{6.1}
\]

again with no aspect-length loss. Consequently a viscous linearized
operator has the additional frequency-\(p\) transport error

\[
                         O(\varepsilon\nu p^2).               \tag{6.2}
\]

It is lower than the inviscid \(O(\varepsilon p)\) term only in the
high-Reynolds window \(\nu p\lesssim1\).

The geometric/Hodge obstruction is therefore removed for compact
fixed-cross-section packets. The unresolved obligations are now more
specific:

1. prove the uniform analytic coalescing-edge scalar pseudomode and its
   frequency envelope (5.1), including every lower coefficient;
2. choose \(\varepsilon,p\), and the edge deficit so that the algebraic
   curvature residual in (5.5) remains small over the required gain time;
3. retain \(\varepsilon\nu p^2\) if the viscous operator, rather than the
   Euler linearization, is transported;
4. keep the packet support in a fixed normalized cross-section. If the
   support radius grows like \(\varepsilon^{-1}\), then \(h\) loses
   uniform ellipticity and this lemma gives no bound.

In particular, (5.5) is a quasimode estimate, not an isolated spectral
point. It closes the aspect-aware Piola bookkeeping but not the
coalescing-WKB theorem or the nonlinear return construction.

## 7. Reproducibility

[The dependency-free exact checker](../checks/aspect_uniform_piola_packet.py)
verifies with rational arithmetic:

* the Jacobian and weighted \(L^2\) density;
* exact Piola divergence preservation;
* every component of the geometric commutator (3.4);
* the pressure commutator (3.5);
* the full decomposition (3.7) on independent first jets; and
* the corrected axial winding inequality (5.7).

The checker does not claim to verify the analytic estimates (4.3),
(5.1), (5.4)--(5.5), or (6.1); those follow from the displayed coefficient
formulas under their stated hypotheses.
