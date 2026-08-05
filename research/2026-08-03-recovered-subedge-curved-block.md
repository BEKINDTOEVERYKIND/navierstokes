# Recovered forced-subedge and curved pressure-block checkpoint

**Date:** 2026-08-03

**Status:** exact local algebra, a source-backed scalar pseudomode reduction,
and a quantitative curvature comparison.  The analytic WKB construction and
the global finite-curvature carrier are self-audited only.  This is not a
Navier--Stokes blow-up theorem.

## 1. Purpose and outcome

This note restores the last spectral checkpoint that was derived after the
live checkout at commit `2e7b481`.  It also removes one ambiguity in that
checkpoint.

The useful conclusions are:

1. the forced velocity--pressure system can be reduced without the singular
   Albritton--Ożański quotient;
2. a fixed subedge spectral parameter produces a scalar analytic
   semiclassical symbol with a genuine pseudospectral bracket;
3. the frozen velocity--pressure pencil has an explicitly invertible
   three-dimensional constraint block and the scalar Rayleigh factor is its
   exact Schur complement;
4. the principal Hodge symbol on a thin torus is an order-
   \(\varepsilon\) perturbation of the straight one; and
5. the original global Hodge-projector comparison is unsupported on the
   long normalized tube and is withdrawn. A safer ambient route transports
   an exactly divergence-free compact straight packet by Piola, compares
   the unprojected momentum equations, and applies the fixed ambient Leray
   projector only to the residual. This route has the formal curvature cost
   \(C\varepsilon p\), conditional on the complete velocity--pressure
   reconstruction and aspect-uniform coefficient bounds.

The last estimate alone is not enough when a fixed subedge gap is propagated
for \(O(j^2)\) time.  The companion note
`2026-08-03-edge-tracking-coalescing-pseudomode.md` removes that exponential
gap by allowing the pseudomode growth to approach the BAS edge.

The principal references are the ring-mode construction of
[Albritton--Ożański](https://arxiv.org/abs/2310.20674), the scalar analytic
pseudomode theorem of
[Dencker--Sjöstrand--Zworski](https://arxiv.org/abs/math/0301242), and the
system extension of [Dencker](https://arxiv.org/abs/0705.4561).

## 2. Forced q-free reconstruction

Use the Fourier convention

\[
       e^{i(\alpha z-n\theta-\omega t)},\qquad
       \alpha=\beta n,
\]

and put

\[
 \Omega={V\over r},\qquad \Gamma=rV,\qquad
 \Lambda=\beta W-\Omega,\qquad
 \gamma=n\Lambda-\omega,\qquad D=1+\beta^2r^2.
\]

For a force \(F=(F_r,F_\theta,F_z)\), the Fourier velocity equations are

\[
\begin{aligned}
 i\gamma u_r-2\Omega u_\theta+P'&=F_r,\\
 i\gamma u_\theta+{\Gamma'\over r}u_r-{in\over r}P&=F_\theta,\\
 i\gamma u_z+W'u_r+i\alpha P&=F_z,\\
 {1\over r}(ru_r)'-{in\over r}u_\theta+i\alpha u_z&=0.
\end{aligned}                                                    \tag{2.1}
\]

The last three equations give the pressure correction

\[
 Q_F={ir(F_\theta-\beta rF_z)\over nD},\qquad P=P_0+Q_F,          \tag{2.2}
\]

where \(P_0\) is the unforced expression in (6.3a) of the locked-pitch
profile note.  Substitution into the radial equation gives the scalar force

\[
 G_F={in^2\over\gamma}\left[
 F_r-Q_F'-{2i\Omega\over\gamma}F_\theta
       +{2n\Omega\over r\gamma}Q_F
 \right].                                                        \tag{2.3}
\]

Equations (2.2)--(2.3) contain only the physical coefficients.  In
particular, they remain meaningful where
\(q_{\rm AO}=-(rV)'/W'\) has a pole.

## 3. Fixed-subedge analytic scalar family

Let \(p=n\), \(s=p^{-1/2}\), and let \(r_*\) be the locked-pitch ring where

\[
 b(r_*)=b_*=\lambda_*^2>0,\qquad
 \Lambda'(r_*)=0,\qquad \Lambda''(r_*)=\kappa>0.
\]

Choose fixed \(\sigma\ne0\) and \(0<y_0<\lambda_*\), and set

\[
 r_{c,p}=r_*+\sigma s,\qquad
 \omega_p=p\Lambda(r_{c,p})+iy_0,\qquad
 Y={r-r_{c,p}\over s}.                                         \tag{3.1}
\]

After the exact q-free scalar equation is divided by its elliptic weight,
it has the analytic semiclassical form

\[
 P_s=-s^2A_s(Y)\partial_Y^2-s^3B_s(Y)\partial_Y
       +1+{b_s(Y)\over\gamma_s(Y)^2}
       +s^2{a_s(Y)\over\gamma_s(Y)}.                            \tag{3.2}
\]

On compact \(Y\)-sets,

\[
 \gamma_s(Y)=\kappa\left(\sigma Y+{Y^2\over2}\right)-iy_0+O(s),
\]

and the limiting symbol is

\[
 \mathfrak p_0(Y,\eta)
 =A_*\eta^2+1+{b_*\over
   [\kappa(\sigma Y+Y^2/2)-iy_0]^2}.                            \tag{3.3}
\]

At \(Y=0\), choose

\[
 \eta_0^2={b_*/y_0^2-1\over A_*}>0.                            \tag{3.4}
\]

Then \(\mathfrak p_0(0,\eta_0)=0\), and, with the convention
\(\{f,g\}=f_\eta g_Y-f_Yg_\eta\),

\[
 \{\Re\mathfrak p_0,\Im\mathfrak p_0\}(0,\eta_0)
 ={4A_*\eta_0b_*\kappa\sigma\over y_0^3}\ne0.                \tag{3.5}
\]

Changing the sign of \(\eta_0\) selects the pseudospectral orientation.
The analytic scalar theorem therefore gives a straight-column packet with
residual

\[
             \|P_s\phi_s\|_2
             \le C p^D e^{-c\sqrt p}\|\phi_s\|_2,              \tag{3.6}
\]

for fixed \(y_0,\sigma\).  Its radial oscillation wavelength is \(p^{-1}\),
while its Gaussian physical width is \(p^{-3/4}\).

This is a pseudomode statement, not an isolated-eigenvalue statement.

## 4. Exact frozen velocity--pressure Schur complement

The scalar factor is not an artifact of eliminating pressure.  Freeze the
straight coefficients and use the high-frequency radial covector \(\eta\).
Put

\[
 A={\Gamma'\over r},\qquad T=W'.
\]

After scaling physical pressure by \(p\), the frozen pencil for
\((u_r,u_\theta,u_z,P)\) is

\[
 \mathcal M=
 \begin{pmatrix}
 i\gamma&-2\Omega&0&i\eta\\
 A&i\gamma&0&-i/r\\
 T&0&i\gamma&i\beta\\
 \eta&-1/r&\beta&0
 \end{pmatrix}.                                                \tag{4.1}
\]

The block formed by rows 2--4 and columns \((u_\theta,u_z,P)\) has

\[
             \det\mathcal C={\gamma D\over r^2}.               \tag{4.2}
\]

It is uniformly invertible whenever \(|\Im\gamma|\ge y_0/2\).  Its exact
radial Schur complement is

\[
 R=i\gamma\left({r^2\over D}\eta^2+1+{b\over\gamma^2}\right)
       -{r^2\over D}\eta\Lambda',                              \tag{4.3}
\]

where

\[
 b=-{2\beta V(T+\beta\Gamma')\over D},\qquad
 \Lambda'=\beta T-\Omega'.                                    \tag{4.4}
\]

Consequently

\[
                    \det\mathcal M={\gamma D\over r^2}R.       \tag{4.5}
\]

At the designed ring \(\Lambda'=0\), the zero of \(R\) is exactly the
zero of (3.3), multiplied by the elliptic factor \(i\gamma\).  At
\(r_{c,p}\), the omitted term is only \(O(s)\).  Thus the characteristic
is simple, the complement stays elliptic, and the scalar bracket controls
the full pencil.

## 5. Thin-torus Hodge symbol

In normalized tubular coordinates the Euclidean metric is

\[
 ds^2=dr^2+r^2d\theta^2+H_\varepsilon^2dz^2,\qquad
 H_\varepsilon=1+\varepsilon r\cos\theta.                       \tag{5.1}
\]

Its volume density is \(J_\varepsilon=rH_\varepsilon\), and

\[
 \operatorname{div}_{g_\varepsilon}u
 ={1\over rH_\varepsilon}\partial_i(rH_\varepsilon u^i),
 \qquad
 \nabla_{g_\varepsilon}q
 =(q_r,r^{-2}q_\theta,H_\varepsilon^{-2}q_z).                   \tag{5.2}
\]

For a nonzero covector \(\xi\), the principal Hodge projector is exactly

\[
 \Pi_\varepsilon(\xi)^i{}_j
 =\delta^i_j-{(g_\varepsilon^{ik}\xi_k)\xi_j
                    \over g_\varepsilon^{k\ell}\xi_k\xi_\ell}.
                                                                    \tag{5.3}
\]

It satisfies

\[
 \Pi_\varepsilon^2=\Pi_\varepsilon,qquad
 \xi_i\Pi_\varepsilon^i{}_j=0,                                 \tag{5.4}
\]

and is \(O(\varepsilon)\)-close to the straight-cylinder projector in
every fixed symbol seminorm on a compact tubular chart.  There is no
frequency power in this principal-symbol comparison.

The Piola transform supplies the corresponding exact divergence-preserving
identification of vector fields.  For a tubular map \(F_\varepsilon\) with
Jacobian \(J_F\),

\[
 u_\varepsilon={1\over J_F}DF_\varepsilon\,u_0\circ F_\varepsilon^{-1}
 \quad\Longrightarrow\quad
 \operatorname{div}u_\varepsilon
 ={1\over J_F}(\operatorname{div}u_0)\circ F_\varepsilon^{-1}.  \tag{5.5}
\]

## 6. The global Hodge comparison is withdrawn; an ambient Piola route

The former version placed global Hodge projectors
\(\mathbb P_\varepsilon\) and \(\mathbb P_0\) on an unspecified "fixed
compact tubular manifold" and asserted

\[
 \|\mathbb P_\varepsilon-\mathbb P_0\|_{L^2\to L^2}
 \le C\varepsilon.                                           \tag{6.1}
\]

That did not prove the needed statement. In the normalized metric (5.1),
the central \(z\)-period is \(2\pi/\varepsilon\). A scalar resolvent argument
on a nominally fixed tube is not automatically uniform on this long
cylinder, and it does not identify either projector with the Leray
projector on the ambient \(\mathbb R^3\) or \(\mathbb T^3\) Clay domain.
Equation (6.1) is therefore withdrawn.

There is a safer comparison that never subtracts changing global
projectors.

1. Reconstruct a compactly supported straight velocity--pressure
   pseudomode \((u_p,\pi_p)\) from the scalar mode using the exact formulas
   in Sections 2--4. This makes \(\operatorname{div}u_p=0\) before bending.
2. Transport \(u_p\) into the physical torus by the Piola map (5.5), which
   preserves divergence exactly, and transport \(\pi_p\) as a scalar.
3. Keep the packet support strictly inside the tubular chart. Extension by
   zero then gives an ambient divergence-free velocity.
4. Compare the **unprojected** straight and curved momentum equations in
   the common tubular coordinates.
5. Apply the one fixed ambient Leray projector to the resulting residual.
   Its \(L^2\) norm is one.

The corrected hollow Gavrilov expansion gives, on the packet support,

\[
 \|g_\varepsilon-g_0\|_{C^2}
 +\|U_\varepsilon-U_0\|_{C^1}\le C\varepsilon.                \tag{6.2}
\]

Consequently the unprojected residual has the conditional estimate

\[
\begin{aligned}
 \|\operatorname{Res}_\varepsilon(u_p,\pi_p)\|_2
 \le{}& Cp^D e^{-c\sqrt p}\|u_p\|_2\\
 &+C\varepsilon
 \big(\|u_p\|_{H^1}+\|\pi_p\|_{H^1}\big).                    \tag{6.3}
\end{aligned}
\]

The q-free reconstruction has elliptic denominators on the packet and hence
only polynomial pressure/velocity loss. Since
\(\|\nabla u_p\|_2\lesssim p\|u_p\|_2\), the target normalized power is

\[
 {\|\mathbb P_{\rm amb}\operatorname{Res}_\varepsilon\|_2
       \over\|u_p\|_2}
 \le \operatorname{poly}(p,\eta^{-1})e^{-c/\hbar_{\rm eff}}
   +\operatorname{poly}(p,\eta^{-1})\,\varepsilon p.          \tag{6.4}
\]

Thus \(m=1\) is the correct **local differential** curvature power, but the
ambient claim remains conditional until the exact reconstructed packet,
Piola normalization, support cutoff, and polynomial \(\pi_p\) bounds are
written in one normed operator lemma. Unlike (6.1), this route has no
aspect-ratio-dependent global Hodge inverse.

## 7. Exact curved-symbol route versus the comparison route

There are now two distinct curvature strategies.

1. **Exact curved symbol.**  Analytically eliminate the uniformly elliptic
   block (4.2), use the implicit-function theorem to continue the simple
   characteristic and its bracket, and apply the system/scalar analytic
   pseudomode theorem to the true curved operator.  This would retain an
   exponentially small residual but still requires a fully written
   two-microlocal toroidal calculus.
2. **Bent straight comparison.** Use the ambient Piola route (6.2)--(6.4).
   This avoids the exact curved normal form, but a fixed gap
   \(\lambda_*-y_0\) would amplify the
   \(\varepsilon p\) term over a \(j^2\) gain interval.  The edge-tracking
   construction in the companion note removes that fixed gap.

The second strategy is now shorter.

## 8. Remote-sideband correction retained for audit

The alternative adjacent-sideband route must use the full scalar symbol
\(\xi^2+q_p\), not the potential \(q_p\) alone.  In the nearest target
sideband the negative pocket has depth \(O(p^{-1/2})\), width \(O(p^{-1})\),
interior wavelength \(O(p^{-3/4})\), and Airy scale \(O(p^{-5/6})\).
An Airy estimate by itself does not prove the needed half-line inverse;
the unresolved object is a coalescing-pocket global resolvent estimate.
Because Sections 3--7 avoid the adjacent sideband altogether, that route is
kept only as a parallel audit target.

## 9. Remaining obligations

This note does **not** establish a Clay-admissible singularity.  The
load-bearing open points are:

1. the complete exact velocity--pressure reconstruction and ambient
   Piola residual lemma described in Section 6;
2. uniform analytic WKB constants in the edge-tracking, coalescing-root
   regime of the companion note;
3. quantitative normalization on the long aspect-ratio-dependent torus,
   including its corrected axial Fourier lattice; and
4. the nonlinear material-phase endpoint inverse with its global viscous
   wake.
