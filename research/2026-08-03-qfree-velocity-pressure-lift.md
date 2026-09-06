# Exact q-free velocity--pressure lift after scalar localization

**Date:** 2026-08-03

**Status:** exact Fourier operator identity and polynomial norm ledger;
curved Piola comparison remains open.

**Scope:** this note repairs the straight velocity reconstruction obligation
in C79/C84. It does not prove the uniform analytic scalar pseudomode or the
finite-curvature ambient estimate.

## 1. Outcome

The localized scalar Rayleigh residual can be lifted to a complete
divergence-free velocity--pressure residual without a Hodge correction.

Use the Fourier convention

\[
 e^{i(\alpha z-n\theta-\omega t)},\qquad
 \alpha=\beta n,
\]

and define

\[
\begin{gathered}
 \Omega=\frac Vr,\qquad \Gamma=rV,\qquad
 \Lambda=\beta W-\Omega,\qquad
 \gamma=n\Lambda-\omega,\qquad D=1+\beta^2r^2,\\
 a=r\frac d{dr}
 \left(\frac{\beta r^2W'-\Gamma'}{rD}\right),\qquad
 b=-\frac{2\beta V(W'+\beta\Gamma')}{D}.                     \tag{1.1}
\end{gathered}
\]

For an arbitrary compactly supported scalar \(H(r)\) in an annulus where
\(|\gamma|\ge y_0>0\), put

\[
 u_r=\frac Hr,                                                \tag{1.2}
\]

\[
 P=-\frac{i\gamma r}{n^2D}H'
   +\frac{i(\beta r^2W'-\Gamma')}{nrD}H,                      \tag{1.3}
\]

and

\[
 u_\theta=\frac{i\Gamma'}{r\gamma}u_r
           +\frac{nP}{r\gamma},\qquad
 u_z=\frac{iW'}{\gamma}u_r-\frac{\alpha P}{\gamma}.           \tag{1.4}
\]

Define the q-free scalar residual

\[
 G=
 \frac d{dr}\left(\frac rD H'\right)
 -\frac{n^2}{r}
 \left(1+\frac{a}{n\gamma}+\frac{b}{\gamma^2}\right)H.        \tag{1.5}
\]

Then the azimuthal momentum, axial momentum, and divergence equations hold
exactly, while the radial residual is

\[
 F_\theta=F_z=0,\qquad
                         F_r=\frac{\gamma}{in^2}G.             \tag{1.6}
\]

Thus scalar localization should be performed **before** reconstruction.
All cutoff errors are included in \(G\), and (1.2)--(1.4) produce a
divergence-free compact straight packet directly.

## 2. Algebraic verification

The forced Fourier system is

\[
\begin{aligned}
 i\gamma u_r-2\Omega u_\theta+P'&=F_r,\\
 i\gamma u_\theta+\frac{\Gamma'}r u_r-\frac{in}rP&=F_\theta,\\
 i\gamma u_z+W'u_r+i\alpha P&=F_z,\\
 \frac1r(ru_r)'-\frac{in}r u_\theta+i\alpha u_z&=0.
                                                               \tag{2.1}
\end{aligned}
\]

Substitution of (1.4) makes the middle two momentum equations identities.
The divergence equation becomes

\[
 \frac{H'}r
 -\frac{n(\beta r^2W'-\Gamma')}{r^3\gamma}H
 -\frac{in^2D}{r^2\gamma}P=0,                                \tag{2.2}
\]

which is exactly (1.3). Differentiating (1.3), substituting it and (1.4)
into the radial equation, and collecting the physical coefficients gives

\[
 \frac{in^2}{\gamma}
 \left(i\gamma u_r-2\Omega u_\theta+P'\right)=G.              \tag{2.3}
\]

No quotient by \(W'\), \(\Gamma'\), or the singular AO auxiliary function
appears. The exact rational/Gaussian-rational checker evaluates (2.1)--(2.3)
for independent background and scalar two-jets.

## 3. Polynomial reconstruction bounds

Fix a compact annulus \(r_0\le r\le r_1\), take \(n\ge1\), and assume

\[
 |\gamma|\ge y_0,\qquad
 \|V,W\|_{C^{m+2}}\le C_m,\qquad
 0<\beta_0\le\beta\le\beta_1,\qquad
 |\omega|\le C_0(1+n).                                      \tag{3.1}
\]

The denominators \(r,D,\gamma\) in (1.2)--(1.4) are then uniformly elliptic.
For every fixed \(m\), these background assumptions give
\(\|\gamma\|_{C^{m+1}}\le C_m(1+n)\). Directly,

\[
 \|P\|_{H^m}
 \le C_m\left(n^{-2}\|\gamma H'\|_{H^m}
                  +n^{-1}\|H\|_{H^m}\right)
 \le C_m n^{-1}\|H\|_{H^{m+1}},                              \tag{3.2}
\]

where the final constant uses the stated linear spectral-frequency bound.
Since \(\gamma'=n\Lambda'\), derivatives of \(\gamma^{-1}\) in (1.4) cost
at most a fixed power of \(n/y_0\). Likewise,

\[
 \|(u_r,u_\theta,u_z)\|_{H^m}
 \le C_m\operatorname{poly}(n,y_0^{-1})
                  \|H\|_{H^{m+1}},                            \tag{3.3}
\]

and

\[
 \|F_r\|_{H^m}
 \le C_m n^{-2}\|\gamma G\|_{H^m}
 \le C_m n^{-1}\|G\|_{H^m}.                                  \tag{3.4}
\]

In the edge-tracking regime \(y_j\to\lambda_*>0\), so
\(|\gamma|\) stays bounded below. Every reconstruction loss is polynomial
in \(n\), \(\eta_j^{-1}\), and fixed Sobolev order. It is absorbed by a
conditional scalar residual \(e^{-c/\hbar_j}\).

These are annular estimates. The packet is cut off strictly inside the
analytic ring neighborhood, so no axis or exterior endpoint analysis is
needed for the pseudomode.

## 4. Consequence for the curved route

The corrected construction order is:

1. build and localize the scalar packet \(H_j\);
2. compute its exact scalar residual \(G_j\);
3. reconstruct \((u_j,P_j)\) by (1.2)--(1.4);
4. regard (1.6) as the complete straight momentum residual;
5. Piola-transport the already divergence-free compact velocity into the
   physical tube; and
6. compare the unprojected curved and straight momentum equations before
   applying the one fixed ambient Leray projector to the residual.

Steps 1 and 5--6 remain open at the uniform analytic/ambient norm level.
This note closes the algebraic straight lift and removes the unsupported
instruction to repair divergence using a global projector on the changing
long tube.

## 5. Reproducibility

[The exact checker](../checks/qfree_velocity_pressure_lift.py) uses several
independent Gaussian-rational jets to verify all four Fourier equations and
the scalar-to-radial residual factor (1.6), with no numerical tolerance.
