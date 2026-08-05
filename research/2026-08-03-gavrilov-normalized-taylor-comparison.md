# Hollow-supported Gavrilov conormal comparison for the curved carrier

**Date:** 2026-08-03
**Status:** corrected after adversarial source audit; local coefficient
comparison repaired, ambient/global Hodge comparison still open
**Scope:** the finite-curvature base-flow comparison required by C79 and
C82. This note does not construct the coalescing WKB amplitudes, a global
periodic pseudomode, or the nonlinear transition.

## 1. Correction

The first version of this note asserted a raw analytic expansion

\[
                         u(\zeta)=L\zeta+O(|\zeta|^2).          \tag{1.1}
\]

That assertion is false. Gavrilov explicitly notes that his unmodulated
field is not smooth on the seed circle. The axial component has a radial
cusp. In local normal coordinates \(\zeta=(X,Z)\), put

\[
 r=|\zeta|,\qquad J(X,Z)=(Z,-X).
\]

In Gavrilov's cylindrical frame \((e_\rho,e_z,e_\phi)\), the source-level
component expansion is

\[
 \begin{aligned}
 p(\zeta)&=\frac{r^2}{2}+\frac34Xr^2+O(r^4),\\
 u_\rho&=Z+\frac12XZ+O(r^3),\\
 u_z&=-X-\frac54X^2-\frac34Z^2+O(r^3),\\
 u_\phi&=\frac r{\sqrt2}-\frac{Xr}{4\sqrt2}+O(r^3).
 \end{aligned}                                                \tag{1.2}
\]

In particular the degree-one part is

\[
 U_1(X,Z)=Ze_\rho-Xe_z+\frac r{\sqrt2}e_\phi.                 \tag{1.3}
\]

The remainder is analytic with uniform derivative bounds on every closed
annulus \(0<c\le r\le C\), but no smooth Taylor expansion through \(r=0\)
is claimed. These formulas follow from Gavrilov's source variables

\[
 \alpha=2r^2+3Xr^2+O(r^4),\qquad p=\frac{\alpha}{4},
 \qquad b=\frac14\sqrt{H(\alpha)},\quad H(\alpha)=4\alpha+O(\alpha^2).
                                                                    \tag{1.4}
\]

The error in (1.1) does **not** destroy the hollow pressure-modulated
carrier. The exact Gavrilov multiplier is identically zero near \(p=0\).
Its normalized support is therefore a fixed annulus, where \(r\) is analytic
and the conormal expansion has ordinary uniform \(C^k\) bounds.

## 2. Correct normalized theorem

Let \(g\) be a fixed smooth or Gevrey multiplier with

\[
 \operatorname{supp}g\subset[q_-,q_+],\qquad 0<q_-<q_+<\infty.
                                                                    \tag{2.1}
\]

Define the exact compact steady Euler field

\[
 W_\varepsilon(\zeta)
 =g\left(\frac{p(\zeta)}{\varepsilon^2}\right)u(\zeta),       \tag{2.2}
\]

and in normalized coordinates \(y=\zeta/\varepsilon\) set

\[
 \widehat U_\varepsilon(y)
 =\frac1\varepsilon W_\varepsilon(\varepsilon y).             \tag{2.3}
\]

Put

\[
 p_2(y)=\frac{|y|^2}{2},\qquad
 U_0(y)=g(p_2(y))
 \left(y_2e_\rho-y_1e_z+\frac{|y|}{\sqrt2}e_\phi\right).      \tag{2.4}
\]

Then, on every fixed normalized tube and for every fixed \(k\),

\[
             \|\widehat U_\varepsilon-U_0\|_{C^k}
             \le C_{k,g}\varepsilon.                          \tag{2.5}
\]

The norm in (2.5) is an ordinary Cartesian \(C^k\) norm. It is legitimate
because both sides vanish on a fixed neighborhood of \(y=0\). The leading
field \(U_0\) is exactly the hollow locked-pitch straight column. In the
cross-sectional polar frame its first two components are \(-r e_\theta\).
Choose the column axial direction as \(e_z^{\rm col}=-e_\phi\). Then

\[
 V(r)=-r g(r^2/2),\qquad
 W(r)=-\frac{r g(r^2/2)}{\sqrt2}=\frac{V(r)}{\sqrt2}.         \tag{2.6}
\]

This fixed orientation change can be absorbed into the helical Fourier
convention. It does not change the strain formula or the full BAS edge.

Thus the original **raw linearization** is refuted, while its needed
**hollow-supported coefficient conclusion** survives in the corrected
form (2.4)--(2.5).

## 3. Proof on the active annulus

On the fixed annulus selected by (2.1), scaling (1.2) gives

\[
 \frac{p(\varepsilon y)}{\varepsilon^2}
 =p_2(y)+\varepsilon\frac34y_1|y|^2+O_{C^k}(\varepsilon^2),
                                                                    \tag{3.1}
\]

and

\[
 \frac{u(\varepsilon y)}{\varepsilon}
 =U_1(y)+\varepsilon U_2(y)+O_{C^k}(\varepsilon^2),           \tag{3.2}
\]

where

\[
\begin{aligned}
 U_2(y)={}&\frac12y_1y_2e_\rho
 -\left(\frac54y_1^2+\frac34y_2^2\right)e_z\\
 &-\frac{y_1|y|}{4\sqrt2}e_\phi.
\end{aligned}                                                 \tag{3.3}
\]

Composition and multiplication imply

\[
\begin{aligned}
 \widehat U_\varepsilon
 &=
 g\left(p_2+\varepsilon\frac34y_1|y|^2
                +O_{C^k}(\varepsilon^2)\right)\\
 &\qquad\cdot
 \left(U_1(y)
                +O_{C^k}(\varepsilon)\right)\\
 &=U_0+O_{C^k}(\varepsilon).                                 \tag{3.4}
\end{aligned}
\]

It remains to justify the regions outside that annulus. Since \(g=0\) on a
neighborhood of \([0,q_-]\), continuity of (3.1) gives a fixed inner
neighborhood on which both the exact and limiting products vanish for all
small \(\varepsilon\). The same argument applies beyond the outer support.
The transition regions are compact subsets of \(r>0\), so (3.4) applies
there. This proves (2.5).

The key distinction is:

* \(u\) itself has only a conormal radial leading term at the seed circle;
* \(g(p/\varepsilon^2)u\) is smooth because the multiplier removes a full
  neighborhood of that circle;
* all normalized derivatives are estimated only after this exact
  multiplication.

## 4. Matching the strain-capped profile

For \(p_2=r^2/2\), the magnitude of the transverse leading term is
\(r|g(r^2/2)|\), while the axial magnitude is smaller by \(1/\sqrt2\).
After fixing the orientation in (2.6), choose

\[
                    g(q)=-\Omega_{\rm cap}(\sqrt{2q}).         \tag{4.1}
\]

The C82 strain-capped profile is hollow, so (4.1) is smooth at the only
place where the square root could be singular. Near the selected spectral
ring it agrees with the analytic log-normal jet; the Gevrey ramp is used
only away from that packet. For every \(\varepsilon>0\), (2.2) remains an
exact compact steady Euler flow.

The tubular metric in normalized coordinates is

\[
 ds^2=dr^2+r^2d\theta^2+
       (1+\varepsilon r\cos\theta)^2dz^2.                     \tag{4.2}
\]

On the fixed support, all fixed derivatives of the metric, orthonormal
frame, and connection differ from the straight cylinder by
\(O(\varepsilon)\). Combining this with (2.5) gives the local coefficient
estimates

\[
 \|g_\varepsilon-g_0\|_{C^2}
 +\|\widehat U_\varepsilon-U_0\|_{C^1}
 \le C\varepsilon,                                           \tag{4.3}
\]

and

\[
 \|\operatorname{sym}_{g_\varepsilon}\nabla\widehat U_\varepsilon
 -\operatorname{sym}_{g_0}\nabla U_0\|_\infty
 \le C\varepsilon.                                           \tag{4.4}
\]

Equation (4.4) is enough for the curved **energy** exponent in C82:
the exact base is divergence free, is supported in this fixed tube, and its
strain maximum moves by at most \(C\varepsilon\).

## 5. What this does not yet prove

The local first-order generator coefficients differ by \(O(\varepsilon)\),
so applying the curved differential expression to a frequency-\(p\) packet
creates the formal local residual

\[
                         O(\varepsilon p).                    \tag{5.1}
\]

However, (4.3) alone does not prove the ambient/global Hodge comparison
claimed in C79. The normalized axial coordinate in (4.2) has period
\(2\pi/\varepsilon\); a Laplacian-resolvent argument on an unspecified
compact tube is not automatically uniform on this long cylinder, nor does
it by itself identify the tube projector with the Leray projector on the
ambient \(\mathbb R^3\) or \(\mathbb T^3\) Clay domain.

Accordingly:

1. the hollow-supported \(C^1\) base-flow comparison is repaired;
2. the C82 curved strain/energy bound follows locally from that comparison;
3. the \(O(\varepsilon p)\) differential residual is the correct formal
   power;
4. converting it into a global divergence-free ambient pseudomode with the
   same normalized \(L^2\) residual remains conditional on an
   aspect-uniform Piola/Hodge or direct solenoidal construction.

The periodic winding correction is recorded separately in C84: because the
central \(z\)-period is \(2\pi/\varepsilon\), its Fourier lattice is
\(\varepsilon\mathbb Z\).

## 6. Reproducibility

[The corrected scaling checker](../checks/gavrilov_normalized_taylor.py)
records:

* the source expansion \(p=r^2/2+(3/4)Xr^2+O(r^4)\);
* the radial axial leading term and its failure to be linear;
* the \(O(\varepsilon)\) annular scaling powers;
* the locked-pitch squared magnitude ratio \(1/2\); and
* the necessity of a hollow support \(q_->0\).
