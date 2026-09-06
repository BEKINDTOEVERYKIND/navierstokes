# Strain-capped locked-pitch compactification and the full PDE semigroup

**Date:** 2026-08-03

**Status:** exact energy identity and algebraic critical-point analysis;
Gevrey compactification by an open-ellipse path self-derived; finite-curvature
coefficient comparison still requires the common-coordinate Taylor ledger.

**Scope:** the locked-pitch Gavrilov profile of C68.  This note supplies an
\(L^2\) semigroup estimate for a slightly redesigned compact modulation.  It
does not construct the nonlinear transition or a Navier--Stokes singularity.

## 1. Outcome

The missing frequency-uniform PDE semigroup bound does not require a
Fourier-block resolvent estimate.

For a straight locked-pitch column
\[
 U=V(r)e_\theta+W(r)e_z,\qquad
 V=r\Omega,\qquad W={V\over\sqrt2},
\]
the operator norm of the symmetric velocity gradient is
\[
 {\cal S}(r)
 ={1\over2}\left[
      (r\Omega')^2+{1\over2}(\Omega+r\Omega')^2
   \right]^{1/2}.                                               \tag{1.1}
\]

At the designed spectral ring,
\[
                         {\cal S}(r_*)=\lambda_*                 \tag{1.2}
\]
exactly.  The unmodified log-normal continuation has a second, larger strain
maximum on its spectrally stable inner flank.  That flank is not needed by
the AO packet.  It can be replaced, before the extra maximum, by a long
hollow Gevrey ramp satisfying
\[
                         {\cal S}(r)<\lambda_*                   \tag{1.3}
\]
throughout the replacement.

The outer cutoff can likewise be chosen below the strict margin.  The
resulting compact Gavrilov pressure modulation preserves the complete jet
at the ring, its isolated phase level, and its unique full BAS edge, while
also obeying
\[
                       \sup_x\|\operatorname{sym}\nabla U(x)\|_{\rm op}
                       =\lambda_* .                              \tag{1.4}
\]

The linearized Euler energy identity then gives, for the **full PDE** and
every Fourier sector simultaneously,
\[
                         \|e^{tL_U}\|_{L^2\to L^2}
                         \le e^{\lambda_*|t|}.                   \tag{1.5}
\]

For a normalized curved torus whose metric and base coefficients are
\(O(\varepsilon)\)-close in \(C^1\), the same argument gives
\[
                         \|e^{tL_{U_\varepsilon}}\|_{2\to2}
                         \le e^{(\lambda_*+C\varepsilon)|t|}.    \tag{1.6}
\]

This supplies the semigroup input in C80 with \(d=0\).  Since
\(\varepsilon_jT_j\to0\) for geometric aspect ratio and polynomial gain
time, the \(O(\varepsilon)\) edge shift does not restore an exponential
Duhamel loss.

## 2. Exact strain formula

Use the orthonormal cylindrical frame.  For
\(U=V e_\theta+W e_z\),
\[
 \nabla U=
 \begin{pmatrix}
 0&-V/r&0\\
 V'&0&0\\
 W'&0&0
 \end{pmatrix}.
\]
Hence
\[
 \operatorname{sym}\nabla U=
 {1\over2}
 \begin{pmatrix}
 0&V'-V/r&W'\\
 V'-V/r&0&0\\
 W'&0&0
 \end{pmatrix}.                                                 \tag{2.1}
\]

Its eigenvalues are \(0,\pm{\cal S}\), where
\[
 {\cal S}={1\over2}\sqrt{(V'-V/r)^2+(W')^2}.                    \tag{2.2}
\]

Writing
\[
 t=\log r,\qquad h={\Omega_t\over\Omega},
\]
gives
\[
 V'-V/r=r\Omega'=\Omega h,\qquad
 W'={\Omega(1+h)\over\sqrt2},
\]
and therefore
\[
 {\cal S}^2
 ={\Omega^2\over8}(3h^2+2h+1).                                 \tag{2.3}
\]

At the selected ring,
\[
 3h_*^2+10h_*+1=0,\qquad
 \lambda_*^2=-h_*\Omega_0^2.
\]
Equation (2.3) gives
\[
 {\cal S}(r_*)^2
 ={\Omega_0^2\over8}(3h_*^2+2h_*+1)
 =-h_*\Omega_0^2
 =\lambda_*^2,                                                  \tag{2.4}
\]
which proves (1.2).

## 3. The extra inner strain maximum

For the uncut log-normal profile,
\[
 h=h_*+c_*t,\qquad
 c_*={8h_*^2\over3h_*+1}={4+28h_*\over9}<0,
\]
and
\[
 \Omega^2
 =\Omega_0^2\exp\left({h^2-h_*^2\over c_*}\right).
\]
Thus
\[
 {\cal S}^2(h)
 ={\Omega_0^2\over8}
  \exp\left({h^2-h_*^2\over c_*}\right)
  (3h^2+2h+1).                                                  \tag{3.1}
\]

The critical points of (3.1) are the zeros of
\[
 P(h)=3h^3+2h^2+(1+3c_*)h+c_* .                                \tag{3.2}
\]
Using \(3h_*^2+10h_*+1=0\), this factors exactly as
\[
 P(h)=(h-h_*)\left[
 3h^2+(2+3h_*)h+{4(h_*+1)\over3}
 \right].                                                       \tag{3.3}
\]

The quadratic factor has negative product.  Moreover,
\[
 Q(0)={4(h_*+1)\over3}<0,\qquad
 Q(-1)={7-5h_*\over3}>0.                                       \tag{3.4}
\]
Hence its negative zero lies in \((-1,0)\), its other zero is positive,
and \(h_*\) is the only critical point on \((-\infty,h_*]\).
Because (3.1) tends to zero as \(h\to-\infty\), the ring is the unique
maximum on the complete outer flank.

On \([h_*,0]\), the only additional critical point is the negative root in
\((-1,0)\), which is a minimum.  Also
\[
 {\cal S}^2(0)
 ={\Omega_0^2\over8}
   \exp\left(-{h_*^2\over c_*}\right)
 <{\Omega_0^2\over2}
 <\lambda_*^2.                                                  \tag{3.5}
\]
Here
\[
 -{h_*^2\over c_*}=-{3h_*+1\over8}
 ={4+\sqrt{22}\over8}<{9\over8},
\]
and the rational Taylor bound \(e^{9/8}<4\) proves the first strict
inequality.  The second follows from \(-h_*>3\).

The remaining positive critical point is the unwanted inner strain maximum.
It occurs after \(h=0\), where the BAS exponent is nonpositive.  We may
therefore replace the profile on that side without touching any spectral
jet used at \(r_*\).

## 4. Strain ellipse and hollow inner extension

Put
\[
 q=\Omega_t.
\]
For the locked-pitch family, (1.1) becomes the norm
\[
 {\cal N}(\Omega,q)
 ={1\over2}\sqrt{q^2+{1\over2}(\Omega+q)^2}.                    \tag{4.1}
\]

The set
\[
 {\cal E}_{\lambda_*}
 =\{(\Omega,q):{\cal N}(\Omega,q)<\lambda_*\}                   \tag{4.2}
\]
is an open convex ellipse containing the origin.  At the join point \(h=0\),
the log-normal jet has
\[
 q=0,\qquad {\cal N}(\Omega,0)<\lambda_*                       \tag{4.3}
\]
by (3.5).  A small interval immediately inside the join also lies strictly
in the ellipse and has \(q/\Omega=h>0\).

Choose a long logarithmic interval to the left.  A smooth monotone path from
\((0,0)\) to the join jet can be made with
\[
 \Omega\le0,\qquad q=\Omega_t\le0,\qquad
 (\Omega,q)\in{\cal E}_{\lambda_*}.                             \tag{4.4}
\]
Indeed, away from a short endpoint-matching layer, take a fixed flat step
of height \(|\Omega(0)|\) over length \(L\).  Its derivative is \(O(L^{-1})\),
so (4.1) converges uniformly to \(|\Omega|/(2\sqrt2)\), whose maximum has
the strict margin (3.5).  In the endpoint layer, follow the original jet
inside the same open ellipse.  Standard \(C^\infty\), or
Gevrey-\(\sigma\) with any \(\sigma>1\), interpolation can be chosen
arbitrarily \(C^1\)-close to this piecewise path and therefore remains in
the ellipse.

Set the extension identically zero to the left of the flat endpoint.
It is a hollow compact inner profile.  Since \(\Omega<0\) and \(q\le0\) on
the ramp,
\[
                             h={q\over\Omega}\ge0,               \tag{4.5}
\]
so the ramp introduces no positive locked-pitch BAS exponent.

The outer log-normal flank already satisfies
\({\cal S}<\lambda_*\) away from the ring and tends to zero.  A far-tail
Gevrey cutoff can be inserted with arbitrarily small \(C^1\) strain cost.
This proves the existence of the strain-capped compactification.

Because the modification is a common locked-pitch pressure modulation and
is unchanged near the ring, all identities and strict local geometry in C68
survive.  The isolated phase level also survives directly.  Write
\[
 z(r)={\beta_*r\over\sqrt2},\qquad \Lambda=\Omega(z-1).
\]
On the part of the new inner ramp where \(\Lambda<0\),
\(1<z\le z_{\rm join}\), while monotonicity gives
\(|\Omega|\le|\Omega_{\rm join}|\).  Hence
\[
 |\Lambda_{\rm new}|
 \le|\Omega_{\rm join}|(z_{\rm join}-1)
 =|\Lambda_{\rm old}(r_{\rm join})|
 <|\Lambda_*|.                                                  \tag{4.6}
\]
Farther inward \(z\le1\), so \(\Lambda_{\rm new}\ge0\).  The outer cutoff
may be taken as the original factor in \([0,1]\).  Thus no new copy of the
negative minimum is introduced.

## 5. Full linearized Euler semigroup

Let \(u\) be a smooth divergence-free solution of the Euler equation
linearized about the steady base \(U\):
\[
 \partial_tu+U\cdot\nabla u+u\cdot\nabla U+\nabla p=0.
\]
Transport and pressure are skew in \(L^2\), so
\[
 {1\over2}{d\over dt}\|u(t)\|_2^2
 =-\int u^T(\operatorname{sym}\nabla U)u\,dx
 \le\lambda_*\|u(t)\|_2^2.                                    \tag{5.1}
\]
Gronwall gives (1.5).  No Fourier decomposition, pressure reconstruction,
or high-frequency parametrix enters this estimate.

On a curved tube, use the Riemannian volume and covariant symmetric
gradient.  Under the normalized \(C^1\) coefficient comparison,
\[
 \|\operatorname{sym}_{g_\varepsilon}\nabla U_\varepsilon
   -\operatorname{sym}_{g_0}\nabla U_0\|_{L^\infty}
 \le C\varepsilon.                                              \tag{5.2}
\]
The covariant energy identity then gives (1.6).

For the edge-tracking carrier,
\[
 (\lambda_*+C\varepsilon_j-y_j)T_j
 =\delta_jT_j+C\varepsilon_jT_j=O(1)+o(1),                     \tag{5.3}
\]
so C80's Duhamel factor stays bounded.

## 6. Remaining caveats

1. The open-ellipse gluing is an existence construction.  A final theorem
   should give one explicit Gevrey modulation and record its derivative
   constants.
2. The actual finite-curvature Gavrilov base still has to be put into the
   common normalized coordinates proving (5.2).  Analytic rescaling gives
   the expected \(O(\varepsilon)\) rate, but the Taylor jets are not yet
   written in this repository.
3. The \(L^2\) semigroup bound closes the spectral Duhamel ledger.  It does
   not supply the all-order nonlinear endpoint inverse or control the global
   viscous wake.
