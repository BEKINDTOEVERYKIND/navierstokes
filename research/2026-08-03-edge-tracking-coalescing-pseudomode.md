# Edge-tracking coalescing pseudomode and material-phase ledger

**Date:** 2026-08-03

**Status:** exact rescaling, uniform local analytic scalar WKB construction,
and exponent ledger self-derived; nonlinear endpoint realization remains
open.

**Scope:** the locked-pitch Gavrilov carrier of C68 and the forced scalar
family restored in `2026-08-03-recovered-subedge-curved-block.md`.  This note
does not claim an Euler eigenvalue, a Navier--Stokes singularity, or closure
of the transition construction.

## 1. Outcome

A fixed subedge pseudomode has an exponentially dangerous gap from the full
BAS edge during a gain interval of length \(T_j\asymp j^2\).  That makes the
geometric curvature residual \(\varepsilon_jp_j\) unusable even though it is
small at each instant.

The gap can be removed without constructing an exact curved pseudomode.
Let the imaginary spectral parameter approach the BAS edge at the inverse
gain scale.  The radial characteristic roots then coalesce, but their exact
normal form has effective semiclassical parameter

\[
                       \hbar_{\rm eff}={s\over\eta_j^3}.          \tag{1.1}
\]

For logarithmic gain \(G_j\asymp j^g\), carrier frequency
\(p_j=j^A\), and \(s=p_j^{-1/2}\), both the analytic WKB residual and the
material-phase coherence condition close precisely when

\[
                              A>3g.                              \tag{1.2}
\]

For the cascade target \(g=2\), this is

\[
                              A>6.                               \tag{1.3}
\]

This is stronger than the earlier Gevrey/edge-deficit condition \(A>4\),
but it removes the need for an exact finite-curvature eigenmode and replaces
the missing curved two-microlocal theorem by a one-dimensional coalescing
normal form.

## 2. Edge tracking

Keep the fixed nonzero off-ring parameter \(\sigma\) in

\[
 r_{c,j}=r_*+\sigma s_j,
 \qquad s_j=p_j^{-1/2}.
\]

Let \(\lambda_*=\sqrt{b_*}\) be the full BAS edge and choose

\[
 y_j=\lambda_*-\delta_j,
 \qquad \delta_j=c_\delta j^{-g},
 \qquad
 \omega_j=p_j\Lambda(r_{c,j})+iy_j.                            \tag{2.1}
\]

At \(Y=0\), the scalar characteristic equation is

\[
 A_*\eta_j^2+1-{b_*\over y_j^2}=0,
\]

so

\[
 \eta_j^2={b_*/y_j^2-1\over A_*}
          ={2\delta_j\over A_*\lambda_*}+O(\delta_j^2),
 \qquad
 \eta_j\asymp j^{-g/2}.                                       \tag{2.2}
\]

The bracket from the fixed-subedge calculation is now nonzero but tends to
zero like \(\eta_j\).  A fixed-bracket invocation of the analytic
pseudomode theorem is therefore not uniform.  The coalescing rescaling below
is necessary.

## 3. Exact coalescing-root rescaling

Write the scalar symbol near \((Y,\eta)=(0,\eta_j)\) as

\[
 p_j(Y,\eta)=A_*\eta^2+V_j(Y)+\text{lower terms},                \tag{3.1}
\]

where

\[
 V_j(0)=-A_*\eta_j^2,
 \qquad
 V_j'(0)=i c_*+O(\delta_j+s_j),
 \qquad
 c_*={2b_*\kappa\sigma\over\lambda_*^3}\ne0.                 \tag{3.2}
\]

The two radial roots are separated by \(2\eta_j\).  Since the potential
changes by their squared separation over a \(Y\)-distance
\(O(\eta_j^2)\), set

\[
                         Y=\eta_j^2X,
 \qquad                    \eta=\eta_j\Xi.                      \tag{3.3}
\]

After division by \(\eta_j^2\), the differential operator has principal
part

\[
 -\left({s_j\over\eta_j^3}\right)^2A_*\partial_X^2
 -A_*+ic_*X+O(\eta_j^2+s_j/\eta_j^2).                          \tag{3.4}
\]

Thus (1.1) is exact.  If \(s_j/\eta_j^3\to0\), the limiting symbol is

\[
                  A_*(\Xi^2-1)+ic_*X,                          \tag{3.5}
\]

whose two characteristics \((X,\Xi)=(0,\pm1)\) have a fixed nonzero
pseudospectral bracket.

The apparently dangerous lower coefficient in (3.4) is harmless:

\[
 {s_j\over\eta_j^2}
 =\eta_j\,{s_j\over\eta_j^3}=o(1).                             \tag{3.6}
\]

The first-derivative term \(s_j^3B\partial_Y\) is smaller still after the
same rescaling.

The analytic WKB recursion therefore uses \(\hbar_{\rm eff}\), not \(s_j\),
as its expansion parameter.  Truncation at
\(N\asymp\hbar_{\rm eff}^{-1}\) gives the quantitative target

\[
 {\|P_{s_j}\phi_j\|_2\over\|\phi_j\|_2}
 \le \operatorname{poly}(s_j^{-1},\eta_j^{-1})
       \exp\left(-c{\eta_j^3\over s_j}\right).                 \tag{3.7}
\]

This is also visible from the local Gaussian.  Its width in \(Y\) is

\[
                         w_Y\asymp\sqrt{s_j\eta_j},              \tag{3.8}
\]

whereas the root-separation neighborhood has radius \(O(\eta_j^2)\).
Their ratio is

\[
 {w_Y\over\eta_j^2}
 \asymp\left({s_j\over\eta_j^3}\right)^{1/2}.                  \tag{3.9}
\]

At the edge of that neighborhood the Gaussian action is
\(c\eta_j^3/s_j\), agreeing with (3.7).

For \(p_j=j^A\) and (2.2),

\[
 {\eta_j^3\over s_j}
 \asymp j^{(A-3g)/2}.                                          \tag{3.10}
\]

It diverges exactly under (1.2).

### 3.1 Uniform analytic pseudomode theorem after dilation

The rescaling does more than identify a target.  It puts the exact operator
inside the standard analytic pseudomode theorem with uniform constants.

Choose the small positive root from the exact **principal** value
\[
 \eta_j^2=
 \frac{b(r_{c,j})/y_j^2-1}{A(r_{c,j})}.
\]
The generally imaginary \(s_j^2a/\gamma\) value is an
\(O(\hbar_j^2\eta_j^4)\) transport coefficient, not part of this eikonal
root.  Use the unitary dilation
\[
 \phi(Y)=\eta_j^{-1}v(X),\qquad Y=\eta_j^2X.                    \tag{3.11}
\]
Then
\[
 Q_j=\eta_j^{-2}T_j^{-1}P_{s_j}T_j
\]
is an analytic semiclassical differential operator on a fixed complex
\(X\)-disc with parameter
\[
                         \hbar_j={s_j\over\eta_j^3}.             \tag{3.12}
\]
The exact fixed-chart calculation in
2026-08-03-uniform-analytic-coalescing-scalar.md proves that its
principal symbol satisfies
\[
 q_j(X,\Xi)
 =A_*(\Xi^2-1)+ic_*X
 +O_{\mathcal A}(\eta_j^2+s_j),                               \tag{3.13}
\]
where \(O_{\mathcal A}\) denotes a uniform analytic-symbol norm on that
fixed disc.  Exact principal centering removes the formerly apparent
\(s_j/\eta_j^2\) constant loss.  The parameter and remaining differential
coefficients have the exact sizes
\[
 s_j=\hbar_j\eta_j^3,\qquad
 {s_j^2\over\eta_j^2}=\hbar_j^2\eta_j^4,\qquad
 {s_j^3\over\eta_j^4}=\hbar_j^3\eta_j^5.                       \tag{3.14}
\]

The same note verifies the uniform analytic seminorms from the unchanged
log-normal complex radial core, writes the explicit phase
\[
 \phi_j(X)=\xi_j\int_0^X
 \sqrt{-W_{0,j}(t)/A_j(t)}\,dt,
\]
and proves \(\operatorname{Im}\phi_j(X)\ge cX^2\) on a fixed real
interval.  Uniform analytic transport and optimal truncation therefore
give
\[
 \|Q_jv_j\|_2\le Ce^{-c/\hbar_j}\|v_j\|_2.                     \tag{3.15}
\]
Returning through (3.11) proves (3.7), with constants uniform in \(j\).
The former conditional analytic-seminorm obligation is thereby closed at
SELF/analytic-proof status.

The reconstruction must be applied **after** localizing the scalar packet.
C93 gives the exact q-free velocity--pressure lift, with
\(F_\theta=F_z=0\), \(F_r=\gamma G/(in^2)\), and no Hodge correction.
C97 then transports the compact divergence-free field by Piola and proves
the aspect-uniform ambient residual estimate.  Fixed Sobolev reconstruction
powers are absorbed by \(e^{-c/\hbar_j}\).

On the normalized thin torus, the central arclength coordinate in
\[
 ds^2=dr^2+r^2d\theta^2+
 (1+\varepsilon_jr\cos\theta)^2dz^2
\]
has period \(2\pi/\varepsilon_j\). Its axial Fourier lattice is therefore
\(\varepsilon_j\mathbb Z\), not \(\mathbb Z\). Choose integers
\[
 p_j=\lceil j^A\rceil,\qquad
 m_j=\operatorname{round}\left({\beta_*p_j\over\varepsilon_j}\right).
\]
Then the physical helical ratio obeys
\[
 \left|{\varepsilon_jm_j\over p_j}-\beta_*\right|
 \le{\varepsilon_j\over2p_j}
 =O(\varepsilon_js_j^2).                                     \tag{3.16}
\]
Without recentering the trapping ring, the resulting analytic-symbol
perturbation is bounded by
\[
 O(\varepsilon_js_j)
 +O\left({\varepsilon_js_j^2\over\eta_j^2}\right)
 =
 O(\varepsilon_j\hbar_j\eta_j^3)
 +O(\varepsilon_j\hbar_j^2\eta_j^4).                         \tag{3.17}
\]
Both terms are below the retained budget when \(\varepsilon_j\) is
geometric. Re-centering at the nearby
\((\varepsilon_jm_j/p_j)\)-dependent phase ring removes the first term.
The winding integer is \(m_j\asymp p_j/\varepsilon_j\), while the physical
axial covector \(\varepsilon_jm_j\) remains \(O(p_j)\).

Finally, adjoin the conjugate Fourier packet.
For a nonzero integer charge, sine--cosine orthogonality makes the real and
imaginary packets a uniformly conditioned real two-dimensional block.
Taking the real part of the complex time-dependent pseudoorbit therefore
gives a real periodic pseudoorbit with comparable norm and the same residual
bound up to a fixed constant.

## 4. Duhamel ledger with no exponential edge loss

Assume the relevant high-frequency Euler semigroup satisfies

\[
 \|S_j(t)\|\le C(1+t)^d e^{\lambda_*t}                         \tag{4.1}
\]

after any finite-dimensional spectral part with faster exact growth has
been separated.  The latter case is beneficial and would supply a true
unstable mode instead of a pseudomode.

For a quasimode with growth \(y_j=\lambda_*-\delta_j\), residual \(r_j\),
and duration \(T_j\asymp G_j/y_j\asymp j^g\), Duhamel gives the relative
bound

\[
 E_j\lesssim r_j(1+T_j)^{d+1}e^{\delta_jT_j}.                   \tag{4.2}
\]

The edge-tracking choice (2.1) makes

\[
                         \delta_jT_j=O(1),                       \tag{4.3}
\]

so there is no \(e^{cj^g}\) penalty.

For the bent finite-curvature packet from the companion note,

\[
 r_j\lesssim
 \operatorname{poly}(j)e^{-c j^{(A-3g)/2}}
 +C\varepsilon_jj^A.                                           \tag{4.4}
\]

If \(\varepsilon_j\) is geometric (or merely faster than every required
inverse polynomial), (4.2) tends to zero and is summable after a harmless
strengthening of the geometric rate.  Thus the curvature comparison
residual is now compatible with long gain.

The later strain-capped compactification C82 supplies (4.1) with \(d=0\)
for a redesigned compact locked-pitch profile, directly from the full
\(L^2\) energy identity.  Its finite-curvature version has exponent
\(\lambda_*+C\varepsilon_j\), which changes (4.3) only by
\(C\varepsilon_jT_j=o(1)\).  That compactification and its curved
coefficient comparison remain self-verified pending cross-audit.

## 5. Material-phase coherence has the same threshold

The real spatial phase of the normal packet is transported by the straight
base flow if its temporal frequency is allowed to depend on radius.  The
static pseudomode instead uses \(p_j\Lambda(r_{c,j})\).  Across the packet,
the accumulated phase discrepancy after time \(T_j\) is

\[
 p_jT_j\,[\Lambda(r)-\Lambda(r_{c,j})].                         \tag{5.1}
\]

Because \(\Lambda'(r_{c,j})=\kappa\sigma s_j+O(s_j^2)\) and the physical
radial width is

\[
 w_r=s_jw_Y\asymp s_j\sqrt{s_j\eta_j},                         \tag{5.2}
\]

the linear variation dominates and

\[
 p_jT_j\sup_{\rm packet}|\Lambda-\Lambda(r_{c,j})|
 \lesssim T_j\sqrt{s_j\eta_j}
 \asymp j^{(3g-A)/4}.                                          \tag{5.3}
\]

This tends to zero exactly when \(A>3g\), the same condition as the WKB
action.  The radial wavevector chirp is weaker.  Hence edge tracking does
not force a second, larger carrier exponent to convert the amplified mode
into a transported material phase.

For \(g=2\), any \(A>6\) simultaneously gives:

\[
 \hbar_{\rm eff}\to0,\qquad
 e^{-c/\hbar_{\rm eff}}\to0,\qquad
 T_j\sqrt{s_j\eta_j}\to0.                                     \tag{5.4}
\]

## 6. Finite-charge compatibility

Treat the packet as a variable WKB phase rather than expanding its amplitude
into curved Fourier sidebands.  Three real phases have fundamental charge
set

\[
                         \{\pm\Phi_1,\pm\Phi_2,\pm\Phi_3\}.
\]

Their quadratic products produce only

\[
 0,\qquad \pm2\Phi_i,\qquad \pm\Phi_i\pm\Phi_k\quad(i<k),       \tag{6.1}
\]

at most nineteen distinct charges.  This count is independent of \(p_j\).
The exact C93 q-free lift preserves the Fourier/WKB charge without a Hodge
correction.  The C97 Piola map then preserves divergence.  Thus the
edge-tracking packet does not by itself create an infinite charge ladder at
quadratic order.

This does not solve the all-order multiphase endpoint problem: repeated
nonlinear corrections generate higher integer combinations, and every
cross-colour charge still has to be routed through the global pressure wake.
It does show that the localized spectral packet is compatible with the
finite quadratic charge chart used at the first transition step.

## 7. Design point and exact window

For the target \(g=2\), choose for example

\[
                       A=8,qquad \delta_j=c_\delta j^{-2}.
\]

Then

\[
 \eta_j\asymp j^{-1},\qquad
 s_j=j^{-4},\qquad
 \hbar_{\rm eff}\asymp j^{-1},\qquad
 w_Y\asymp j^{-5/2},                                            \tag{7.1}
\]

and the material-phase discrepancy is \(O(j^{-1/2})\).  The carrier
frequency \(p_j=j^8\) remains polynomial, while the WKB residual is
\(e^{-cj}\) up to polynomial factors.  Taking any \(A>8\) strengthens both
margins; the sharp ledger threshold is \(A>6\).

These thresholds serve different purposes.  Long-gain fidelity and
summability need \(A>6\).  If the residual itself is retained as a terminal
force and must beat every fixed-order loss \(e^{C_Nj}\), direct domination
needs \(A>8\).  A convenient \(e^{-cj^2}\) residual needs \(A\ge10\).

## 8. What is now load-bearing

The **local scalar** carrier component is now closed at
SELF/analytic-proof status by the exact fixed-chart construction above and
its companion checker.  C82 supplies the full \(L^2\) semigroup edge, C93
the exact straight velocity--pressure lift, and C97 the aspect-uniform
Piola/ambient comparison.  The broader carrier chain remains subject to
independent audit and, more importantly, to the nonlinear material-phase
endpoint and wake construction.

The constructive front is now the localized nonlinear endpoint problem:
either launch the one-phase bath on a bounded terminal pulse and carry its
charged plus zero-charge wakes, or fall back to the all-order three-phase
material-transition inverse.  Both routes must retain the global second-
and higher-jet wake.
