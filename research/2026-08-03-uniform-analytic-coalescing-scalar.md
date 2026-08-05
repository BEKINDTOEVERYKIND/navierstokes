# Uniform analytic coalescing-edge scalar pseudomode

**Date:** 2026-08-03

**Status:** exact fixed-chart operator and uniform one-dimensional analytic
WKB construction, self-verified.  The nonlinear transition and the
Navier--Stokes singularity problem remain open.

**Scope:** this note closes the local scalar obligation left open in C80,
C84, and C97 for the analytic core of the strain-capped locked-pitch
profile.  It uses the exact q-free lift C93 and the aspect-uniform Piola
transport C97 after scalar localization.  It does not claim an eigenvalue,
a nonlinear instability theorem, or a Clay solution.

Sections 2--4 are a local scalar theorem under the displayed fixed analytic
core.  The ambient curved ledger in Section 5 additionally uses the
normalized curved-base comparison and semigroup hypotheses stated in C82
and C97.  Those inputs are self-derived elsewhere in the repository but are
not re-proved or independently audited here.

## 1. Verdict

The coalescing scaling is valid, and its constants can be made uniform.
The point that was missing from the earlier power ledger is that the
strain-capped profile is **identically equal to the analytic log-normal
profile on a fixed neighborhood of the spectral ring**.  The Gevrey
cutoffs occur outside that neighborhood.  All coefficients used by the
packet therefore extend holomorphically to one fixed complex radial disc.

Let

\[
 s=p^{-1/2},\qquad y=\lambda_*-\delta,\qquad
 r_c=r_*+\sigma s,\qquad \sigma\ne0,
\]

and assume

\[
             \delta\longrightarrow0,\qquad
             h:=\frac{s}{\eta^3}\longrightarrow0,             \tag{1.1}
\]

where the **exact** positive number \(\eta\) is defined in (2.5) below.
Then there is a compactly supported scalar packet \(H_{s,\delta}\) in the
fixed analytic ring neighborhood such that

\[
 \frac{\|P_sH_{s,\delta}\|_2}{\|H_{s,\delta}\|_2}
 \le C\eta^2\exp\!\left(-\frac c h\right)
 =C\eta^2\exp\!\left(-c\frac{\eta^3}{s}\right).              \tag{1.2}
\]

Shrinking \(c\) absorbs every fixed polynomial reconstruction or Sobolev
factor.  The constants are independent of \(s,\delta\).

For

\[
 p_j=j^A,\qquad \delta_j\asymp j^{-g},\qquad T_j\asymp j^g,
\]

one has

\[
 \eta_j\asymp j^{-g/2},\qquad
 h_j\asymp j^{-(A-3g)/2}.                                    \tag{1.3}
\]

Thus \(A>3g\) is the exact scalar WKB gate.  At \(g=2,A=8\),
the analytic residual is \(e^{-cj}\), not \(e^{-cj^2}\).  A Gaussian
\(e^{-cj^2}\) target requires \(A\ge10\).  The latter is unnecessary for
summability: every strict \(A>6\) gives a summable stretched exponential.

There is a separate endpoint-flatness gate.  If this residual itself is
retained as a terminal force and must beat derivative losses
\(e^{C_Nj}\) for every fixed order \(N\), its action power must be strictly
larger than one:

\[
                 \frac{A-3g}{2}>1,\qquad A>3g+2.              \tag{1.4}
\]

For \(g=2\), long-gain fidelity needs \(A>6\), all-order terminal flatness
by this direct domination needs \(A>8\), and the convenient
\(e^{-cj^2}\) target needs \(A\ge10\).  These are three different
requirements.

## 2. Exact scalar operator before and after dilation

Use

\[
 D(r)=1+\beta^2r^2,\qquad
 A(r)=\frac{r^2}{D(r)},\qquad
 B(r)=r\frac d{dr}\left(\frac r{D(r)}\right),                \tag{2.1}
\]

and the q-free coefficients

\[
 a=r\frac d{dr}
 \left(\frac{\beta r^2W'-\Gamma'}{rD}\right),\qquad
 b=-\frac{2\beta V(W'+\beta\Gamma')}{D}.                    \tag{2.2}
\]

For \(n=p=s^{-2}\), set

\[
 \omega=p\Lambda(r_c)+iy,\qquad
 \gamma_s(Y)=s^{-2}\{\Lambda(r_c+sY)-\Lambda(r_c)\}-iy.     \tag{2.3}
\]

Multiplying the exact scalar residual by \(-r/n^2\) gives

\[
\begin{aligned}
 P_s={}&-s^2A(r_c+sY)\partial_Y^2
        -s^3B(r_c+sY)\partial_Y\\
 &+1+\frac{b(r_c+sY)}{\gamma_s(Y)^2}
       +s^2\frac{a(r_c+sY)}{\gamma_s(Y)}.                    \tag{2.4}
\end{aligned}
\]

The correct exact centering uses the principal potential only:

\[
 \boxed{\quad
 \eta^2=\frac{b(r_c)/y^2-1}{A(r_c)}.
 \quad}                                                       \tag{2.5}
\]

Since \(b(r_*)=\lambda_*^2\), \(b'(r_*)=0\), and
\(A(r_*)=A_*>0\),

\[
 \eta^2=\frac{2\delta}{A_*\lambda_*}
       +O(\delta^2+s\delta+s^2).                              \tag{2.6}
\]

Condition (1.1) implies \(s=o(\delta^{3/2})\), so (2.6) gives
\(\eta\asymp\delta^{1/2}\) and positivity for all small parameters.

Use the unitary dilation

\[
 (T_\eta v)(Y)=\eta^{-1}v(X),\qquad Y=\eta^2X,
\]

and put \(Q_{s,\delta}=\eta^{-2}T_\eta^{-1}P_sT_\eta\).  With

\[
 r(X)=r_c+s\eta^2X,
 \quad
 \gamma(X)=s^{-2}\{\Lambda(r(X))-\Lambda(r_c)\}-iy,         \tag{2.7}
\]

the exact fixed-chart operator is

\[
\boxed{
\begin{aligned}
 Q_{s,\delta}={}&-h^2A(r(X))\partial_X^2
  -h^3\eta^5B(r(X))\partial_X\\
 &+\eta^{-2}\left(1+\frac{b(r(X))}{\gamma(X)^2}\right)
  +h^2\eta^4\frac{a(r(X))}{\gamma(X)}.
                                                               \tag{2.8}
\end{aligned}}
\]

This also resolves a small ambiguity in the earlier note.  At \(X=0\),
the last term in (2.8) is generally imaginary, so the *complete*
differential expression does not vanish at \((0,\pm1)\).  It is an
\(O(h^2\eta^4)\) lower term.  Equation (2.5) makes
\((0,\pm1)\) an exact characteristic of the principal analytic symbol,
which is the quantity entering the eikonal equation.

## 3. Uniform analytic coefficient domain

Fix a complex radial disc

\[
                       |r-r_*|<\rho                           \tag{3.1}
\]

inside the unchanged log-normal core.  On this disc
\(A,B,a,b,\Lambda\) are holomorphic, real on the real axis, and have fixed
Cauchy bounds.  The zeros of \(D\) lie outside a smaller such disc.

For every fixed \(R\), (2.7) maps \(|X|<R\) into (3.1) when the parameters
are small.  Taylor expansion gives, uniformly on that complex \(X\)-disc,

\[
 \gamma(X)=-iy
 +\eta^2\frac{\Lambda'(r_c)}sX
 +\frac12\eta^4\Lambda''(r_c)X^2
 +O_{\mathcal A}(s\eta^6),                                  \tag{3.2}
\]

where

\[
             \frac{\Lambda'(r_c)}s=\kappa\sigma+O(s).        \tag{3.3}
\]

Consequently \(|\gamma(X)|\ge\lambda_*/2\) on \(|X|<R\) for
all sufficiently small parameters.  Every reciprocal in (2.8) has a
uniform analytic seminorm there.

Define

\[
 W_0(X)=\eta^{-2}
 \left(1+\frac{b(r(X))}{\gamma(X)^2}\right).                 \tag{3.4}
\]

The exact constant and linear Taylor coefficients are

\[
 W_0(0)=-A(r_c),                                              \tag{3.5}
\]

\[
 W_0'(0)=d_s+ic_s,\qquad
 d_s=-\frac{s b'(r_c)}{y^2},\qquad
 c_s=\frac{2b(r_c)\Lambda'(r_c)}{s y^3}.                     \tag{3.6}
\]

Because \(b'(r_*)=\Lambda'(r_*)=0\),

\[
 d_s=O(s^2),\qquad
 c_s\longrightarrow c_*=\frac{2b_*\kappa\sigma}{\lambda_*^3}\ne0. \tag{3.7}
\]

Moreover,

\[
 W_0(X)=-A(r_c)+(d_s+ic_s)X
                   +\eta^2X^2R_{s,\delta}(X),                \tag{3.8}
\]

where \(R_{s,\delta}\) is bounded in every analytic seminorm on a
slightly smaller fixed disc.  Also

\[
 A(r(X))=A(r_c)+O_{\mathcal A}(s\eta^2).                     \tag{3.9}
\]

Thus the principal symbols

\[
 q_{s,\delta}^{(0)}(X,\Xi)=A(r(X))\Xi^2+W_0(X)              \tag{3.10}
\]

converge in a uniform analytic-symbol topology to

\[
                   A_*(\Xi^2-1)+ic_*X.                       \tag{3.11}
\]

Exact centering removes the apparent \(s/\eta^2\) loss from the constant
coefficient.  It is not necessary to assume a uniform analytic estimate as
an additional hypothesis: it follows from (3.1)--(3.9).

## 4. Explicit complex phase and analytic WKB

Put

\[
 K(X)=-\frac{W_0(X)}{A(r(X))},\qquad K(0)=1.                 \tag{4.1}
\]

On a sufficiently small fixed complex disc, let \(K^{1/2}\) be the branch
equal to one at zero and choose

\[
 \xi_0=-\operatorname{sign}(c_s),\qquad
 \phi(X)=\xi_0\int_0^XK(t)^{1/2}\,dt.                       \tag{4.2}
\]

This is an explicit analytic solution of the eikonal equation

\[
                 A(r(X))\phi'(X)^2+W_0(X)=0.                \tag{4.3}
\]

Its Gaussian sign is uniform.  Indeed,

\[
 \operatorname{Im}\phi''(0)
 =-\xi_0\frac{c_s}{2A(r_c)}
 =\frac{|c_s|}{2A(r_c)}>0.                                   \tag{4.4}
\]

After decreasing a fixed \(R>0\),

\[
             \operatorname{Im}\phi(X)\ge c_0X^2
             \quad (X\in\mathbb R,\ |X|\le R)              \tag{4.5}
\]

with \(c_0>0\) independent of \(s,\delta\).  Hence the packet has
\(X\)-width \(h^{1/2}\), \(Y\)-width
\(\eta^2h^{1/2}=\sqrt{s\eta}\), and physical radial width
\(s\sqrt{s\eta}\).

For completeness, the analytic transport construction can be written
directly.  Conjugating (2.8) by \(e^{i\phi/h}\), and using (4.3), gives

\[
 e^{-i\phi/h}Q_{s,\delta}e^{i\phi/h}
 =-ih\mathcal L+h^2\mathcal M_0+h^3\mathcal M_1,             \tag{4.6}
\]

where

\[
\begin{aligned}
 \mathcal L&=A(r(X))\{2\phi'\partial_X+\phi''\},\\
 \mathcal M_0&=-A(r(X))\partial_X^2
 +\eta^4\frac{a(r(X))}{\gamma(X)}
 -i\eta^5B(r(X))\phi',\\
 \mathcal M_1&=-\eta^5B(r(X))\partial_X.                    \tag{4.7}
\end{aligned}
\]

Let \(a_{-1}=0\), take \(a_0=(\phi')^{-1/2}\), and for \(k\ge1\)
solve

\[
 \mathcal La_k=-i\{\mathcal M_0a_{k-1}
                         +\mathcal M_1a_{k-2}\},\qquad a_k(0)=0. \tag{4.8}
\]

The inverse used here is not abstract.  If \(f\) is analytic, the solution
of \(\mathcal La=f\), \(a(0)=0\), is

\[
 a(X)=(\phi'(X))^{-1/2}
 \int_0^X
 \frac{f(t)}
 {2A(r(t))(\phi'(t))^{1/2}}\,dt.                            \tag{4.8a}
\]

Indeed,

\[
 2A(\phi')^{1/2}\partial_X\{(\phi')^{1/2}a\}
 =A(2\phi'a'+\phi''a)=\mathcal La.
\]

All coefficients in (4.7), \((A\phi')^{-1}\), and their analytic norms
are uniformly bounded.  On nested discs
\(D_{\rho'}\Subset D_\rho\), Cauchy's inequality gives
\[
 \|\partial_X^mf\|_{\rho'}
 \le m!(\rho-\rho')^{-m}\|f\|_\rho.
\]
Equation (4.8a) recovers one of the two derivatives lost through
\(\mathcal M_0\).  Induction in (4.8), with a linearly decreasing disc
radius, therefore gives the uniform analytic majorant

\[
                       \|a_k\|\le C_1^{k+1}k!.                \tag{4.9}
\]

Truncate \(\sum h^ka_k\) at \(N=\lfloor c_1/h\rfloor\).
Stirling's bound and (4.9) make the interior error \(O(e^{-c_2/h})\).
Choose a fixed cutoff supported in \((-R,R)\), equal to one on
\((-R/3,R/3)\), with derivatives supported where \(|X|\ge R/3\).
By (4.5), its commutator with \(Q_{s,\delta}\) gains
\(e^{-c_0R^2/(9h)}\).  Taking the truncation constant \(c_1\) small enough
makes this Gaussian dominate the uniform exponential bound on the
truncated amplitude.  The packet norm is bounded below by a constant times
\(h^{1/4}\); changing the final exponent absorbs this polynomial
normalization.  This proves

\[
 \frac{\|Q_{s,\delta}v_{s,\delta}\|_2}
      {\|v_{s,\delta}\|_2}
 \le C e^{-c/h},                                              \tag{4.10}
\]

and unitarity of \(T_\eta\) proves (1.2).

This direct one-dimensional construction is the uniform version of the
analytic pseudomode mechanism; no compactness assertion about unspecified
seminorms is being used.

### 4.1 The frequency envelope required by C97

The physical radial coordinate is

\[
                 r-r_c=s\eta^2X.
\]

The three possible radial derivative scales are all bounded by the carrier
frequency \(p=s^{-2}\).  The oscillatory phase has covector

\[
 \frac1{s\eta^2h}=\frac{\eta}{s^2}=p\eta.                    \tag{4.11}
\]

The fixed-\(X\) cutoff has covector

\[
 \frac1{s\eta^2}=p\frac{s}{\eta^2}=p\,h\eta,                 \tag{4.12}
\]

and the Gaussian envelope, whose physical width is
\(s\sqrt{s\eta}\), has covector

\[
 \frac1{s\sqrt{s\eta}}
 =p\sqrt{\frac{s}{\eta}}
 =p\sqrt{h\eta^2}.                                            \tag{4.13}
\]

All three ratios to \(p\) are uniformly bounded and tend to zero under
(1.1).  The angular and axial Fourier covectors are already \(O(p)\).
Consequently, for \(0\le k\le2\),

\[
                 \|\partial_r^kH\|_2\le C_kp^k\|H\|_2.       \tag{4.14}
\]

Since \(\gamma\) is bounded above and below on the packet, the exact C93
formulas then give

\[
 \|\nabla u\|_2+\|\nabla P\|_2
 \le Cp\|u\|_2+C\|F\|_2.                                     \tag{4.15}
\]

Thus the frequency-envelope hypothesis in C97 is a consequence of the
constructed scalar packet, not an additional assumption.

## 5. Exact velocity lift, curvature, viscosity, and gain time

Let \(G\) be the scalar residual in C93.  Comparing its definition with
(2.4) gives the exact identity

\[
                         G=-\frac{n^2}{r}P_sH.                \tag{5.1}
\]

The q-free velocity--pressure lift therefore has

\[
 F_\theta=F_z=0,\qquad
 F_r=\frac{\gamma}{in^2}G
             =-\frac{\gamma}{ir}P_sH.                        \tag{5.2}
\]

There is no local or global Hodge correction.  Scalar localization is
followed by the exact C93 reconstruction.  For the ambient step, assume the
normalized curved-base \(W^{1,\infty}\) comparison in C97 (4.1), the packet
frequency envelope in C97 (5.1)/(5.4), and the curved semigroup edge in C82
(5.2).  The frequency envelope is verified directly below for the present
packet; the base comparison and semigroup hypotheses remain the prior
self-derived inputs.  Under those hypotheses, the divergence-preserving
Piola map of C97 supplies the aspect-uniform curved residual

\[
 r_j\le \operatorname{poly}(p_j,\eta_j^{-1})e^{-c/h_j}
       +C\varepsilon_jp_j
       +C\mu_jp_j^2
       +C\varepsilon_j\mu_jp_j^2,                             \tag{5.3}
\]

where \(\mu_j\) is the viscosity normalized at stage \(j\).  The third
term is the straight viscous defect of an Euler pseudomode; the fourth is
the curved Laplacian commutator.  The corrected axial winding is

\[
 m_j=\operatorname{round}(\beta_*p_j/\varepsilon_j),
\]

or, more cleanly, the scalar construction is carried out at the exact
periodic pitch \(\beta_j=\varepsilon_jm_j/p_j\).  The pitch error is below
the curvature budget.

Under the cited C82 comparison, the strain estimate gives a semigroup edge
\(\lambda_*+O(\varepsilon_j)\).  Since
\(\delta_jT_j=O(1)\), Duhamel closes under the transparent sufficient
conditions

\[
\boxed{
 A>3g,\qquad
 \varepsilon_jj^{A+g}\longrightarrow0,\qquad
 \mu_jj^{2A+g}\longrightarrow0.}                             \tag{5.4}
\]

For stage-by-stage absolute summability, replace the last two limits by

\[
 \sum_j\varepsilon_jj^{A+g}<\infty,
 \qquad
 \sum_j\mu_jj^{2A+g}<\infty.                                \tag{5.5}
\]

Geometric aspect and Reynolds schedules satisfy both for every fixed
polynomial \(p_j\).  At the convenient \(g=2,A=8\) point, the two
algebraic requirements are

\[
                   \varepsilon_jj^{10}\to0,
                   \qquad \mu_jj^{18}\to0.                   \tag{5.6}
\]

At \(A=10\), which upgrades the analytic residual to \(e^{-cj^2}\), they
become \(\varepsilon_jj^{12}\to0\) and \(\mu_jj^{22}\to0\).

The scale schedule already used by the constructive cascade is a concrete
specialization.  Let

\[
 \ell_j=\ell_0q^j\quad(0<q<1),\qquad
 \varepsilon_j=\ell_j^\beta,\qquad
 \mu_j=\nu\ell_j^{\gamma_{\rm sc}-1}p_j^{-\gamma_{\rm sc}},
                                                                    \tag{5.7}
\]

with \(1<\gamma_{\rm sc}<3/2\).  Then the three gain-integrated algebraic
errors are

\[
\begin{aligned}
 T_j\varepsilon_jp_j
   &\asymp \ell_j^\beta j^{A+g},\\
 T_j\mu_jp_j^2
   &\asymp \nu\ell_j^{\gamma_{\rm sc}-1}
               j^{A(2-\gamma_{\rm sc})+g},\\
 T_j\varepsilon_j\mu_jp_j^2
   &\asymp \nu\ell_j^{\beta+\gamma_{\rm sc}-1}
               j^{A(2-\gamma_{\rm sc})+g}.                  \tag{5.8}
\end{aligned}
\]

Every line is summable for \(\beta>0\) and
\(\gamma_{\rm sc}>1\), since a positive geometric power beats every
polynomial.  The independent two-viscosity/third-jet wake gate imposes the
stronger upper restriction.  The common nonempty window is

\[
                    0<\beta<\frac{\gamma_{\rm sc}-1}{2}.      \tag{5.9}
\]

Thus the spectral packet does not shrink the existing geometric wake
window.

## 6. Corrections to the earlier checkpoint

This calculation resolves four inconsistencies or stale statements in the
edge-tracking note.

1. Its Section 3.1 correctly called the uniform analytic seminorm estimate
   open, while Section 8 prematurely called the internal scalar WKB chain
   closed.  Sections 2--4 above now supply that missing estimate and
   construction.  The **full nonlinear carrier chain is still not closed**.
2. The exact root is (2.5), for the principal potential.  The
   \(s^2a/\gamma\) value at the center is a transport-level \(O(h^2)\)
   coefficient, not part of the eikonal root.
3. \(g=2,A=8\) gives \(e^{-cj}\).  Long-gain summability only needs
   \(A>6\); direct all-order terminal flatness against \(e^{C_Nj}\) needs
   \(A>8\); and the stronger \(e^{-cj^2}\) rate needs \(A\ge10\).
4. The phrase "local Hodge correction preserves the phase" is obsolete.
   C93 gives the exact straight divergence-free velocity--pressure lift,
   and C97 transports it to the ambient curved tube by Piola without a
   changing Hodge projector.

## 7. Remaining boundary

The local scalar theorem is no longer the spectral bottleneck.  The ambient
claim still inherits C82/C97's normalized curved-base hypotheses and needs
independent audit.  The remaining prize-level obligations are:

* the amplified real packet must enter a material-phase endpoint chart;
* charged interactions and the global pressure/viscous wake must be carried
  through all orders;
* the stage-scaled viscosity schedule in (5.4) must be realized by one
  physical smooth finite-energy solution, not assigned independently; and
* the exact endpoint inverse and infinite return map remain open.

The calculation above supplies a genuine exponentially accurate local
carrier.  It does not by itself imply finite-time blow-up.

## 8. Reproducibility

The [companion dependency-free checker](../checks/uniform_analytic_coalescing_scalar.py)
verifies with exact rational and Gaussian-rational arithmetic:

* the exact centered constant and linear potential coefficients;
* the nonzero limiting bracket and Gaussian phase sign;
* every dilation power in (2.8);
* the \(A>3g\) action threshold and the \(A=8\) versus \(A=10\) rates; and
* the combined curvature/viscosity/gain exponents in (5.4) and the
  geometric cascade window (5.7)--(5.9).

It checks the algebra and exponent ledger, not the analytic majorant proof
in (4.8)--(4.10), the actual profile's fixed complex core, the cutoff
commutator and frequency-envelope estimates, or the C82/C97 ambient
comparison hypotheses.
