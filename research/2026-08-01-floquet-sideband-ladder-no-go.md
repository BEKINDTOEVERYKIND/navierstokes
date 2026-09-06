# Floquet sideband ladder: exact capture obstruction and the envelope escape

**Date:** 2026-08-01
**Status:** leading WKB/Bloch algebra; self-derived and not cross-audited
**Scope:** local straight-column symbol.  This note does not construct a
Navier--Stokes singularity or a nonlinear localized transition.

## 1. Decision

A real low child necessarily couples a resonant high-frequency carrier to
an infinite sideband ladder.  The ladder has a useful stable Bloch sector,
but it does **not** possess a nonzero finite-support or square-summable
packet which simultaneously

1. remains captured after the child turns on,
2. writes a nonzero low stress, and
3. keeps every positive-positive high sum purely pressure.

The obstruction is exact for the leading translation-invariant symbol.
It invalidates the finite central-doublet closure.  The only surviving
version is an approximate, spatially localized envelope with a damped
high-sum corrector and a full charge-zero child profile.  Compact smooth
localization then creates a separate core--buffer problem.

## 2. Resonant ladder

Use an orthonormal frame \((e_r,e_h,e_t)\), put

\[
 q=Qe_r,\qquad
 K_n=\Lambda e_t+\left(n+\frac12\right)q,
 \qquad \epsilon={Q\over\Lambda},\qquad n\in\mathbb Z,
 \tag{2.1}
\]

and define the exactly transverse vector

\[
 p_n=e_r-\left(n+\frac12\right)\epsilon e_t,
 \qquad K_n\cdot p_n=0.                              \tag{2.2}
\]

Write the positive carrier amplitudes as

\[
                 U_n=a_np_n+b_ne_h.                   \tag{2.3}
\]

At the central resonant covector, let the parent-column amplitude matrix
have unstable exponent \(g>0\) and unstable eigenvector
\(p+He_h\), with \(H<0\).  To leading order for
\(|n|\epsilon\ll1\), the parent and a real child
\(-2Z\sin(Qr)e_h\) give

\[
 \begin{aligned}
  \dot a_n&={g\over H}b_n-\kappa a_n,\\
  \dot b_n&=gH a_n+QZ(a_{n-1}+a_{n+1})-\kappa b_n,
 \end{aligned}                                        \tag{2.4}
\]

where \(\kappa=\nu\Lambda^2\) is the leading carrier heat rate.
The two shifts have the same sign.  Reality is precisely what supplies
both of them.

The exact heat rate is

\[
 \kappa_n=\nu\left(\Lambda^2+left(n+\frac12\right)^2Q^2\right).
 \tag{2.5}
\]

It adds a common discrete confining operator to the two polarizations and
does not restore finite support.

## 3. No finite-support invariant packet

Suppose \((a_n,b_n)\) has finite support and \(Zg\ne0\).  At an extreme
index \(N\), if \(a_N\ne0\), equation (2.4) immediately creates
\(b_{N+1}\).  If \(a_N=0\) but \(b_N\ne0\), the parent first creates
\(a_N\), after which the child creates \(b_{N+1}\).  Applying the same
argument at the other edge proves

\[
 \boxed{\text{the zero sequence is the only finite-support invariant
 packet of (2.4).}}                                    \tag{3.1}
\]

Thus every exact real-child closure contains an infinite sideband tail.

## 4. Bloch symbol and the unstable opposite sector

For

\[
 A(\theta)=\sum_na_ne^{-in\theta},\qquad
 B(\theta)=\sum_nb_ne^{-in\theta},                     \tag{4.1}
\]

the translation-invariant ladder is multiplication by

\[
 M(\theta)=
 \begin{pmatrix}
 -\kappa & g/H\\
 gH+2QZ\cos\theta&-\kappa
 \end{pmatrix}.                                        \tag{4.2}
\]

Its eigenvalues are

\[
 \boxed{\lambda_\pm(\theta)
 =-\kappa\pm
 \sqrt{g\left(g+{2QZ\over H}\cos\theta\right)}.}     \tag{4.3}
\]

Because \(H<0\), the child weakens the hyperbolic parent near
\(\theta=0\) but strengthens it near \(\theta=\pi\).  In particular,

\[
 \lambda_+(\pi)
 =-\kappa+\sqrt{g\left(g+{2QZ\over |H|}\right)}.
 \tag{4.4}
\]

Whenever the unperturbed carrier grows, \(g>\kappa\), (4.4) is positive.
No real child can stabilize the whole \(\ell^2\) ladder.

There is nevertheless a legitimate invariant subspace of the *linear*
ladder: take \(A,B\in L^2(\mathbb T)\) supported in a small arc around
\(\theta=0\).  Once

\[
        g\left(g-{2QZ\over|H|}\cos\theta\right)<\kappa^2
                                                               \tag{4.5}
\]

on that arc, the packet is spectrally quenched.  A compactly supported
smooth Bloch profile gives a Schwartz sideband sequence.  Its low
\(q\)-stress is generically nonzero; up to normalization, its
\(e_h\)-component is

\[
       2Q\operatorname{Re}\int_{\mathbb T}
            e^{i\theta}A\overline B\,d\theta.           \tag{4.6}
\]

This establishes capture only for the linear ladder, not closure of the
quadratic Navier--Stokes system.

## 5. Exact high-sum pressure obstruction

The positive-positive interaction of \(K_n\) and \(K_m\) lands at

\[
             K_n+K_m=2\Lambda e_t+(n+m+1)q.             \tag{5.1}
\]

Since

\[
       U_n\cdot K_m=(m-n)Qa_n,                          \tag{5.2}
\]

the \(e_h\)-component of the force at
\(2\Lambda e_t+jq\) is, apart from the common factor \(iQ\),

\[
       C_j=\sum_{n+m=j-1}(m-n)a_nb_m.                   \tag{5.3}
\]

With the generating functions
\(A(z)=\sum_na_nz^n\), \(B(z)=\sum_nb_nz^n\), all these
coefficients are encoded by

\[
          \sum_j C_jz^{j-1}=A(z)zB'(z)-zA'(z)B(z).      \tag{5.4}
\]

The direction \(e_h\) is perpendicular to every wavevector in (5.1), so
it cannot be pressure.  Consequently pressure-only high sums require

\[
                       AB'-A'B=0.                       \tag{5.5}
\]

For an \(H^1\) Bloch packet, (5.5) says

\[
                       B(\theta,t)=\rho(t)A(\theta,t)   \tag{5.6}
\]

on every connected component where \(A\ne0\).  This is only a necessary
condition: finite-\(\epsilon\) in-plane terms impose further high-sum
constraints.  They are \(O(\epsilon^2)\) for a fixed envelope and
optimistically \(O((\epsilon/\eta)^2)\) for an envelope of width \(\eta\).

Substitute (5.6) into (4.2).  The common heat term cancels, leaving

\[
      \rho'+{g\over H}\rho^2
      =gH+2QZ\cos\theta.                                \tag{5.7}
\]

The left side is independent of \(\theta\).  If \(Z\ne0\), the right
side can be constant on the support of \(A\) only when that support is
contained in a level set of \(\cos\theta\).  Such level sets have measure
zero.  An \(L^2\) function supported there vanishes.  Therefore

\[
 \boxed{\text{An active real child admits no nonzero captured
 \(\ell^2\) packet with pressure-only high sums.}}       \tag{5.8}
\]

In particular, no such packet can retain a nonzero \(q\)-stress.

This is the exact obstruction missed by a central-pair truncation.

## 6. Allowing the full charge-zero child profile

One apparent escape from (5.7) is to retain all low harmonics.  In the
physical Bloch coordinate \(y\), let the child be \(z(y,t)e_h\).  The
leading envelope equations have the form

\[
 \begin{aligned}
 a_t&={g\over H}b,\\
 b_t&=(gH-z_y)a,\\
 z_t&=-2\partial_y\operatorname{Re}(a\overline b).
 \end{aligned}                                         \tag{6.1}
\]

If high sums remain pressure, write \(b=\rho(t)a\).  The first two
equations force the child strain \(s=-z_y\) to be spatially constant
wherever \(a\ne0\).  For a real common polarization, the last equation
then gives

\[
                 s_t=2\rho(t)\partial_y^2|a|^2.         \tag{6.2}
\]

The spatial shape of \(a\) is fixed up to a scalar under the first
equation.  Hence \(|a_0|^2\) must have constant second derivative on
every component of its support.  It is quadratic there.  A nonzero
quadratic cannot meet a zero exterior with all derivatives vanishing.
Thus the inviscid leading one-envelope system has no nonzero
\(C_c^\infty\) solution satisfying the exact pressure and uniform-strain
conditions.

This does not rule out a multi-colour construction.  It identifies its
job: realize a quadratic weighted intensity on the next-child core while
pushing taper defects into controlled collars and wakes.

## 7. The approximate envelope window

The obstruction is exact, but it leaves a nonempty asymptotic window.
Let a Bloch packet have width \(\eta\ll1\) near \(\theta=0\), and let the
required logarithmic gain be \(G\) over time \(T\simeq G/g\).  Then

* \(\cos\theta=1+O(\eta^2)\) on the packet;
* polarization drift is \(O(G\eta^2)\);
* the Wronskian high-sum defect is of relative size \(O(G\eta)\);
* WKB incompressibility errors are small when \(\epsilon/\eta\ll1\).

Choose the auxiliary carrier by

\[
                  \nu\Lambda^2=\vartheta g,
                  \qquad {1\over4}<\vartheta<1.         \tag{7.1}
\]

Then the carrier still has positive net exponent
\((1-\vartheta)g\), the child heat ratio is
\(\nu Q^2/g=\vartheta\epsilon^2\), and the homogeneous \(2\Lambda\)
high-sum band has negative leading exponent
\((1-4\vartheta)g\).  Heat spreads the Bloch packet by
\(O(\epsilon\sqrt G)\).  The hierarchy

\[
          \epsilon\sqrt G\ll\eta\ll {1\over G},
          \qquad\text{equivalently}\qquad
          \epsilon G^{3/2}\ll1,                         \tag{7.2}
\]

keeps diffusion, WKB error, and the high-sum defect perturbative.
Leakage from \(\theta=0\) to the unstable sector near \(\pi\) is then of
Gaussian size

\[
                     \exp\!\left[-{c\over\epsilon^2G}\right].
                                                               \tag{7.3}
\]

Equations (7.1)--(7.3) are a consistency window, not a nonlinear theorem.
The missing estimate must control all child harmonics, the damped
\(2\Lambda\) corrector, curvature, pressure, and their feedback for the
whole gain interval.

## 8. Palasek-volume audit of the forced \(2\Lambda\) corrector

The homogeneous damping in Section 7 is not enough.  The forced response
must be measured against the small projection onto the desired child.

Use the lagged-intermittency notation

\[
 q=N^b,\qquad g=N^\beta,
 \qquad \alpha=1+{3\over2b},\qquad 2b<\beta<\alpha,      \tag{8.1}
\]

and the intermediate-carrier, child, and parent volumes

\[
 \begin{aligned}
 V_a&=N^{-\{2b(\alpha-1)-(b-1)\beta\}},\\
 {V_a\over V_c}&=N^{\beta(b-1)},\\
 {V_a\over V_p}
 &=N^{-(b-1)(2\alpha-\beta-2)}.
 \end{aligned}                                         \tag{8.2}
\]

Because \(2b(\alpha-1)=3\), one has \(V_c=N^{-3}\).
Put

\[
 f={V_c\over V_a}=N^{-\beta(b-1)},\qquad
 h={V_a\over V_p}=N^{-w},\qquad
 w=(b-1)\left({3\over b}-\beta\right).                  \tag{8.3}
\]

Let \({\cal A}\) be the \(L^2\) amplitude of the intermediate carrier.
Its approximately uniform low stress has pointwise size
\(q{\cal A}^2/V_a\).  Projection onto an \(L^2\)-normalized child
supported in \(V_c\) gives

\[
 F_c\asymp {q\sqrt{V_c}\over V_a}{\cal A}^2.            \tag{8.4}
\]

This is the desired trilinear coefficient.  If the nonpressure fraction of
the positive-positive interaction is \(\delta_2\), its full high-band norm
is instead

\[
 F_2\asymp {\delta_2 q\over\sqrt{V_a}}{\cal A}^2,
 \qquad
 {F_2\over F_c}\asymp\delta_2 f^{-1/2}.                 \tag{8.5}
\]

Let the damped \(2\Lambda\) band have gap
\(\gamma_2g\), where optimistically
\(\gamma_2=4\vartheta-1>0\).  The forced response and the intended child
therefore obey

\[
 {\|Y_2\|_2\over\|X_c\|_2}
 \asymp {\delta_2\over\gamma_2}f^{-1/2}.                \tag{8.6}
\]

For Palasek's terminal amplitude
\(\|X_c\|_2=q^{\beta-\alpha}\), the identities in (8.1)--(8.4) give

\[
 {\|Y_2\|_2^2\over{\cal A}^2}
 \asymp {\delta_2^2\over\gamma_2^2f},
 \qquad
 {\text{high-band dissipation rate}
  \over\text{parent work rate}}
 \asymp {\delta_2^2\over f}.                            \tag{8.7}
\]

Thus the necessary condition is not merely \(\delta_2=o(1)\), but

\[
                 \boxed{\delta_2=o(\sqrt f)
                  =o\!\left(N^{-\beta(b-1)/2}\right).}  \tag{8.8}
\]

A Bloch packet occupying fraction \(\eta\) of a slow carrier torus needs
container volume \(V_a/\eta\).  Fitting it inside the parent forces

\[
                         \eta\ge {V_a\over V_p}=h.       \tag{8.9}
\]

Generically the child-induced Wronskian defect is at least of effective
order \(\eta\); even the optimistic second-order envelope ledger also
contains \((\epsilon/\eta)^2\).  But (8.8)--(8.9) have no common window:

\[
 h\ll\sqrt f
 \quad\Longleftrightarrow\quad
 (b-1)\left({3\over b}-\beta\right)
 >{\beta(b-1)\over2}
 \quad\Longleftrightarrow\quad
 \beta<{2\over b}.                                     \tag{8.10}
\]

This contradicts \(\beta>2b\) for every \(b>1\).  At the smallest
geometrically possible packet width, the unwanted dissipation exceeds the
available parent work by

\[
 {h^2\over f}
 =N^{\,2(b-1)(3\beta/2-3/b)}\longrightarrow\infty.      \tag{8.11}
\]

Heat diffusion makes the restriction stronger.  With
\(\Lambda^2\asymp g/\nu\),

\[
 \epsilon={q\over\Lambda}=N^{b-\beta/2},\qquad
 G=N^{\beta(b-1)/b},\qquad
 \epsilon\sqrt G=N^{b-\beta/(2b)}.                     \tag{8.12}
\]

Hence a localized stable-sector packet also requires
\(\beta>2b^2\) and
\(\eta\gg N^{b-\beta/(2b)}\).  These conditions cannot repair (8.10).
Equation (8.12) uses the literal Palasek dormant gain.  If a different
architecture preloads a seed needing only \(G\asymp(\log N)^2\), the heat
lower bound is much weaker, but the volume/work contradiction (8.10)--(8.11)
is unchanged because it does not involve \(G\).

Finally, the low stress on \(V_a\setminus V_c\) has no small
\(\delta_2\).  Its orthogonal wake is larger than the desired child
projection by

\[
                         f^{-1/2}
             =N^{\beta(b-1)/2}.                         \tag{8.13}
\]

Its frequency is only \(q\), whose normalized heat rate is
\(O(\epsilon^2g)\), so the spectral filter cannot remove it.  Exact
pressure/cancellation in the buffer, or a genuinely invariant localization
graph concentrating the stress on \(V_c\), is therefore necessary before
the spectral-filter route has a Palasek energy ledger.

## 9. Global symmetry boundary

If one keeps a single cylindrical phase \(m\theta+kz\), all radial
sidebands and all its generated harmonics remain helically symmetric.
Mahalov, Titi, and Leibovich proved global regularity of the corresponding
Navier--Stokes helical invariant subspace.  Therefore a globally periodic
or infinite-column realization confined to this ladder cannot itself be a
blowup solution.  The ladder may only be a local amplifier.  Axial
localization, changing phase directions, or rotating the helical axis
between stages is essential.

This symmetry observation does not exclude localized Gavrilov bubbles:
their cutoffs and multi-axis assembly break the global helical symmetry.

## 10. Strategic consequence

The finite doublet is closed neither algebraically nor dynamically.  A
single exact \(\ell^2\) replacement is also impossible.  Further numerical
optimization of a finite colour set cannot repair this.

The prize-level survivor is more specific:

1. a narrow, Schwartz sideband envelope in the stable Bloch sector;
2. a full localized child profile rather than only \(\pm q\);
3. at least a multi-phase/core--buffer stress chart;
4. a damped high-sum invariant graph; and
5. stage-to-stage breaking of the helical symmetry.

No GPU experiment addresses the missing uniform estimates.  They are an
analytic realization problem.

## References

* A. Mahalov, E. S. Titi, and S. Leibovich, *Invariant helical subspaces
  for the Navier--Stokes equations*, Arch. Rational Mech. Anal. 112
  (1990), 193--222, https://doi.org/10.1007/BF00381234.
* See also the repository notes
  `2026-07-29-two-colour-endpoint-rank.md`,
  `2026-07-29-global-scaled-return-obstruction.md`, and
  `2026-07-29-gavrilov-active-transition-ledger.md` for the charge grading,
  global-return obstruction, and localized transition target.
