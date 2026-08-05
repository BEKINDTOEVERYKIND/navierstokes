# Embedded-quasimode nonlinear-instability gate

**Date:** 2026-08-01
**Status:** new conditional lemma; self-derived, not cross-audited
**Purpose:** decide whether the unfinished thin-torus/AO calculation really
needs an exact finite-frequency eigenvalue.

This note is deliberately narrower than a Navier--Stokes construction.  It
shows that a sufficiently accurate, sufficiently fast growing Euler
**quasimode** can imply nonlinear instability even when its growth rate lies
inside the full bicharacteristic essential growth edge.  It does not turn that
instability into a recurrent cascade.

**Update (2026-08-02).**  The later C54/C63/C65 calculations recover the
straight-column AO asymptotic and match its selected Batchelor ring to the
full principal-symbol BAS edge.  The missing data are no longer those
principal constants.  They are the frequency-uniform PDE semigroup bound,
the real-mode lower bound, a compact base-profile realization, and the
normalized finite-curvature residual in the same norm.

## 1. The recovered endpoint and the remaining data

The uncommitted endpoint of the preceding calculation was reported as follows.

* A high-frequency thin-torus mode, with charge of the form
  \(\eta_M=\alpha/M\), had a formal curvature defect
  \(O(M^{-1})+e^{-cM}\).
* The local bounded-sideband/Floquet complement had been checked.
* The resulting AO growth multiplier lay just **inside** the full
  bicharacteristic essential edge.  Thus an isolated-eigenvalue spectral-gap
  theorem could not be applied.
* The stated missing theorem was a finite-frequency transverse Weber/global
  resolvent estimate.

Later artifacts now give the AO eigenvalue asymptotic and the full BAS edge
for a certified straight Batchelor profile.  They do not give the norm in
which a compatible compact-ring curvature defect is \(O(M^{-1})\), and C67
shows that the standard Gavrilov seed is not that Batchelor base profile.
Consequently no curved-ring application is asserted here.  Sections 3--6
give the remaining data that must be proved.

## 2. What the Friedlander--Strauss--Vishik theorem actually needs

For

\[
  \partial_t w=Lw+N(w),\qquad X\hookrightarrow Z,
\]

Friedlander--Strauss--Vishik (FSV) assume, locally in \(X\),

\[
 \|N(w)\|_Z\le C_0
 \|w\|_X^{1-\eta}\|w\|_Z^{1+\eta},
 \qquad 0<\eta\le 1.                                      \tag{2.1}
\]

Writing

\[
 \Lambda=\lim_{t\to\infty}t^{-1}\log\|e^{tL}\|_Z,
\]

their Theorem 2.2 asks for one orbit with two-sided growth
\(e^{\lambda t}\), where

\[
                 (1+\eta)\lambda>\Lambda.                 \tag{2.2}
\]

It does **not** say that the growing mode must be spectrally isolated.  Hence
an exact smooth eigenmode embedded in essential spectrum is already enough if
(2.2) holds.  The next lemma replaces the exact mode by a quasimode.

Primary sources:

* S. Friedlander, W. Strauss and M. Vishik,
  [*Nonlinear instability in an ideal fluid*](https://www.numdam.org/article/AIHPC_1997__14_2_187_0.pdf),
  Theorem 2.2 and Lemma 3.2.
* R. Shvydkoy,
  [*The essential spectrum of advective equations*](https://arxiv.org/abs/math-ph/0412019),
  for the bicharacteristic-amplitude description of the essential spectrum.

## 3. Quasimode extension

### Proposition 3.1 (conditional abstract instability lemma)

Let \(S(t)=e^{tL}\) be a strongly continuous semigroup on \(Z\).  Assume
local well-posedness in \(X\), (2.1), and, for some \(\Gamma>0\),

\[
                    \|S(t)\|_{Z\to Z}\le C e^{\Gamma t}.
                                                                  \tag{3.1}
\]

Suppose there are real quasimodes \(\phi_j\in D(L)\cap X\) and real numbers
\(0<\lambda_j<\Gamma\) such that

\[
 \|\phi_j\|_Z=1,\qquad
 K_j:=\|\phi_j\|_X\longrightarrow\infty,\qquad
 \|(L-\lambda_j)\phi_j\|_Z\le\varepsilon_j,              \tag{3.2}
\]

and \(\lambda_j\to\lambda\in(0,\Gamma)\).  Put

\[
 b_j={\Gamma-\lambda_j\over\lambda_j}.
\]

If

\[
 \inf_{j\gg1}\big((1+\eta)\lambda_j-\Gamma\big)>0,
 \qquad
 \varepsilon_j K_j^{b_j}\longrightarrow0,                \tag{3.3}
\]

then the zero solution is nonlinearly unstable in \(X\).

For a semigroup with growth bound \(\Lambda\), one may take
\(\Gamma=\Lambda+\delta\) for any fixed \(\delta>0\).  Thus the strict
inequalities can be tested with arbitrarily small growth-bound loss.

### Proof ledger

Let \(r_j=(L-\lambda_j)\phi_j\).  Variation of constants gives

\[
 S(t)\phi_j-e^{\lambda_jt}\phi_j
   =-\int_0^t S(t-s)e^{\lambda_js}r_j\,ds,
\]

and therefore, uniformly for large \(j\),

\[
 \|S(t)\phi_j-e^{\lambda_jt}\phi_j\|_Z
       \le C\varepsilon_j e^{\Gamma t}.                   \tag{3.4}
\]

Set \(P_j=\varepsilon_jK_j^{b_j}\).  By (3.3) one can choose
\(d_j\downarrow0\) so slowly that \(P_jd_j^{-b_j}\to0\), and then take

\[
                  a_j={d_j\over K_j}.                     \tag{3.5}
\]

The initial \(X\)-norm \(a_jK_j=d_j\) tends to zero.  For a fixed small
\(\theta>0\), define the exit-scale time

\[
       T_j={1\over\lambda_j}\log{\theta\over a_j}.
                                                                  \tag{3.6}
\]

The principal linear term at \(T_j\) has \(Z\)-norm \(\theta\), while
(3.4)--(3.6) give

\[
 \begin{aligned}
 a_j\|S(T_j)\phi_j-e^{\lambda_jT_j}\phi_j\|_Z
 &\le C\theta\varepsilon_j
       \left({\theta\over a_j}\right)^{b_j}\\
 &\le C_\theta\theta P_jd_j^{-b_j}=o(\theta).             \tag{3.7}
 \end{aligned}
\]

Assume for contradiction that all these small data remain in the local
\(X\)-stability ball.  Under the usual continuity bootstrap
\(\|w_j(t)\|_Z\lesssim a_je^{\lambda_jt}\), Duhamel's formula and (2.1)
give, using the first condition in (3.3),

\[
 \begin{aligned}
 \left\|\int_0^tS(t-s)N(w_j(s))\,ds\right\|_Z
 &\lesssim \rho_0^{1-\eta}a_j^{1+\eta}
   \int_0^t e^{\Gamma(t-s)}e^{(1+\eta)\lambda_js}\,ds\\
 &\lesssim \rho_0^{1-\eta}a_j^{1+\eta}
             e^{(1+\eta)\lambda_jt}.                     \tag{3.8}
 \end{aligned}
\]

At \(t=T_j\), (3.8) is at most
\(C\rho_0^{1-\eta}\theta^{1+\eta}\).  Choose \(\theta\) once so that
this is below \(\theta/4\), then take \(j\) large enough that (3.7) is
below \(\theta/4\).  It follows that
\(\|w_j(T_j)\|_Z\ge\theta/2\), hence also that its \(X\)-norm has a fixed
positive lower bound.  The standard first-exit alternative closes the
bootstrap: exiting the \(X\)-stability ball earlier is already the desired
instability.  This contradicts nonlinear stability for data whose
\(X\)-norm tends to zero.  \(\square\)

The real-quasimode hypothesis avoids an irrelevant complexification issue.
A complex version needs a uniform lower-growth statement for a real or
imaginary part at the chosen exit times.

## 4. A sharper Euler interpolation exponent

For \(n\)-dimensional Euler, take \(X=H^s\), \(Z=L^2\), and put

\[
                         \beta={n\over2}+1.                \tag{4.1}
\]

FSV's printed Lemma 3.2 uses

\[
              \eta_{\rm FSV}={1\over2}-{n+2\over4s}
                    ={1\over2}-{\beta\over2s}.             \tag{4.2}
\]

There is a standard sharper distribution of derivatives.  Choose
\(a,b\ge0\) with \(a+b>n/2\).  Sobolev multiplication and interpolation
give

\[
\begin{aligned}
 \|\mathbb P(w\cdot\nabla w)\|_2
 &\lesssim \|w\|_{H^a}\|w\|_{H^{b+1}}\\
 &\lesssim
 \|w\|_{H^s}^{(a+b+1)/s}
 \|w\|_2^{2-(a+b+1)/s}.                                  \tag{4.3}
\end{aligned}
\]

Consequently (2.1) holds for every

\[
                       0<\eta<1-{\beta\over s}.            \tag{4.4}
\]

The endpoint in (4.4) is not claimed.  The arbitrarily small loss in
\(a+b>n/2\) is harmless for all strict inequalities below.  Equation (4.3)
is an analytic input, not something certified by the arithmetic checker.

## 5. The polynomial-frequency gate

Assume a normalized frequency-\(M\) quasimode has

\[
 K_M\lesssim M^s,\qquad \varepsilon_M\lesssim M^{-q},
 \qquad {\lambda_M\over\Lambda}\longrightarrow\chi\in(0,1).
                                                                  \tag{5.1}
\]

Letting the semigroup loss \(\delta\downarrow0\), the two conditions in
(3.3), with the sharper estimate (4.4), become

\[
 s>{\beta\chi\over2\chi-1},
 \qquad
 s<{q\chi\over1-\chi}.                                   \tag{5.2}
\]

There exists such an \(s\) if and only if

\[
                  \boxed{\displaystyle
                  \chi>{\beta+q\over\beta+2q}}.           \tag{5.3}
\]

In three dimensions \(\beta=5/2\).  For comparison, using only FSV's
printed exponent (4.2) yields the more conservative gate

\[
                  \chi>{\beta+2q\over\beta+3q}.            \tag{5.4}
\]

| residual order \(q\) | sharper product gate (5.3) | printed-FSV gate (5.4) |
|---:|---:|---:|
| 1 | \(7/9\) | \(9/11\) |
| 2 | \(9/13\) | \(13/17\) |
| 3 | \(11/17\) | \(17/23\) |
| 4 | \(13/21\) | \(21/29\) |
| \(q\to\infty\) | \(1/2\) | \(2/3\) |

For example, if \(q=1\) and \(\chi=4/5\), then \(s=7/2\) works:
the quasimode leakage exponent is
\(s(\chi^{-1}-1)=7/8<1\), while one may choose
\(\eta=15/56<1-(5/2)/(7/2)=2/7\), giving
\((1+\eta)\chi=71/70>1\).

Equivalently, for a measured ratio \(\chi>1/2\), the required residual
order is

\[
                    q>{\beta(1-\chi)\over2\chi-1}.         \tag{5.5}
\]

Thus higher-order curvature WKB has a precise target; “take more terms” is
not an unquantified instruction.

### 5.1 The \(L^p\) lever

The \(L^2/H^s\) choice is not forced.  FSV explicitly note that their Euler
argument extends to

\[
                 Z=L^p,\qquad X=W^{s,p},qquad
                 1<p<\infty,\quad s>{n\over p}+1.         \tag{5.6}
\]

For fixed finite \(p\), the Sobolev product calculation in (4.3) becomes

\[
 \|w\cdot\nabla w\|_{L^p}
 \lesssim \|w\|_{W^{a,p}}\|w\|_{W^{b+1,p}},
 \qquad a+b>{n\over p}.                                  \tag{5.7}
\]

Thus (4.4) and (5.3) hold with

\[
             \beta_p={n\over p}+1,qquad
             \boxed{\displaystyle
             \chi_p>{\beta_p+q\over\beta_p+2q}}.          \tag{5.8}
\]

Here \(\chi_p=\lambda/\Lambda_p\) must use the **full \(L^p\) semigroup
growth bound**.  It is not legitimate to insert the \(L^2\) edge without
proving that the two edges agree.

For \(n=3\) and a first-order residual:

| fixed phase space | \(\beta_p\) | required \(\chi_p\) |
|---|---:|---:|
| \(p=2\) | \(5/2\) | \(7/9\) |
| \(p=3\) | \(2\) | \(3/4\) |
| \(p=6\) | \(3/2\) | \(5/7\) |
| \(p\to\infty\) | \(1\) | tends to \(2/3\) |

The limit row is only an infimum over fixed finite \(p\); no boundedness of
the Leray projection on \(L^\infty\) is assumed.  Also, if the thin-torus
support shrinks with \(M\), both the normalized \(L^p\) residual and the
\(W^{s,p}\) cost must be recomputed.  The \(L^p\) lever is useful only if
those powers and \(\Lambda_p\) do not erase its gain.

## 6. AO decision gate

The unfinished AO branch can bypass the exact transverse-Weber eigenvalue
theorem if all four items below are proved in a common normalization.

1. **Growth ratio:** compute
   \(\chi=\liminf\lambda_M/\Lambda\), where \(\Lambda\) is the growth bound
   on the same \(L^2\) space, not merely a restricted charge sector.
2. **Residual norm:** prove
   \(\|(L-\lambda_M)\phi_M\|_2\lesssim M^{-q}\) after
   \(\|\phi_M\|_2=1\).  A pointwise or formal curvature remainder is not
   enough.
3. **Strong norm cost:** prove
   \(\|\phi_M\|_{H^s}\lesssim M^s\), including cutoffs, pressure recovery,
   torus geometry and all bounded sidebands.
4. **Domain/reality:** put the mode in \(D(L)\), keep it divergence-free,
   and supply either a real quasimode or the real-part lower-growth variant.

If the reported \(O(M^{-1})\) remainder satisfies item 2 and
\(\chi_2>7/9\), Proposition 3.1 closes a nonlinear Euler-instability lemma
without any isolated eigenvalue.  More generally, a fixed finite \(p\) can
lower that threshold toward \(2/3\), provided the full \(L^p\) edge and all
normalized powers are established.  If \(1/2<\chi_p\) lies below the
first-order threshold, the \(p\)-analogue of (5.5) states the minimum WKB
order to build.  If \(\chi_p\le1/2\) for every admissible phase space, this
particular mixed-norm argument cannot work, even with superalgebraic
residual.

## 7. What this does and does not buy

This gate may remove one exact spectral bottleneck.  It does **not** solve
the downstream problem.  A Millennium-level construction would still have
to:

* convert instability into a controlled finite-time endpoint map rather
  than mere departure from a steady state;
* join that map to the verified disjoint Gavrilov DSS wake;
* replace the refuted one-carrier positive-stress target by a genuinely
  multi-colour transition;
* retain and sum the nonlocal pressure/velocity wake forced by the second
  viscous endpoint jet; and
* obtain an exact smooth Navier--Stokes solution satisfying one of the Clay
  alternatives, not a forced model or a formal cascade.

The result is therefore a conditional bridge lemma and a quantitative
research gate, not evidence that singularity formation has been proved.
