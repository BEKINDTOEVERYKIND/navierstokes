# Factorial one-cell schedule: recovered claims C126--C131

**Date:** 2026-08-05

**Status:** exact scale algebra and summability, plus explicitly conditional
localization and exit-chart gates.  This note does **not** construct a
localized Navier--Stokes stage map.  In particular, it does not promote the
finite A2 ladder to an invariant subsystem of the PDE.

This is a conservative reconstruction of the missing C126--C131 artifact
layer.  It uses only the one-cell geometry already selected by the program.
The purpose is to freeze one schedule on which energy, time, viscosity,
separated pressure, and retained-wake budgets all converge, so that the
remaining question cannot be hidden inside a change of scales.

## 1. C126: exact monochromatic compact localization is impossible

Let \(u\in C_c^\infty(\mathbb R^3)\) and suppose

\[
             -\Delta u=\kappa^2u,\qquad \kappa>0.                 \tag{1.1}
\]

Taking Fourier transforms gives

\[
        (|\xi|^2-\kappa^2)\widehat u(\xi)=0.                     \tag{1.2}
\]

Thus \(\widehat u\) vanishes on the nonempty open set
\(\{|\xi|\ne\kappa\}\).  Because compact support makes \(\widehat u\)
entire, analytic continuation gives \(\widehat u\equiv0\), hence \(u=0\).
The same conclusion applies componentwise to a compactly supported Beltrami
eigenfield.

Consequently an exact eigenshell pump cannot also be compactly localized.
Every one-cell theorem must retain a collar defect, use a non-monochromatic
packet, or prove an exact nonlinear cancellation.  Treating localization as
free would remove the actual obstruction.

This statement is **PROVED**.

## 2. C127: the factorial schedule

For stages \(j\ge1\), put \(n=j+1\) and define

\[
 \ell_j=(j!)^{-8},\qquad
 a_j=(j!)^{10},\qquad
 \kappa_j=\ell_j^{-1}=(j!)^8.                                \tag{2.1}
\]

The exact one-step frequency and amplitude ratios are

\[
 q_j={\kappa_{j+1}\over\kappa_j}=n^8,
 \qquad
 g_j={a_{j+1}\over a_j}=n^{10}.                              \tag{2.2}
\]

Introduce the dormant-seed fraction and the ideal three-dimensional
concentration factor

\[
 b_j=n^{-2},\qquad F_j=q_j^{3/2}=n^{12}.                     \tag{2.3}
\]

Then the handoff identity is exact:

\[
                      b_jF_j=g_j.                            \tag{2.4}
\]

The natural energy, one-inertial-time scale, Reynolds number, normalized
one-inertial-time heat action, and corresponding dissipation scale are

\[
\begin{aligned}
 E_j&=a_j^2\ell_j^3=(j!)^{-4},\\
 \tau_j&={\ell_j\over a_j}=(j!)^{-18},\\
 {\rm Re}_j&={a_j\ell_j\over\nu}=\nu^{-1}(j!)^2,\\
 \mu_j&={\nu\kappa_j\over a_j}=\nu(j!)^{-2},\\
 D_j&=E_j\mu_j=\nu(j!)^{-6}.
                                                               \tag{2.5}
\end{aligned}
\]

Therefore

\[
 \sum_jE_j+\sum_j\tau_j+\sum_j\mu_j+\sum_jD_j<\infty,
\qquad {\rm Re}_j\longrightarrow\infty.                    \tag{2.6}
\]

The decaying-pump and routing blocks use a logarithmic number of inertial
times.  Since \(\log q_j=8\log n\), the larger time, heat, and dissipation
budgets are

\[
 \widetilde\tau_j=8(\log n)\tau_j,\qquad
 \widetilde\mu_j=8(\log n)\mu_j,\qquad
 \widetilde D_j=8(\log n)D_j.                               \tag{2.7}
\]

They are summable as well (for example \(\log n\le n\), followed by C130).
Thus the logarithmic gain time required later is charged explicitly rather
than hidden in the nominal inertial-time scale.

Also

\[
 {E_{j+1}\over E_j}=n^{-4},\qquad
 {\tau_{j+1}\over\tau_j}=n^{-18},\qquad
 {\mu_{j+1}\over\mu_j}=n^{-2}.                              \tag{2.8}
\]

These are arithmetic identities, not asymptotic comparisons.  The gradient
scale \(a_j/\ell_j=(j!)^{18}\) diverges while the logarithmic-window time
budget and the sum of the nominal packet energies remain finite.  This is
the scalar singularity schedule that a full stage map would iterate.

This statement is **PROVED AT THE SCALAR-LEDGER LEVEL**.

## 3. C128: exact curl cutoff, with a nonlocal conditional gate

Let \(W\) be a smooth Beltrami field on a neighborhood containing the
cutoff support (for example, a global Beltrami field),
\(\nabla\times W=\kappa W\), and let

\[
 \chi_M(x)=\chi(\kappa x/M),
\]

where \(\chi=1\) on the unit core and is compactly supported in a fixed
larger ball.  The field

\[
 U_M=\nabla\times\left({\chi_M W\over\kappa}\right)
     =\chi_MW+{\nabla\chi_M\over\kappa}\times W              \tag{3.1}
\]

is exactly divergence free, compactly supported, and agrees with \(W\) on
the core.  If the carrier-normalized derivatives of \(W\) are uniformly
bounded, every fixed carrier-normalized derivative of the collar correction
is \(O(M^{-1})\); each additional derivative that lands on the cutoff gives
another factor \(M^{-1}\).

Equation (3.1) is useful but is not a localized Euler or Navier--Stokes
solution.  The cutoff stress produces a Leray-pressure field everywhere,
and C126 shows that no choice of a literal compact monochromatic cutoff can
remove all of it.  The usable one-source/one-target hypothesis is:

\[
 \boxed{
 \text{before endpoint-chart inversion, the separated collar contribution}
 \ \le C{\log q_j\over M_j}.}                                \tag{3.2}
\]

Only this displayed estimate is assumed below.  Uniformity under the full
nonlinear evolution, repeated sources, and carried wakes is **OPEN**.

Thus the curl identity is **PROVED** and the PDE estimate (3.2) is
**CONDITIONAL**.

## 4. C129: an exponential localization window

Choose the separation/collar width

\[
                         M_j=n^{7/2}.                         \tag{4.1}
\]

If the already-constructed local packet has the displayed Gaussian tail
\(Ce^{-c|y|^2}\) in carrier coordinates, truncating at \(M_j\) gives

\[
                    \varepsilon_j^{\rm loc}
                    \le Ce^{-cM_j^2}=Ce^{-cn^7}.             \tag{4.2}
\]

Every fixed polynomial derivative or chart factor is absorbed:

\[
             \sum_j n^C e^{-cn^7}<\infty
             \quad\text{for every fixed }C,c>0.             \tag{4.3}
\]

This does not assert that the full A2 pump has the required spatial tail;
it records the precise implication if the local packet construction supplies
one.  The tail implication (4.2)--(4.3) is **PROVED CONDITIONALLY ON THE
DISPLAYED LOCAL DECAY**.

## 5. C130: factorial decay absorbs every fixed polynomial gain

For \(\alpha>0\) and fixed \(C\), the ratio test gives

\[
 { (j+1)^C/((j+1)!)^\alpha \over j^C/(j!)^\alpha}
 =\left(1+{1\over j}\right)^C(j+1)^{-\alpha}\longrightarrow0.
                                                                    \tag{5.1}
\]

Hence

\[
             \sum_j{j^C\over(j!)^\alpha}<\infty.             \tag{5.2}
\]

In particular, any *fixed* polynomial conditioning loss can be placed on
the energy, viscosity, or dissipation quantities in (2.5) without changing
summability.  This statement does not absorb exponential-in-\(j^2\) endpoint
losses or an uncontrolled number of derivative losses.

This statement is **PROVED**.

## 6. C131: the shrinking exit chart

Let \(\mathcal E_j\) denote the finite-dimensional exit map after all
infinite-dimensional PDE estimates have been established.  Suppose, on a
fixed control ball, that \(D\mathcal E_j(0)\) has a right inverse \(R_j\)
with

\[
                       \|R_j\|\le L_j,\qquad L_j\le n^2,     \tag{6.1}
\]

and that the usual quantitative inverse-function smallness condition is
uniform after multiplication by \(L_j\).  Then a fixed control ball covers
an output ball of radius \(c/L_j\ge cn^{-2}\).  The dormant seed exponent in
(2.3) is exactly the allowed exponent, since

\[
                         L_jb_j\le1.                          \tag{6.2}
\]

Equation (6.2) fixes only the power of \(n\), not the radius constant.  A
unit-coefficient target \(b_j=n^{-2}\) is covered only if the quantitative
radius constant is large enough (or after a fixed normalization of the seed
and focus).  This is the reason for the \(n^{-2}\) seed scale and the allowed
\(j^2\) chart loss.  It is not a proof that \(\mathcal E_j\) exists, that its
derivative is onto, that its nonlinear Lipschitz constant is uniform, or
that the required radius constant is available.  Those are parts of the
one-cell closure theorem.

This statement is an **EXACT CONDITIONAL INVERSE-MAP GATE**.

## 7. The complete summable scalar ledger

Under the one-source/one-target estimate (3.2), applying the worst allowed
exit-chart loss gives

\[
 \varepsilon_j^{\rm press}
 \le C L_j{\log q_j\over M_j}
 \le 8C{\log n\over n^{3/2}}.                               \tag{7.1}
\]

The integral test shows that its sum converges.  Prescribe the retained wake
increment by

\[
                    \varepsilon_j^{\rm wake}\le Cn^{-4}.     \tag{7.2}
\]

Here \(\varepsilon_j^{\rm press}\) is the size of the *control correction*
after applying the right inverse: the pre-chart output defect is
\(O((\log n)n^{-7/2})=o(n^{-2})\), so it lies inside the conditional output
ball from C131, and the chart is supposed to cancel it.  It is not an
allowed terminal active-coordinate error and cannot be substituted for the
\(n^{-6}\) pre-chart BAFL target.  Exact cancellation and the absence of new
leakage from that control remain conditional on the endpoint map.

Then, at the scalar level,

\[
 \sum_j\left(
 E_j+\widetilde\tau_j+\widetilde\mu_j+\widetilde D_j+
 \varepsilon_j^{\rm loc}+\varepsilon_j^{\rm press}+
 \varepsilon_j^{\rm wake}
 \right)<\infty.                                            \tag{7.3}
\]

No mode-leakage term appears in (7.3), because no such estimate has been
proved.  It is not licensed to set that term to zero.  The next note names
the exact backward-weighted estimate needed to add it.

## Claim boundary

* C126, C127, and C130 are exact.
* C128 contains an exact divergence-free cutoff identity and a clearly
  isolated conditional pressure-tail estimate.
* C129 is an exact implication from a stated packet-tail hypothesis.
* C131 is a quantitative conditional exit-chart gate.
* None of C126--C131 proves an unforced localized Navier--Stokes stage.

The accompanying dependency-free checker verifies every exponent and
summability comparison used above.  A passing checker is not evidence for
the conditional PDE hypotheses.
