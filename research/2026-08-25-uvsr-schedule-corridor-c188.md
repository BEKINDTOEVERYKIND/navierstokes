# C188: the UVSR demand corridor and the epsilon-regularity destruction test

**Date:** 2026-08-25

**Status:** exact optimization of two explicitly declared scalar schedule
classes; exact viscosity-coordinate law; explicit smooth-point
epsilon-functional bounds; source-scoped CKN/Tao audit; adversarial
same-session re-derivation recorded in
[`audit/AUDIT-2026-08-25-c188.md`](../audit/AUDIT-2026-08-25-c188.md);
no UVSR profile or singular solution constructed

**Checker:**
[checks/uvsr_schedule_corridor_c188.py](../checks/uvsr_schedule_corridor_c188.py)

## 0. Verdict and claim boundary

The pre-registered corridor is **nonempty**, but two corrections are
load-bearing.

First, C135's raw ledger multiplier \(q^{3/2}\) is not an accident of the
factorial exponents.  Under its equal normalized \(L^2\)-shape convention,
it is forced whenever the same leading seed energy is transferred from a
parent volume into a child volume smaller by \(q^{-3}\).  Uniformly
comparable seed/child shape constants preserve the exponent but insert an
explicit constant ratio; see (1.3).  What the schedule can minimize is the
**net amplitude gain**

\[
             g=bq^{3/2}=q^\gamma.
\]

Second, Caffarelli--Kohn--Nirenberg epsilon regularity does not require
every smooth intermediate packet to remain above one epsilon threshold.
Its contrapositive lower bound applies to shrinking backward cylinders
centered at an actual singular point.  The local energy inequality also
contains an indefinite pressure/advective boundary flux, so it gives no
schedule-only local upper wall.

There was no repository definition of “the entire admissible schedule
space”: C127 deliberately froze one schedule.  C188 therefore declares
the admissible demand class explicitly.  It first optimizes arbitrary
one-step sequences under the inherited per-return high-Re rule and, for
the legacy branch, the displayed C176 worst-case envelope.  It then uses
generalized factorial schedules to witness compatibility with the landed
summability ledgers.  No class is changed after the test.

1. In the **bounded-profile legacy-collar envelope class**, the landed
   C176 upper majorant \(C(1+\Lambda)q^{-1}J^{7/2}\), with
   \(J\le C_J\log q\) and \(\sup_j\Lambda_j<\infty\), is required to be
   smaller than the inherited \(b^3\) active tolerance without a new
   cancellation.  The exact polynomial corridor for this sufficient
   worst-case certificate is

   \[
                         \boxed{\frac76<\gamma<\frac32}.      \tag{0.1}
   \]

   Its lower endpoint is a sharp infimum within this declared envelope and
   is not attained.  This is not a lower bound on the actual collar
   response, since C176 proved no matching lower estimate.  The current
   C127 exponent \(5/4\) is admissible but not envelope-optimal.
2. In the **direct-UVSR class**, the collar is a component of the full
   renormalization residual and is certified or cancelled there rather
   than requiring \(q^{-1}J^{7/2}/b^3\to0\) separately.  The bare
   finite-energy/high-Re corridor is

   \[
                         \boxed{1<\gamma<\frac32}.             \tag{0.2}
   \]

   Again the lower endpoint is an infimum, not an attained pure power.

No additional numerical power restriction follows from CKN or Tao using
only the landed scalar ledger.  The published thresholds and exponents are
not numerical, and their theorems do not imply a per-stage amplitude
condition absent singular-center, scale-occupancy, pressure, and wake
hypotheses.  Consequently the immediate
UVSR search specification is (0.1) as long as the bounded-profile C176
worst-case majorant is paid separately.  Replacing it by (0.2) requires a
certified smaller collar response or full-residual cancellation, not a
relabeling of the old upper bound.

There is also an exact fixed-point consequence.  The normalized viscosity
obeys \(\mu_{j+1}=(q_j/g_j)\mu_j\).  Hence a fixed supercritical scaling
\(g>q\) cannot have an autonomous fixed point at \(\mu>0\).  The actual
factorial schedule is nonautonomous in \(n=j+1\), so its honest object is a
sequence of stage maps, or the augmented skew-product
\((X,\mu,n)\mapsto(\mathcal R_{n,\mu}X,(q_n/g_n)\mu,n+1)\).  Although
\(\mu_j\to0\), an autonomous inviscid limiting operator and fixed point
would additionally require a proved normalized limit
\(\mathcal R_{n,0}\to\mathcal R_{\infty,0}\); C188 does not prove one.

## 1. Schedule-independent identities

Let

\[
 q_j={\ell_j\over\ell_{j+1}}>1,\qquad
 g_j={a_{j+1}\over a_j}>1,\qquad
 E_j=a_j^2\ell_j^3,\qquad
 {\rm Re}_j={a_j\ell_j\over\nu}.                          \tag{1.1}
\]

Then, with no power-law assumption,

\[
 {E_{j+1}\over E_j}={g_j^2\over q_j^3},\qquad
 {{\rm Re}_{j+1}\over{\rm Re}_j}={g_j\over q_j},\qquad
 {\tau_{j+1}\over\tau_j}={1\over g_jq_j},\quad
 {\mu_{j+1}\over\mu_j}={q_j\over g_j},\quad
 \tau_j={\ell_j\over a_j},\quad
 \mu_j={\nu\over a_j\ell_j}.                              \tag{1.2}
\]

A dormant seed of relative amplitude \(b_j\) carries ledger energy
\(b_j^2E_j\).  Under the equal normalized \(L^2\)-shape convention, if
that same energy is the child energy, then

\[
 b_j^2={g_j^2\over q_j^3},\qquad
 b_j={g_j\over q_j^{3/2}},\qquad
\boxed{F_j={g_j\over b_j}=q_j^{3/2}}.                    \tag{1.3}
\]

More generally, if the seed and child normalized profiles have squared
\(L^2\) constants \(\kappa_{s,j}\) and \(\kappa_{c,j}\), equality of their
actual energies gives
\[
 F_j=q_j^{3/2}\left({\kappa_{s,j}\over\kappa_{c,j}}\right)^{1/2}.
                                                               \tag{1.3a}
\]
Thus \(q^{3/2}\) is an exact ledger multiplier, and remains the physical
power when \(0<\kappa_-\le\kappa_{s,j},\kappa_{c,j}\le\kappa_+<\infty\).
An exact UVSR profile must certify these constants rather than silently
identify two changing shapes.

Thus changing the schedule cannot lower the equal-shape ledger multiplier.
It changes \(b_j\) and hence \(g_j\).  Under that normalization, a proposal
with \(F_j<q_j^{3/2}\) uses only a fraction of the seed energy; the
untransferred energy must then be dissipated, exported, or retained in the
wake and charged inside the complete UVSR state.

If a stage renormalization is autonomous and fixes a positive normalized
viscosity, (1.2) forces \(g_j=q_j\).  That is the critical boundary
\(\gamma=1\), where the Reynolds number does not increase.  Every
high-Re schedule with \(\gamma>1\) instead has

\[
                {\mu_{j+1}\over\mu_j}=q_j^{1-\gamma}<1.   \tag{1.4}
\]

Thus a fixed-scaling model may seek an inviscid-face fixed point with a
proved viscous stable/trapping direction.  For the actual factorial
schedule, however, any future stage map must have the nonautonomous
augmented coordinate form
\[
 (X,\mu,n)\longmapsto
 \left(\mathcal R_{n,\mu}X,{q_n\over g_n}\mu,n+1\right). \tag{1.5}
\]
The \((\mu,n)\) coordinates have no finite-\(n\) fixed point.  The state
space and \(\mathcal R_{n,\mu}\) are not yet defined; their construction is
open.  An autonomous boundary fixed point is conditional on a separately
certified normalized limiting map.

## 2. The generalized factorial ledger

Under the equal-shape ledger convention in (1.3), before imposing any
factorial or power ansatz, define
\[
 g_j=q_j\rho_j,\qquad b_j=\rho_jq_j^{-1/2},\qquad
 1<\rho_j<q_j^{1/2}.                                     \tag{2.0}
\]
The lower inequality is exactly the inherited per-return high-Re rule; the
upper inequality is exactly the dormant-seed rule \(b_j<1\).  Therefore the
polynomial demand floor over every one-step sequence in this direct class
is \(1\) when \(q_j\to\infty\) and
\(\log\rho_j/\log q_j\to0\).  A fixed \(\rho_j=\rho_*>1\) realizes that
polynomial order, although the normalized pure monomial convention
\(g=q^\gamma\) requires \(\gamma>1\).  This
statement is local to the per-return rule; allowing Reynolds-decreasing
stages offset by later stages would define a different architecture and is
not admitted by the landed stage ledger.

The exact additional global energy condition for an arbitrary sequence is
\[
 E_{j_0}\sum_{j=j_0}^{\infty}
       \prod_{k=j_0}^{j-1}b_k^2<\infty,                  \tag{2.0a}
\]
with the empty product equal to one.  Equations (2.0) and (3.2a)--(3.2c)
are one-step demand identities; the explicit schedules below separately
verify (2.0a) and the inherited feedback sum.

Put \(n=j+1\), and take fixed exponents \(Q,S,P>0\):

\[
 \ell_j=(j!)^{-Q},\qquad a_j=(j!)^P,\qquad
 q_j=n^Q,\qquad b_j=n^{-S}.                               \tag{2.1}
\]

The same-energy identity (1.3) is equivalent to

\[
       P={3Q\over2}-S,\qquad
       g_j=n^P=q_j^\gamma,\qquad
       \gamma={P\over Q}={3\over2}-{S\over Q}.            \tag{2.2}
\]

Every scalar ledger is then exact:

\[
\begin{aligned}
 E_j&=(j!)^{-2S},\\
 \tau_j&=(j!)^{-(5Q/2-S)},\\
 {\rm Re}_j&=\nu^{-1}(j!)^{Q/2-S},\\
 \mu_j&={\nu\ell_j^{-1}\over a_j}
       =\nu(j!)^{S-Q/2},\\
 D_j&=E_j\mu_j=\nu(j!)^{-(S+Q/2)}.                       \tag{2.3}
\end{aligned}
\]

For the explicitly certified class \(S\ge2\), the parent-frequency
feedback has the numerical bound

\[
 \sum_{n=2}^\infty b_j\log q_j
 =Q\sum_{n=2}^\infty {\log n\over n^S}
 \le 2Q
 \quad\text{for }S\ge2.                                  \tag{2.4}
\]

Here \(\log n\le n^{1/2}\) and
\(\sum_{n=2}^\infty n^{-3/2}\le2\).  The Reynolds number grows precisely
when

\[
                         Q>2S,                             \tag{2.5}
\]

which is equivalent to \(\gamma>1\).  For \(S>0\), energy is factorially
summable and \(\gamma<3/2\).  This proves the direct-UVSR corridor (0.2).

The general-sequence boundary in (2.0) has an exact scalar witness, not
only an approaching pure-power family.  For \(n=j+1\), take
\[
 \ell_j=(j!)^{-4},\qquad a_j=2^j(j!)^4,qquad
 q_j=n^4,qquad g_j=2n^4=2q_j,qquad b_j=2n^{-2}.          \tag{2.6}
\]
For \(n\ge2\), \(0<b_j<1\), and the exact ledgers are
\[
 E_j=4^j(j!)^{-4},\quad
 \tau_j=2^{-j}(j!)^{-8},\quad
 {\rm Re}_j=\nu^{-1}2^j,\quad
 \mu_j=\nu2^{-j},\quad
 D_j=\nu2^j(j!)^{-4}.                                    \tag{2.7}
\]
Thus every return doubles Reynolds while the gain has polynomial order
exactly one.  The feedback sum has the explicit bound
\[
 \sum_{n=2}^{\infty}b_j\log q_j
 =8\sum_{n=2}^{\infty}{\log n\over n^2}
 \le8\sum_{n=2}^{\infty}n^{-3/2}\le16.                  \tag{2.8}
\]
Moreover \(E_{j+1}/E_j=4/(j+1)^4\le1/4\) for \(j\ge1\), so
\[
             \sum_{j=1}^{\infty}E_j\le{4\over3}E_1={16\over3}. \tag{2.9}
\]
This is a scalar boundary witness; changing \(Q=8\) still requires new
shell geometry if it is used physically.

No statement here promotes these schedules to a PDE stage.  They are the
core scalar demand ledgers on which a UVSR profile would be tested.

## 3. Exact optimization with the landed collar loss

C176's largest displayed bounded-profile collar majorant before the active
chart has the form

\[
 \varepsilon_j^{\rm col}
 \le C_{\rm col}J_j^{7/2}q_j^{-1},\qquad
 J_j\le C_J\log q_j,                                     \tag{3.1}
\]

where \(C_{\rm col},C_J\ge1\) are fixed constants.  Here a uniform bound
on C176's profile factor \(1+\Lambda_j\) is absorbed into
\(C_{\rm col}\).  If \(\Lambda_j\) is not uniformly bounded, it must stay
in the comparison and the exponent window below does not follow by itself.
Requiring this specific inherited upper majorant to enter the \(b_j^3\)
active tolerance with a vanishing relative margin gives

\[
 {\varepsilon_j^{\rm col}\over b_j^3}
 \le C_{\rm col}(C_JQ\log n)^{7/2}n^{\,3S-Q}
 \longrightarrow0.                                      \tag{3.2}
\]

This optimization does not require a power ansatz.  Set
\[
 K_*=C_{\rm col}C_J^{7/2},\qquad
 \omega_j={b_jq_j^{1/3}\over
                 K_*^{1/3}(\log q_j)^{7/6}}.             \tag{3.2a}
\]
Then the declared upper-envelope ratio is exactly
\[
 {K_*q_j^{-1}(\log q_j)^{7/2}\over b_j^3}
 ={1\over\omega_j^3}.                                    \tag{3.2b}
\]
Consequently this worst-case certificate has a vanishing relative margin
exactly when \(\omega_j\to\infty\), and its equal-shape all-sequence net
demand is
\[
 g_j=K_*^{1/3}q_j^{7/6}(\log q_j)^{7/6}\omega_j.          \tag{3.2c}
\]
The energy ledger additionally requires \(b_j<1\) and (2.0a), so not every
rapidly growing \(\omega_j\) is admissible.  Choices satisfying
\[
 \omega_j\to\infty,qquad
 {\log\omega_j\over\log q_j}\to0                         \tag{3.2d}
\]
approach the \(7/6\) polynomial floor.  The explicit slowly divergent
choice below is admissible, and the factorial calculation checks that the
other inherited scalar ledgers do not raise the floor.

For fixed power exponents, this displayed worst-case certificate holds
exactly when

\[
                           Q>3S.                           \tag{3.3}
\]

At equality the positive logarithmic power diverges; below equality the
polynomial factor diverges.  Combining (2.2) and (3.3),

\[
 0<{S\over Q}<{1\over3}\qquad\Longleftrightarrow\qquad
 {7\over6}<\gamma<{3\over2}.                             \tag{3.4}
\]

This is the sharp corridor for the declared sufficient legacy-collar
envelope (0.1), not a necessity theorem for the physical collar.

### 3.1 An exact boundary-approaching integer family

For every integer \(m\ge1\), set

\[
 Q_m=6m+2,\qquad S_m=2m,\qquad P_m=7m+3.                  \tag{3.5}
\]

Then

\[
\begin{gathered}
 P_m={3Q_m\over2}-S_m,\qquad Q_m-3S_m=2,\\
 \gamma_m={7m+3\over6m+2}
          ={7\over6}+{1\over9m+3}\downarrow{7\over6}.   \tag{3.6}
\end{gathered}
\]

The corresponding exact core ledgers are

\[
\begin{aligned}
 \ell_j&=(j!)^{-(6m+2)},& a_j&=(j!)^{7m+3},\\
 q_j&=n^{6m+2},& b_j&=n^{-2m},\\
 F_j&=n^{9m+3}=q_j^{3/2},& g_j&=n^{7m+3}=q_j^{\gamma_m},\\
 E_j&=(j!)^{-4m},& \tau_j&=(j!)^{-(13m+5)},\\
 {\rm Re}_j&=\nu^{-1}(j!)^{m+1},&
 \mu_j&=\nu(j!)^{-(m+1)},\\
 D_j&=\nu(j!)^{-(5m+1)}.&&                              \tag{3.7}
\end{aligned}
\]

One compatible scalar chart/separation choice is
\[
                   L_j=n^2,\qquad M_j=n^{S_m+2}.          \tag{3.8}
\]
It gives \(L_jb_j=n^{2-2m}\le1\),
\[
 {\log q_j\over M_j}\le {Q_m\over n}\,b_j
 \quad(\log n\le n),
\]
and
\(L_j(\log q_j)/M_j=Q_m(\log n)n^{-2m}\), whose sum is at
most \(2Q_m\) by (2.4).  This prevents the family from silently retaining
C127's frozen \(M=n^{7/2}\), which fails for \(m\ge2\).  It remains a
scalar choice, not a shell or PDE realization.

This family has an explicit collar margin, not just exponent notation.
The elementary inequality

\[
             \log n\le {3\over2}n^{1/4}\quad(n\ge1)       \tag{3.9}
\]

follows because the maximum of \((\log x)x^{-1/4}\) is \(4/e<3/2\),
and \(e>65/24>8/3\).  Consequently

\[
 {\varepsilon_j^{\rm col}\over b_j^3}
 \le K_m n^{-9/8},\qquad
 K_m=C_{\rm col}\left({3C_JQ_m\over2}\right)^{7/2}.     \tag{3.10}
\]

For the threshold explicit in the declared constants

\[
 N_m=\left\lceil
 C_{\rm col}^{\,8}\left({3C_JQ_m\over2}\right)^{28}
 \right\rceil,                                           \tag{3.11}
\]

every \(n\ge\max(2,N_m)\) satisfies

\[
             \varepsilon_j^{\rm col}\le n^{-1}b_j^3.    \tag{3.12}
\]

The landed C127 schedule is exactly \(m=1\):
\((Q,S,P)=(8,2,10)\) and \(\gamma=5/4\).  The first strictly lower-demand
member of the particular family (3.5) is \(m=2\):

\[
 (Q,S,P)=(14,4,17),\qquad \gamma={17\over14}<{5\over4}.  \tag{3.13}
\]

Increasing \(m\) lowers the exponent but makes the one-step frequency
ratio more severe.  C188 certifies the family; it does not assert that
the \(m=2\) geometry has a smaller full UVSR residual.

### 3.2 The unattained boundary and its logarithmic price

The strictness in (3.3) cannot be removed for a pure power.  A
boundary-approaching non-power seed makes the price visible.  For \(q>e\),
choose

\[
 b(q)=K_*^{1/3}q^{-1/3}(\log q)^{7/6}\log\log q.          \tag{3.14}
\]

Then exactly

\[
 {K_*q^{-1}(\log q)^{7/2}\over b(q)^3}
 ={1\over(\log\log q)^3}\longrightarrow0,               \tag{3.15}
\]

while the net same-energy demand is

\[
 g(q)=b(q)q^{3/2}
 =K_*^{1/3}q^{7/6}(\log q)^{7/6}\log\log q.              \tag{3.16}
\]

Thus \(7/6\) is the exact polynomial infimum for saturation of the
declared worst-case envelope, with this sufficient slowly-diverging margin.
It is not a physical necessity: C176 supplies neither a lower bound on
\(J\) nor a matching lower bound on the actual collar response.

### 3.3 Preserve the landed C180 shell while lowering demand

The family (3.5) is a scalar sharpness family.  C180 has only proved its
factorial-shell arithmetic for \(Q=8\); changing \(Q\) would require a new
shell theorem.  The demand can be lowered without making that change.
Keep

\[
 q=n^8,\qquad b=n^{-S},\qquad
 g=n^{12-S}=q^{\,3/2-S/8}.                               \tag{3.17}
\]

For every \(2\le S<8/3\), the parent feedback obeys (2.4), Reynolds grows,
and the C176 collar has the strict power margin

\[
 {q^{-1}(\log q)^{7/2}\over b^3}
 =8^{7/2}(\log n)^{7/2}n^{\,3S-8}\longrightarrow0.       \tag{3.18}
\]

For example, the exact choice \(S=5/2\) keeps the C180 frequency schedule
and gives

\[
\begin{gathered}
 a_j=(j!)^{19/2},\qquad b_j=n^{-5/2},\qquad
 g_j=n^{19/2}=q_j^{19/16},\\
 E_j=(j!)^{-5},\quad \tau_j=(j!)^{-35/2},\quad
 {\rm Re}_j=\nu^{-1}(j!)^{3/2},\quad
 \mu_j=\nu(j!)^{-3/2},\quad D_j=\nu(j!)^{-13/2}.          \tag{3.19}
\end{gathered}
\]

The inherited \(L_j\le n^2\) chart obeys
\(L_jb_j\le n^{-1/2}\), the retained-wake scale is
\(b_j^2=n^{-5}\), and the active tolerance is
\(b_j^3=n^{-15/2}\).  Equation (3.18) is then the explicit comparison
\[
 {q^{-1}(\log q)^{7/2}\over b^3}
 =8^{7/2}(\log n)^{7/2}n^{-1/2}\longrightarrow0.         \tag{3.20}
\]

Indeed, \(\log n\le6n^{1/14}\) because the maximum of
\((\log x)x^{-1/14}\) is \(14/e<6\).  Therefore every
\(n\ge48^{28}\) satisfies the fully numerical bound

\[
 {q^{-1}(\log q)^{7/2}\over b^3}\le n^{-1/8}.             \tag{3.21}
\]

This seed change requires an explicit respecification of C161's frozen
normalization.  Keep its preparation size \(\varepsilon=n^{-28}\), but
replace \(H=n^{26}\) by
\[
 H=n^{51/2},\qquad c_{\rm seed}={\varepsilon\over q}=n^{-36},
 \qquad c_0=Hc_{\rm seed}=n^{-21/2}={b\over q}.           \tag{3.22}
\]
After the normalized star, \(c_1=c_0/\sqrt q=n^{-29/2}\), so
\[
 qc_0=b,\qquad \sqrt{q^3}\,c_1=b,\qquad q^3c_1=n^{19/2}=g. \tag{3.23}
\]
The split count must also be integral.  Set
\[
 J_{\rm split}=\lceil n^{5/2}\rceil,\qquad
 \theta=J_{\rm split}^{-1}.                              \tag{3.24}
\]
Since \(n^{5/2}\le J_{\rm split}\le n^{5/2}+1\le2n^{5/2}\),
\[
 {b\over2}\le\theta\le b,\quad J_{\rm split}\theta=1,
 \quad J_{\rm split}b\theta^2={b\over J_{\rm split}}\le b^2,
 \quad J_{\rm split}b\theta^3={b\over J_{\rm split}^2}\le b^3. \tag{3.25}
\]
With C127's old \(L=n^2,M=n^{7/2}\), the remaining scalar comparisons
also survive.  Besides \(Lb=n^{-1/2}\), (3.9) gives
\[
 {\log q\over M}\le12n^{-3/4}b,
 \qquad {L\log q\over M}\le12n^{-5/4},
 \qquad \sum_{n=2}^{\infty}{L\log q\over M}\le48.        \tag{3.26}
\]

This lowers the net exponent from \(5/4\) to \(19/16\) without altering
C180's \(q=n^8\) arithmetic, but it does alter C161's seed and split
normalization exactly as in (3.22)--(3.25).  The logarithmically corrected
boundary choice (3.14), also with \(q=n^8\), approaches the envelope's
\(7/6\) polynomial infimum.  These are scalar admissibility statements;
all physical UVSR residual and chart constants remain to be certified.

## 4. What epsilon regularity actually tests

Normalize a physical viscosity \(\nu>0\) by \(s=\nu t\), \(v=u/\nu\), and
\(\pi=p/\nu^2\).  For a backward physical cylinder

\[
 Q_r^\nu(z_0)=B_r(x_0)\times(t_0-r^2/\nu,t_0),           \tag{4.1}
\]

the standard dimensionless quantities include

\[
\begin{aligned}
 A_\nu(r;z_0)
 &= {1\over\nu^2r}\mathop{\rm ess\,sup}_{t_0-r^2/\nu<s<t_0}
       \int_{B_r(x_0)}|u(x,s)|^2\,dx,\\
 C_\nu(r;z_0)
 &= {1\over\nu^2r^2}\int_{Q_r^\nu(z_0)}|u|^3,\\
 D_\nu(r;z_0)
 &= {1\over\nu^2r^2}\int_{Q_r^\nu(z_0)}
       |p-(p)_{B_r}|^{3/2},\\
 \mathcal E_\nu(r;z_0)
 &= {1\over\nu r}\int_{Q_r^\nu(z_0)}|\nabla u|^2.         \tag{4.2}
\end{aligned}
\]

The original CKN theorem proves partial regularity for suitable weak
solutions.  In one standard CKN--Lin velocity--pressure formulation there
is an existential universal \(\varepsilon_{CD}>0\) such that small
\(C_\nu(r)+D_\nu(r)\) gives regularity on a smaller cylinder.  Its
contrapositive says that at a singular point

\[
       C_\nu(r;z_*)+D_\nu(r;z_*)\ge\varepsilon_{CD}         \tag{4.3a}
\]

for every sufficiently small \(r\).  A dissipation formulation is weaker
in scale: it gives
\(\limsup_{r\downarrow0}\mathcal E_\nu(r;z_*)\ge
\varepsilon_{\mathcal E}\), hence bad arbitrarily small radii, not every
small radius.  The threshold depends on the chosen functional; it is
universal but is not supplied as a certified decimal in these sources. See
[Caffarelli--Kohn--Nirenberg](https://doi.org/10.1002/cpa.3160350604) and
[Lin](https://doi.org/10.1002/(SICI)1097-0312(199803)51:3%3C241::AID-CPA2%3E3.0.CO;2-A)
for the velocity--pressure framework.  As a distinct modern example,
[Albritton--Barker--Prange, Theorem B](https://arxiv.org/abs/2211.16188)
gives a one-scale \(L^4\) criterion in a first-time-singularity setting;
it is not being identified with (4.3a).

Both contrapositives are centered and terminal.  Neither applies to an
arbitrary packet center at a smooth time \(t_j<T_*\).

### 4.1 Explicit smooth-point destruction test

Suppose \(u,p\) are smooth on a neighborhood of
\(\overline{Q_R^\nu(z_0)}\), and after fixing one pressure gauge let

\[
 M_u=\sup_{Q_R^\nu}|u|,\qquad
 M_p=\sup_{Q_R^\nu}|p|,\qquad
 M_\nabla=\sup_{Q_R^\nu}|\nabla u|.                       \tag{4.3}
\]

For every \(0<r<R\), direct integration gives the explicit bounds

\[
\begin{aligned}
 A_\nu(r;z_0)
 &\le {4\pi\over3}\left({M_u\over\nu}\right)^2r^2,\\
 C_\nu(r;z_0)+D_\nu(r;z_0)
 &\le {4\pi\over3}
 \left[\left({M_u\over\nu}\right)^3
 +2^{3/2}\left({M_p\over\nu^2}\right)^{3/2}\right]r^3,\\
 \mathcal E_\nu(r;z_0)
 &\le {4\pi\over3}\left({M_\nabla\over\nu}\right)^2r^4.
                                                               \tag{4.4}
\end{aligned}
\]

Every quantity in (4.4) falls below every fixed positive epsilon at small
enough \(r\).  This does not mean viscosity has erased the packet; it
certifies the smoothness already assumed.  Hence the proposed rule
“every intermediate scale must stay above \(\varepsilon_0\) at its smooth
stage center” fails its own destruction test.  This does not control a
natural radius \(r\asymp\ell_j\) in a backward cylinder centered at the
eventual singular point: such a cylinder includes later times, and the
smooth bounds in (4.4) are not uniform in \(j\).  Terminal-center tracking
and scale occupancy therefore remain open.

For comparison, instantaneous concentration theorems require additional
center, radius, and sometimes type-I hypotheses; see
[Kang--Miura--Tsai, Theorem 1.6](https://arxiv.org/abs/2006.13145).
Their existential constants likewise do not provide a numerical UVSR
window.

### 4.2 The local energy inequality has no schedule-only upper wall

For a nonnegative cutoff \(\phi\), a suitable unforced solution obeys

\[
\begin{aligned}
 \int |u(t)|^2\phi(t)+2\nu\int_0^t\!\int|\nabla u|^2\phi
 \le{}&\int |u(0)|^2\phi(0)
 +\int_0^t\!\int |u|^2(\partial_s\phi+\nu\Delta\phi)\\
 &+\int_0^t\!\int (|u|^2+2p)u\cdot\nabla\phi.             \tag{4.5}
\end{aligned}
\]

The last term has no sign and includes inward pressure/advective flux.
For a packet of amplitude \(a\), radius \(\ell\), and one turnover duration
\(\tau=\ell/a\), its advective flux scale is exactly

\[
        a^3\ell^{-1}\ell^3\tau=a^2\ell^3,                \tag{4.6}
\]

the packet energy itself.  Global energy \(E_0=\sup_t\int|u(t)|^2\) gives
only

\[
                  A_\nu(r;z_0)\le {E_0\over\nu^2r},       \tag{4.7}
\]

which worsens as \(r\downarrow0\).  A genuine local upper corridor would
need certified initial local energy, pressure oscillation, boundary flux,
time occupancy, and wake bounds.  None follows from (4.5) alone.

The same caveat applies to a numerical approximate UVSR profile: the
unforced epsilon criteria apply to an exact suitable solution.  A forced
or residual-stable version must include the full residual explicitly.

## 5. Tao's critical-norm rate adds no fixed power

For viscosity one on \(\mathbb R^3\), Tao proves that a classical unforced
solution with finite blow-up time \(T_*\) satisfies

\[
 \limsup_{t\uparrow T_*}
 {\|u(t)\|_{L^3(\mathbb R^3)}
  \over
  \left(\log\log\log{1\over T_*-t}\right)^c}=\infty       \tag{5.1}
\]

for an unspecified absolute \(c>0\); see
[Tao, *Quantitative bounds for critically bounded solutions to the
Navier--Stokes equations*](https://arxiv.org/abs/1908.04958).  The
regularity constants have triple-exponential dependence with untracked
absolute exponents.  Tao notes that one may take an internal large constant
\(C_0=10^5\), but this does not turn \(c\) or the final theorem into a
certified numerical threshold.

For physical viscosity \(\nu\), the normalization used in Section 4 gives
the corresponding statement

\[
 \limsup_{t\uparrow T_*}
 {\,\|u(t)\|_{L^3(\mathbb R^3)}/\nu
  \over
  \left(\log\log\log{1\over \nu(T_*-t)}\right)^c}
 =\infty,                                                \tag{5.1a}
\]

after fixing the same dimensionless time convention.  Fixed \(\nu>0\)
does not change the schedule exponent.

For a fixed nondegenerate packet shape, the active contribution scales as

\[
             \|u_j\|_3=c_{\rm sh}a_j\ell_j
             =c_{\rm sh}(j!)^{Q/2-S}.                     \tag{5.2}
\]

Every fixed-power schedule with \(Q>2S\) grows factorially.  Along the
nominal turnover tail defined below, it dominates Tao's triple-logarithmic
scale.  In the legacy family (3.5), (5.2) is
\(c_{\rm sh}(j!)^{m+1}\).  This compatibility calculation supplies no
positive numerical gap above the high-Re boundary \(\gamma=1\).

More explicitly, set \(s=5Q/2-S>0\) and define the **nominal**
one-turnover tail
\(\Delta_j^{\rm nom}=\sum_{k=j}^\infty(k!)^{-s}\).  Since
\(\tau_{k+1}/\tau_k=(k+1)^{-s}\le2^{-s}\),

\[
 (j!)^{-s}\le\Delta_j^{\rm nom}
 \le{(j!)^{-s}\over1-2^{-s}}.                            \tag{5.3}
\]

Consequently, for every fixed \(c>0\),

\[
 { (j!)^{Q/2-S}
  \over[\log\log\log(1/\Delta_j^{\rm nom})]^c}
  \longrightarrow\infty.                                \tag{5.4}
\]

The limit is elementary; it does not supply a threshold index because the
source theorem does not specify \(c\).

The landed stages provide only an upper charge by a fixed multiple of
\((\log q_k)\tau_k\), which is summable in the certified \(S\ge2\) class.
That upper bound does not give a lower bound on the actual remaining time,
so C188 does not replace \(\Delta_j^{\rm nom}\) by the physical blow-up
tail.  A two-sided stage-duration comparison is part of the open
stage-time-coverage problem.

This comparison is only a compatibility check.  Tao's theorem is global,
unforced, is stated on \(\mathbb R^3\), and holds along a time subsequence
of an actual blowing-up solution.  A retained wake may carry its \(L^3\)
norm, and the distinguished times need not be stage endpoints.  A
per-return inequality would require a torus/domain transfer, a
residual-stable forced extension for approximate profiles, stage-time
coverage, and two-sided active-plus-wake \(L^3\) shape bounds.  C188 proves
none of those additional hypotheses.

## 6. Corridor verdict and next test

The requested CKN/Tao corridor is not empty at the level of the landed
scalar ledger.  The strongest numerical statement justified there is the
bounded-profile legacy-collar envelope window (0.1), with the exact family
(3.5) and the explicit majorant threshold (3.11).  No further numerical
power restriction follows from the cited regularity theory without the
open singular-center, time-occupancy, pressure, and wake hypotheses.  The
smooth-center calculation does refute the proposed universal
stagewise-\(\varepsilon_0\) rule.

Accordingly:

* UVSR is no longer hard-coded to the C127 net exponent \(5/4\).  While the
  bounded-profile C176 worst-case majorant is paid separately, its net
  same-energy target is any \(\gamma\in(7/6,3/2)\), with the sufficient
  boundary correction (3.16).  The fixed-shell point \(\gamma=19/16\) in
  (3.19) is one convenient exact lower-demand specification that does not
  require a new C180 shell theorem; it does require the explicit C161
  respecification (3.22)--(3.25).
* The equal-shape ledger focus remains \(q^{3/2}\); the physical formula is
  (1.3a).  Calling a smaller power a schedule optimization would require
  either nonuniform shape constants or an explicit energy disposal/wake
  obligation.
* The broader direct specification \(g=q\rho\), \(\rho=bq^{1/2}>1\),
  becomes admissible only after a direct UVSR full-residual certificate
  removes the separate C176 comparison.  Its normalized pure-power window
  is \((1,3/2)\).
* The factorial stage operator is the nonautonomous augmented map (1.5).
  A fixed positive normalized viscosity is incompatible with a fixed
  supercritical scaling, while an autonomous inviscid limiting map remains
  an additional open convergence problem.
* No profile search is justified by an alleged explicit CKN epsilon
  number.  The next pre-registered computation remains the physical PPRG
  two-episode witness test.

## Claim boundary

* **PROVED:** (1.2)--(1.4) and the necessary nonautonomous coordinate form
  in (1.5), the generalized factorial ledger (2.2)--(2.5),
  the bounded-profile worst-case-envelope infimum \(7/6\), the direct
  infimum \(1\), the exact families (3.5) and (3.17), their scalar ledgers,
  the polylogarithmic envelope identity (3.14)--(3.16), the fixed-scaling
  positive-viscosity incompatibility, and the smooth-point bounds (4.4).
* **SOURCE-SCOPED:** the stated CKN/Lin epsilon-regularity implication and
  Tao's \(\mathbb R^3\) subsequential \(L^3\) lower rate.
* **REFUTED:** using epsilon regularity to demand a uniform lower threshold
  at every smooth intermediate stage center.
* **OPEN:** a structured state and well-defined \(\mathcal R_{n,\mu}\), a
  normalized autonomous limit of the factorial stage maps, a UVSR profile,
  a residual-stable forced epsilon criterion,
  torus transfer of Tao's quantitative theorem, singular-center tracking,
  pressure/flux/wake bounds, the PPRG witness, the three-loss viscous
  bridge, nonlinear trapping, and a Navier--Stokes singularity.
