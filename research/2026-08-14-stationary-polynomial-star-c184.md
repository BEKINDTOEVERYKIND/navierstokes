# C184: the stationary polynomial row has an exponential collateral face

**Date:** 2026-08-14

**Status:** exact Fourier-range characterization, diagonal interpolation,
selected-row tangent Gram, zero-loop jet, outer-face collateral lower bound,
and initial viscous-drift ledger; no finite-frequency packet, acceptable
reservoir, physical focus, or one-cell stage

**Checker:**
[checks/stationary_polynomial_star_c184.py](../checks/stationary_polynomial_star_c184.py)

## 0. Claim boundary

C183 shows that every stationary inviscid scalar

\[
                         \Theta=H(f),                     \tag{0.1}
\]

on the C152 planar pump admits an exact common return frame after one
covector retuning.  This note asks whether the polynomial freedom in
\(H\) can realize C179's planar translation star without leaving that
frame.

The answer has two distinct parts.

1. **The selected tangent row is algebraically possible.**  If
   \(\deg H\le q\), the coefficients on the diagonal frequencies
   \(r(r_1+r_2)\) form triangular coordinates.  A fixed positive fraction
   of the radii \(r\asymp q\) can be prescribed at the normalized
   \(q^{-3/2}\) scale.  On every C176 source sub-sheet pruned into C179's
   compact separated cone, their full-polarization tangent Gram is
   uniformly positive and bounded.  Thus this route does not die by rank.
2. **The complete reservoir is exponentially too large.**  A nonzero
   outer diagonal coefficient \(\tau_q\) forces, with no possible
   lower-degree cancellation, an off-line binomial face whose
   coefficient-\(\ell^2\) size is at least

   \[
       \boxed{
       |\tau_q|{(9/4)^q-1\over\sqrt q}.}                  \tag{0.2}
   \]

   For the normalized star scale \(|\tau_q|\asymp bq^{-3/2}\), this is
   \(\gtrsim b(9/4)^q/q^2\), exponentially above every polynomial
   reservoir/wake budget.

Moreover, fixing \(H(0)=0\) and \(H'(0)=c\) makes the nonlinear part
\(H(f)-cf\) vanish to second order on the zero loop.  It contributes
neither velocity nor strain to the central principal Kelvin cocycle.  Its
bright Fourier row is therefore a finite-frequency decomposition of a
globally huge passive scalar, not an additional principal Floquet gain.

The exact obstruction is called **stationary-polynomial collateral
explosion (SPCE)**.  It rules out the natural degree-\(q\), line-shift,
polynomial-norm completion of PPRG.  It does not rule out a nonpolynomial
profile, a palette without an outer vertex coefficient, deliberate use of
the collateral modes, or a finite-frequency mechanism whose reservoir
budget is paid by some other construction.

## 1. Exact coefficient space

Use phase coordinates

\[
 a=r_1\cdot x,\qquad b=r_2\cdot x,qquad d=(1,1),
\]

and write a phase-lattice frequency as \((j,k)\in\mathbb Z^2\).  The C152
Hamiltonian is

\[
 f(a,b)=\cos a+\cos b+{4\over5}\cos(a+b).                \tag{1.1}
\]

Its six nonzero Fourier coefficients form the finitely supported measure

\[
 \mu_{\pm(1,0)}=\mu_{\pm(0,1)}={1\over2},\qquad
 \mu_{\pm(1,1)}={2\over5}.                               \tag{1.2}
\]

Let \(C_n=\mu^{*n}\), with \(C_0=\delta_0\).  For

\[
                         H(t)=\sum_{n=0}^q h_nt^n,         \tag{1.3}
\]

the coefficient vector of \(H(f)\) is exactly

\[
                    \boxed{\widehat{H(f)}=\sum_{n=0}^q h_nC_n.} \tag{1.4}
\]

The map \(H\mapsto H(f)\) is injective: the continuous nonconstant
function \(f\) has an interval as its range, so a polynomial vanishing on
that range is zero.  Hence

\[
 \mathcal V_q=\operatorname{span}_{\mathbb R}\{C_0,\ldots,C_q\}
 \quad\hbox{has dimension }q+1.                           \tag{1.5}
\]

The conditions in C183 are simply

\[
                         h_0=0,\qquad h_1=c.               \tag{1.6}
\]

Thus the admissible coefficient vectors form the affine space

\[
                   cC_1+\operatorname{span}\{C_2,\ldots,C_q\}, \tag{1.7}
\]

of dimension \(q-1\).  Equivalently,

\[
                         H(t)-ct=t^2S(t).                  \tag{1.8}
\]

All coefficients are real and obey

\[
 \widehat{H(f)}_{j,k}=\widehat{H(f)}_{k,j}
 =\widehat{H(f)}_{-j,-k}=\widehat{H(f)}_{-k,-j}.          \tag{1.9}
\]

The support lies in the hexagon

\[
                 \max\{|j|,|k|,|j-k|\}\le q.             \tag{1.10}
\]

These are strong coefficient correlations: there are \(O(q^2)\) available
Fourier sites but only \(q-1\) free real parameters.

## 2. Diagonal interpolation is triangular

Put \(\alpha=2/5\).  Reaching \((r,r)\) in fewer than \(r\) steps is
impossible, while the only \(r\)-step path is \(r\) copies of \((1,1)\).
Therefore

\[
                    C_n(r,r)=0\quad(n<r),\qquad
                    C_r(r,r)=\alpha^r.                    \tag{2.1}
\]

It follows that

\[
 (h_2,\ldots,h_q)\longmapsto
 (\widehat{H(f)}_{(2,2)},\ldots,
  \widehat{H(f)}_{(q,q)})                                 \tag{2.2}
\]

is an upper-triangular isomorphism.  More locally, choose
\(r_0\ge2\), set \(h_2=\cdots=h_{r_0-1}=0\), and prescribe arbitrary real
numbers \(\tau_r\) for \(r_0\le r\le q\).  Back substitution uniquely
chooses \(h_{r_0},\ldots,h_q\) so that

\[
                       \widehat{H(f)}_{(r,r)}=\tau_r.      \tag{2.3}
\]

The fixed value \(h_1=c\) does not enter these high coefficients.
Reality automatically supplies the negative frequencies.

Consequently the line palette

\[
 \mathcal G_q=\{\pm r(r_1+r_2):\lfloor q/2\rfloor<r\le q\} \tag{2.4}
\]

has \(q+O(1)\) signed modes whose coefficients can be prescribed at the
first-row level.  Constant-fraction subbands of (2.4) may be used to avoid
C179's equal-radius determinant surface.

### 2.1 The selected tangent Gram really is a star

Let \(n\) be the unit vertical vector, \(g\perp n\), \(m=n\cdot p\ne0\),
and \(a\in p^\perp\).  Omitting the common factor \(-i\), one scalar
coefficient \(\tau_g\) produces the exact tangent edge

\[
 E_{p,g}a=\tau_gP_{p+g}\{ma+(a\cdot g)n\}.                \tag{2.5}
\]

On a fixed compact set satisfying C179's separated-cone conditions,
with \(|p|,|g|\asymp q\), the determinant formula in C179 and compactness
give

\[
             c q|a|\le |E_{p,g}a|/|\tau_g|
                    \le Cq|a|.                            \tag{2.6}
\]

If \(Q=|{\cal G}|\asymp q\) and

\[
                         |\tau_g|={\lambda\over q\sqrt Q}, \tag{2.7}
\]

then the selected collective row has

\[
 \boxed{
 c^2\lambda^2I_{p^\perp}
       \le\sum_{g\in{\cal G}}E_{p,g}^*E_{p,g}
       \le C^2\lambda^2I_{p^\perp}.}                    \tag{2.8}
\]

This is the normalized full-polarization tangent Gram.  The checker also
gives an exact rational example where the reality-paired Gram eigenvalues
lie between \(7/3\) and \(4\).  Thus SPCE is not a disguised rank defect.

Equation (2.8) is only the row selected from (1.4).  The complete
\(H(f)\) contains every collateral coefficient forced by the same
polynomial, and those edges cannot be deleted from the Navier--Stokes
evolution.

There is also no hidden cancellation inside the instantaneous sourcewise
Gram.  For one input frequency \(p\), different shifts \(g\) land at
orthogonal Fourier outputs \(p+g\), so the complete tangent Gram is

\[
       K_p=\sum_g|\widehat{H(f)}_g|^2A_{p,g}^*A_{p,g}.    \tag{2.9}
\]

Every collateral summand is positive semidefinite.  In fact, for a
degree-\(q\) polynomial, the forced face contains \(g_0=qr_1\) and
\(g_1=qr_1+r_2\), with distinct squared lengths \(2q^2\) and
\(2(q^2-q+1)\).  A C176 source has nonzero vertical charge, so
\(p+g_i\ne0\); C179's determinant formula shows at most one of these two
edges can lie on the equal-radius rank-loss surface.  The other is
invertible on both polarizations.  Thus the collateral is dynamically
visible already in the tangent row, though a later time-ordered propagator
could in principle recycle it.

## 3. Exact exponential collateral face

Let \(H\) have degree at most \(q\), and suppose its outer diagonal
coefficient is

\[
                         \tau_q=\widehat{H(f)}_{(q,q)}\ne0. \tag{3.1}
\]

Here and below \(q\ge2\).  In particular, the fixed linear term \(cf\)
has no coefficient on the face with first phase coordinate \(q\), so it
cannot cancel any face coefficient used below.

Only \(h_qf^q\) reaches the face whose first phase coordinate is \(q\).
To end at \((q,s)\), \(0\le s\le q\), every step must be either
\((1,0)\), of weight \(\beta=1/2\), or \((1,1)\), of weight
\(\alpha=2/5\).  Hence, exactly,

\[
 \widehat{H(f)}_{(q,s)}
       =h_q{q\choose s}\beta^{q-s}\alpha^s,
 \qquad
 \tau_q=h_q\alpha^q.                                    \tag{3.2}
\]

Every coefficient in (3.2) has the same sign.  The off-diagonal members
\(0\le s<q\) therefore obey

\[
 \begin{aligned}
 \sum_{s=0}^{q-1}|\widehat{H(f)}_{(q,s)}|
   &=|\tau_q|\left\{\left({\alpha+\beta\over\alpha}\right)^q-1\right\}\\
   &=|\tau_q|\{(9/4)^q-1\}.                              \tag{3.3}
 \end{aligned}
\]

There are \(q\) terms.  Cauchy--Schwarz and normalized Parseval prove

\[
 \boxed{
 \|H(f)-cf\|_2
 \ge\left(\sum_{s=0}^{q-1}
       |\widehat{H(f)}_{(q,s)}|^2\right)^{1/2}
 \ge |\tau_q|{(9/4)^q-1\over\sqrt q}.}                  \tag{3.4}
\]

No coefficient \(h_n\), \(n<q\), can cancel this face.  The same bound
holds for \(L^\infty\), since \(\|F\|_\infty\ge\|F\|_2\) on the
normalized torus.

At the normalized star scale \(Q\asymp q\) and
\(|\tau_q|\asymp b/(q\sqrt Q)\), (3.4) becomes

\[
                         \boxed{
 \|H(f)-cf\|_2\gtrsim {b\over q^2}(9/4)^q.}              \tag{3.5}
\]

Thus no polynomial reservoir budget \(bq^A\), with fixed \(A\), can
carry this selected star.  Conversely, imposing such a budget makes the
outer selected coefficient, and hence any selected row with comparable
coefficients, exponentially smaller than (2.7).

An exact line-supported repair is impossible even without the quantitative
bound.  If a nonconstant \(H(f)\) had Fourier support in one rank-one
subgroup, it would depend on one linear phase.  Differentiation in the
orthogonal direction would give \(H'(f)w\cdot\nabla f=0\).  On an open
set \(H'(f)\ne0\), so analyticity would force \(w\cdot\nabla f\equiv0\),
contradicting the two noncollinear frequencies in (1.1).

## 4. The high polynomial freedom has zero principal jet

Let

\[
                         R(f)=H(f)-cf=f^2S(f).             \tag{4.1}
\]

On the zero loop \(f=0\),

\[
                         R(f)=0,\qquad \nabla R(f)=0.      \tag{4.2}
\]

Therefore adding \(R(f)n\) changes neither the trajectory nor the velocity
gradient along the exact C159 orbit.  The covector resonance and the
principal Kelvin connection are exactly those of the linear scalar
\(cf\).  The differential from the \(q-1\) high polynomial parameters to
the central velocity/strain jet is zero; only varying \(c\) supplies one
jet direction.

This zero jet is not a degeneracy of the global coefficient coordinates.
For parameters \(x_2,\ldots,x_q\), Parseval gives the global Gram

\[
 G_{rs}=\langle f^r,f^s\rangle_{L^2}
       =\int f^{r+s},\qquad 2\le r,s\le q.                \tag{4.3}
\]

It is positive definite because
\(x^TGx=\int|\sum_{r=2}^qx_rf^r|^2\), and the range of the nonconstant
continuous function \(f\) contains an interval.  Thus the admissible
fixed-\(c\) family has a full-rank global \(L^2\) tangent Gram but an
identically zero first-jet Gram on \(f=0\).

This reconciles (2.8) with C183.  Equation (2.8) isolates first-generation
global Fourier daughters.  The exact Lagrangian gauge sees their complete
sum as multiplication by \(R(f)\), which has zero first jet on the central
orbit and comes with the exponentially larger collateral face (3.4).
Counting the selected daughters while discarding that face would violate
the exact passive evolution.

## 5. Exact initial viscous-drift ledger

In parent-rescaled variables, let the passive scalar solve

\[
 \partial_t\Theta+A(t)v_h\cdot\nabla\Theta=\mu\Delta\Theta,
 \qquad \Theta(0)=H(f).                                  \tag{5.1}
\]

Since \(v_h\cdot\nabla H(f)=0\),

\[
                         \partial_t\Theta(0)=\mu\Delta H(f). \tag{5.2}
\]

The physical wavevector represented by \((q,s)\) is
\(qr_1+sr_2\), whose squared norm is

\[
                  2(q^2+s^2-qs)\ge{3\over2}q^2.          \tag{5.3}
\]

Combining (3.4), (5.2), and (5.3) gives the exact initial-rate bound

\[
 \boxed{
 \|\partial_t\Theta(0)\|_2
 \ge {3\over2}\mu q^2|\tau_q|
          {(9/4)^q-1\over\sqrt q}.}                     \tag{5.4}
\]

This uses only mutually orthogonal Fourier modes on the outer face.  The
term \(cf\) has no such mode because \(q\ge2\), and modes off that face
cannot cancel them in \(L^2\).

For \(|\tau_q|\asymp bq^{-3/2}\), this is
\(\gtrsim\mu b\{(9/4)^q-1\}\).  On C127's schedule

\[
 q=(j+1)^8,\qquad b=(j+1)^{-2},\qquad
 \mu=\nu(j!)^{-2},                                      \tag{5.5}
\]

the right side diverges super-exponentially relative to the factorial
denominator.  Thus factorial viscosity does not make this exponentially
large stationary polynomial into a slowly drifting reservoir.  Equation
(5.4) is an initial-rate statement, not a lower bound for the solution
difference at a later fixed time.

## 6. What remains

C184 closes one tempting use of C183's common-return family:

> A degree-\(q\) stationary polynomial whose outer vertex-line
> coefficients form a normalized \(q\)-star necessarily carries an
> exponentially larger off-line reservoir and an exponentially large
> normalized initial viscous drift.  It cannot meet the polynomial
> C161/C176/C179 reservoir and wake ledgers.

The surviving same-geometry questions are narrower:

1. can a nonpolynomial \(H(f)\) or an exactly evolved passive scalar
   concentrate a useful translation palette without SPCE;
2. can the collateral face itself be used coherently, rather than treated
   as leakage, while retaining C176 injectivity and full-polarization
   control; or
3. can a non-returning PPRG profile exploit C183's explicit affine/quadratic
   covector drift and still pay finite-frequency, localization, heat,
   C125, RIGM, and BAFL?

No new geometry, one-cell stage, or singularity is claimed.

## 7. Verification boundary

The dependency-free checker constructs the Laurent convolution powers of
\(f\) with exact rational arithmetic; verifies the hexagonal support,
symmetries, triangular diagonal interpolation, and outer-face binomial
formula; checks the exact selected-row Gram in a separated model cone and
the positive sourcewise nature of collateral contributions;
proves the collateral \(\ell^1/\ell^2\) ledger and Laplacian lower bound on
finite cases; and verifies the factorial-versus-\((9/4)^q\) schedule
comparison.  It does not prove a localized finite-frequency parametrix,
an approximate nonpolynomial no-go, or any one-cell obligation.
