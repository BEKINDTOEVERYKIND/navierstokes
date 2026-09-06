# C200: smooth finite-horizon A2 polarization filters

**Date:** 2026-09-06

**Status:** explicit-constant principal-symbol theorem on the C195 finite-horizon tube. No invariant splitting, infinite-horizon limit, finite-frequency PDE solution, energy balance, retained band, viscosity, nonlinear return, or singularity is asserted here.

**Checker:** [a2_smooth_finite_horizon_filter_c200.py](../checks/a2_smooth_finite_horizon_filter_c200.py).

## 1. Data, normalization, and result

Use C195's orbit, period `T<76/25`, and any positive integer block count `R` for which every block is in its certified tube. Put `H=RT`, `L=1+H`. Work in a common initial-position lift `y`. Let `K=|k_{0,*}|`, so `5<K<6`, and write the transported physical covector as

\[
 k(t,y,\xi)=K D\Phi_t(y)^{-T}\xi . \tag{1}
\]

The real normalized initial covector `xi` is near `k_{0,*}/K`. The domain consists of **all real `(y,xi)` whose entire R-block trajectory satisfies C195's tube**; it may in particular be obtained by C195's sufficient initial-radius bound. Derivative assertions hold on its interior and extend to compact subdomains by continuity. They do not assume an invariant open set for an infinite-time dynamics.

Let `E_j=[E_1(k_j),E_2(k_j)]` be C195's physical frame at `jT`, and let `M_j` map frame coefficients from `(j-1)T` to `jT`, including the changing frame. C195 gives, with `I=[137/1000,1/5]` and `J=diag(1,-1)`, for both `M_j` and `J M_j^{-1}J`:

\[
 f_M(s)=\frac{M_{21}+M_{22}s}{M_{11}+M_{12}s}\in I,
 \quad g_M(s)=M_{11}+M_{12}s>3000,
 \quad 0<f_M'(s)<\frac1{4500000}. \tag{2}
\]

The interval images are strictly interior by C195. Its determinant bound is `0<det M<2`. Fix `s_*=17/100`. Define forward slopes `s_0=s_*`, `s_j=f_{M_j}(s_{j-1})`, and backward slopes `r_R=s_*`, `r_{j-1}=f_{J M_j^{-1}J}(r_j)`. These are finite-horizon selections, not canonical bundles.

There are real transverse Kelvin solutions `b_+(j;y,xi)` and `b_-(j;y,xi)` satisfying

\[
 |b_+(R)|=|b_-(0)|=1,\qquad
 |b_+(j)|\le4\,3000^{-(R-j)},\qquad
 |b_-(j)|\le4\,3000^{-j}. \tag{3}
\]

Their forward and backward lines have the C195 angle and projector bounds `sin(angle)>=520/569` and `||P||<=569/520` at every comparison section. Thus the same tube supplies a **terminally normalized expanding family and an initially normalized contracting family**.

Set

\[
 D_*=4\cdot10^{40}L,\qquad
 d_+(j)=2999^{-(R-j)},\quad d_-(j)=2999^{-j}. \tag{4}
\]

For either sign, every nonnegative integer `r` (in particular `r<=4` required for C197 synthesis), and all comparison sections,

\[
 \boxed{\|D_\xi^r b_\pm(j)\|\le30\,r!D_*^r d_\pm(j)}, \tag{5}
\]

\[
 \boxed{\|D_yD_\xi^r b_\pm(j)\|
 \le10^{45}(R+1)L^4\,r!D_*^r d_\pm(j)}. \tag{6}
\]

A second spatial derivative is also controlled:

\[
 \boxed{\|D_y^2D_\xi^r b_\pm(j)\|
 \le10^{90}(R+1)^2L^8\,r!D_*^r d_\pm(j)}. \tag{6a}
\]

These are Euclidean multilinear operator norms. In particular all entrance and endpoint symbol jets have **polynomial horizon costs**; the expanding entrance and contracting endpoint also retain exponential smallness. No factor `exp(cH)` is hidden in these constants.

## 2. One-block estimates on complex covector neighborhoods

Only the covector is complexified. The base `y` and its flow stay real, so no complexified long-time flow estimate is used. C195 and C193 imply throughout the real tube

\[
 3<|k|<14;\qquad 3<|P_nk_j|<4,\quad |m|<4,\quad |k_j|<6
 \quad\hbox{at comparison sections}. \tag{7}
\]

At those sections `||E_j||<4`, `||E_j^+||<2`. For a complex perturbation `zeta` of the normalized initial covector, C194 gives

\[
 |\Delta k(t)|\le1296L|\zeta|=:\eta. \tag{8}
\]

Take

\[
 \delta=10^{-40}L^{-1},\qquad |\zeta|<\delta . \tag{9}
\]

All matrix transposes in the analytic extension are bilinear transposes, not conjugate transposes; norms remain Hermitian Euclidean norms. For `eta<=10^{-6}`, direct quotient subtraction gives

\[
 |k\cdot k|\ge8,\quad \|P(k)-P(k_{\rm real})\|\le100\eta,
 \quad\|P(k)\|\le2,\quad \|DP(k)\|\le110. \tag{10}
\]

For example, `|Delta(k dot k)|<=29 eta`; the projector difference is bounded by `29 eta/8+196*29 eta/72<100 eta`. At a section, `D=p dot p` has modulus at least 8, and

\[
 E=[n-mk/(k\cdot k),k\times n],\quad
 E^+=\begin{bmatrix}((k\cdot k)n-mk)^T/D\\(k\times n)^T/D\end{bmatrix}.
 \tag{11}
\]

Quotient differentiation or subtraction gives, with slack,

\[
 \|\Delta E\|,\|\Delta E^+\|\le100\eta,
 \quad\|D_kE\|,\|D_kE^+\|\le100,
 \quad\|E\|\le5,\quad\|E^+\|\le3. \tag{12}
\]

C197's ambient Kelvin reflection identity gives real propagator norm `exp(6T)<3^19<2*10^9`. Its complex generator differs from the real one by at most `1200 eta<1`, so both forward and inverse complex propagators are bounded by `exp(7T)<3^22<10^11`. Duhamel gives their difference at most `10^15 eta`. Therefore, for the actual coefficient block and its reflected inverse,

\[
 \|M_{\rm real}\|<10^{11},\quad \|M_{\rm complex}\|<2\cdot10^{12},
 \quad \|\Delta M\|\le10^{20}L|\zeta|\le10^{-20}. \tag{13}
\]

The same estimate for an inverse block uses the actual backwards ambient propagator and swapped endpoint frames; no ill-conditioned matrix inversion estimate is needed.

For the spatial derivative, C194/C197 give

\[
 \|D_yk\|\le K|\xi|J_1^2J_2\le2\cdot10^{11}L^4,
 \qquad\|D_yA(\Phi_t(y))\|\le1620L. \tag{14}
\]

Here the normalized initial covector has `|xi|<2`, also on the complex neighborhoods. Using (10), the Kelvin generator derivative is at most `3*10^14 L^4`. Thus the complex propagator's spatial derivative is at most `2*10^26 L^4`. Differentiating the two endpoint frames in its coefficient representation gives, for both kinds of block,

\[
 \boxed{\|D_yM\|\le10^{28}L^4}. \tag{15}
\]

This derivative is with respect to the **original** initial position. Hence (15) already charges transport through all earlier blocks.

## 3. A holomorphic graph transform with an explicit disk

Put `sigma=10^{-14}` and let `U_sigma={s in C:dist(s,I)<sigma}`. Fix a real parameter point in the tube and compare each complex block to the corresponding real block there. For `s in U_sigma`, choose `s_real in I` with `|s-s_real|<sigma`. Equations (2) and (13) give

\[
 |M_{11}+M_{12}s|>3000-10^{11}\sigma-2\cdot10^{-20}>2999. \tag{16}
\]

For a fixed real block the exact difference formula is

\[
 f_M(s)-f_M(s_{\rm real})
 =\frac{\det M\,(s-s_{\rm real})}
 {(M_{11}+M_{12}s)(M_{11}+M_{12}s_{\rm real})}. \tag{17}
\]

Thus this difference is smaller than `sigma/10^6`. Changing the matrix by (13) adds less than `10^-23`: use the quotient identity, `|s|<1/4`, and `|f_{M_real}(s)|<1/4`. Consequently every complex graph map sends `U_sigma` into itself. Also `|det M_complex|<3`, because the determinant perturbation is at most `2*10^11*10^-20+10^-40<1`, and

\[
 |\partial_s f_{M_{\rm complex}}|<3/2999^2<10^{-6}. \tag{18}
\]

Both the forward and the backward finite iterates are holomorphic on (9), remain in `U_sigma`, and have modulus below 1/4, independently of `R`.

At fixed slope, the parameter derivative of the quotient is bounded by `||D_yM||`: its numerator is at most `(1+|s|)(1+|f_M(s)|)||D_yM||`, and its denominator exceeds 2999. Differentiating the finite graph recurrence and summing the geometric series in (18) gives

\[
 \|D_ys_j\|,\|D_yr_j\|\le2\cdot10^{28}L^4. \tag{19}
\]

This is the precise reason smooth finite-horizon selectors do not acquire an exponential horizon cost.

## 4. Exact normalized formulas and their jets

Define

\[
 a_j^+=E_j(1,s_j)^T,\quad a_j^-=E_j(1,-r_j)^T,
 \quad q_j^\pm=|a_j^\pm| \quad\hbox{on real parameters},
\]

and let `g_i=g_{M_i}(s_{i-1})` and `h_i=g_{J M_i^{-1}J}(r_i)`. The required Kelvin solutions are exactly

\[
 b_+(j)=\frac{a_j^+}{q_R^+\prod_{i=j+1}^R g_i},\qquad
 b_-(j)=\frac{a_j^-}{q_0^-\prod_{i=1}^j h_i}. \tag{20}
\]

Empty products equal one. The coefficient identities defining `s_j,r_j` verify propagation by `M_j` directly. Orthogonality of the real frame and (7) give

\[
 \frac12<q_j^\pm<2. \tag{21}
\]

Equations (2), (20), and (21) prove (3), with no constancy of an individual off-ray Floquet multiplier assumed.

On each complex covector neighborhood extend `q_j^pm` by the holomorphic square root of `a_j^pm dot a_j^pm` that is positive on real parameters. This square root exists: the displacement between corresponding complex and real graph iterates obeys `Delta_j<=10^-6 Delta_previous+10^-23`, starting from zero, in either direction. Thus every displacement is less than `2*10^-23<sigma/10^5`; (12) then gives `|a_complex-a_real|<10^-16`. Thus

\[
 |q_j^\pm|>\frac13,\quad |a_j^\pm|<3,
 \quad\|D_ya_j^\pm\|<10^{30}L^4,
 \quad\|D_yq_j^\pm\|<10^{31}L^4. \tag{22}
\]

The last two estimates follow from (12), (14), and (19), followed by differentiating `(q_j)^2=a_j dot a_j`; the complex square root is used only to bound jets of the real Euclidean normalization.

Every complex gain in (20) has modulus at least 2999. Its spatial derivative is at most `10^41 L^4`, by (13), (15), and (19). Therefore its logarithmic derivative has modulus at most `10^38 L^4`. Logarithm branches need not be chosen: this assertion is the quotient `D_yg/g`. Product differentiation of (20), using at most `R` factors, gives

\[
 |b_\pm(j)|\le30d_\pm(j),\qquad
 \|D_yb_\pm(j)\|\le10^{45}(R+1)L^4d_\pm(j) \tag{23}
\]

throughout each complex covector ball.

For a holomorphic function bounded by `B` on a complex Euclidean ball of radius `delta`, its real r-th multilinear derivative is bounded by `B r!(4/delta)^r`. To see the constant, inscribe a Cartesian polydisk of radius `delta/2`, apply the multivariate Cauchy formula, and sum directional coordinate factors using `||v||_1<=sqrt(3)||v||_2<2||v||_2`. Apply this once to (23) and once to its spatial derivative. This proves (5)--(6).

At the two physically unit sections there is no product of gains in the normalization. The estimates in (22) therefore also give the useful sharper statements

\[
 \|D_\xi^r b_+(R)\|,\|D_\xi^r b_-(0)\|\le10r!D_*^r,
 \quad
 \|D_yD_\xi^r b_+(R)\|,\|D_yD_\xi^r b_-(0)\|
 \le10^{33}L^4r!D_*^r. \tag{23a}
\]

## 5. Second spatial derivatives

This section makes C194's checked auxiliary third-flow-jet ledger load-bearing for (6a). Differentiating its exact intrinsic-frame identity `D Psi_t=Q_t S Q_0^{-1}` twice, the same quotient calculations give

\[
 \|D^2Q\|\le341,\quad\|D^2Q^{-1}\|\le7926,
 \quad\|D^2a_0\|\le9396.
\]

For transparency the last bound uses the following numerator/denominator ledger for `a_0=n_0/h^2`:

\[
 |n_0|\le65/2,\quad |Dn_0|\le123,\quad |D^2n_0|\le479,
 \quad |D(h^{-2})|\le8,\quad |D^2(h^{-2})|\le222,
\]

so `D²a_0 <=(4/9)*479+2*8*123+(65/2)*222<9396`. The two frame derivatives use `D²(g/h)<=327`, `D²(3Jg)<=14`, and `2*15*14²+6*341=7926`. Here is the full six-channel bound before expansion. Put

\[
 w=1+4t,\quad j_1=18w,\quad
 j_2=756w^2+16524(t+2t^2)+672w,
\]
\[
 \alpha_1=918(t+2t^2),\qquad
 \alpha_2(t)=\int_0^t[9396j_1(s)^2+51j_2(s)]\,ds.
\]

The product-rule channels, using `||Q||<8`, `||Q^-1||<3`, and C194's first frame derivatives, give

\[
 J_3(t)\le9\{3[341j_1^2+14j_2]w+24\alpha_2+63408w
       +84j_1\alpha_1+2352j_1w+1344\alpha_1\}.
\]

The factor nine converts the phase third derivative to physical coordinates. Expanding this bound gives

\[
 J_3(t)\le4474548+749731896t+3071113056t^2+4031918208t^3
 \le4031918208(1+t)^3. \tag{23b}
\]

The checker reconstructs the unexpanded six-channel polynomial as well as checking the final coefficients. No second WKB corrector is inferred from this flow bound.

Differentiating the inverse matrix `F^{-T}` twice now gives

\[
 \|D_y^2k\|\le12[2J_1^3J_2^2+J_1^2J_3]\le10^{21}L^7,
 \qquad\|D_y^2A(\Phi_t(y))\|\le3\cdot10^6L^2. \tag{23c}
\]

On the same complex covector disks, quotient differentiation gives `||D_k²P||<=1000` and `||D_k²E||,||D_k²E^+||<=10^5`. For instance `P=(kk^T)/(k dot k)` gives the sum
`2/8+2*30²/64+225*2/64+2*225*30²/512<1000`.
The product rule then supplies

\[
 \|D_y^2\mathcal L\|\le10^{28}L^8,\quad
 \|D_y^2H_{\rm block}\|\le2\cdot10^{41}L^8,\quad
 \|D_y^2M\|\le10^{45}L^8, \tag{23d}
\]

where `mathcal L=(-I+2P)A`. The middle bound is
`exp(7T)[T sup||D_y²mathcal L||+T²(sup||D_ymathcal L||)²]`.
The frame second derivatives are at most `10^28L^8`; the six channels of `D_y²(E_end^+H E_start)` give the last estimate.

The graph second derivative must be charged separately from its contraction. Write `D=M11+M12 s`, `N=M21+M22 s`, and `f=N/D`. At fixed `s`,

\[
 f_{yy}=\frac{N_{yy}-fD_{yy}-2f_yD_y}{D},\qquad
 f_{ys}=\frac{N_{ys}-f_sD_y-fD_{ys}-f_yD_s}{D},\qquad
 f_{ss}=-\frac{2M_{12}\det M}{D^3}. \tag{23e}
\]

For two different spatial directions the two symmetric terms replace the displayed factor two. Using (13), (15)--(19), and (23d), these formulas imply

\[
 \|f_{yy}\|\le10^{57}L^8,\quad
 \|f_{ys}\|\le10^{42}L^4,\quad |f_{ss}|\le10^3.
\]

Differentiating the graph recurrence a second time gives the inhomogeneous bound

\[
 10^{57}L^8+2\cdot10^{42}L^4(2\cdot10^{28}L^4)
 +10^3(2\cdot10^{28}L^4)^2<10^{71}L^8.
\]

The geometric series with ratio below `10^-6` therefore proves
`||D_y²s_j||,||D_y²r_j||<=10^73L^8`.
Consequently

\[
 \|D_y^2a_j^\pm\|\le10^{74}L^8,\quad
 \|D_y^2q_j^\pm\|\le10^{76}L^8,\quad
 \|D_y^2g_i\|,\|D_y^2h_i\|\le10^{86}L^8.
\]

The second logarithmic derivative of a gain is at most `10^84L^8`, by `g_yy/g-(g_y/g)^2`. For either full denominator in (20), the first and second logarithmic derivatives are therefore at most
`10^39(R+1)L^4` and `10^85(R+1)L^8`. Two product derivatives of (20) now give the complex-disk bound `10^90(R+1)^2L^8 d_pm(j)`. The same Cauchy estimate used for (5) proves (6a).

## 6. Clock and usage boundary

The harmless analytic replacement 3000 by 2999 retains the C193 clock: `2999>exp(8)` (the checker verifies a rational upper bound for `exp(8)`). Hence for

\[
 R=\left\lceil\frac38\log n\right\rceil+1,\qquad q=n^8,
\]

\[
 2999^{-R}<e^{-8}n^{-3}=e^{-8}q^{-3/8}. \tag{24}
\]

The derivative factors in (5)--(6) then cost explicit powers of `1+log q`; they do not remove the seed's `q^{-3/8}` factor. A logarithmically shrinking certified tube is therefore compatible with smooth polarization selection. One must still charge any cutoff and Fourier-parameter chart derivatives when using these symbols in a synthesized packet.

C197's existing theorem starts from a particular projected constant vector. Its numerical constants cannot simply be reused for (20): the new initial symbol derivatives in (5)--(6) must be inserted into the Kelvin and residual estimates. This note supplies those derivatives and the finite-horizon endpoint normalization, but does not itself perform that finite-frequency substitution, spatial/frequency synthesis, band retention, or physical energy balance.
