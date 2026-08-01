# Affine-core viscous capture and the chirped-cell residual

## Status

This note isolates two statements that can be checked without any localization or
asymptotic argument.

1.  The spatially affine core model has a genuinely captured viscous orbit.  For
    a fixed viscosity-to-growth ratio one obtains a quantitative nonzero child,
    only quadratic parent depletion, and eventual decay of both carrier
    polarizations.
2.  The nonlinear relative-phase cell writes its low ridge exactly, but its
    common high sum is not exactly pressure unless the phase slope is constant.
    The surviving high residual is nevertheless quadratically small in
    `phase scale / carrier frequency`.

Neither result localizes the affine core, transports it through a curved column,
or closes an infinite Navier--Stokes cascade.  Those remain separate PDE gaps.

## 1. Closed affine-core ODE

Assume the exact core reduction

\[
 \dot P=4dxy,\qquad
 \dot x=-dPy-\kappa x,\qquad
 \dot y=(-dP+cZ)x-\kappa y,\qquad
 \dot Z=-2cxy.                                      \tag{1}
\]

Here `P` is the parent strain amplitude, `x,y` are the two full carrier
polarizations, and `Z` is the affine child.  The scalar truncation `y=-x` is
not invariant after `Z` becomes nonzero; retaining both `x,y` is essential.

There is an exact affine relation

\[
 Z+\frac{c}{2d}(P-P_0)=0.                            \tag{2}
\]

Put

\[
 r=\frac{c}{2d},\quad \rho=r^{-2},\quad
 \tau=dP_0t,\quad \alpha=\frac{\kappa}{dP_0},
\]
\[
 p=r^2\left(1-\frac P{P_0}\right),\qquad
 U=\frac r{P_0}(x+y),\qquad
 V=\frac r{P_0}(x-y).                               \tag{3}
\]

Then (1) is exactly

\[
\begin{aligned}
 U'&=((1+\rho)p-1-\alpha)U+pV,\\
 V'&=-pU+(1-(1+\rho)p-\alpha)V,\\
 p'&=V^2-U^2 .                                      \tag{4}
\end{aligned}
\]

The equilibrium at the origin has eigenvalues `-1-alpha`, `1-alpha`, and
`0`.  For `0<alpha<1`, choose the positive branch of its one-dimensional
unstable manifold.  As `tau -> -infinity`,

\[
 V\sim e^{(1-\alpha)\tau},\qquad
 p=\frac{V^2}{2(1-\alpha)}+O(V^4),\qquad
 U=\frac{V^3}{4(1-\alpha)(2-\alpha)}+O(V^5).        \tag{5}
\]

Thus this orbit is the vanishing-seed limit; changing the tiny seed only
translates time.

## 2. Exact confinement and convergence

Let

\[
 J_\rho=U^2+V^2+(1+\rho)p^2-2p .                   \tag{6}
\]

Direct differentiation in (4) gives

\[
 J_\rho'=-2\alpha(U^2+V^2).                        \tag{7}
\]

Set `delta=1-alpha` and

\[
 I_U(\tau)=\int_{-\infty}^{\tau}U(s)^2\,ds .
\]

Since `p'=V^2-U^2`, integration of (7) from the zero equilibrium gives the
stronger pointwise identity

\[
 U^2+V^2+(1+\rho)p^2-2\delta p+4\alpha I_U=0.      \tag{8}
\]

The nontrivial branch can never return to `p=0`; (8) would force the whole
past orbit to vanish.  Consequently

\[
 0<p(\tau)\le \frac{2\delta}{1+\rho}               \tag{9}
\]

for every finite `tau`.  This proves global boundedness.  Equation (7) gives
`U,V in L^2`; (4) gives uniform continuity, hence `U,V -> 0`.  Moreover
`p' in L^1`, so `p -> p_infinity`.  At the endpoint,

\[
 p_\infty\bigl(2\delta-(1+\rho)p_\infty\bigr)
       =4\alpha I_U(\infty).                        \tag{10}
\]

Nonzero dissipation already implies `0<p_infinity<2 delta/(1+rho)`.  The next
estimate makes the lower bound explicit.

## 3. Quantitative capture

Rewrite the first equation of (4) as

\[
 U'+a(\tau)U=pV,\qquad
 a=1+\alpha-(1+\rho)p .                             \tag{11}
\]

By (9),

\[
 a(\tau)\ge a_0:=2-3\delta .                        \tag{12}
\]

Assume `delta<2/3`.  Multiply (11) by `U`, integrate over the complete orbit,
and use `U(-infinity)=U(+infinity)=0`:

\[
 a_0 I_U
 \le \int p|UV|
 \le \frac{2\delta}{1+\rho}\sqrt{I_UI_V}.
\]

Therefore

\[
 I_U\le q I_V,\qquad
 q=\left[\frac{2\delta}{(1+\rho)(2-3\delta)}\right]^2.       \tag{13}
\]

Also `I_V-I_U=p_infinity`.  If `q<1`, then

\[
 I_U\le\frac q{1-q}p_\infty .                      \tag{14}
\]

Substitution in (10) proves

\[
 \boxed{
 p_\infty\ge
 \frac{2\delta-4\alpha q/(1-q)}{1+\rho}}
                                                               \tag{15}
\]

whenever the numerator is positive.  Thus

\[
 \frac{Z_\infty}{P_0}=\frac{p_\infty}{r},\qquad
 1-\frac{P_\infty}{P_0}=\frac{p_\infty}{r^2}.       \tag{16}
\]

This is the desired weak-depletion geometry: a child of order `r^-1` costs
only order `r^-2` of the parent.

A convenient fully quantitative choice is

\[
 \alpha=0.9,\qquad r\ge10.
\]

Then (15) gives

\[
 0.1489<p_\infty<0.1981.                            \tag{17}
\]

For comparison, direct integration at `rho=0` gives
`p_infinity approximately 0.18643`; the numerical value is not used in the
proof.

## 4. Terminal carrier damping

At a terminal point `(0,0,p)`, the carrier matrix in (4) has eigenvalues

\[
 \lambda_\pm=-\alpha\pm\sqrt{D_\rho(p)},\qquad
 D_\rho(p)=1-2(1+\rho)p+\rho(2+\rho)p^2.            \tag{18}
\]

For `alpha=0.9`, `r>=10`, and the interval (17), `D_rho` is decreasing in
`p` and

\[
 \sqrt{D_\rho(p_\infty)}<0.839<0.9.                 \tag{19}
\]

Both polarizations therefore decay with a certified dimensionless rate at
least `0.061`.  Initial growth has rate `1-alpha=0.1`.  A seed `e^-G` needs
dimensionless plateau length about `10G`; after capture, carrier cleanup occurs
on the same `O(1/(delta dP0))` scale.

For `delta -> 0`, the slow rescaling

\[
 p=\delta P,\qquad V=\delta W,\qquad U=\delta^2R,
 \qquad T=\delta\tau
\]

has the leading system

\[
 R=PW/2,\qquad P_T=W^2,\qquad W_T=(1-P)W.           \tag{20}
\]

Its vanishing-seed separatrix satisfies `W^2=P(2-P)` and ends at `P=2`.
Thus `p_infinity=2 delta+O(delta^2)`.  The first correction can also be
computed.  At `rho=0`, the fast variable is `R=PW/2+O(delta)`, while

\[
 W_T=(1-P)W-\delta PR+O(\delta^2),\qquad
 P_T=W^2+O(\delta^2).
\]

For `H_slow=W^2+P^2-2P`, this gives

\[
 (H_{\rm slow})_T=-\delta P^2W^2+O(\delta^2).
\]

Along the leading separatrix,
`integral P^2 W^2 dT=integral_0^2 P^2 dP=8/3`.  Comparing the terminal
value of `H_slow` yields the formal expansion

\[
 p_\infty=2\delta-\frac43\delta^2+O(\delta^3).
\]

This explains the numerical endpoint; (15), rather than the formal
expansion, is the rigorous estimate.

## 4A. Robust capture, affine defect, and the seed qualification

The exact affine relation (2) must not be silently assumed after localization.
Introduce the independently normalized child variable

\[
 q=\frac{rZ}{P_0}.
\]

Allowing additive normalized errors, the four-dimensional system is

\[
\begin{aligned}
 U'&=(q+\rho p-1-\alpha)U+qV+e_U,\\
 V'&=-qU+(1-\alpha-\rho p-q)V+e_V,\\
 p'&=V^2-U^2+e_p,\\
 q'&=V^2-U^2+e_q.
\end{aligned}                                                   \tag{R1}
\]

The affine defect `a=q-p` obeys the exact identity

\[
 a'=e_q-e_p.                                                   \tag{R2}
\]

If the physical `P,Z` equations have errors `E_P,E_Z`, then

\[
 e_p=-\frac{r^2E_P}{dP_0^2},\qquad
 e_q=\frac{rE_Z}{dP_0^2},\qquad
 a'=\frac{r(E_Z+rE_P)}{dP_0^2}.                               \tag{R3}
\]

Thus a parent-equation error is amplified by `r^2` in the invariant-breaking
norm.  A residual estimate that controls only the unweighted physical error
can miss the dangerous component.

There is a uniform robust-capture statement, but it starts at a fixed incoming
section rather than at an exponentially small seed.  For `alpha=0.9` and
`0<=rho<=0.01`, choose a compact incoming section of the exact unstable orbit
before its nonlinear turn.  Ordinary `L1` continuous dependence on the fixed
transition segment gives constants `eta_0,T_*>0`, independent of `rho`, such
that data within `eta_0` of that section and

\[
 \int_{\tau_{\rm in}}^{\tau_{\rm in}+T_*}|e|\,d\tau\le\eta_0,
 \qquad
 |a(\tau_{\rm in})|+\int|e_q-e_p|\,d\tau\le\eta_0             \tag{R4}
\]

reach a terminal section with small carrier radius
`R=(U^2+V^2)^(1/2)` and

\[
 0.13\le p,q\le0.22.                                          \tag{R5}
\]

This step is a finite-time ODE theorem: take a fixed tubular neighborhood of
the exact compact transition, a Lipschitz bound `L_*` for its vector field,
and a positive margin `m` to the boundary of (21e).  Gronwall gives the claim
whenever the incoming displacement plus `||e||_1` is less than
`m exp(-L_*T_*)`.  No asymptotic PDE input is hidden here.

On the terminal box, the carrier has a direct Lyapunov estimate:

\[
 \frac12(R^2)'
 =(q+\rho p-1.9)U^2+(0.1-\rho p-q)V^2+Ue_U+Ve_V
 \le-0.03R^2+R|e_c|,                                         \tag{R6}
\]

where `e_c=(e_U,e_V)`.  Consequently

\[
 R(\tau)\le e^{-0.03(\tau-\tau_c)}R_c
 +\int_{\tau_c}^{\tau}e^{-0.03(\tau-s)}|e_c(s)|\,ds.          \tag{R7}
\]

Moreover,

\[
 |\Delta p|,|\Delta q|
 \le\int_{\tau_c}^{\tau}R^2\,ds+\|e_p\|_1+\|e_q\|_1.        \tag{R8}
\]

Small fixed tail norms therefore retain the terminal box, and there

\[
 \boxed{\frac{0.12}{r}\le\frac ZP
 =\frac{q}{r(1-\rho p)}\le\frac{0.23}{r}}.                    \tag{R9}
\]

For a cleanup interval of length `G/0.03`, (21g) yields

\[
 R_{\rm out}\le R_ce^{-G}
 +\int e^{-0.03(\tau_{\rm out}-s)}|e_c(s)|\,ds.               \tag{R10}
\]

The carrier therefore exits exponentially small only if the residual is
switched off before cleanup or its backward-weighted tail is exponentially
small.  An unweighted `o(1)` residual gives only an `o(1)` exit.

There is an equally important qualification during the gain phase.  No
theorem allowing arbitrary additive `o(1)` error can preserve an `e^{-G}`
seed: an additive impulse of mass `e^{-G}` in the unstable polarization can
cancel it exactly.  Any PDE transition must instead provide at least one of:

1. a theorem beginning at the fixed incoming section;
2. multiplicative carrier error `e_c=B(\tau)(U,V)` with
   `integral ||B|| << 1`; or
3. a backward-weighted unstable forcing smaller than the seed,
   \[
   \left|\int_{\tau_-}^{\tau_{\rm in}}
   e^{-0.1(s-\tau_-)}\ell_u\cdot e_c(s)\,ds\right|
   \le\theta e^{-G},\qquad\theta<1.                           \tag{R11}
   \]

Under one of these seed-survival hypotheses, gain takes `10G+O(1)`
dimensionless time and cleanup takes at most `G/0.03`.  The full capture still
fits in a `C G` action window.  This is the robustness target that a localized
PDE block must meet; a bare small-residual statement is insufficient.

## 5. Exact low ridge from a nonlinear relative phase

Fix an orthonormal frame `(r,t,h)`, write `s=r dot x`, and set

\[
 w=r+Hh,\qquad
 \phi_\sigma=\Lambda t\cdot x+\sigma\psi(s)/2,
 \qquad
 a_\sigma=w-\sigma\frac{\psi'(s)}{2\Lambda}t,
 \quad \sigma=\pm1.                                \tag{21}
\]

Each complex wave `u_sigma=A_sigma a_sigma exp(i phi_sigma)` is exactly
divergence-free.  Indeed `div(a_sigma)=0` and
`a_sigma dot grad(phi_sigma)=0`.

The ordered low interaction of `u_+` with the conjugate of `u_-`, plus its
transpose, is

\[
 2iA_+\overline{A_-}\psi'(s)(r+Hh)e^{i\psi(s)}.
\]

Its radial part is a gradient.  With the Euler sign, the projected complex
force is exactly

\[
 -2iH A_+\overline{A_-}\psi'e^{i\psi}h.             \tag{22}
\]

Taking `A_+ conjugate(A_-)=ic`, adding the real conjugate, and choosing

\[
 \sin\psi(s)=\frac1{4Hc}\int_0^s f(q)\,dq          \tag{23}
\]

writes the real ridge `f(s)h` exactly.  In particular `f(s)=gamma s` is
obtained from `sin(psi)=gamma s^2/(8Hc)` on any core where the right side
stays uniformly away from `+/-1`.

## 6. Correction: the common high sum is not exactly pressure

Let `delta_phase=psi'/(2 Lambda)`.  A direct ordered-pair calculation gives

\[
 (u_+\cdot\nabla)u_-+(u_-\cdot\nabla)u_+
 =-iA_+A_-\frac{(\psi')^2}{\Lambda}
       t\,e^{i2\Lambda t\cdot x}.                  \tag{24}
\]

The polarization-derivative terms do cancel, but the coefficient in (24)
depends on the transverse variable `s`.  It is a gradient only when
`psi'` is constant.

For an exact formula put `k=2 Lambda`,

\[
 F(s)=-iA_+A_-\frac{(\psi')^2}{\Lambda},\qquad
 L_k=(\partial_s^2-k^2)^{-1}.
\]

On the whole line (or mode by mode in a periodic core),

\[
 \mathbb P\{F(s)t e^{iky}\}
 =\{-ik(L_kF)'r+(L_kF)''t\}e^{iky}.                 \tag{25}
\]

For transverse frequencies much smaller than `k`,

\[
 \mathbb P N_{+,-}^{\rm high}
 =A_+A_-\frac{\psi'\psi''}{\Lambda^2}
       r\,e^{i2\Lambda y}
 +O\!\left(\frac{|A_+A_-|Q^4}{\Lambda^3}
             +\frac{|A_+A_-|Q^5}{\Lambda^4}\right),           \tag{26}
\]

where `|partial^m psi| <= C_m Q^m`.  The first displayed remainder is the
longitudinal component; the second is the next radial term.

The self interaction is

\[
 (u_\sigma\cdot\nabla)u_\sigma
 =-\sigma A_\sigma^2\frac{\psi''}{2\Lambda}
 t\,e^{i(2\Lambda y+\sigma\psi)}.                  \tag{27}
\]

Applying (25) with
`F_sigma=-sigma A_sigma^2 psi'' exp(i sigma psi)/(2 Lambda)` gives the leading
radial residual

\[
 \mathbb P N_{\sigma,\sigma}
 =\frac{A_\sigma^2}{4\Lambda^2}
   \bigl(\psi'\psi''-i\sigma\psi'''\bigr)
   r\,e^{i(2\Lambda y+\sigma\psi)}
 +O\!\left(\frac{|A_\sigma|^2Q^4}{\Lambda^3}
             +\frac{|A_\sigma|^2Q^5}{\Lambda^4}\right).      \tag{28}
\]

Thus a formula retaining only `psi' psi''` is incomplete unless
`psi'''=0` on the exact core.

There is also a lower-frequency term which is absent from a complex-wave-only
calculation.  Pairing each wave with its real conjugate gives

\[
 N_{\rm zero}
 =\bigl(|A_-|^2-|A_+|^2\bigr)\frac{\psi''}{\Lambda}
 t.                                                        \tag{29}
\]

It cancels for equal carrier magnitudes.  The reciprocal child feedback,
however, does not preserve `A_-=-A_+`: in (4), it creates the stable component
`U`.  If `A_+=x,A_-=y`, then

\[
 y^2-x^2=-\frac{P_0^2}{r^2}UV,\qquad
 xy=-\frac{P_0^2}{4r^2}p'.                              \tag{30}
\]

Using (13)--(14),

\[
 \int |UV|\le\frac{\sqrt q}{1-q}p_\infty .             \tag{31}
\]

Consequently the time-integrated zero-mode defect, relative to the desired
cross beat, is bounded by

\[
 \frac{4\sqrt q}{1-q}\frac Q\Lambda .                  \tag{32}
\]

For `alpha=0.9,r>=10` the constant before `Q/Lambda` is below `0.48`.  This
is one order larger than the high-band chirp residual.  If `psi` is exactly
quadratic, `psi''` is constant and (29) is merely a spatially uniform
Galilean acceleration; its strain vanishes.  The exact arcsine profile in
(23) has nonconstant `psi''`, so its induced strain must be budgeted at
order `Q/Lambda`.

A useful compromise is therefore

\[
 \psi(s)=a s^2,\qquad \chi=|a|R^2\ll1.              \tag{33}
\]

On `|s|<=R`, its written ridge is

\[
 8Hca\,s\cos(as^2)=\gamma s\,[1+O(\chi^2)],
 \qquad \gamma=8Hca.                                \tag{34}
\]

It loses exact affinity but removes the larger zero-mode strain defect.
Taking `G chi^2 -> 0` makes the nonlinear child-gradient error perturbative
over gain action `G`.  At fixed `gamma,R`, the required carrier product is
`c=gamma R^2/(8H chi)`, so choosing, for example, `chi=G^-1` costs only a
polylogarithmic factor while giving `G chi^2=G^-1`.  A rigorous PDE use still
needs a robustness estimate for (4), but the terminal gap in (19) provides a
fixed margin for such an estimate.

If each carrier has size `A`, the desired low force has size `A^2 Q`, while
(26) and (28) have size

\[
 A^2\frac{Q^3}{\Lambda^2}
   =A^2Q\left(\frac Q\Lambda\right)^2.              \tag{35}
\]

In the power-law ledger `Q=N^b`, `Lambda=N^(beta/2)`, this relative error is

\[
 (Q/\Lambda)^2=N^{-(\beta-2b)}.                    \tag{36}
\]

It remains negligible over a gain interval `G` provided
`G N^{-(beta-2b)} -> 0` and the `2 Lambda` bands do not have their own
resonant instability.  If `nu Lambda^2` is comparable to the parent growth
rate, those bands instead see about four times the viscous damping.

## 7. What this does and does not resolve

The exact affine core has no phase ladder: `U_child=Z(t)s h` has constant
gradient, `U_child dot grad` annihilates carriers independent of `h`, and
`u_carrier dot grad(U_child)` stays at the same carrier phase.  Its reciprocal
feedback is precisely the closed two-polarization system (1).

The remaining obstruction is localization.  A cutoff makes the child gradient
nonconstant in the buffer, introduces new sidebands, and creates Euler and
viscous seam errors.  A Bogovskii correction restores divergence but does not
make those errors small.  Sequential pulses additionally require redesigning
the next carrier at the newly created strain scale.  None of these claims is
settled by the finite-dimensional capture lemma.

## 8. Audit of the global-periodic-ridge shortcut

Working on the periodic Clay formulation removes endcaps, and a global ridge
`u=A h f(Qs)` is indeed an exact steady Euler shear.  It does not, however,
meet the finite-energy scale ledger.  A next-stage strain of size `Q^beta`
requires

\[
 A\asymp Q^{\beta-1}.
\]

For a nonintermittent ridge on a fixed torus this is also its `L2` size, which
diverges.  A straight exact shear can be intermittent in at most the two
directions perpendicular to `h`; it must remain invariant along `h`.  If its
effective transverse bandwidth is `Q`, the two-dimensional Bernstein estimate
gives

\[
 \|\nabla u\|_\infty\lesssim Q^2\|u\|_2,
 \qquad
 \|u\|_2\gtrsim Q^{\beta-2}.                         \tag{37}
\]

The BAS/viscous window under consideration has `beta>2`, so even an optimally
thin straight tube has diverging energy.  Making the tube thinner merely
changes the effective bandwidth and reproduces (37).  Full three-dimensional
intermittency is therefore essential; it requires a curved/compact Euler
geometry or a cutoff with a genuinely controlled seam.  The global periodic
ridge is a useful local model but not a finite-energy cascade workaround.

The cleaner next *cell* lemma is instead a phase-localized quadratic chirp:

* `psi(s)=a s^2` on the child core;
* `psi` is smoothly constant outside a slightly larger interval;
* no amplitude cutoff is made in the `s` direction.

It combines the `O(chi^2)` affine accuracy in (34)--(35) with a purely
Galilean zero-mode defect on the exact inner core.  All `psi'''` and
nonconstant zero-mode strain are confined to the phase-transition buffer.
This is analytically cleaner than the exact arcsine profile.  By itself it
still produces a slab-like child and leaves the carrier nonintermittent in the
other directions, so it is not yet a finite-energy construction.  Its proper
use would be inside an independently localized BAS/Gavrilov block, followed by
a quantitative buffer estimate.
