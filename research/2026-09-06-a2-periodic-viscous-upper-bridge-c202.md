# C202: a periodic multi-beam upper bridge through the actual viscous pump

**Date:** 2026-09-06

**Status:** explicit periodic linearized-Navier--Stokes comparison for the
compact C197 family, with a nonempty budget on C191's factorial schedule.
No expanding aperture, retained-band lower bound, physical normalization,
nonlinear stage return, UVSR, or singularity is asserted.

**Checker:** [a2_periodic_viscous_upper_bridge_c202.py](../checks/a2_periodic_viscous_upper_bridge_c202.py).

## 1. Periodicity does not require a new localization estimate

Use the torus `T3=(R/(2 pi Z))^3`, with **unnormalized** Lebesgue L2 and
H1 norms on a fundamental cube. The A2 field of C194 is smooth and
periodic, satisfies `Delta U=-2U`, and obeys the global bounds

\[
 \|DU\|_{\rm op}\le6,\qquad
 \|D^2U\|_{\rm mult}\le15/2.                         \tag{1}
\]

Fix a lift of `y0` with `f(y0)=0`, a real smooth profile chi supported in
the unit ball, a unit real polarization e, `0<hbar<=eps<=1/40`, and
integer carriers

\[
 p\in p_0+d\mathbb Z^3,\quad p_0\in\mathbb Z^3,
 \quad d\in\mathbb N,\quad d\varepsilon\ge1,\quad
 \xi_p=\hbar p\in\xi_c+[-1/4,1/4]^3,
 \quad |\xi_c|=1.                                      \tag{2}
\]

Use precisely C197's symbols B, c and phases on the lifted flow:
`B(0,y,xi)=(I-P(xi))e`, `c=-k cross b/|k|^2`,
`b=chi_eps(Y_t(x)) B(t,Y_t(x),xi)`, and

\[
 u_{\rm app}(t,x)=\sum_p a_p\frac{\hbar}{i}
       \operatorname{curl}\big(e^{ip\cdot Y_t(x)}c_p(t,x)\big).
                                                               \tag{3}
\]

The initial compact patch is periodized over `2 pi Z^3`. The lifted flow
obeys `Phi_t(x+2 pi m)=Phi_t(x)+2 pi m`; the same holds for its inverse.
Thus `exp(i p dot Y_t)` is periodic for every integer p. All symbol
coefficients are equivariant under the same translations.

The projection of the transported initial patch remains injective. Indeed,
`Phi_t(x)=Phi_t(y)+2 pi m` implies `x=y+2 pi m`, and two distinct points of
the initial ball of radius eps cannot differ by a nonzero `2 pi m`.
Distinct translated supports therefore never overlap. Periodizing the
compact R3 approximation and its local residual preserves each L2 norm
exactly. Its exact-curl form also gives zero torus mean and divergence.
No globally single-valued real phase is claimed or needed.

Let uE be the exact torus linearized-Euler solution with datum
`u0=u_app(0)`. The residual identity and torus Leray energy estimate give
the same constants as C197:

\[
 \boxed{\|u_E(t)-u_{\rm app}(t)\|_2\le
 \hbar e^{6t}\|a\|_{\ell^2}
 \left[10^{94}(1+t)^{37}\|\chi\|_\infty
  +10^{80}\varepsilon^{-1}(1+t)^{31}
       \|\nabla\chi\|_\infty\right].}                 \tag{4}
\]

This argument periodizes the explicit approximation and residual, then
solves the global torus pressure problem. It does not periodize the exact
R3 Euler solution, which need not have compact support.

## 2. A count-independent initial H1 estimate

Put

\[
 C_0=\|\chi\|_\infty,\quad C_1=\|\nabla\chi\|_\infty,
 \quad C_2=\|D^2\chi\|_{F,\infty},\qquad
 A_2=\max(C_0,C_1,C_2).                                \tag{5}
\]

Here F is the Frobenius norm of the Hessian. The symbol `A_2` in (5) is
only a profile constant, unrelated to the name of the background field.
For initial symbols write

\[
 B_p=(I-P(\xi_p))e,\qquad c_p=-\frac{\xi_p\times e}{|\xi_p|^2}.
                                                               \tag{6}
\]

The parameter cube ensures `1/2<|xi_p|<2`, and hence

\[
 |B_p|\le1,\quad |c_p|\le2,\quad
 |p|\,|B_p|\le2\hbar^{-1},\quad
 |p|\,|c_p|\le\hbar^{-1}.                            \tag{7}
\]

The carrier-cell proof of C197(15) applies to vector and tensor
coefficients, as well as to multiplication by a vector- or tensor-valued
envelope. It gives constant 27 for each of the following product-rule
terms, independently of the carrier count. At time zero,

\[
 u_0=\chi_\varepsilon\sum_pa_pB_pe^{ip\cdot x}
  +\frac{\hbar}{i}\nabla\chi_\varepsilon\times
      \sum_pa_pc_pe^{ip\cdot x},                     \tag{8}
\]

up to an irrelevant constant phase for the patch center. Consequently

\[
 \|u_0\|_2\le27(C_0+2\hbar\varepsilon^{-1}C_1)\|a\|_2,
                                                               \tag{9}
\]

\[
 \|\nabla u_0\|_{2,F}\le27\left[
  2\hbar^{-1}C_0+2\varepsilon^{-1}C_1
    +2\hbar\varepsilon^{-2}C_2\right]\|a\|_2.       \tag{10}
\]

In (10) the four differentiated products cost respectively
`eps^-1 C1`, `2 hbar^-1 C0`, `2 hbar eps^-2 C2`, and `eps^-1 C1`.
The last coefficient uses the sharper last bound in (7). Combining (9)
and (10), using `sqrt(x^2+y^2)<=x+y` and `hbar<=eps<=1`, yields

\[
 \boxed{\|u_0\|_{H^1}\le
 27[3\hbar^{-1}C_0+4\varepsilon^{-1}C_1
       +2\hbar\varepsilon^{-2}C_2]\|a\|_2
 \le243\hbar^{-1}A_2\|a\|_2.}                       \tag{11}
\]

The weaker H1 growth estimate below is chosen so that no second or third
transported symbol derivatives are required.

## 3. The full PDE viscous comparison

For an arbitrary smooth mean-zero divergence-free Euler perturbation w,
let `X=||w||2`, `Y=||grad w||2,F`. Differentiating the linearized equation
and integrating against w and its derivatives gives

\[
 \frac12\frac d{dt}\|w\|_{H^1}^2
 \le6X^2+12Y^2+13XY\le19(X^2+Y^2).                  \tag{12}
\]

There are two gradient-gradient terms, each bounded by `6Y^2`: one from
the stretching term and one from differentiation of transport. The
coefficient derivative is at most
`sqrt(3)*(15/2)XY<13XY`. Pressure terms vanish against the
divergence-free derivatives. The final bound follows from
`13XY<=(13/2)(X^2+Y^2)`. Thus

\[
                \|u_E(t)\|_{H^1}\le e^{19t}\|u_0\|_{H^1}. \tag{13}
\]

Let mu>0 be the actual normalized viscosity. The heat-decaying background

\[
 U_\mu(s)=\alpha(s)U,\qquad \alpha(s)=e^{-2\mu s},
 \qquad \tau(s)=\frac{1-e^{-2\mu s}}{2\mu}           \tag{14}
\]

is an exact unforced Navier--Stokes solution, since the A2 field is
Beltrami with `Delta U=-2U`. Let u_mu solve the genuine full torus
linearized Navier--Stokes equation around (14), with datum u0:

\[
 \partial_su_\mu+\alpha(s)(U\cdot\nabla u_\mu+Au_\mu)
       -\mu\Delta u_\mu+\nabla p_\mu=0,
 \qquad\operatorname{div}u_\mu=0.                   \tag{15}
\]

The inviscid comparison `w(s)=uE(tau(s))` has the same transport and
stretching coefficients. With `z=u_mu-w`, its energy satisfies

\[
 \frac d{ds}\|z\|_2^2
  \le12\alpha(s)\|z\|_2^2+\mu\|\nabla w\|_2^2,
 \qquad z(0)=0.                                      \tag{16}
\]

Indeed the diffusion discrepancy pairs as
`-mu <grad w,grad z>`; Young's inequality consumes half the positive
viscous dissipation and gives (16). Applying (13) and the integrating
factor `exp(-12 tau)` gives the explicit full-PDE comparison

\[
 \boxed{\|u_\mu(s)-u_E(\tau(s))\|_2
 \le\sqrt{\mu s}\,e^{19\tau(s)}\|u_0\|_{H^1}
 \le243\hbar^{-1}A_2\|a\|_2\sqrt{\mu s}\,e^{19\tau(s)}.}
                                                               \tag{17}
\]

The only exponential loss here is `e^{19 tau}`; the viscosity factor is
`sqrt(mu s)`. This is an additive upper error, not a multiplicative lower
damping factor for a certified growing packet. Equations (4) and (17)
give the actual periodic viscous-to-parametrix bound by the triangle
inequality at `t=tau(s)`, with no dependence on the number of carriers.
Taking twice the real part gives a real solution and costs at most two
relative to the positive-carrier coefficient norm in all displayed upper
bounds. No positive/negative output orthogonality is used.

## 4. A concrete budget on the actual factorial stage scale

Use C191's fixed shell schedule and C197's two distinct upper-error tests:

\[
 q=n^8,\quad \mu=\nu((n-1)!)^{-3/2},\quad\hbar=q^{-1},
 \quad t\le t_*:=\frac{57}{400}\log q+\frac{152}{25}. \tag{18}
\]

For the stable-width test take `eps=q^-1/12`, `a=a_ref`; for the
expanding-width test take `eps=q^-1/3`, `a=q^-3/8 a_ref`. These labels
describe the inherited parameter choices, not proved stable or expanding
bundles. Assume `A2 ||a_ref||2>0` and the completely explicit condition

\[
                     n\ge\max(10^{1250},\nu).        \tag{19}
\]

For every integer n>=4, at least `(n-1)/2` factorial factors are at least
`n/2>=sqrt(n)`. Therefore

\[
 (n-1)!\ge n^{(n-1)/4},\qquad
 \mu\le n^{1-3(n-1)/8}\le n^{-800}=q^{-100},          \tag{20}
\]

where the final inequality already holds for n>=2137. In particular
`q>=10^10000`, `t_*<=log q`, and

\[
 2\mu t_*\le2q^{-100}\log q\le q^{-99}<1/100.        \tag{21}
\]

There is a physical pump time s with `tau(s)=t`, and the same clock
estimate as C191 gives

\[
 t\le s=-\frac{\log(1-2\mu t)}{2\mu}
  \le\frac{100}{99}t\le2\log q\le q,\qquad
 \alpha(s)\ge99/100.                                \tag{22}
\]

Now `19*(57/400)=1083/400`, and
`exp(19*152/25)<3^116<10^56`. The real-completed viscous discrepancy
from (17), divided by `A2 ||a_ref||2`, obeys the single explicit chain

\[
 \begin{split}
 \frac{\|u_{\mu,\mathbb R}(s)-u_{E,\mathbb R}(t)\|_2}
      {A_2\|a_{\rm ref}\|_2}
 &\le486q\sqrt{\mu s}\,e^{19t}\\
 &<486\cdot10^{56}
       q^{1-50+1/2+1083/400}\\
 &=486\cdot10^{56}q^{-18317/400}
  <10^{59}q^{-45}\le10^{-449941}.                    \tag{23}
 \end{split}
\]

C197(25) and (4) give the other error at most `4*10^-319` under the same
normalization, since `max(C0,C1)<=A2`. Consequently

\[
 \boxed{\|u_{\mu,\mathbb R}(s)-u_{{\rm app},\mathbb R}(t)\|_2
 <\frac1{100}A_2\|a_{\rm ref}\|_2.}                  \tag{24}
\]

The normalized physical stage time is s; C191's physical time is
`((n-1)!)^-35/2 s`. This uses the actual heat-decaying A2 pump rather
than replacing its clock by a stationary viscous coefficient.

The carrier conditions have explicit instances: for `hbar=q^-1`, let
p0 be the coordinatewise nearest integer to `q xi_c`, choose
`d=ceil(eps^-1)`, and take any finite subset of

\[
 p=p_0+dr,\qquad r\in\mathbb Z^3,\qquad
 |r_j|\le\left\lfloor\frac q{8d}\right\rfloor.        \tag{25}
\]

Then `d eps>=1` and
`|p/q-xi_c|_infty<=1/(2q)+1/8<=1/4` for q>=4. Thus the positive
family is nonempty and periodic, and real completion is defined exactly.
For the expanding-width choice and q>=64, the number in (25) is at least
`q^2/4096`: `d<=2q^1/3`, `q/(8d)>=1`, and
`2 floor(x)+1>=x` give this lower bound. This count is availability of
carriers in the upper-estimate cube, not availability of growing carriers.

## 5. Smooth initial polarizations and the C200 selectors

The bridge also applies to spatially varying smooth transverse initial
symbols `B0(y,xi)`, including C200's finite-horizon selectors. The needed
parameter domain can be local. Replace C197's coordinate cube by

\[
 \xi(\theta)=\xi_c+\mathsf A\sin\theta,
 \qquad\|\mathsf A\|_{\rm op}\le1/4,\qquad
 1/2<|\xi(\theta)|<2.                                  \tag{26}
\]

Assume this entire image and the initial support lie in the symbol's
domain, and that every chosen carrier has a theta label. No extension
outside that domain is required. The coordinate derivatives of the sine
map are columns of `A` times sine or cosine; mixed derivatives vanish.
Thus the same product-rule count as C197(12) gives

\[
 \|(1-\Delta_\theta)^2S(\xi(\theta))\|
 \le216\sum_{r=0}^4\frac{\sup\|D_\xi^rS\|}{r!}.       \tag{27}
\]

The synthesis constant is again `262440`. In particular an anisotropic
small parameter parallelepiped does not need an unproved global
polarization extension.

For j=0,1,2 define the finite initial seminorms

\[
 M_j=\sum_{r=0}^4\frac1{r!}
  \sup\|D_y^jD_\xi^rB_0(y,\xi)\|_{\rm mult}.         \tag{28}
\]

The suprema are on the initial support and full parameter image. If a
neighborhood is required to define derivatives, the same stated bounds
are assumed there. Spatial derivatives have multilinear operator norm.
Transversality `xi dot B0=0` is retained.

Equation (4) holds multiplied by `M0+M1`. To verify this directly with
C197's exact recurrences, the initial r-th xi derivative is at most
`r! M0<=C_r h^r M0`. Its propagated amplitude is therefore bounded by
`M0` times C197's E-majorant. The initial y derivative propagates by the
same homogeneous system, contributing at most
`M1 e^{6t}h^r E_r(1+t)^{3r}` before the final Eulerian factor `J1`.
The inhomogeneous mixed derivative is bounded by `M0` times the old
G-majorant. The finite rational recurrences satisfy `G_r>=E_r` for
0<=r<=4, as the checker verifies. Hence the whole spatial derivative
is bounded by `M0+M1` times C197's spatial majorant, with unchanged time
powers, and the residual follows by linearity in B and its derivative.

For the H1 estimate put `W=J_4(xi/|xi|^2)=1654` as an upper bound on
the parameter image: explicitly

\[
 W=\sum_{r=0}^4\frac{V_r2^{r+1}}{r!}
   =2+12+56+272+1312=1654.                            \tag{29}
\]

Then `J_4(D_y^j c0)<=W Mj`. Decomposing the exact initial curl as

\[
 u_0=\sum_pa_pe^{ip\cdot y}\left[
  \chi_\varepsilon B_0+\frac{\hbar}{i}
  (\nabla\chi_\varepsilon\times c_0
        +\chi_\varepsilon\operatorname{curl}_y c_0)\right]
                                                               \tag{30}
\]

and applying (27) to each product gives

\[
 \boxed{\|u_0\|_{H^1}
 \le10^{10}\hbar^{-1}A_2(M_0+M_1+M_2)\|a\|_2.}     \tag{31}
\]

Here is the complete constant ledger. The `M0` terms in L2 plus
gradient are bounded, after factoring out `hbar^-1 A2`, by `5+5W`.
The `M1` terms are bounded by `2+(29/2)W`, and the `M2` terms by `5W`.
These follow by differentiating the three summands in (30): the phase
derivative uses `J_4(xi)<=3`; curl uses `sqrt(6)<5/2`; converting a
single operator spatial derivative to Frobenius costs less than two;
and the derivative of curl costs less than five for its full gradient.
The envelope Hessian already has Frobenius norm in (5). All three
coefficients are at most 23985, and
`262440*23985<10^10`. The bounds `hbar<=eps<=1` absorb each remaining
envelope factor. These are tensor-symbol applications of (27), so the
constant is independent of carrier count.

Combining (31) with (17) replaces the constant `243` there by
`10^10(M0+M1+M2)`. This is an explicit extension to actual initial
symbols, rather than an assumption that the selected polarization has
the fixed-e jets of C197.

For a concrete substitution, [C200](2026-09-06-a2-smooth-finite-horizon-filter-c200.md)
gives its smooth initial finite-horizon polarizations, on its stated
joint base/covector domain, the bounds

\[
 \begin{split}
 \|D_\xi^rB_0\|&\le30r!\mathcal D^r d,\\
 \|D_yD_\xi^rB_0\|&\le10^{45}(N+1)L^4r!\mathcal D^r d,\\
 \|D_y^2D_\xi^rB_0\|&\le10^{90}(N+1)^2L^8r!\mathcal D^r d,
 \end{split}\qquad
 L=1+NT,\quad\mathcal D=4\cdot10^{40}L.              \tag{32}
\]

For the expanding selector `d=2999^-N`; for the stable selector `d=1`.
The fixed factor K in C200's covector `K D Phi_t^-T xi` does not change
the Kelvin amplitude equation, which has degree zero in the covector;
the physically normalized B0 and its normalized-xi derivatives are
therefore the same input symbols here.
Since T>3, `N+1<=L`. Summing (32) for r=0,...,4 yields

\[
 M_0+M_1\le2\cdot10^{208}L^9d,\qquad
 M_0+M_1+M_2\le2\cdot10^{253}L^{14}d.                 \tag{33}
\]

Suppose the chosen carrier image (26) and initial support actually lie
in that C200 domain. On the clock (18), put `x=log q`, so `L<=x`.
Take `a=a_ref` in both cases, `eps=q^-1/12` for the stable selector and
`eps=q^-1/3` for the expanding selector. For the expanding case require
`2999^N>=q^{3/8}`; this places the small initial factor in the selector,
so no second factor is inserted in the coefficient sequence. The real
Euler error then satisfies the explicit common bound

\[
 \frac{\|u_{E,\mathbb R}-u_{{\rm app},\mathbb R}\|_2}
      {A_2\|a_{\rm ref}\|_2}
 \le8\cdot10^{320}x^{46}q^{-37/600}
 \le8\cdot10^{-66}.                                  \tag{34}
\]

The first coefficient includes the two profile terms, real completion,
`exp(6t)<=10^18q^{171/200}`, and (33). For the second, at
`q0=10^10000` one has `x0<10^5`, `x0^46<10^230`, and
`q0^-37/600<10^-616`. The function `x^46 exp(-37x/600)` is decreasing
for `x>=20000`, since `46*600/37<20000`.

Under (19)--(22), the real viscous discrepancy from (31) is bounded by

\[
 \frac{\|u_{\mu,\mathbb R}-u_{E,\mathbb R}\|_2}
      {A_2\|a_{\rm ref}\|_2}
 \le4\cdot10^{319}x^{14}q^{-45}
 \le4\cdot10^{-449611}.                              \tag{35}
\]

The first inequality follows from the stronger intermediate power
`q^-18317/400` and `d<=1`. The second again follows by monotonicity
and `x0^14<10^70`. Thus (34)+(35) is less than `1/100` for the actual
periodic viscous evolution of those selected initial symbols, whenever
the carrier/support inclusion stated above holds. That inclusion must
be checked for the same proposed endpoint geometry; (32) is not a
statement that every carrier in the fixed broad cube of (2) grows.

### 5.1 A fixed smooth profile with a numerical derivative bound

The profile constant need not remain uninstantiated. Put `g(t)=exp(-1/t)`
for t>0 and g(t)=0 for t<=0, and define

\[
 H(t)=\frac{g(t)}{g(t)+g(1-t)},\qquad
 \eta(x)=H(2-4|x|),\qquad
 \chi(x)=\eta(x_1)\eta(x_2)\eta(x_3).                 \tag{36}
\]

The denominator never vanishes. H is smoothly zero for t<=0 and one
for t>=1. The absolute value in eta causes no nonsmoothness at zero
because eta is constant on a neighborhood there. This chi has the
plateau and support required by C201 and takes values in [0,1].

For 0<t<1, put u=1/t. The elementary maxima of `u^m exp(-u)` give
`|g'|<=1` and `|g''|<=8`; the latter uses
`4^4/e^4+2*3^3/e^3<8`, with `e>8/3`. One of t,1-t is at least 1/2,
so the denominator D obeys `D>=e^-2>1/9`, while `|D'|<=2` and
`|D''|<=16`. Two quotient derivatives therefore give

\[
 |H'|\le9+162=171,\qquad
 |H''|\le72+324+1296+5832=7524.                       \tag{37}
\]

It follows that `|eta'|<=684` and `|eta''|<=120384`. The product profile
satisfies

\[
 \|\nabla\chi\|_\infty<1368,\qquad
 \|D^2\chi\|_{F,\infty}
 \le3\cdot120384+6\cdot684^2<4\cdot10^6,
 \qquad \boxed{A_2\le4\cdot10^6.}                    \tag{38}
\]

Thus the absolute real error from (34)--(35), when
`||a_ref||2<=1`, is less than `10^-58` for this one fixed profile. This
last bound is an explicit physical L2 error; it still does not assert
that a growing or concentrated exact solution has been proved without
the endpoint hypotheses discussed above.

## 6. Boundary and verification

This result removes the periodicity and viscous **upper-error** obstacles
for the specified compact C197 family. It supplies no lower comparison
between `A2 ||a_ref||2` and physical entrance energy, and no coherent
endpoint, output Fourier retention, or growing common aperture. The
stage threshold is explicit and intentionally enormous; it is not an
instantiated approximate nonlinear profile. No part of the C185
operator-norm citation restriction or the terminal singular-center
obligation is discharged here.

The dependency-free checker verifies the product-rule constants,
H1-energy slack, integrating-factor exponents, factorial cutoff,
carrier arithmetic, and every rational power and integer comparison in
the finite-stage chain. It does not substitute arithmetic assertions for
the PDE integration-by-parts arguments displayed above.
