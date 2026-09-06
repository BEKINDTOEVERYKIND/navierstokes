# C197: uniform compact multi-beam Kelvin/WKB upper estimates

**Date:** 2026-09-06

**Status:** proved explicit-constant upper-error theorem on C194's invariant annulus. It proves no expanding bundle across a
fixed aperture, retained output band, viscosity estimate, or nonlinear return.

**Checker:** [a2_covector_jet_synthesis_c197.py](../checks/a2_covector_jet_synthesis_c197.py).

All spatial and covector derivatives below use Euclidean multilinear operator
norms. The ambient field is C194's smooth periodic A2 field on R3, with
`D=6`, `A1<=15/2`, and invariant annulus `|f|<=1/10`. Its flow obeys

\[
 J_1(t)\le54(1+4t)\le216(1+t),\qquad
 J_2(t)\le286992(1+t)^2.                                      \tag{1}
\]

## 1. The full ambient Kelvin propagator retains exponent six

Let `k` solve the covector equation, `q=k/|k|`, and let `H(t,s)` be the
full three by three propagator for the Kelvin amplitude equation. Its
coefficient matrix has the exact factorization

\[
 L=(-I+2qq^T)A,\qquad (-I+2qq^T)^T(-I+2qq^T)=I.             \tag{2}
\]

Thus `||L||=||A||<=6`, even on nontransverse inputs, and the full ambient
energy inequality gives

\[
 \boxed{\|H(t,s)\|\le e^{6(t-s)}}\qquad(t\ge s).             \tag{3}
\]

No transversality is needed for a derivative forcing term. In particular,
repeated covector differentiation need not replace the time exponent by
`18` or introduce a loss once per orbit period. The orthogonal reflection
factorization is valid for real covectors; complex amplitudes are allowed.

## 2. Explicit covector and mixed spatial jet recurrences

Let `y` be the initial position in the annulus, `F=D_y Phi_t(y)`,
`k=F^{-T}xi`, and `kappa<=|xi|<=R`. Set

\[
 P(k)=kk^T/|k|^2,\quad v(k)=k/|k|^2,\quad
 L=-A+2PA,\quad B(0,y,\xi)=(I-P(\xi))e,\quad |e|=1.          \tag{5}
\]

`e` is a fixed real vector. `B` solves the Kelvin equation and is transverse. For `r>=0`, the
following constants bound the r-th Frechet derivatives of the rational
functions:

\[
 \begin{array}{c|rrrrrr}
 r&0&1&2&3&4&5\\\hline
 C_r&1&4&20&144&1392&16800\\
 V_r&1&3&14&102&984&11880
 \end{array}                                                \tag{6}
\]

Specifically `||D^r P(k)||<=C_r|k|^{-r}` and
`||D^r v(k)||<=V_r|k|^{-r-1}`. These follow by differentiating
`|k|^2P=kk^T` and `|k|^2v=k`: for `r>=3`, both recurrences are
`Z_r=2rZ_{r-1}+r(r-1)Z_{r-2}`, with the explicitly listed initial values.

Write `z=1+t`, `h=46656/kappa`. For `0<=r<=4`, define rational constants

\[
 \begin{split}
 E_0&=1,\\
 E_r&=\left[C_r+12\sum_{j=1}^r
 {r\choose j}\frac{C_jE_{r-j}}{3r-j+1}\right] \quad(r\ge1),\\
 F_r&=(\mathbf1_{r=0}+2C_r)1620
   +12(C_{r+1}+rC_r)216^3\,286992\,R/\kappa,\\
 G_r&=\left[12\sum_{j=1}^r
 {r\choose j}\frac{C_jG_{r-j}}{3r-j+7}
 +\sum_{j=0}^r {r\choose j}\frac{F_jE_{r-j}}{3r-j+6}\right].
 \end{split}                                                \tag{7}
\]

Then the actual Kelvin symbol satisfies

\[
 \boxed{\|D_\xi^r B\|\le e^{6t}h^rE_rz^{3r}},\qquad
 \boxed{\|D_xD_\xi^r B\|\le216e^{6t}h^rG_rz^{3r+7}}.          \tag{8}
\]

Here the first bound is along `x=Phi_t(y)`, and the second is the
Eulerian spatial derivative. To verify (7), put `eta=J1^2/kappa`.
For `r>=1`, `||Dxi^r L||<=12C_r eta^r`. In Lagrangian coordinates,
`||D_y F^{-T}||<=J1^2 J2`; consequently

\[
 \|D_yD_\xi^r L\|\le h^rF_r z^{2r+5}.                       \tag{9}
\]

The first term in `F_r` bounds the derivative of `A(Phi_t(y))`; its
coefficient `1620` is `(15/2)*216`. The second follows from

\[
 \|D_yD_\xi^rP\|
 \le(C_{r+1}+rC_r)\eta^r J_1^3J_2R/\kappa.                  \tag{10}
\]

Apply (3) in Duhamel's formula. Integrating `z^{3r-j}` produces the
denominator `3r-j+1` in `E`; integrating the two mixed-derivative
channels produces `3r-j+7` and `3r-j+6` in `G`. All resulting powers are
bounded by `z^{3r}` and `z^{3r+6}`, respectively. Passing from `y` to `x`
costs one `J1<=216z`. The initial y derivative is zero. The initial xi
derivatives cost `C_r kappa^{-r}<=C_r h^r`.

## 3. Uniform synthesis without a beam-count loss

Take a unit vector `xi_c`, `r0=1/4`, and

\[
 \xi(\theta)=\xi_c+r_0(\sin\theta_1,\sin\theta_2,\sin\theta_3).
                                                               \tag{11}
\]

This parameter torus stays within `1/2<|xi|<2`. Every xi in the closed
cube `xi_c+[-1/4,1/4]^3` has a theta label; choose arcsines coordinatewise.
Thus no unproved parameter extension is required. For any symbol `s`,
the chain rule gives

\[
 \|(1-\Delta_\theta)^2s(\xi(\theta))\|
 \le S_0+9r_0S_1+33r_0^2S_2+30r_0^3S_3+9r_0^4S_4
 \le216\sum_{j=0}^4S_j/j!,                                \tag{12}
\]

where `S_j=sup ||Dxi^j s||`. The same bound holds for matrix and tensor
symbols with their operator norms.

Fourier integration by parts on the theta torus yields

\[
 \sum_{\ell\in\mathbb Z^3}\|\widehat s_\ell\|_\infty
 \le45\|(1-\Delta_\theta)^2s\|_\infty.                     \tag{13}
\]

Indeed a max-norm shell n contains `24n^2+2` lattice points and
`sum n^{-2}<5/3`, `sum n^{-4}<5/3`, giving
`1+(24+2)*5/3=133/3<45`. Each theta Fourier term is multiplication by
`s_hat_ell(y)` followed by a scalar unit-modulus carrier multiplier.

For a torus envelope eta with Fourier support `[-L,L]^3` and carriers
`p in p0+d Z^3`, `p0 in Z^3`, integer `d>2L`, Parseval gives exactly

\[
 \left\|\eta\sum_pa_p e^{i\ell\cdot\theta_p}e^{ip\cdot y}
 \right\|_2=\|\eta\|_2\|a\|_{\ell^2}.                     \tag{14}
\]

For R3 let `chi` be real and in `Cc^infty((-1,1)^3)` and
`chi_eps(y)=eps^{-3/2}chi((y-y0)/eps)`. If `d eps>=1`, periodizing
`|chi_eps|^2` over the carrier cell gives

\[
 \left\|\chi_\varepsilon\sum_pa_p
 e^{i\ell\cdot\theta_p}e^{ip\cdot y}\right\|_2
 \le27\|\chi\|_\infty\|a\|_{\ell^2}.                     \tag{15}
\]

The unsimplified constant is `(2+2pi/(d eps))^{3/2}`. Each cell sees at
most `(d eps/pi+1)^3` support translates; `pi<22/7` makes (15) valid.
The same estimate for a derivative of chi has an additional eps^{-1}.
Volume-preserving Lagrangian pullback changes neither norm. Combining
(12)--(15) gives the fully uniform constant `262440=216*45*27` for the
compact case. No factor depending on the number of carriers occurs.

## 4. Explicit uniform first-order residual majorant

For a symbol T, use the submultiplicative jet norm

\[
 \mathcal J(T)=\sum_{r=0}^4\frac1{r!}\sup\|D_\xi^rT\|.       \tag{16}
\]

Define the following polynomials with nonnegative rational coefficients,
using (6)--(7), `kappa=1/2`, and `R=2`:

\[
 \begin{split}
 P(z)&=\sum_{r=0}^4 C_r h^r z^{2r}/r!,\\
 P_x(z)&=\frac{216\,286992R}{\kappa}
   \sum_{r=0}^4(C_{r+1}+rC_r)h^rz^{2r+3}/r!,\\
 V(z)&=\frac{216}{\kappa}
   \sum_{r=0}^4V_rh^rz^{2r+1}/r!,\\
 V_x(z)&=\frac{216^2\,286992R}{\kappa^2}
   \sum_{r=0}^4(V_{r+1}+rV_r)h^rz^{2r+4}/r!,\\
 B(z)&=\sum_{r=0}^4h^rE_rz^{3r}/r!,\\
 B_x(z)&=216\sum_{r=0}^4h^rG_rz^{3r+7}/r!.
 \end{split}                                                \tag{17}
\]

These bound the respective jet norms, with the factor `e^{6t}` removed
only from the two B polynomials. The derivative bounds for P and v follow
by differentiating `k=DY_t(x)^Txi`, using `||Dx k||<=J2 R` and (6).
In particular the second term of each product derivative is absorbed
using `kappa/(J1^2 R)<=1`.

Let

\[
 \begin{split}
 R_0(z)&=[(225/2)P+90P_x+105/2]VB
       +[90P+72](V_xB+VB_x),\\
 R_1(z)&=[90P+72]V(216z)B,\\
 K_0&=262440R_0(1),\qquad K_1=787320R_1(1).
 \end{split}                                                \tag{18}
\]

The exact-arithmetic checker also verifies the convenient integer bounds
`K0<=10^94` and `K1<=10^80`. These conservative constants can be used
directly in (21).

For clarity the residual algebra behind (18) is explicit. C194 has
`c=-v cross b`, `pi=2i v^TAb`, and

\[
 h_c=D_tc=-2\operatorname{tr}(AP)(v\times b)
            +(A^Tv)\times b+v\times Ab.                    \tag{19}
\]

Use `|tr(AP)|<=3D||P||`, `||curl c||<=sqrt6||Dx c||<5/2||Dx c||`,
and bound the curl commutator and `A curl c` together by
`5D||Dx c||`. Product differentiation gives

\[
 \begin{split}
 \mathcal J(D_td+Ad+\nabla\pi)
 \le{}&[15A_1\mathcal J(P)+15D\mathcal J(D_xP)+7A_1]
         \mathcal J(v)\mathcal J(b)\\
 &+[15D\mathcal J(P)+12D]
       [\mathcal J(D_xv)\mathcal J(b)
        +\mathcal J(v)\mathcal J(D_xb)].
 \end{split}                                                \tag{20}
\]

For `b=chi_eps(y)B`, separate the residual into chi_eps times a symbol
and the three first derivatives of chi_eps times symbols. The former is
bounded by `e^{6t}R0(z)`. Each latter symbol is bounded by
`e^{6t}R1(z)`, since `||Dx y||<=216z`. The degrees of R0 and R1 are at
most 36 and 30. This explains both powers and the factor three in K1.

For finite carriers `p=xi_p/hbar in p0+d Z^3` with all xi_p in the cube
(11), define each C194 phase `phi_p=xi_p dot Y_t(x)` and the exact curl
approximation, then sum with scalar coefficients a_p. Assume
`d eps>=1`, `eps<=1/40`, and chi is real with support in the unit ball, and `f(y0)=0`, so
the support stays in the invariant annulus. Let u be the exact global
linearized Euler solution with the same initial datum as the sum. Then

\[
 \boxed{
 \|u(t)-u_{\rm app}(t)\|_2
 \le\hbar e^{6t}\|a\|_{\ell^2}
 \left[K_0(1+t)^{37}\|\chi\|_\infty
 +K_1\varepsilon^{-1}(1+t)^{31}\|\nabla\chi\|_\infty\right].}
                                                               \tag{21}
\]

The Leray energy estimate costs `e^{6(t-s)}` and the residual costs
`e^{6s}`, giving precisely `e^{6t}`. The constants in (18) are explicit
rational numbers, computable without floating point. Their size is not
hidden: they are extremely conservative, but independent of hbar, eps,
the carrier count, and the chosen coefficients. Equation (21) is complex
linear. Taking `2 Re(u_app,+)` and the corresponding exact solution gives
a real solution and costs at most exactly twice the displayed bound in
terms of `||a_plus||_l2`, or a factor `sqrt(2)` relative to the combined
positive-and-negative coefficient norm. No orthogonality of transported
positive and negative output bands is assumed.

## 5. Clock consequence and boundary

For `t<=19g/50 log q+152/25`, `hbar` a fixed constant times `q^{-1}`,
and fixed kappa/R/parameter cube, (21) retains the exact exponential
factor `q^{57g/25}`. The new costs are powers of `1+log q` and the
displayed constants, not a power of the number of beams. At `g=3/8`,
the two C196 first-order powers are again `37/600` and `14/75` after
the stated stable/expanding normalization, provided the putative same
solution supplies those signals. This is a sufficient upper-error
statement; finite q must satisfy the explicit inequality obtained from
(21), and no signal lower bound follows from it.

### 5.1 An explicit nonempty finite-q upper-error budget

The constants above give a concrete, deliberately conservative threshold.
Let `Achi=max(||chi||_infty,||grad chi||_infty)`, assume
`Achi ||a_ref||_l2>0`, take `hbar<=q^{-1}`, and retain all of the geometric and common-lattice assumptions of (21).
Suppose

\[
 q\ge10^{10000},\qquad
 0\le t\le\frac{57}{400}\log q+\frac{152}{25},              \tag{22}
\]

where log is natural. For the stable-width test set
`eps=q^{-1/12}` and `a=a_ref`. For the expanding-width test set
`eps=q^{-1/3}` and `a=q^{-3/8}a_ref`. These are two separate upper-error
comparisons. In both cases (21) implies, for a positive carrier box,

\[
 \|u-u_{\rm app}\|_2
 \le2\cdot10^{112}(\log q)^{37}q^{-37/600}
       A_\chi\|a_{\rm ref}\|_{\ell^2}.                     \tag{23}
\]

The real completion `2 Re(u_app,+)` consequently satisfies

\[
 \boxed{\|u_{\mathbb R}-u_{{\rm app},\mathbb R}\|_2
 <\frac1{100}A_\chi\|a_{\rm ref}\|_{\ell^2}.}               \tag{24}
\]

Here is the entire numerical comparison. Put `x=log q`. Equation (22)
gives `1+t<=x` and `exp(6t)<=10^{18}q^{171/200}`; the latter follows
from `exp(912/25)<3^{37}<10^{18}`. The stable test's two decay powers
are `29/200` and `37/600`. The expanding test's powers, including its
specified coefficient scaling, are `13/25` and `14/75`. Each is at least
`37/600`, and `K0,K1<=10^{94}`. This proves (23).

The function `x^{37}exp(-37x/600)` is decreasing for `x>=600`. At
`q0=10^{10000}`, `20000<log q0<30000`, so

\[
 \begin{split}
 4\cdot10^{112}(\log q)^{37}q^{-37/600}
 &\le4\cdot10^{112}(\log q_0)^{37}q_0^{-37/600}\\
 &<4\cdot10^{112}10^{185}10^{-616}
 =4\cdot10^{-319}<\frac1{100}.                              \tag{25}
 \end{split}
\]

The factor four includes the real-completion cost. The denominator in
(24) is a declared coefficient/profile scale. It is neither the exact
initial physical L2 norm nor a certified growing signal. Thus (24)
certifies a nonempty upper-error budget without claiming concentration,
fixed-energy growth, or a successful endpoint.

The theorem does not put a fixed-aperture family inside C195's narrow
expanding tube, produce C196's periodic endpoint under the flow, retain
an output Fourier band, or treat viscosity. The compact and torus
synthesis statements have distinct envelope hypotheses; (21) uses the
compact R3 hypothesis only.
