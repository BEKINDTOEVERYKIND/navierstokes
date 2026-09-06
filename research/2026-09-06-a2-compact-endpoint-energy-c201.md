# C201: compact A2 endpoint geometry and physical energy normalization

**Date:** 2026-09-06

**Status:** exact return geometry, integer-lattice count, compact vector Riesz estimates, and exact two-solution energy balancing. The endpoint-energy estimates use C200's declared smooth-symbol bounds. The final conditional lemma explains precisely how independently verified PDE error bounds transfer these estimates to exact solutions. No retained Fourier band, viscous evolution, nonlinear return, or singularity is proved here.

**Checker:** [a2_compact_endpoint_energy_c201.py](../checks/a2_compact_endpoint_energy_c201.py).

All L2 norms use Lebesgue measure, on R3 or on a fundamental 2pi-periodic torus cell. A compact patch and its torus image have equal integrals. In the torus version all carriers below are integer vectors and the phase exponential is globally defined by flow pullback.

## 1. The reference return has rank-one shear

Let H=RT, L=1+H, and F=D Phi_H(y0) at the C159 returning reference orbit. C194's exact factorization is

\[
 F=P^{-1}Q(y_H)S(H)Q(y_0)^{-1}P,
 \quad S-I=(\alpha,0,-\sqrt2H)^T(0,1,0).
\]

On f=0 the axial velocity vanishes, while the planar reference returns modulo the physical torus. The frame Q depends on its returning planar gradient. Thus Q(y_H)=Q(y_0), and conjugation gives

\[
 F=I+a\otimes b,\quad a\cdot b=0,\quad
 F^{-1}=I-a\otimes b,\quad \det F=1,
 \quad \|a\otimes b\|\le\|F\|+1\le217L.                 \tag{1}
\]

Choose an orthonormal endpoint frame e1,e2,e3 with e3 parallel to a; if a=0, choose any frame. Put c=10^-100 and

\[
 w_1=w_2=c/L,\qquad w_3=c/L^2,\qquad
 A=F^T[e_1w_1,e_2w_2,e_3w_3].                           \tag{2}
\]

For v in the unit Euclidean ball, only its third coordinate contributes to a dot sum(ej wj vj). Consequently

\[
 \|A\|\le c/L+217Lc/L^2\le218c/L<220c/L.                \tag{3}
\]

The same elementary bound for a parameter cube is
`|A z| <=220 sqrt(3)c/L` when `|z_j|<=1`. Notice that this statement describes a box of **central endpoint covectors**. A nonlinear transported packet is not thereby proved to have Fourier support in that box.

## 2. Count a full three-dimensional integer grid

Let xi_c=k0,*/|k0,*|. The returning central covector obeys F^-T xi_c=xi_c. Fix

\[
 \varepsilon_+=q^{-1/3},\quad \varepsilon_-=q^{-1/12},\quad
 d=\lceil8\pi q^{1/3}\rceil,\quad h=d/q.                  \tag{4}
\]

Use all p in d Z3 for which xi_p=p/q belongs to xi_c+A[-1,1]^3; call this positive set S and write M=|S|. Its negative is disjoint from S, since (3) places S in the ball of radius 1/4 about the unit xi_c. Both signs belong to the **same** lattice d Z3. A generic coset p0+d Z3 would not justify the same signed-lattice Riesz proof.

Since pi<22/7,

\[
 d\le27q^{1/3},\qquad d\varepsilon_\pm\ge8\pi.          \tag{5}
\]

Suppose

\[
 11664\,c^{-1}L^3q^{-2/3}\le1.                          \tag{6}
\]

Every point in xi_c+A[-1/2,1/2]^3 is within Euclidean distance sqrt(3)h/2 of a nearest grid point. For every row of A^-1,

\[
 \|\operatorname{row}_j A^{-1}\|
 \le\|F^{-T}\|/w_j\le216c^{-1}L^3.
\]

Equations (5)-(6) place that grid point in xi_c+A[-1,1]^3. The grid cells around the selected points therefore cover the half box. Its volume is det(A)=w1 w2 w3, because det(F)=1. Hence

\[
 Mh^3\ge w_1w_2w_3,\qquad
 M\ge\frac{c^3q^2}{27^3L^4},\qquad
 \boxed{\varepsilon_+^{-3/2}\sqrt M
       \ge\frac{c^{3/2}}{162}\frac{q^{3/2}}{L^2}.}        \tag{7}
\]

The irrational orientation and shear have not been rounded away. Condition (6) pays for them. For the C193 clock and q>=10^10000, one has L<=log q; (6) follows by checking the decreasing function x^3 exp(-2x/3) at x=10000 log 10, using 20000<x<30000. The same threshold makes all later displayed small-correction requirements nonempty.

## 3. Compact scalar and vector Riesz bounds

Choose the following explicit smooth plateau. Let `theta(t)=exp(-1/t)` for t>0 and zero otherwise, let `S(t)=theta(t)/(theta(t)+theta(1-t))`, let `rho(s)=S(2-4|s|)`, and set `chi(y)=rho(y1)rho(y2)rho(y3)`. Then 0<=chi<=1, chi=1 on [-1/4,1/4]^3, and its support is contained in [-1/2,1/2]^3, strictly inside the unit ball. Smoothness at s=0 follows because rho is constant near zero. On (0,1), `theta'<=4/e^2<1` and `theta(t)+theta(1-t)>=e^-2>1/9`, so `|S'|<9`, `|rho'|<=36`, and `||grad chi||_infty<=36 sqrt(3)<72`. Set chi_eps(y)=eps^-3/2 chi((y-y0)/eps), and retain the explicit derivative constant `Achi=max(1,||grad chi||_infty)<=72`.

For any finite vector coefficients z_p on d Z3 and d eps>=8pi,

\[
 \boxed{\frac18\|z\|_{\ell^2}
 \le\left\|\chi_\varepsilon\sum_pz_p e^{ip\cdot(y-y_0)}\right\|_2
 \le27\|z\|_{\ell^2}.}                                 \tag{8}
\]

For the lower bound periodize |chi_eps|^2 over a cube of side 2pi/d. In each coordinate the plateau has at least floor(d eps/(4pi)) >=d eps/(8pi) representatives. The periodized weight is thus at least d^3/(8pi)^3. Parseval on that cell contributes (2pi/d)^3, giving squared lower bound 1/64. The upper bound is C197's explicit periodization estimate. The same argument is valid on the torus because d is an integer, the compact envelope stays in one initial coordinate patch, and the torus contains exactly d^3 such cells.

Now let a(y,xi) be any real unit symbol whose covector derivatives through order four satisfy

\[
 \|D_\xi^r a\|\le30 r!D^r,\qquad D=4\cdot10^{40}L,       \tag{9}
\]

on the full parameter image xi=xi_c+A sin(theta). This is an enlarged smooth-parameter box, not just the discrete set S; the admissible tube hypothesis must hold on it. C200 supplies (9) at the expanding endpoint and at the contracting entrance once that hypothesis has been checked.

Put eta=D||A||. By (3), eta<=10^43 c=10^-57. The difference
`s(y,theta)=a(y,xi_c+A sin(theta))-a(y,xi_c)` obeys

\[
 \|(1-\Delta_\theta)^2s\|
 \le30[\sqrt3\eta+9\eta+66\eta^2+180\eta^3+216\eta^4]
 \le14190\eta.                                         \tag{10}
\]

The coefficient list follows by differentiating the three coordinate sines, exactly as in C197; replacing each sine derivative vector by its bound ||A|| preserves that derivation. The zeroth term uses the mean-value theorem along the parameter cube, and eta<=1 bounds its higher powers.

Assign arcsine labels to the positive xi_p; assign the same label and the same real symbol to -p. For arbitrary scalar coefficients on the signed grid, C197's theta-Fourier estimate and (8) give

\[
 \left\|\chi_\varepsilon\sum_pz_p
    [a(y,xi_p)-a(y,xi_c)]e^{ip\cdot(y-y_0)}\right\|_2
 \le\delta_a\|z\|_{\ell^2},\quad
 \delta_a=17240850\eta<2\cdot10^{-50}.                  \tag{11}
\]

For negative p the notation a(y,xi_p) means the parity extension a(y,-xi_p). This is consistent with the Kelvin equation, which depends on covector direction modulo sign. The proof uses scalar unit-modulus multipliers on the entire signed lattice, so there is no reality-completion triangle loss.

Multiplication by the unit vector a(y,xi_c) preserves a scalar pointwise modulus. Consequently (8)-(11) prove the physical vector bound

\[
 \boxed{(1/8-\delta_a)\|z\|_{\ell^2}
 \le\left\|\chi_\varepsilon\sum_pz_pa(y,xi_p)e^{ip\cdot(y-y_0)}\right\|_2
 \le(27+\delta_a)\|z\|_{\ell^2}.}                       \tag{12}
\]

This is a lower bound for the actual physical principal field, not a declaration that coefficient norm equals velocity norm. Volume-preserving flow pullback preserves both bounds.

There is also a direct uniform scalar projection at the expanding endpoint. Set `e_*=a(y0,xi_c)`, a real unit vector. On the envelope support, (9), (14), and the mean-value theorem give

\[
 |a(y,xi_p)-e_*|
 \le10^{45}(R+1)L^4\varepsilon_++60\eta.                \tag{12a}
\]

If this displayed quantity is at most 1/2, then `e_* dot a(y,xi_p)>=1/2` simultaneously for every positive carrier and every support point. The same parity convention covers the negative carriers. This proves a common physical projection directly from the normalized symbol, without confusing a coefficient-cone slope with a Euclidean angle. The condition holds on the declared clock for q>=10^10000.

For a symbol of size at most 4 sigma and covector jets bounded by the right side of (9) times sigma, the same argument gives the simpler upper bound `(108+delta_a)sigma ||z||`. We will allow the convenient larger bound `811 sigma ||z||`; it also covers a zeroth-order bound 30 sigma.

## 4. Exact real curls and their physical correction

For a transverse real symbol a_p at a comparison section write k_p=D Phi_t^-T xi_p and c_p=-k_p cross a_p/|k_p|^2. The factor K used in C200 to describe its reference covector does not affect the Kelvin amplitude; this k_p is the gradient of the normalized phase. With phase phi_p=xi_p dot(Y_t-y0), use

\[
 V_{\rm app}=\frac1{iq}\operatorname{curl}
 \sum_{p\in S\cup(-S)}z_p\chi_\varepsilon(Y_t)
            c_p e^{iq\phi_p},\qquad z_p=(2M)^{-1/2}.      \tag{13}
\]

Set c_-p=-c_p and z_-p=z_p. Equation (13) is exactly real, periodic in the torus version, and divergence free. Its principal term is the field in (12) with combined coefficient norm one. Equivalently it is `2 Re` of the positive sum with coefficient `1/sqrt(2M)` per carrier. If an endpoint-test lemma is stated with positive coefficient `1/sqrt(M)`, its linear signal here is multiplied by `1/sqrt(2)>=1/2`; that factor must be retained when applying such a lemma.

Here is a deliberately conservative explicit correction bound. Suppose, in addition to (9),

\[
 \|D_yD_\xi^r a\|
 \le10^{45}(R+1)L^4 r!D^r                              \tag{14}
\]

as in C200. Define

\[
 \mathcal E_{\rm curl}(q,\varepsilon)
 =10^{250}(R+1)L^{20}A_\chi/(q\varepsilon).             \tag{15}
\]

Then the L2 norm of the exact-curl correction in (13) is at most E_curl. The same conclusion multiplied by sigma holds when (9),(14) and the symbol itself carry sigma.

For reproducibility, use J(T)=sum_{r=0}^4 ||Dxi^r T||/r!. With kappa=1/2, R_cov=2, and C194's J1,J2, C197's explicit rational derivative arrays give

\[
 J(v)\le10^{25}L^9,\quad J(D_xv)\le10^{35}L^{12},\quad
 J(a)\le10^{165}L^4,\quad
 J(D_xa)\le10^{211}(R+1)L^9,
\]

where v=k/|k|^2. Thus

\[
 J(c)\le10^{190}L^{13},\qquad
 J(D_xc)\le10^{237}(R+1)L^{18}.                         \tag{16}
\]

Product differentiation, `||curl||<=5/2||D_x||`, the three envelope derivatives, and the theta-synthesis constant 262440 bound the correction by

\[
 \frac{262440}{q}
 [\tfrac52 10^{237}(R+1)L^{18}
  +648\,10^{190}L^{14}\varepsilon^{-1}A_\chi],
\]

which is below (15). These loose powers and constants charge every derivative; they do not depend on M. They use ||A||<=1, so no small parameter has been presumed for this correction calculation.

In particular the terminally normalized expanding multibeam and the initially normalized contracting **single beam** satisfy

\[
 1/8-\delta_a-E_{\rm curl}\le\|V_{\rm app}\|_2
       \le27+\delta_a+E_{\rm curl}.                     \tag{17}
\]

For the single beam use any positive lattice carrier from S, its negative, and coefficient 1/sqrt(2). In the two small branches, with sigma=2999^-R, the upper bound is `(811+E_curl)sigma`. Thus one stable beam suffices to carry physical entrance energy; it need not be a second M-beam construction.

## 5. Exact energy balancing for two genuine PDE solutions

This lemma does not identify an approximate field with an exact solution. Let U,V be any two real solutions of the same linear PDE on [0,H]. Assume their **actual** endpoint norms satisfy

\[
 1/10\le\|U(H)\|_2,\|V(0)\|_2\le29,\qquad
 \|U(0)\|_2,\|V(H)\|_2\le10^{-6}.                      \tag{18}
\]

The inequalities follow from (17) only after the specific PDE approximation errors have also been bounded. Define

\[
 A=\|U(H)\|_2^2-\|U(0)\|_2^2>0,\quad
 B=\|V(0)\|_2^2-\|V(H)\|_2^2>0,
\]

\[
 C=\Re\langle U(H),V(H)\rangle
       -\Re\langle U(0),V(0)\rangle,\qquad
 \lambda=\frac{C+\sqrt{C^2+AB}}{B}>0.                   \tag{19}
\]

Then Z=U+lambda V is a solution and obeys the exact identity

\[
 \|Z(0)\|_2=\|Z(H)\|_2=:N.                            \tag{20}
\]

There is no orthogonality premise. Expanding the difference of squared norms gives A+2lambda C-lambda^2 B=0. Triangle inequalities in (20) give

\[
 \frac{\|U(H)\|_2-\|U(0)\|_2}{\|V(0)\|_2+\|V(H)\|_2}
 \le\lambda\le
 \frac{\|U(H)\|_2+\|U(0)\|_2}{\|V(0)\|_2-\|V(H)\|_2}.
\]

Substituting (18) yields

\[
 \boxed{1/300<\lambda<300,\qquad 1/20<N<30.}             \tag{21}
\]

For any desired b>0, `(b/N)Z` therefore has physical L2 norm **exactly b** at both endpoints. This does not assert no intermediate overshoot; a later stage budget must charge the proved trajectory bound. Nor does it create a concentration lower bound by itself: that lower bound must be established for Z, including the stable term and the PDE error, before division by N.

## 6. Scope

C201 proves that the C180-shaped central endpoint box contains enough separated, genuinely three-dimensional periodic integer carriers to attain the required L^-2 phase-space tax; that C200's normalized compact symbols have physical energy comparable to their declared coefficients with explicit constants; and that exact energy equality is available for actual growing/contracting solutions once (18) is verified. It supplies no Fourier-support theorem for the transported field. Its role is to remove the carrier arithmetic, coefficient-normalization, and exact endpoint-energy ambiguities from the next same-solution endpoint calculation.
