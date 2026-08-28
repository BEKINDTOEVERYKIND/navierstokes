# C194: an exactly solenoidal pressure-resolved \(A_2\) WKB beam has an explicit first-order error

**Date:** 2026-08-28

**Status:** exact curl/pressure identity and explicit local \(L^2\)
linearized-Euler error on the \(A_2\) invariant annulus; the action--angle
and residual chains are displayed and their arithmetic is checked; no
periodic integer carrier, off-ray Floquet bundle, stable
ballast, retained band, fixed-energy concentration, viscosity, nonlinear
stage, UVSR, or singularity claim

**Checker:**
[A2 curl/WKB bridge](../checks/a2_curl_wkb_bridge_c194.py)

## 0. The proved bridge and its boundary

Let

\[
 U=N\times\nabla f-\sqrt2 fN,\qquad
 f=\cos a+\cos b+\frac45\cos(a+b),\qquad A=DU.          \tag{0.1}
\]

Work on \(\mathbb R^3\), with Lebesgue \(L^2\), taking the periodic \(A_2\)
field (0.1) as a smooth bounded coefficient field.  This fixes a genuine
global linearized-Euler pressure problem rather than treating a local chart
as a standalone PDE.  Suppose the transported support lies in the invariant
annulus

\[
                              |f|\le\frac1{10}.          \tag{0.2}
\]

This condition persists because \(U\cdot\nabla f=0\).  Let
\(\phi\) solve \(D_t\phi=0\), put \(k=\nabla\phi\), and let a transverse
amplitude \(b\) solve

\[
 D_tk=-A^Tk,\qquad
 D_tb=-Ab+2k\frac{k^TAb}{|k|^2},\qquad k\cdot b=0.      \tag{0.3}
\]

For \(0<\hbar\le1\), define

\[
 c=-\frac{k\times b}{|k|^2},\qquad
 d=\frac1i\operatorname {curl}c,\qquad
 \pi=2i\frac{k^TAb}{|k|^2},                            \tag{0.4}
\]

and

\[
 v_{\rm app}=\frac{\hbar}{i}\operatorname {curl}
                  \left(e^{i\phi/\hbar}c\right)
             =e^{i\phi/\hbar}(b+\hbar d),\qquad
 p_{\rm app}=\hbar e^{i\phi/\hbar}\pi.                 \tag{0.5}
\]

Then \(v_{\rm app}\) is exactly divergence free, not merely divergence
free to leading order, and direct product-rule cancellation gives

\[
 (D_t+A)v_{\rm app}+\nabla p_{\rm app}
  =\hbar e^{i\phi/\hbar}
       \left(D_td+Ad+\nabla\pi\right).                  \tag{0.6}
\]

Fix \(x_0\) on the zero level, a unit covector \(|k_0|=1\),
\(b_*\in k_0^\perp\), and
\(\chi\in C_c^\infty(\mathbb R^3)\) supported in the unit ball.  For
\(0<\varepsilon\le1/40\), set

\[
 \phi(0,x)=k_0\cdot(x-x_0),\qquad
 \chi_\varepsilon(x)=\varepsilon^{-3/2}
       \chi\!\left(\frac{x-x_0}{\varepsilon}\right),\qquad
 b(0,x)=b_*\chi_\varepsilon(x).                        \tag{0.7}
\]

Indeed
\(\|\nabla_xf\|_\infty\le(14/5)\sqrt2<4\), so the unit-ball support,
\(f(x_0)=0\), and \(\varepsilon\le1/40\) place the initial support in
(0.2).  Invariance then keeps it there.
Here \(x_0,k_0,b_*,\chi\) are fixed independently of the carrier parameter
introduced in Section 3.  Let \(v\) be the exact \(\mathbb R^3\)
linearized-Euler solution with \(v(0)=v_{\rm app}(0)\).  Then

\[
\begin{aligned}
 \|v(t)-v_{\rm app}(t)\|_2
 \le \frac{\hbar e^{6t}|b_*|}{|k_0|}\bigg[&
 4{,}199{,}040\,\varepsilon^{-1}(1+t)^3
                  \|\nabla\chi\|_2\\
 &+2{,}898{,}006{,}000{,}000{,}000\,(1+t)^7
                  \|\chi\|_2\bigg].                   \tag{0.8}
\end{aligned}
\]

Thus C192's proposed first-order exponent is genuinely attained:

\[
                  \boxed{\Gamma=6<\frac{350}{57}}       \tag{0.9}
\]

with exact margin \(8/57\).  This is a positive pressure-resolving PDE
estimate, not just a bound for \(DU\).

It is not yet the C192 same-witness theorem.  In particular, (0.8) does not
prove that a three-dimensional off-ray packet lies in C192's expanding
bundle, nor that the broad C193 ballast lies in the covariant stable bundle.
It also does not make the generally nonintegral C159 carrier periodic on the
torus.  Accordingly no \(q^{3/8}\) endpoint growth is asserted here.

## 1. The exact \(A_2\) annulus and flow constants

Writing \(p=(a+b)/2\), \(X=\cos^2p\), and \(E=f\), direct elimination of
\(\cos((a-b)/2)\) gives

\[
\begin{aligned}
 50X\left(|\nabla_{(a,b)}f|^2-\frac32\right)
  ={}&-128X^3+164X^2-43X+16\\
    &+(25-50X)E^2+40E.                                \tag{1.1}
\end{aligned}
\]

On \(0<X\le1\), the first four terms are at least \(9\); on
\(|E|\le1/10\), the last two are at least \(-17/4\).  Therefore

\[
       \boxed{|\nabla_{(a,b)}f|^2\ge\frac{319}{200}}.   \tag{1.2}
\]

The case \(X=0\) is disjoint from the annulus because then \(E=-4/5\).
If
\(\Phi_t\) is the physical flow of \(U\), define, for \(r=1,2,3\),

\[
 J_r(t)=\sup_{|\tau|\le t}
 \max\{\|D^r\Phi_\tau\|_{L^\infty(\mathcal A)},
        \|D^r\Phi_\tau^{-1}\|_{L^\infty(\mathcal A)}\},
 \qquad \mathcal A=\{|f|\le1/10\}.                     \tag{1.3}
\]

Here is the exact action--angle derivation; the checker verifies its
constant ledger.  Put

\[
 y=Px=(a,b,z),\qquad z={N\cdot x\over3},\qquad
 P=\begin{pmatrix}1&-1&0\\0&1&-1\\1/3&1/3&1/3\end{pmatrix}.
 \tag{1.3a}
\]

The eigenvalues of \(PP^T\) are \(3,1,1/3\), so
\(\|P\|=\|P^{-1}\|=\sqrt3\).  With

\[
 g=\nabla_{(a,b)}f,\quad H=D^2_{(a,b)}f,\quad h=|g|^2,
 \quad J=\begin{pmatrix}0&-1\\1&0\end{pmatrix},
\]

the physical flow becomes

\[
             \dot y=(3Jg,-\sqrt2 f).                  \tag{1.3b}
\]

On the annulus, (1.2) gives \(h>3/2\).  The elementary identity

\[
 h=2\sin^2p\left(\cos q+{8\over5}\cos p\right)^2
       +2\cos^2p\sin^2q,qquad q={a-b\over2},
\]

and
\(r(1-r^2)\le2/5\), \(r^2(1-r^2)\le1/4\) for
\(0\le r\le1\), give

\[
 h\le2+{64\over25}+{32\over25}={146\over25}<6.
 \tag{1.3c}
\]

For a unit vector \(u\), the modal formula gives
\[
 |u^THu|\le u_1^2+u_2^2+{4\over5}(u_1+u_2)^2
       \le1+{8\over5}={13\over5}.
\]
For unit \(u,v,w\), the two coordinate modes cost at most \(2\), while
the diagonal mode costs at most
\((4/5)(\sqrt2)^3<12/5\).  Thus the three Fourier modes give

\[
             \|H\|_{\rm op}\le{13\over5},\qquad
             \|DH\|_{\rm mult}\le{22\over5}.         \tag{1.3d}
\]

Define the intrinsic frame and its scalar shear by

\[
 e={g\over h},\qquad
 Q(y)=\big[(3Jg,0),(e,0),e_z\big],qquad
 a_0(y)={g^THg-(Jg)^TH(Jg)\over h^2},                 \tag{1.3e}
\]

\[
 \alpha(t,y_0)=\int_0^t a_0(y_s)\,ds,qquad
 S(t,y_0)=
 \begin{pmatrix}
 1&\alpha&0\\0&1&0\\0&-\sqrt2t&1
 \end{pmatrix}.                                      \tag{1.3f}
\]

Since \(df(e)=1\), direct differentiation gives

\[
 D(3Jg)e-De(3Jg)=a_0(3Jg).
\]

Consequently the full three-dimensional variation has the exact factorization

\[
                    D_y\Psi_t(y_0)=Q(y_t)S(t,y_0)Q(y_0)^{-1}.
                                                            \tag{1.3g}
\]

Indeed the first and third frame columns return to their corresponding
columns, while the second maps to
\((e_t,0)+\alpha(3Jg_t,0)-\sqrt2t e_z\).  No undefined angle coordinate
is being used here.

Orthogonality of the columns of \(Q\), (1.3c)--(1.3d), and
\(D(Q^{-1})=-Q^{-1}(DQ)Q^{-1}\) give

\[
 \|Q\|\le3\sqrt6,\quad \|Q^{-1}\|\le\sqrt6,\quad
 \|DQ\|\le8+6=14,\quad \|D(Q^{-1})\|\le84.           \tag{1.3h}
\]

Moreover

\[
 |a_0|\le {2\|H\|\over h}<{7\over2},
\]

and, writing the numerator of \(a_0\) as \(n_0\),
\(|n_0|\le2\|H\|h\),
\(|Dn_0|\le4\|H\|^2\sqrt h+2\|DH\|h\), and
\(|Dh|\le2\|H\|\sqrt h\).  Since
\(h^{3/2}>(3/2)(6/5)=9/5\),

\[
 \|Da_0\|\le {12(13/5)^2\over9/5}
                  +{2(22/5)\over3/2}
             ={764\over15}<51.                        \tag{1.3i}
\]

Thus \(\|S^{\pm1}\|\le1+4t\), and (1.3g) first gives the
phase-coordinate bounds

\[
 \|D_y\Psi_{\pm t}\|\le18(1+4t),\qquad
 \|D\alpha\|\le918(t+2t^2).                          \tag{1.3j}
\]

Differentiating \(Q_tSQ_0^{-1}\), the three product-rule terms are bounded
respectively by

\[
 756(1+4t)^2,\qquad16524(t+2t^2),\qquad672(1+4t).
 \tag{1.3k}
\]

Finally \(\Phi_t=P^{-1}\Psi_tP\), so first derivatives pay the exact
factor \(3\), while second derivatives pay
\(3\sqrt3<6\).  This proves, for positive and negative time, the following
physical flow-jet bounds:

\[
 J_1(t)\le54(1+4t)\le216(1+t),                         \tag{1.4}
\]

\[
\begin{aligned}
 J_2(t)&\le6\left[756(1+4t)^2
              +16524(t+2t^2)+672(1+4t)\right]\\
       &=8568+151560t+270864t^2
        \le286992(1+t)^2.                              \tag{1.5}
\end{aligned}
\]

For future second-order work it also verifies the auxiliary, presently
non-load-bearing ledger

\[
                   J_3(t)\le4{,}031{,}918{,}208(1+t)^3.\tag{1.6}
\]

C193 supplies the exact physical coefficient bounds

\[
 \|A\|_{\rm op}\le6,\qquad
 \|DA\|_{\rm mult}\le3\sqrt6,\qquad
 \|D^2A\|_{\rm mult}\le9.                              \tag{1.7}
\]

Vectors use Euclidean norm, matrices induced Euclidean operator norm, and
higher derivatives multilinear operator norm.  Curl estimates use the
Frobenius norm and explicitly pay
\(|\operatorname {curl}c|\le\sqrt2|Dc|_F\).  Since the \(A_2\) field is
axially invariant and \(AN=0\), the associated rank/tensor conversions are
already included in the numerical coefficients of (0.8); no uncharged
coordinatewise norm is substituted for (1.7).

## 2. Curl and pressure cancellation

Because \(k\cdot b=0\), (0.4) has

\[
                         k\times c=b.                   \tag{2.1}
\]

Expanding the curl in (0.5) therefore gives its second equality.  At order
\(\hbar^0\), equations (0.3) and

\[
 i k\pi=-2k\frac{k^TAb}{|k|^2}                         \tag{2.2}
\]

cancel the full vector equation.  The unexpanded curl makes divergence
zero exactly.  What remains is precisely (0.6), with no Leray remainder and
no omitted pressure channel.

For reference, if

\[
 D=\|A\|_{L^\infty,{\rm op}},\quad
 A_1=\|DA\|_{L^\infty,{\rm mult}},\quad
 K_- =\inf|k|,\quad K_1=\|Dk\|_{L^\infty,{\rm op}},
 \quad B_0=\|b\|_2,\quad B_1=\|Db\|_{L^2,F},          \tag{2.3}
\]

then the complete product-rule chain is as follows.  For every spatial
direction \(\xi\), with \(r=|k|\),

\[
 Dc[\xi]=-{Dk[\xi]\times b+k\times Db[\xi]\over r^2}
       +{2(k\cdot Dk[\xi])(k\times b)\over r^4}.       \tag{2.4}
\]

If \(K_{1,F}=\|Dk\|_{L^\infty,F}\), this gives

\[
 C_1:=\|Dc\|_2\le {B_1\over K_-}
                    +{3K_{1,F}\over K_-^2}B_0.         \tag{2.5}
\]

Direct use of (0.3) gives, with \(\beta=k^TAk\),

\[
 h_c:=D_tc=-{2\beta\over r^4}k\times b
       +{(A^Tk)\times b+k\times Ab\over r^2}.          \tag{2.6}
\]

Thus, if
\[
 A_{1,F}=\sup_x\left(\sum_j\|\partial_jA(x)\|_{\rm op}^2\right)^{1/2},
\]

\[
 \|Dh_c\|_2\le H_0B_1+H_1B_0,\qquad
 H_0={4D\over K_-},\qquad
 H_1={4A_{1,F}\over K_-}+{20DK_{1,F}\over K_-^2}.     \tag{2.7}
\]

The constants are not schematic: differentiating the three summands in
(2.6) gives the \((Db,DA,Dk)\) coefficient triples
\((2,2,14),(1,1,3),(1,1,3)\).

The material derivative does not commute with curl.  In components,

\[
 D_t(\operatorname {curl}c)_i
  =(\operatorname {curl}h_c)_i
       -\epsilon_{ijk}A_{\ell j}\partial_\ell c_k.     \tag{2.8}
\]

Here \(|\operatorname {curl}z|\le\sqrt2|Dz|_F\), and
\(AN=0\) gives \(|A|_F\le\sqrt2D\).  The commutator in (2.8) and the
term \(A\operatorname {curl}c\) are therefore each bounded by
\(\sqrt2D|Dc|_F\).  Finally, direct differentiation of \(\pi\) gives

\[
 \|\nabla\pi\|_2\le {2D\over K_-}B_1
   +\left({2A_1\over K_-}
          +{6DK_1\over K_-^2}\right)B_0.              \tag{2.9}
\]

The \(DA\) and \(Dk\) channels in (2.9) are scalar directional linear
functionals; their gradient norm is therefore controlled directly by the
multilinear/operator norms \(A_1,K_1\), with no Frobenius conversion.
The coefficient field and transported phase gradient are axially
invariant.  Hence \(DA[N]=Dk[N]=0\), and the Frobenius norm of either
derivative map is at most \(\sqrt2\) times its multilinear/operator norm.
Applying this
conversion term by term in (2.5)--(2.9), rather than silently identifying
the norms, yields the sharp-interface estimate

\[
\boxed{
 \|D_td+Ad+\nabla\pi\|_2
 \le (6\sqrt2+2){D\over K_-}B_1
      +\left({10A_1\over K_-}
              +{58DK_1\over K_-^2}\right)B_0.}        \tag{2.10}
\]

In particular \(6\sqrt2+2<25/2\), so the wholly rational version is

\[
 \|D_td+Ad+\nabla\pi\|_2
 \le {25D\over2K_-}B_1
      +\left({10A_1\over K_-}
              +{58DK_1\over K_-^2}\right)B_0.         \tag{2.11}
\]

The unrounded first coefficient in (2.10) is retained in the next
substitution; that is what pays the remaining Frobenius conversion without
changing the headline constants.

### 2.1 Transported phase and amplitude derivatives

Let \(Y_t=\Phi_{-t}\).  Factor the full amplitude as

\[
 b(t,x)=\chi_\varepsilon(Y_t(x))\,\widetilde b(t,x),
 \qquad k(t,x)=DY_t(x)^Tk_0,                           \tag{2.12}
\]

where \(\widetilde b(0,x)=b_*\).  From (1.3),

\[
 K_-^{-1}\le {J_1\over|k_0|},\qquad
 K_1\le |k_0|J_2.                                     \tag{2.13}
\]

Put \(q=k/|k|\) and \(P_q=I-q\otimes q\).  Then
\(\|Dq\|_{\rm op}\le J_1J_2\), while transversality gives

\[
                         |\widetilde b(t)|\le e^{6t}|b_*|. \tag{2.14}
\]

For completeness, take a material direction \(h_s\), put
\(r_q=D_{h_s}q\), \(z=D_{h_s}\widetilde b\), and
\(w=P_qz\).  Differentiating the Kelvin equation and projecting gives the
exact covariant equation

\[
 \nabla_t^\perp w=-P_qAw-P_q(DA[h_s])\widetilde b
   +2r_q(q\cdot A\widetilde b)
   +(r_q\cdot\widetilde b)(P_qAq+\dot q),
 \qquad \dot q=-P_qA^Tq.                              \tag{2.15}
\]

The last two terms in (2.15) cost at most
\(4D|r_q||\widetilde b|\).  Back-propagating a unit final direction costs
\(J_1(t)\), and the transverse homogeneous equation costs only
\(e^{6(t-s)}\).  Therefore, with

\[
 \mathcal D_b(t)=\int_0^t
       \left[A_1J_1(s)+4D J_1(s)J_2(s)\right]ds,       \tag{2.16}
\]

one has

\[
 \|P_qD\widetilde b\|_{\rm op}
 \le e^{6t}|b_*|J_1(t)\mathcal D_b(t).                \tag{2.17}
\]

There is also an essential algebraic longitudinal derivative:

\[
 D_h\widetilde b=P_qD_h\widetilde b
       -q\big((D_hq)\cdot\widetilde b\big).            \tag{2.18}
\]

It follows that

\[
 \|D\widetilde b\|_{\rm op}
 \le e^{6t}|b_*|J_1(t)\big(\mathcal D_b(t)+J_2(t)\big).\tag{2.19}
\]

The derivative in (2.19) has the axial direction in its kernel, so its
Frobenius norm is at most \(\sqrt2\) times its operator norm.  Volume
preservation and (0.7) now give

\[
 B_0\le e^{6t}|b_*|\|\chi\|_2,                        \tag{2.20}
\]

\[
 B_1\le e^{6t}|b_*|\left[
 J_1\varepsilon^{-1}\|\nabla\chi\|_2
 +\sqrt2J_1(\mathcal D_b+J_2)\|\chi\|_2\right].      \tag{2.21}
\]

Substituting (2.13), (2.20), and (2.21) into the unrounded estimate
(2.10), and using \(D=6\), \(A_1=3\sqrt6\), gives

\[
\begin{aligned}
 \|D_td+Ad+\nabla\pi\|_2
 \le {e^{6t}|b_*|\over|k_0|}\bigg[&
 90J_1^2\varepsilon^{-1}\|\nabla\chi\|_2\\
 &+\left(90J_1^2\mathcal D_b+30\sqrt6J_1
             +462J_1^3J_2\right)\|\chi\|_2\bigg].   \tag{2.22}
\end{aligned}
\]

Every conversion in (2.22) has explicit slack:

\[
 36\sqrt2+12<66<90,\qquad
 (36\sqrt2+12)\sqrt2=72+12\sqrt2<90,\qquad
 90+58\cdot6=438<462.                                 \tag{2.23}
\]

The last inequality is precisely where the longitudinal term (2.18) is
charged; it is not discarded.  The replacement of \(J_1^2J_2\) by
\(J_1^3J_2\) in (2.22) is harmless because \(J_1(t)\ge1\).

Using (1.4)--(1.5), \(3\sqrt6<15/2\), and \(D=6\), (2.16) yields

\[
 \mathcal D_b(t)\le371{,}942{,}442(1+t)^4.             \tag{2.24}
\]

Consequently the gradient coefficient in (2.22) is
\(90\cdot216^2=4{,}199{,}040\).  The three amplitude coefficients are

\[
\begin{gathered}
 4{,}199{,}040\cdot371{,}942{,}442
       =1{,}561{,}801{,}191{,}655{,}680,\\
 30\cdot{5\over2}\cdot216=16{,}200,\\
 462\cdot216^3\cdot286{,}992
       =1{,}336{,}204{,}776{,}259{,}584,
\end{gathered}                                        \tag{2.25}
\]

whose sum is
\(2{,}898{,}005{,}967{,}931{,}464
 <2{,}898{,}006{,}000{,}000{,}000\).  Finally the error solves a forced
linearized-Euler equation.  Its energy inequality costs
\(e^{6(t-s)}\), while (2.22) costs \(e^{6s}\); these exponentials multiply
to exactly \(e^{6t}\).  Bounding the two time integrals by
\((1+t)^3\) and \((1+t)^7\) proves (0.8).

## 3. The C192 clock arithmetic

Set the physical carrier parameter to \(\hbar=Q^{-1}\).  C193's
three-dimensional concentration scale has spatial width
\(\varepsilon=Q^{-1/4}\), so

\[
 \hbar=Q^{-1},\qquad \varepsilon=Q^{-1/4},\qquad
 \hbar\varepsilon^{-1}=Q^{-3/4}.                       \tag{3.1}
\]

The explicit support condition above holds for this choice when
\(Q\ge40^4\); for the narrower \(Q^{-1/2}\) test below it holds when
\(Q\ge1600\).

On C192's raw clock

\[
 t<\frac{57}{400}\log Q+\frac{76}{25},                 \tag{3.2}
\]

one has

\[
 e^{6t}<e^{456/25}Q^{171/200}.                         \tag{3.3}
\]

After division by a *hypothetical same-packet* \(Q^{3/8}\) signal, (0.8)
therefore gives the explicit arithmetic comparison

\[
\begin{aligned}
 \frac{\|v-v_{\rm app}\|_2}
      {( |b_*|/|k_0|)Q^{3/8}}
 <e^{456/25}\bigg[&4{,}199{,}040(1+t)^3
       \|\nabla\chi\|_2 Q^{-27/100}\\
 &+2{,}898{,}006{,}000{,}000{,}000(1+t)^7
       \|\chi\|_2 Q^{-13/25}\bigg].                   \tag{3.4}
\end{aligned}
\]

The first exponent is exactly

\[
 \frac34+\frac38-\frac{171}{200}=\frac{27}{100},       \tag{3.5}
\]

and the second is \(13/25\).  If one deliberately narrows the envelope to
\(\varepsilon=Q^{-1/2}\), then
\(\hbar\varepsilon^{-1}=Q^{-1/2}\) and the first exponent becomes

\[
 \frac12+\frac38-\frac{171}{200}=\frac1{50},           \tag{3.6}
\]

which is exactly C192's stress-test threshold.  Thus (0.8) is compatible
both with the \(Q^{-1/4}\) concentration width and with C192's narrower
test.  The qualifier in italics remains load-bearing: C194 supplies only
the error half, not the off-ray lower-growth half needed to divide by that
signal.

C193's filter pays one additional return.  Its clock is

\[
 t<\frac{57}{400}\log Q+\frac{152}{25},                \tag{3.7}
\]

so a future composite use must replace \(e^{456/25}\) in (3.4) by
\(e^{912/25}\).  This changes no power exponent.  The raw-clock constant
must not be quoted as the completed C193-filter constant.

## 4. What remains open

C194 does not certify:

1. an integer-carrier periodic phase or reality completion;
2. a quantitative expanding/stable Floquet splitting on the required
   three-dimensional shrinking cap;
3. initial or endpoint Fourier projection and transported-band retention;
4. the C193 fixed-energy \(L^\infty/L^2\) endpoint for the same exact
   solution;
5. viscous damping, active extraction, depletion, wake, C125, UVSR, or a
   Navier--Stokes singularity.

An exact second-corrector telescoping recurrence is available algebraically,
but (1.6) alone does not bound its forced stable component or its complete
residual.  No \(O(Q^{-1})\) second-order bridge is claimed.

## 5. Verification boundary

The dependency-free checker verifies the annulus identity, all rational
flow-jet and residual constants, the curl/pressure signs, and the C192 power
margins.  It does not integrate the C159 ray or an off-ray polarization
column, solve a periodic PDE, certify a band projector, or estimate a
second-corrector residual.
