# C204: one physically normalized viscous A2 retained-band endpoint

**Date:** 2026-09-06

**Status:** finite-frequency, real, periodic, exactly solenoidal, actual
linearized-Navier--Stokes concentration theorem on the factorial stage
scales. The entrance and terminal physical energies agree exactly. This
is a linear perturbation theorem, not an unforced nonlinear scale return.

**Checker:** [a2_physical_viscous_retained_endpoint_c204.py](../checks/a2_physical_viscous_retained_endpoint_c204.py).

The proof uses C200--C203 on the same selected family. It does not use
C185's abstract operator norm as a packet multiplier. Norms on the
side-\(2\pi\) torus use unnormalized Lebesgue measure.

## 1. The actual PDE statement

Fix viscosity \(\nu>0\), an integer

\[
 n\ge\max(10^{12500},\nu),\quad q=n^8,\quad
 \mu=\nu((n-1)!)^{-3/2},\quad
 R=\left\lceil\frac3{64}\log q\right\rceil+1,
 \quad H=RT,\quad L=1+H,\quad c=10^{-100}.             \tag{1}
\]

Here \(T\) is the certified C159 period, \(3<T<76/25\). Let \(U\) be the
A2 Beltrami field of C194, and put

\[
 U_\mu(s)=e^{-2\mu s}U,\quad
 \tau(s)=\frac{1-e^{-2\mu s}}{2\mu},\quad
 S=-\frac{\log(1-2\mu H)}{2\mu}.                     \tag{2}
\]

The clock is well defined and

\[
 \tau(S)=H,\quad H\le S\le\frac{100}{99}H,\qquad
 e^{-2\mu S}\ge\frac{99}{100}.                       \tag{3}
\]

There exists a nonzero smooth real mean-zero divergence-free solution

\[
 \partial_s v+e^{-2\mu s}(U\cdot\nabla v+DU\,v)
       -\mu\Delta v+\nabla p=0                       \tag{4}
\]

for which, for any prescribed \(b>0\),

\[
 \boxed{\|v(0)\|_2=\|v(S)\|_2=b.}                    \tag{5}
\]

Let \(\xi_*=k_{0,*}/|k_{0,*}|\), let \(k_c\) be a nearest integer
vector to \(q\xi_*\), and choose the orthonormal frame of C201(1)--(3)
for the reference return shear. Define

\[
 \mathcal B=\left\{k\in\mathbb Z^3:
 |e_j\cdot(k-k_c)|\le\frac{8cq}{L}\ (j=1,2),\quad
 |e_3\cdot(k-k_c)|\le\frac{8cq}{L^2}\right\},\qquad
 \mathcal B_{\mathbb R}=\mathcal B\cup(-\mathcal B).
                                                               \tag{6}
\]

The exact solution, projected into this **actual finite Fourier band**,
satisfies

\[
 \boxed{\|P_{\mathcal B_{\mathbb R}}v(S)\|_\infty
       \ge10^{-157}\,b\,\frac{q^{3/2}}{L^2}.}         \tag{7}
\]

Thus the \(J^{-2}\) concentration order identified in C196 is attained, with

\[
 L<1+\frac{57}{400}\log q+\frac{152}{25}\le\log q.     \tag{8}
\]

No claim of exact band support for the whole solution is made: (7)
concerns its explicitly retained component. The amplitude constant is
small and the sufficient scale is very large; both are part of the
statement, not suppressed asymptotic qualifications.
The same lower bound holds for \(\|v(S)\|_\infty\), and the construction
also gives a true same-energy concentration gain:

\[
 \boxed{
 \frac{\|v(S)\|_\infty/\|v(S)\|_2}
      {\|v(0)\|_\infty/\|v(0)\|_2}
 >\frac{25}{24}q^{3/8}.}                              \tag{8a}
\]

## 2. One concrete profile and one common admissible tube

For completeness a fixed smooth cutoff with numerical derivative bounds
is available. Set \(f_+(s)=e^{-1/s}\) for \(s>0\) and zero otherwise, and

\[
 \theta(s)=\frac{f_+(1-s)}{f_+(1-s)+f_+(s)},\qquad
 \chi(y)=\prod_{j=1}^3\theta(4|y_j|-1).               \tag{9}
\]

The quotient is defined for all real \(s\); it equals one for \(s\le0\)
and zero for \(s\ge1\). The apparent absolute-value corners occur where
the function is constant. Thus \(\chi\) is smooth, lies in \([0,1]\),
equals one on \([-1/4,1/4]^3\), and is supported in

\[
 [-1/2,1/2]^3\subset B(0,1),\qquad
 \|\nabla\chi\|_\infty<10^3,\quad
 \|D^2\chi\|_{F,\infty}<10^6.                        \tag{10}
\]

On \(0\le s\le1\), the denominator is at least \(e^{-2}>1/8\).
The bounds \(|f_+'|\le1\), \(|f_+''|\le10\), applied to the quotient,
give \(|\theta'|\le128\), \(|\theta''|\le5712\).
Scaling each factor by four bounds the coordinate first derivatives by
512 and all second derivatives by 262144. Hence the gradient is below
1024 and its sharper \(512\sqrt3<10^3\) bound proves (10); the Hessian
Frobenius norm is below \(3\cdot262144<10^6\).

Put \(\varepsilon_+=q^{-1/3}\), \(\varepsilon_-=q^{-1/12}\),

\[
 d=\lceil8\pi q^{1/3}\rceil,\quad
 A=F^T[e_1c/L,e_2c/L,e_3c/L^2],\quad
 \mathcal S=\{p\in d\mathbb Z^3:\ p/q\in\xi_*+A[-1,1]^3\},
 \quad M=\#\mathcal S,                                \tag{11}
\]

where \(F=D\Phi_H(y_0)\) at the reference orbit. C201 gives

\[
 \|A\|\le220c/L,\quad
 M\ge\frac{c^3q^2}{27^3L^4},\quad
 \varepsilon_+^{-3/2}\sqrt M
       \ge\frac{c^{3/2}}{162}\frac{q^{3/2}}{L^2}.       \tag{12}
\]

We verify the joint C200 domain, including the entire sine parameter
image in C202(26), rather than only the discrete carriers. Allow initial
positions with

\[
 |y-y_0|\le\rho_0:=10^{-25}L^{-3},\qquad
 \xi=\xi_*+A\sin\vartheta.                            \tag{13}
\]

The covector displacement is at most \(400c/L\). Rescaling to C195's
reference covector multiplies it by \(K=|k_{0,*}|<6\). The separated
base/covector comparison, obtained directly from C195(3.1)--\(3.2\), is

\[
 |\Delta X(t)|\le216L\rho_0,\qquad
 |\Delta k(t)|\le434000000L^3\rho_0
                       +216L\,6|\Delta\xi|
       <5\cdot10^{-17},\qquad 0\le t\le H.            \tag{14}
\]

This uses the actual \(J_1\) cost for an initial covector change, rather
than assigning it the \(J_1J_2\) base-position cost. The C195 phase,
level, vertical-charge, horizontal-charge, and gamma comparison rows
then lie respectively inside

\[
 10^{-7},\quad10^{-11},\quad10^{-7},\quad10^{-5},\quad10^{-5}.
                                                               \tag{15}
\]

For example \( |\Delta\gamma|\le2|\Delta k|+100|\Delta X|<10^{-13}\),


\( |\Delta f|\le4\rho_0\), and the initial horizontal charge costs at
most \(200\rho_0+42|\Delta\xi|\). The latter follows by separately
varying the base and the full covector in \(k\cdot(N\times\nabla f)\).
At the threshold (1), \(\varepsilon_-<\rho_0/4\); both compact
envelopes and neighborhoods needed for differentiation therefore lie in
this same certified tube. Homogeneity of the Kelvin equation makes the
factor \(K\) in (14) immaterial for the physical amplitude or the
normalized phase used in the exact curl.

## 3. Exact solutions and every upper error

On every positive carrier use C200's terminally unit-normalized expanding
symbol \(b_+\), and its parity extension on the negative carrier.
Use the envelope \(\chi_{\varepsilon_+}\) and coefficient

\[
                    z_p=(2M)^{-1/2}\quad(p\in\pm\mathcal S).
\]

Let \(u_+\) be the exact solution (4) with C201's exact-curl datum.
Choose any one positive carrier from \(\mathcal S\), its negative,
coefficients \(1/\sqrt2\), and C200's initially unit-normalized
contracting symbol \(b_-\), with envelope \(\chi_{\varepsilon_-}\).
Let \(u_-\) be its exact solution of the same PDE. These are two finite
superpositions used as initial data for a full infinite-dimensional PDE,
not truncated evolution equations.

Write \(x=\log q\). C202(34)--(35), with the actual C200 initial jets,
the profile bound \(A_{\chi,2}\le10^6\), and positive coefficient norm
at most one, gives for either solution, throughout the stage,

\[
 E_{\rm PDE}\le8\cdot10^{326}x^{46}q^{-37/600}
                 +4\cdot10^{325}x^{14}q^{-45}.          \tag{16}
\]

It includes the actual heat-decay clock and viscous discrepancy. No
second \(q^{-3/8}\) factor is inserted in the expanding coefficients;
that small factor is already inside the selected initial symbol.
C201's exact-curl corrections at the comparison sections obey

\[
 E_{\rm curl}\le10^{256}x^{21}q^{-2/3},                \tag{17}
\]

using \(R+1\le L\le x\) and the worst of the two envelope widths.
The two small principal branches are bounded by \(811\,2999^{-R}\),
and

\[
                     2999^{-R}<q^{-3/8}.              \tag{18}
\]

Here is a convenient explicit common finite-scale enclosure. For

\[
 q\ge10^{100000},\quad x_0=100000\log10\in(200000,300000),
\]

the decreasing functions \(x^a e^{-\beta x}\) in (16)--(18) give

\[
 E_{\rm PDE}<10^{-5500},\qquad
 E_{\rm curl}<10^{-65000},\qquad
 812q^{-3/8}<10^{-37496}.                              \tag{19}
\]

For instance the largest term in (16) is less than
\(8\cdot10^{326+276-6166}<10^{-5563}\).
Consequently the physical Riesz bounds of C201 imply the **actual**
norm inequalities

\[
 \frac1{10}\le\|u_+(S)\|_2,\|u_-(0)\|_2\le29,
 \qquad \|u_+(0)\|_2,\|u_-(S)\|_2<10^{-6}.            \tag{20}
\]

Apply C201's exact quadratic energy balancing to these two solutions.
It gives \(1/300<\lambda<300\) and

\[
 Z=u_++\lambda u_-,\qquad
 N=\|Z(0)\|_2=\|Z(S)\|_2,\qquad1/20<N<30.             \tag{21}
\]

The promised solution is \(v=(b/N)Z\). This proves (5) by an exact
identity involving the actual solutions, without identifying a
coefficient scale with physical energy.

## 4. The retained-band lower bound for that same solution

Let \(x_*=\Phi_H(y_0)\), choosing its lift for the local argument. Use
C203 with

\[
 w_1^{\rm test}=w_2^{\rm test}=2cq/L,\qquad
 w_3^{\rm test}=2cq/L^2,\qquad
 V_{\rm test}=8c^3q^3/L^4,
 \quad r_0=10^{-30}L^{-4}.                            \tag{22}
\]

All widths exceed one. The central endpoint gradients are
\(\eta_p=F^{-T}p\). Equation (11), and the nearest-integer center
error of at most \(\sqrt3/2\), put them in the inner test box. The
negative affine partners lie outside the outer positive test box:
their distance from \(k_c\) exceeds \(q\), whereas its outer box has
Euclidean radius below \(16cq\).

The initial phases are centered at \(y_0\), so their terminal phases at

\[
 x_*\quad\hbox{are }p\cdot(Y_H(x_*)-y_0)=0.
\]

On \(|x-x_*|\le r_0\), first \(|f(x)|\le4r_0<1/10\), so C194's
inverse-flow derivative estimates apply on this entire local ball.
The inverse flow then stays inside (13), because
\(216Lr_0<\rho_0/2\). C200's terminal expanding amplitudes are unit
vectors throughout this domain. Set \(e=b_+(R;y_0,\xi_*)\), also a
real unit vector. Their covector derivative bound and (13) give

\[
 |b_+(R;y_0,\xi_p)-e|
 \le30(4\cdot10^{40}L)(400c/L)<10^{-50},
 \qquad e\cdot b_+(R;y_0,\xi_p)\ge\frac12.            \tag{23}
\]

Thus C203 has \(\kappa=1/2\), \(B_0=B_g=1\). After converting the
C200 initial-position derivative to an Eulerian derivative, its
remaining constants may be taken to be

\[
 B_1\le10^{48}L^6,\qquad
 D\le10^{49}L^6q^{1/3},\qquad C\le10^6qL^2.           \tag{24}
\]

Equations (10), (14), and (24) give these bounds directly; no
pointwise error estimate for the exact PDE is presumed. Substitution
in C203(10) yields

\[
 \mathcal E\le10^{271}x^{12}q^{-2/3}<10^{-66000}<1/8.
                                                               \tag{25}
\]

Let \(A_{\rm peak}=\varepsilon_+^{-3/2}\sqrt M\). C201's count and
(22) imply

\[
                  A_{\rm peak}/\sqrt{V_{\rm test}}>1/500.       \tag{26}
\]

Apply C203 to \(\sqrt2 Z\); this accounts for the signed coefficient
normalization \(1/\sqrt{2M}\). Its distance from the real leading
expanding field includes the curl correction, the two PDE errors, and
the entire contracting solution. By (16)--(19) and \(\lambda<300\),

\[
 E\le2\left[E_{\rm curl}+E_{\rm PDE}
       +300\big(812q^{-3/8}+E_{\rm PDE}\big)\right]
 <10^{-8}<\frac{(1/2)A_{\rm peak}}{256\sqrt{V_{\rm test}}}.
                                                               \tag{27}
\]

C203 therefore proves a retained-band lower bound for this exact
solution, including its stable ballast:

\[
 \begin{split}
 \|P_{\mathcal B_{\mathbb R}}v(S)\|_\infty
 &\ge\frac{b}{30}\frac{A_{\rm peak}}{432}\\
 &\ge\frac{10^{-150}}{2099520}
              b\frac{q^{3/2}}{L^2}
 >10^{-157}b\frac{q^{3/2}}{L^2}.                       \tag{28}
 \end{split}
\]

This is the complete inequality chain from the same physical solution
to (7). The retained-volume tax has not been removed or hidden.

## 5. Entrance and full-trajectory physical budgets

The same solution obeys the additional explicit estimates

\[
 \|v(0)\|_\infty\le20000bq^{9/8},\qquad
 \sup_{0\le s\le S}\|v(s)\|_2<10^{16}b,                 \tag{29}
\]

\[
 \int_0^S\|v(s)\|_2^2\,ds<4\cdot10^{33}b^2,\qquad
 \mu\int_0^S\|\nabla v(s)\|_2^2\,ds<2\cdot10^{34}b^2.  \tag{30}
\]

These bounds charge the linear perturbation's full trajectory; they do
not assert that its norm never exceeds its endpoint norm.

For the entrance bound, the same mesh argument as C201 puts every
occupied grid cell inside \(\xi_*+A[-3/2,3/2]^3\). Since
\(d>24q^{1/3}\), this yields

\[
 M(d/q)^3\le27\det A,\quad
 M\le\frac{c^3q^2}{512L^4},\quad
 \sqrt{2M}\le\frac{c^{3/2}q}{16L^2}.                  \tag{31}
\]

At time zero \(k=\xi\), so \(|c_p|\le2|b_p|\) and
\(|\operatorname{curl}c_p|\le5\|D_yb_p\|\). With
\(\sigma=2999^{-R}\), C200 and the explicit envelope give

\[
 \|u_+(0)\|_\infty
 \le\sqrt{2M}\varepsilon_+^{-3/2}\sigma
 \left[4+\frac{8000}{q\varepsilon_+}
             +\frac{5\cdot10^{45}(R+1)L^4}{q}\right]
 \le q^{9/8},                                        \tag{32}
\]

\[
 \|u_-(0)\|_\infty
 \le\sqrt2\varepsilon_-^{-3/2}
 \left[1+\frac{2000}{q\varepsilon_-}
             +\frac{5\cdot10^{45}(R+1)L^4}{q}\right]
 \le3q^{1/8}.                                        \tag{33}
\]

Each correction bracket beyond its leading constant is at most one at
the threshold (1); to obtain the last bound in (33), its correction is
in fact below \(1/10\). Thus (21) gives
\(\|v(0)\|_\infty\le20b(q^{9/8}+900q^{1/8})
<20000bq^{9/8}\).

For all action-section times \(s_j\) defined by \(\tau(s_j)=jT\),
C200--C202 and C201's nonunit-symbol synthesis bound give

\[
 \|u_+(s_j)\|_2\le812\rho^{R-j}+E_{\rm PDE},\qquad
 \|u_-(s_j)\|_2\le812\rho^j+E_{\rm PDE},\quad
 \rho=1/2999.                                        \tag{34}
\]

Here the same small upper bound \(E_{\rm PDE}\) holds for every section,
and the exact-curl correction carries the corresponding factor
\(\rho^{R-j}\) or \(\rho^j\). In particular each norm is below 812:
the sharper principal bound \(811\rho^k\) and corrections smaller than
one prove that assertion even when \(k=0\).
The actual PDE energy inequality between sections is
\(\|v(s)\|_2\le e^{6(\tau(s)-jT)}\|v(s_j)\|_2\).
Since \(e^{6T}<3^{19}<2\cdot10^9\), it proves the supremum bound in
(29) with the exact majorant

\[
 2\cdot10^9\cdot20\cdot301\cdot812\,b
 =9776480000000000\,b<10^{16}b.                        \tag{35}
\]

Keeping the geometric factors gives the stronger action estimate.
For \(jT\le t\le(j+1)T\), (21) and (34) imply

\[
 \|v(s(t))\|_2\le4\cdot10^{10}b
       [812\rho^{R-j}+243600\rho^j+301E_{\rm PDE}].
\]

Both sums of squared geometric factors are below two, and
\(R E_{\rm PDE}^2<1\): from (16),
\[
 R E_{\rm PDE}^2
 \le81\cdot10^{652}x^{93}q^{-37/300}<1.
\]
Using \(T<4\) and \((a+b+c)^2\le3(a^2+b^2+c^2)\) gives

\[
 \int_0^H\|v(s(t))\|_2^2\,dt
 \le(4\cdot10^{10})^2\,12
       [2\cdot812^2+2\cdot243600^2+301^2]b^2
 <3\cdot10^{33}b^2.                                  \tag{36}
\]

The factor \(ds/dt\le100/99\) proves the first bound (30).
Finally the exact energy identity and (5) give

\[
 \mu\int_0^S\|\nabla v\|_2^2\,ds
 =-\int_0^S e^{-2\mu s}\langle DU\,v,v\rangle\,ds
 \le6\int_0^H\|v(s(t))\|_2^2\,dt<2\cdot10^{34}b^2,
\]

which proves the other bound (30). No nonlinear cancellation is being
credited by this linearized energy calculation.

### 5.1 A concentration quotient without a logarithmic loss

Keep the unrounded endpoint inequality from (28),
\(\|v(S)\|_\infty\ge(b/N)A_{\rm peak}/432\), using C203's
unprojected version. The extra return in (1) gives
\[
 \sigma=2999^{-R}\le\frac{q^{-3/8}}{2999},
\]
since \(2999>e^8\). At the threshold (1), the two correction terms
in (32) sum to less than \(1/100\). Thus
\[
 \|u_+(0)\|_\infty
 \le\frac{401}{100}\sqrt2\,\sigma A_{\rm peak}
 <\frac{401}{70\cdot2999}A_{\rm peak}q^{-3/8}.          \tag{37}
\]
The stable entrance satisfies
\[
 \lambda\|u_-(0)\|_\infty\le900q^{1/8}
 \le\frac1{10000}A_{\rm peak}q^{-3/8},                 \tag{38}
\]
because (12) reduces this last inequality to
\(q\ge1458\cdot10^{156}L^2\), which holds at (1) and thereafter.
The exact rational comparison
\[
 \frac{401}{70\cdot2999}+\frac1{10000}<\frac1{450}
\]
therefore proves
\[
 \|v(0)\|_\infty<\frac bN\frac{A_{\rm peak}q^{-3/8}}{450},
 \qquad
 \frac{\|v(S)\|_\infty}{\|v(0)\|_\infty}
 >\frac{450}{432}q^{3/8}=\frac{25}{24}q^{3/8}.         \tag{39}
\]
Equation (5) makes this exactly the concentration quotient (8a).
The identical ratio lower bound holds with the retained endpoint
norm in the numerator. This calculation uses the shared peak scale;
it does not divide unrelated upper and lower estimates and cancel
their logarithmic factors without justification.

## 6. What this advances, and what remains

C204 completes the finite-frequency periodic viscous version of the
C193 filter, with physical endpoint normalization, and proves the
C196 log-taxed absolute retained concentration. The filter transfers
entrance energy carried chiefly by a contracting polarization to
endpoint energy carried by a concentrating expanding polarization.
It is not scalar multiplication of a fixed-energy seed.

The equation in (4) is linearized. Adding \(v\) to its pump leaves the
quadratic defect \(P(v\cdot\nabla v)\); this note neither cancels nor
certifies that defect within the nonlinear stage budget. It also does
not prove a nonlinear return of the retained band and wake, the C125
relative residual, or a repeating UVSR profile. Equations (29)--(30)
bound the linear perturbation's trajectory; they do not bound the
quadratic defect or pump depletion. A subsequent nonlinear use must
control those terms, the pressure and wake, and retain the standing
terminal singular-center obligation.

The constructed band has C196's logarithmic width orders, centered at
the C159 covector. Its identification with C180's prescribed carrier
center, frame, splitter output, and exit chart is not proved here.
Neither is C182's stronger entrance hypothesis \(\|v(0)\|_\infty\le bq\);
the proved entrance bound is (29).

The argument is constructive at the level of explicitly defined finite
carrier sums and full PDE solutions. The smooth compact envelopes have
Fourier tails; “finite-frequency” here means a finite specified carrier
scale, not finite Fourier support of the entire initial datum.
Its astronomical sufficient scale
has not been numerically instantiated. No numerical profile residual,
Navier--Stokes singularity, or Clay result is claimed.
