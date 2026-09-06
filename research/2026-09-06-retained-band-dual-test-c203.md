# C203: an explicit retained-band test for a nonlinear-phase packet

**Date:** 2026-09-06

**Status:** explicit Fourier projection lower-bound lemma. This is a
conversion estimate; its application to the A2 solution must verify every
packet and error hypothesis below.

**Checker:** [retained_band_dual_test_c203.py](../checks/retained_band_dual_test_c203.py).

Use unnormalized Lebesgue measure on the torus of side \(2\pi\). All
Euclidean local coordinates below are centered at the proposed endpoint.
Let \(e_1,e_2,e_3\) be an orthonormal frame, \(w_j\ge1\),

\[
 V=w_1w_2w_3,\qquad w_*=\min_jw_j,\qquad k_c\in\mathbb Z^3.
\]

The retained positive band is the set of integer modes

\[
 \mathcal B=\{k:\ |e_j\cdot(k-k_c)|\le4w_j,
                         \ j=1,2,3\}.                 \tag{1}
\]

No rationality of the frame is required. Let \(P_\mathcal B\) be the
orthogonal Fourier projection onto exactly this finite set.

## 1. An actual finite-band dual function

Put \(s(z)=(\sin z/z)^4\), with its removable value at zero, and

\[
 \psi(x)=\sqrt V\prod_{j=1}^3s(w_j e_j\cdot x),\qquad
 g(x)=e^{ik_c\cdot x}\sum_{m\in\mathbb Z^3}\psi(x+2\pi m).
                                                               \tag{2}
\]

The sum defines a nonnegative envelope in \(L^1\cap L^2\); equivalently
it is its finite Fourier series. The transform convention is
\(\widehat h(\eta)=\int_{\mathbb R^3}h(x)e^{-i\eta\cdot x}\,dx\).
The identity

\[
 \widehat s(\theta)=\pi\left(\frac23-\frac{\theta^2}4
                       +\frac{|\theta|^3}{16}\right)
                       \quad (|\theta|\le2)             \tag{3}
\]

follows by convolving four copies of the transform of \(\sin z/z\).
Its full support is \([-4,4]\), it is nonnegative, and
\(\widehat s(\theta)>1\) on \(|\theta|\le1\). Thus

\[
 \widehat\psi(\eta)\ge V^{-1/2}
       \quad\hbox{if }|e_j\cdot\eta|\le w_j,
 \qquad \operatorname{supp}\widehat g\subset\mathcal B. \tag{4}
\]

In particular, this is a genuine retained-band test and is not a
compactly supported function silently treated as bandlimited.

The elementary pointwise bound \(s(z)\le\min(1,|z|^{-4})\) gives

\[
 \int s<3,\qquad \int |z|s\le2,\qquad\int z^2s<3.
\]

Consequently

\[
 \|g\|_1<\frac{27}{\sqrt V},\quad
 \int |x|\psi(x)\,dx\le\frac{54}{w_*\sqrt V},\quad
 \int |x|^2\psi(x)\,dx\le\frac{81}{w_*^2\sqrt V}.
                                                               \tag{5}
\]

Also

\[
                         \|g\|_2<64.                   \tag{6}
\]

Indeed the lattice-point cube argument gives
\(\#\mathcal B\le8\prod_j(4w_j+\sqrt3/2)\le1000V\).
Each Fourier coefficient of \(g\) has magnitude at most
\(27/\((2\pi)^3\sqrt V\)\). Parseval gives
\(\|g\|_2^2\le1000\,27^2/\(2\pi\)^3<3375<64^2\), using
\(\pi>3\). These estimates also justify the periodization in (2).

## 2. Hypotheses on one approximate endpoint

Let \(M\ge1\), \(0<\varepsilon\le1\), and suppose a real leading
endpoint field has the form

\[
 v_{\rm lead}(x)=2\operatorname{Re}\left[
       M^{-1/2}\sum_{p=1}^M e^{i\phi_p(x)} b_p(x)\right]. \tag{7}
\]

It is a well-defined periodic field, even though individual real phases
need only be defined on a local lift. The vectors \(b_p\) are real. Fix a
real unit vector \(e\), a radius \(0<r_0\le1\), and positive constants
\(\kappa,B_0,D,C\). On the local ball \(|x|\le r_0\), suppose

\[
 \begin{gathered}
 \phi_p(0)=0,\quad \eta_p=\nabla\phi_p(0),\quad
 |\phi_p(x)-\eta_p\cdot x|\le C|x|^2,\\
 b_p(0)=\varepsilon^{-3/2}a_p,\quad
 e\cdot a_p\ge\kappa,\quad |a_p|\le B_0,\\
 |b_p(x)-\varepsilon^{-3/2}a_p|
                 \le\varepsilon^{-3/2}D|x|.             \tag{8}
 \end{gathered}
\]

Globally, assume explicitly
\(\|v_{\rm lead}\|_\infty\le2\varepsilon^{-3/2}\sqrt M B_g\),
where \(B_g\) is stated explicitly. The affine comparison built from
the vectors \(a_p\) has norm at most
\(2\varepsilon^{-3/2}\sqrt M B_0\). Individual local phases need not
extend globally; these two distinct sum bounds suffice.
Assume for every carrier

\[
 |e_j\cdot(\eta_p-k_c)|\le w_j,
 \qquad -\eta_p-k_c\notin\prod_j[-4w_j,4w_j]_{e_j}.       \tag{9}
\]

The second condition prevents the negative-frequency affine reality
partner from entering the test band; no orthogonality of nonlinear
transported bands is assumed.

Define the explicit error number

\[
 \mathcal E=\frac{54D}{w_*}
   +\frac{81}{w_*^2}
       \left[B_0C+\frac{B_g+B_0}{r_0^2}\right].          \tag{10}
\]

If a real exact endpoint \(u\) obeys

\[
 \|u-v_{\rm lead}\|_2\le E,\qquad
 \mathcal E\le\frac\kappa4,\qquad
 E\le\frac{\kappa\varepsilon^{-3/2}\sqrt M}{256\sqrt V}, \tag{11}
\]

then the physically retained endpoint satisfies

\[
 \boxed{\|P_\mathcal B u\|_\infty
       \ge\frac\kappa{108}\varepsilon^{-3/2}\sqrt M.}   \tag{12}
\]

For the real reality-complete projection \(P_{\mathcal B\cup(-\mathcal B)}\),
the same lower bound holds: pair it with the same \(g\), since this larger
projection also acts identically on \(g\).
The identical lower bound holds for \(\|u\|_\infty\) itself, by applying
Hölder directly to \(\langle u,g\rangle\). This is not an assertion that
an orthogonal Fourier projection is an \(L^\infty\) contraction.

## 3. Proof and the costs that enter an A2 application

Let \(A=\varepsilon^{-3/2}\sqrt M\). Unfold the periodized test (2) to
\(\mathbb R^3\). The affine positive-frequency comparison contributes
at least \(\kappa A/\sqrt V\), by (4) and (8). Its affine negative
partner contributes zero, by (9). On the local ball use
\(|e^{ia}-e^{ib}|\le|a-b|\), (8), and the first two moments (5).
Outside the ball use the global sum bound, the corresponding affine
bound, and \(\mathbf1_{|x|>r_0}\le|x|^2/r_0^2\). The two reality
partners together cost at most \(2A\mathcal E/\sqrt V\). Therefore

\[
 \operatorname{Re}\langle e\cdot v_{\rm lead},g\rangle
       \ge\frac A{\sqrt V}(\kappa-2\mathcal E).          \tag{13}
\]

Equation (6) bounds the exact-solution replacement error by \(64E\).
As \(P_\mathcal B g=g\), Hölder and (5) now give

\[
 \frac{27}{\sqrt V}\|P_\mathcal B u\|_\infty
 \ge\frac A{\sqrt V}(\kappa-2\mathcal E)-64E
 \ge\frac{\kappa A}{4\sqrt V}.
\]

This proves (12). It uses \(L^2\) error only, yet certifies an absolute
point norm of an actual finite Fourier projection.

For a C194-type phase
\(\phi_p(x)=q\xi_p\cdot(Y_t(x)-y_0)\),
\(|\xi_p|\le R\), and a transported compact envelope with
\(\chi(0)=1\), the concrete substitutions are

\[
 C=\frac{qRJ_2(t)}2,\qquad
 D=B_1+\|\nabla\chi\|_\infty B_0J_1(t)\varepsilon^{-1},
 \qquad B_g=\|\chi\|_\infty B_0,                       \tag{14}
\]

provided \(B_0,B_1\) genuinely bound the chosen terminal amplitude and
its spatial derivative on the stated sets. The error \(E\) must include
the exact-curl correction as well as the full PDE approximation and
viscous errors. They cannot be dropped because (7) is only the leading
field. If \(w_*\) is a fixed constant times \(q/(1+t)^2\), (10) charges
the actual curvature and envelope losses directly. No unspecified
Fourier-tail theorem is being used.

This lemma alone constructs neither the growing family nor its physical
energy normalization. Those hypotheses must be established for the same
solution before (12) receives endpoint credit.
