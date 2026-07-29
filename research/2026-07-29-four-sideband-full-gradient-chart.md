# Four common-carrier sidebands give an exact full-gradient chart

Date: 2026-07-29

## Result and claim boundary

The preferred arbitrary-polarization one-carrier construction is not
restricted to a rank-five symmetric endpoint.

With one fixed parent

\[
 k=Ke_3,\qquad a=e_1,
\]

the four low directions

\[
\begin{aligned}
q_1&=(-45,-36,20),&
q_2&=(-4,-5,9),\\
q_3&=(1,1,1),&
q_4&=(1,2,3)
\end{aligned}
\]

and two independently controlled partner polarizations for each direction
give an exact rank-eight chart of \(\mathfrak{sl}(3)\) for every positive
integer \(K\).  All waves retain the same fast direction, so the
all-generation cancellation of the nominal \(K\) factor remains exact.

In particular, the chart synthesizes

\[
 S_{5/4}
 =\operatorname{diag}\left(-\frac54,\frac94,-1\right).
\]

The affine Euler flow \(U=S_{5/4}x\) multiplies a carrier parallel to
\(e_3\) by \(e^t\) and a Kelvin polarization parallel to \(e_1\) by
\(e^{5t/4}\).  It is therefore the exact local logarithmic deformation
required by the \(\gamma=5/4\) carrier/amplitude ledger.

This removes the former full-gradient algebraic defect without introducing
a second fast direction.  It does **not** localize the affine flow, create
an endpoint return, cancel pressure multipoles, or prove a Navier--Stokes
singularity.

All identities are checked over the rationals in
[`checks/four_sideband_full_gradient.py`](../checks/four_sideband_full_gradient.py).

---

## 1. Exact child planes

For a low direction

\[
 q=(\alpha,\beta,\delta),
\qquad r=q-Ke_3,
\]

and a partner amplitude \(b\perp r\), the exact matched child is

\[
 {\cal L}_{q,K}b
 =
 P_q\left[\alpha b+(b\cdot q)e_1\right].
\tag{1.1}
\]

Use the domain basis

\[
 d_1=(-\beta,\alpha,0),
\qquad
 d_2=(K-\delta,0,\alpha)
\]

of \(r^\perp\), and the child basis

\[
 c_1=(-\beta,\alpha,0),
\qquad
 c_2=(-\delta,0,\alpha)
\tag{1.2}
\]

of \(q^\perp\).  The exact block is

\[
 [{\cal L}_{q,K}]
 =
 \begin{pmatrix}
 \alpha&-2\alpha\beta K/|q|^2\\
 0&\alpha(|q|^2-2K\delta)/|q|^2
 \end{pmatrix},
\tag{1.3}
\]

so

\[
 \det[{\cal L}_{q,K}]
 =
 \frac{\alpha^2(|q|^2-2K\delta)}{|q|^2}.
\tag{1.4}
\]

For the four directions above the exceptional carrier values are

\[
 \frac{3721}{40},\qquad
 \frac{61}{9},\qquad
 \frac32,\qquad
 \frac73.
\tag{1.5}
\]

None is an integer.  Every child plane is therefore exactly controllable
at every positive integer carrier.

---

## 2. Exact rank eight

For each \(q_\alpha\), use the two child columns in (1.2).  A
divergence-free sine mode

\[
 c\sin(q\cdot y)
\]

has gradient \(c\otimes q\) at \(y=0\).  Express trace-free matrices in
the eight coordinates

\[
 (M_{11},M_{12},M_{13},M_{21},
   M_{22},M_{23},M_{31},M_{32});
\tag{2.1}
\]

\(M_{33}=-M_{11}-M_{22}\).

Put the eight columns

\[
 c_{\alpha,s}\otimes q_\alpha,
\qquad
\alpha=1,\ldots,4,\quad s=1,2,
\tag{2.2}
\]

into an \(8\times8\) matrix in the order displayed above.  Exact Gaussian
elimination gives

\[
 \boxed{\det{\cal G}=15\,451\,090\,200.}
\tag{2.3}
\]

Hence

\[
 \operatorname{span}
 \{c_{\alpha,s}\otimes q_\alpha\}
 =\mathfrak{sl}(3).
\tag{2.4}
\]

Composing (2.3) with the four exact child blocks gives

\[
 \boxed{
 \det{\cal G}_K
 =
 15\,451\,090\,200
 \prod_{\alpha=1}^4
 \frac{
 (q_{\alpha,1})^2\,
 \big(|q_\alpha|^2-2Kq_{\alpha,3}\big)}
 {|q_\alpha|^2}.
 }
\tag{2.5}
\]

By (1.5), (2.5) never vanishes at a positive integer \(K\).

---

## 3. Exact Kelvin target and a periodic entry-strain germ

For \(\gamma=5/4\), the coefficients in the ordered bases

\[
 (c_{1,1},c_{1,2},c_{2,1},c_{2,2},
   c_{3,1},c_{3,2},c_{4,1},c_{4,2})
\]

are

\[
\begin{aligned}
x=\bigg(&-\frac{169}{73367},
-\frac{17027}{9537710},
\frac{9753}{62440},
\frac{193829}{2435160},\\
&\frac{9693}{6580},
\frac{144689}{42770},
\frac{145119}{209620},
-\frac{163724}{157215}\bigg).
\end{aligned}
\tag{3.1}
\]

Let

\[
 C_\alpha
 =x_{2\alpha-1}c_{\alpha,1}
  +x_{2\alpha}c_{\alpha,2}.
\tag{3.2}
\]

Then, exactly,

\[
 \sum_{\alpha=1}^4C_\alpha\otimes q_\alpha
 =
 \operatorname{diag}\left(-\frac54,\frac94,-1\right).
\tag{3.3}
\]

Thus

\[
 U_{\rm germ}(y)
 =\sum_{\alpha=1}^4
 C_\alpha\sin(q_\alpha\cdot y)
\tag{3.4}
\]

is periodic, real analytic, divergence free, and satisfies

\[
 U_{\rm germ}(0)=0,\qquad
 \nabla U_{\rm germ}(0)=S_{5/4}.
\tag{3.5}
\]

The exact affine field \(U_{\rm aff}=S_{5/4}x\) solves Euler with

\[
 p_{\rm aff}(x)
 =-\frac12x\cdot S_{5/4}^2x.
\tag{3.6}
\]

For a Kelvin wave with \(k(0)\parallel e_3\) and \(a(0)\parallel e_1\),

\[
\begin{aligned}
 k(t)&=e^{-S_{5/4}^Tt}k(0)=e^tk(0),\\
 a'(t)&=-S_{5/4}a(t)=\frac54a(t),
\end{aligned}
\tag{3.7}
\]

because the Kelvin pressure correction vanishes for these eigendirections.
Therefore

\[
 |k(t)|=e^t|k(0)|,\qquad
 |a(t)|=e^{5t/4}|a(0)|.
\tag{3.8}
\]

This is exactly the local carrier/amplitude scaling in the
\(\gamma=5/4\) cascade.

The periodic germ (3.4) agrees with the affine field only to first order
at one point.  Promoting it to a scale-uniform analytic core, preserving
the Kelvin deformation across the packet, and returning the entire field
remain nonlinear analytic tasks.

---

## 4. The common-carrier null structure survives

Every generated charged wave still has physical wave vector

\[
 p=Kh\,e_3+\xi,
\qquad h\ne0,
\tag{4.1}
\]

where \(\xi\) is a slow integer combination of the four \(q_\alpha\)'s.
For transverse amplitudes

\[
 A\cdot(Kh\,e_3+\xi)=0,
\qquad
 B\cdot(Kg\,e_3+\eta)=0,
\]

one has

\[
\begin{aligned}
 A\cdot(Kg\,e_3+\eta)
 &=A\cdot\left(\eta-\frac gh\xi\right),\\
 B\cdot(Kh\,e_3+\xi)
 &=B\cdot\left(\xi-\frac hg\eta\right).
\end{aligned}
\tag{4.2}
\]

There is no \(K\) on the right.  Adding \(q_4\) enlarges the slow lattice
but does not add a fast direction, so it cannot recreate the former
two-colour \(O(K)\) chain.

There is one bookkeeping qualification.  Four integer vectors in
\(\mathbb Z^3\) cannot be freely labelled.  Here

\[
 -15q_1+47q_2-669q_3+182q_4=0.
\tag{4.3}
\]

A hierarchy should therefore be indexed by the physical slow vector
\(\xi\), or by the quotient of \(\mathbb Z^4\) by (4.3), rather than by
four free slow charges.  This first displayed relation has total
coefficient size \(913\), so it does not affect the low-order chart.  At
higher order it merges interaction paths but does not alter (4.2):
all representatives have the same \(\xi\), while \(h\) remains the only
coefficient of \(K\).

The remaining gates are therefore analytic and global:

1. the sideband chart differentiates low acceleration, so a from-rest
   material deformation has the \(T^2/2\) time factor;
2. the affine/periodic germ must recur as an incoming low strain;
3. active one-carrier covariance is not pressure-quadrupole dark;
4. terminal charged modes must enter the next packet or a genuine wake;
5. viscosity and both cutoff seams require the all-order endpoint
   hierarchy.
