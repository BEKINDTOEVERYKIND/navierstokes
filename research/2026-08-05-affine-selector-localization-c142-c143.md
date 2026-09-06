# C142--C143: exact affine child/wake selector and its localization obstruction

**Date:** 2026-08-05
**Status:** exact global-affine/Kelvin algebra and exact scale identities;
finite-energy unforced localization and the one-cell stage remain open
**Checker:**
[checks/affine_selector_localization_c142_c143.py](../checks/affine_selector_localization_c142_c143.py)

## 0. Claim boundary

This note tests an active-focus mechanism inside the existing C114--C141
one-cell geometry. It introduces no new cascade architecture. The
question is whether a common incompressible affine deformation can amplify
the intended \(A_2\) child by

\[
 h=q^{3/2}
 \tag{0.1}
\]

while leaving the named C117 wakes bounded and then be localized inside one
parent cell.

For global Kelvin modes the answer is favorable in velocity amplitude: an
explicit affine map amplifies the child by \(\asymp h\) and both named wakes
by only \(O(1)\). It also permits viscosity to damp the wake carriers while
leaving the child carrier essentially undamped.

The localization audit gives one exact obstruction and one sharp failure
of the generic cutoff estimate. A genuinely
material three-dimensional child envelope is stretched into a covector of
size \(hqK\), and its viscous action is

\[
 \asymp q^2\log n,
 \tag{0.2}
\]

not \(o(1)\). A stationary compact affine core avoids that material heat
cost, but the direct product estimate for its parent--collar exposure has
integrated, backward-focused scale \(O(\log h)\), independent of how
slowly the strain is run. This estimate does not close the focus step. It
is not a lower bound for every cutoff: algebraic or Lagrangian cancellation
could make the actual interaction smaller.

The surviving affine possibility is narrower: a co-moving
Lagrangian-frame completion must cancel the leading parent transport,
avoid the stretched three-dimensional envelope, and make the integrated
wake-to-active collar block \(O(b)\), rather than its generic
\(O(\log h)\) size. No such completion is proved here.

## 1. Existing \(A_2\) geometry

Use

\[
 \begin{aligned}
 N&=(1,1,1),& k_c&=(1,0,-1),\\
 e_1&=(2,-1,-1),&e_2&=(1,1,-2).
 \end{aligned}
 \tag{1.1}
\]

These are the C116 child and the two C117 wakes. Put

\[
 u={k_c\over\sqrt2},\qquad
 v={e_2-e_1\over\sqrt6},\qquad
 w={N\over\sqrt3}.
 \tag{1.2}
\]

Then \((u,v,w)\) is an oriented orthonormal basis and

\[
 k_c=\sqrt2u,\qquad
 e_{1,2}={3\sqrt2\over2}u\mp{\sqrt6\over2}v,
 \qquad N=\sqrt3w.
 \tag{1.3}
\]

In particular,

\[
 e_1+e_2=3k_c.
 \tag{1.4}
\]

## 2. C142: the exact global affine selector

Let \(h\ge1\), and define a symmetric volume-preserving covector map \(G\)
by

\[
 Gu=u,\qquad
 Gv=sv+dw,\qquad Gw=dv+sw,
 \tag{2.1}
\]

where

\[
 s={h+h^{-1}\over2},\qquad
 d={h-h^{-1}\over2}.
 \tag{2.2}
\]

The eigenvalues of \(G\) are \(1,h,h^{-1}\). Let \(F=G^{-T}=G^{-1}\).
Choose a smooth scalar \(\rho(t)\) and replace \(h\) in (2.1) by
\(e^{\rho(t)}\). Then

\[
 S(t)=F'(t)F(t)^{-1}
 \tag{2.3}
\]

is symmetric and trace free. Indeed, if \(B\) is the fixed symmetric map
with \(Bu=0\), \(Bv=w\), and \(Bw=v\), then
\(G(t)=e^{\rho(t)B}\), \(F(t)=e^{-\rho(t)B}\), and
\(S(t)=-\rho'(t)B\). The whole-space field

\[
 U(t,x)=S(t)x,\qquad
 p(t,x)=-{1\over2}x\mathbin\cdot(S'(t)+S(t)^2)x
 \tag{2.4}
\]

is an exact unforced Euler and Navier--Stokes solution because
\(\nabla\cdot U=0\), \(\Delta U=0\), and \(S'+S^2\) is symmetric. It is
nonperiodic and has infinite energy; (2.4) is an algebraic model, not the
desired stage.

For an inviscid Kelvin mode with initial wavevector \(k\) and polarization
\(N\perp k\), the endpoint wavevector is \(Gk\). The exact velocity
amplitude is

\[
 {F(k\times N)\times Gk\over |Gk|^2}
   =P_{Gk}GN.
 \tag{2.5}
\]

The equality follows from \(\det F=1\):

\[
 F(k\times N)=Gk\times GN.
 \tag{2.6}
\]

Define

\[
 A(k)=|P_{Gk}GN|,\qquad
 D={h^2+h^{-2}\over2}=s^2+d^2.
 \tag{2.7}
\]

Direct substitution of (1.3) into (2.1) gives

\[
 \begin{aligned}
 |Gk_c|^2&=2,& A(k_c)^2&=3D,\\
 |Ge_i|^2&={3\over2}(D+3),&
 A(e_i)^2&={3(3D+1)\over D+3},\qquad i=1,2.
 \end{aligned}
 \tag{2.8}
\]

Relative to \(|N|^2=3\), the child velocity gain squared is \(D\), whereas
each wake gain squared is

\[
 {3D+1\over D+3}\in[1,3).
 \tag{2.9}
\]

Thus the child gain is \(h/\sqrt2+o(h)\), while both wake velocity gains
remain bounded. This is genuine selector algebra, not a small-angle
approximation.

### 2.1 The derivative and conditioning costs

Put

\[
 X_c=|Gk_c|A(k_c),\qquad X_i=|Ge_i|A(e_i).
 \tag{2.10}
\]

Equation (2.8) yields the exact identity

\[
 X_i^2={9\over4}X_c^2+{9\over2}.
 \tag{2.11}
\]

More generally, for every invertible common deformation \(G\), not just
(2.1),

\[
 \boxed{
 |Gk_c|A(k_c)
 \le { |Ge_1|A(e_1)+|Ge_2|A(e_2)\over3}.}
 \tag{2.12}
\]

Indeed, \(|Gk|A(k)=|Gk\times GN|\), so (2.12) is the triangle inequality
applied to (1.4). A common affine map therefore cannot amplify the
child in a one-derivative norm while suppressing both wake derivatives.
This is not by itself an LCE no-go: the physical wake coefficient begins
at \(b^2\), while the child begins at \(b\), and exact Kelvin labels do not
interact in the linear global model.

The full transverse child map is also badly conditioned. Since
\(k_c^\perp=\operatorname{span}\{v,w\}\), its singular values are exactly

\[
 h,\qquad h^{-1}.
 \tag{2.13}
\]

Hence its inverse norm is \(h\) and its condition number is \(h^2\). At
\(h=q^{3/2}=n^{12}\), a chart controlling the whole transverse plane would
pay \(n^{12}\), not the C131 allowance \(n^2\). This does not rule out a
single prescribed scalar ray; it rules out inferring a uniformly
conditioned full chart from a determinant calculation.

## 3. Viscous selection of the three carriers

Take the constant-strain path \(\rho(t)=\lambda t\) on \(0\le t\le T\),
with \(h=e^{\lambda T}\). Viscosity multiplies each Kelvin amplitude by

\[
 \exp\left(-\nu\int_0^T|G(t)k|^2dt\right).
 \tag{3.1}
\]

For the child and either named wake, the exact frequency actions are

\[
 I_c=2T,
 \qquad
 I_w={9T\over2}+{3(h^2-h^{-2})\over8\lambda}.
 \tag{3.2}
\]

Set \(q=n^8\), \(h=q^{3/2}=n^{12}\), and choose, for fixed \(L>0\),

\[
 \lambda={3\nu K^2(h^2-h^{-2})\over8L\log n},
 \qquad T={\log h\over\lambda}.
 \tag{3.3}
\]

Then the second term in \(\nu K^2I_w\) is exactly \(L\log n\), so the
wake amplitude multiplier is at most \(n^{-L}\). The child exponent is

\[
 2\nu K^2T
 ={16L\log n\log h\over3(h^2-h^{-2})}=o(1).
 \tag{3.4}
\]

At the factorial frequency \(K=(j!)^8\),

\[
 T_j={32L(\log n)^2
 \over\nu(j!)^{16}(n^{24}-n^{-24})},
 \tag{3.5}
\]

so \(\sum_jT_j<\infty\). But the strain action is
\(\lambda_jT_j=\log h=12\log n\), not a summable small perturbation.

## 4. C143: material-envelope and collar obstructions

### 4.1 A transported three-dimensional envelope is heat-killed

The carrier calculation (3.2)--(3.4) does not control a localized
three-dimensional envelope. A material envelope covector component
initially of child scale \(qK\) in the expanding covector direction becomes
\(qKe^{\lambda t}\). Its exact viscous action is

\[
 \begin{aligned}
 \nu q^2K^2\int_0^T e^{2\lambda t}dt
 &= {\nu q^2K^2(h^2-1)\over2\lambda}\\
 &= {4Lq^2\log n\,(h^2-1)
       \over3(h^2-h^{-2})}
 \sim {4L\over3}q^2\log n.
 \end{aligned}
 \tag{4.1}
\]

It diverges. Launching that covector at only \(qK/h\) avoids (4.1), but
its initial physical width is

\[
 {h\over qK}={q^{1/2}\over K},
 \tag{4.2}
\]

which is larger than the parent scale \(K^{-1}\). Thus the tuned
viscosity selector is not compatible with an ordinary material
three-dimensional child envelope inside one parent cell.

This is a statement about a common materially transported envelope. It
does not rule out an Eulerian source, an unfolded anisotropic
construction, or a specially slaved collar.

### 4.2 Exact compact affine core and its scale ledger

The alternative is to localize the affine velocity itself. Let

\[
 \ell=K^{-1},\qquad r={\ell\over q},\qquad
 \mu={\nu K\over a},\qquad
 \Lambda={\lambda\over aK}
 ={3\mu(h^2-h^{-2})\over8L\log n}.
 \tag{4.3}
\]

For a constant symmetric trace-free \(S\), put

\[
 {\cal A}(x)={Sx\times x\over3},\qquad
 U_c=\nabla\times\bigl(\chi(x/r){\cal A}(x)\bigr),
 \tag{4.4}
\]

where \(\chi=1\) on the inner core and is supported in the doubled core.
The identity

\[
 \nabla\times{\cal A}=Sx
 \tag{4.5}
\]

shows that \(U_c\) is exactly divergence free, compactly supported, and
equals the affine field on the inner core.

With \(p_0=-x\cdot S^2x/2\) and \(p_c=\chi p_0\), standard product-rule
counting gives, for \(1\le p\le\infty\) and \(m\ge0\),

\[
 \begin{aligned}
 \|\nabla^mU_c\|_p&\lesssim
       \lambda r^{1-m+3/p},\\
 \|\nabla^mR_E\|_p&\lesssim
       \lambda^2r^{1-m+3/p},\\
 \|\nabla^mR_\nu\|_p&\lesssim
       \nu\lambda r^{-1-m+3/p},
 \end{aligned}
 \tag{4.6}
\]

where \(R_E=(U_c\cdot\nabla)U_c+\nabla p_c\) and
\(R_\nu=-\nu\Delta U_c\). Both residuals vanish on the inner core and
outside the doubled core, so they are supported in the collar. The
implicit constants depend on finitely many fixed cutoff seminorms, not on
\(q,K,\lambda\), or \(a\).

Normalize time by \(\tau=(aK)^{-1}\) and velocities by the parent
\(L^p\) scale \(a\ell^{3/p}\). Over \(T=\log h/\lambda\), (4.6) gives

\[
 \begin{aligned}
 \varepsilon_{E,p}
 &=\Lambda q^{-1-3/p}\log h,\\
 \varepsilon_{\nu,p}
 &=\mu q^{1-3/p}\log h.
 \end{aligned}
 \tag{4.7}
\]

For \(p=2\), after substituting (3.3),

\[
 \varepsilon_{E,2}
 ={9\over2L}\mu\left(q^{1/2}-q^{-11/2}\right)
 \sim {9\over2L}\mu q^{1/2},
 \qquad
 \varepsilon_{\nu,2}=12\mu q^{-1/2}\log n.
 \tag{4.8}
\]

Even after multiplication by the full backward gain \(h=q^{3/2}\), both
are eventually \(o(n^{-6})\) because
\(\mu_j=\nu(j!)^{-2}\). The heat-crossing number of the collar is

\[
 {\nu T\over r^2}
 ={8L\log n\log h\,q^2\over3(h^2-h^{-2})}
 \sim {32L(\log n)^2\over q}=o(1).
 \tag{4.9}
\]

Thus absolute self-Euler and viscous sizes are not the fatal terms.
For completeness, the same ledger gives

\[
 {\|U_c\|_2^2\over a^2\ell^3}\lesssim\Lambda^2q^{-5},
 \qquad
 {\nu\int_0^T\|\nabla U_c\|_2^2dt\over a^2\ell^3}
 \lesssim {9\over2L}\mu^2(1-q^{-6}).
 \tag{4.10}
\]

These are again absolute costs. Relative to the affine core itself,

\[
 {T\|R_E\|_p\over\|U_c\|_p}\lesssim\lambda T=\log h,
 \tag{4.11}
\]

so factorial smallness relative to the parent does not establish fidelity
of the localized selector.

### 4.3 The strain-speed-independent parent/collar exposure

Let \(V\) be an \(O(a)\) parent field. For a stationary cutoff,
\(V\cdot\nabla U_c\) is \(O(a\lambda)\) on the collar. Its integrated
size relative to the parent \(L^p\) scale is

\[
 {T\|V\cdot\nabla U_c\|_p\over a\ell^{3/p}}
 \lesssim q^{-3/p}\log h.
 \tag{4.12}
\]

The right-hand side is independent of \(\lambda\), and it is the natural
size when the collar interaction is nondegenerate. In \(L^2\), the
backward focus \(h=q^{3/2}\) turns this bound into

\[
 hq^{-3/2}\log h=\log h.
 \tag{4.13}
\]

Slowing the strain cannot improve this product estimate. A genuinely
co-moving Lagrangian cutoff could cancel the leading material-transport
term. After following the parent velocity at the core center, a Lipschitz
parent remainder gains \(Kr=q^{-1}\), leaving the favorable scale

\[
 q^{-1}\log h=n^{-8}\,12\log n=o(n^{-6}).
 \tag{4.14}
\]

Equation (4.14) is only the scale available after that cancellation; the
co-moving construction and its Leray/pressure control are open.

There is an equivalent wake formulation. Unless it is algebraically
canceled, the natural collar mixing coefficient is \(O(\lambda)\), whose
time integral is \(O(\log h)\). The corresponding direct bound on the
retained C140 wake \(b^2=n^{-4}\) is

\[
 O(b^2\log h),
 \tag{4.15}
\]

which does not imply the pre-chart \(b^3=n^{-6}\) allowance: the ratio of
the natural scales is
\[
 {b^2\log h\over b^3}=12n^2\log n\longrightarrow\infty.
\]
This is a failure of the direct bound, not a lower bound on the actual
coupling. LCE would require the integrated off-diagonal wake-to-active
block to be \(O(b)\), an improvement by \(b/\log h\) over generic
localization.

## 5. Exact surviving target

The global selector settles the carrier algebra but not the stage. The
affine route survives only if one proves a **Lagrangian collar-slaving
estimate** with all of the following properties on the existing \(A_2\)
cell:

1. leading parent transport is canceled in a co-moving frame;
2. the localized child is not carried through the heat-killed envelope
   direction (4.1);
3. the global pressure/cutoff wake remains in the retained channel; and
4. the backward-weighted wake-to-active collar operator is \(O(b)\).

Without those structural cancellations, (4.1) is fatal to the material
envelope and the stationary-cutoff estimates (4.13)--(4.15) do not close.
With them, the scalar residual powers in (4.8), (4.9), and (4.14) fit the
C127 schedule. This note proves neither the required cancellations nor a
universal stationary-cutoff no-go, and it does not claim a
Navier--Stokes stage map or blow-up.
