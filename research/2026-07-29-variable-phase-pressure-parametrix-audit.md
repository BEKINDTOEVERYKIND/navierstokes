# Variable-phase pressure audit for the single-carrier construction

Date: 2026-07-29

## Status and claim boundary

This note audits the pressure/Leray gap left in Section 8.2 of
`2026-07-29-forward-multiphase-parametrix.md`.

The result is positive, but conditional in an important and precise sense.
On nonzero angle charges, the variable-phase pressure has an explicit
semiclassical expansion whose \(K^{-n}\) coefficient costs at most \(n\)
slow derivatives.  Consequently Gevrey-\(s\) slow data produce
\(C^n(n!)^s\), not \(C^n(n!)^{2s}\), coefficient growth.  In particular,
Gevrey two remains compatible with the published polynomial-carrier
condition

\[
K_j=j^A,\qquad A>4.
\]

There is, however, no inverse on the space of *all* angle-mean-zero
profiles.  The extended operator has exact fast gauge modes

\[
H(\vartheta-K\phi(x)).
\]

The correct statement is a slow-envelope (equivalently, semiclassically
microlocal) parametrix.  Its exponentially small divergence error can then
be removed by the ordinary physical pressure solve.  A future theorem must
state this restriction; mean zero in \(\vartheta\) by itself is
insufficient.

This closes the formal pressure derivative count.  It does **not** construct
the nonlinear material-phase hierarchy, prove uniform Gevrey bounds for its
low flow, or route the nonlocal pressure tail through the annular wake.

---

## 1. Exact operator and formal inverse

Let

\[
{\cal D}_i=\partial_i+K\xi_i(x)\partial_\vartheta,
\qquad
\xi=\nabla\phi,
\qquad
a=|\xi|^2\ge c_0^2>0.
\]

For a scalar profile,

\[
\Delta_K:={\cal D}_i{\cal D}_i
=K^2A_2+KA_1+A_0,
\tag{1.1}
\]

where

\[
\begin{aligned}
A_2&=a\,\partial_\vartheta^2,\\
A_1&=\left(2\xi\cdot\nabla+\operatorname{div}\xi\right)
       \partial_\vartheta,\\
A_0&=\Delta.
\end{aligned}
\tag{1.2}
\]

For a vector profile \(F\),

\[
\operatorname{div}_K F=KB_1F+B_0F,
\qquad
B_1F=\xi\cdot\partial_\vartheta F,
\qquad
B_0F=\operatorname{div}F.
\tag{1.3}
\]

There is an exact conjugation behind these formulas.  Define the angle
shear

\[
(S_KU)(x,\vartheta)
=U(x,\vartheta+K\phi(x)).
\tag{1.C1}
\]

On the torus, \(\phi\) is circle-valued and \(K\) is an integer, so this
is globally well-defined; locally the identity needs no topological
assumption.

Then

\[
\partial_iS_K=S_K{\cal D}_i,
\qquad
\Delta_xS_K=S_K\Delta_K,
\qquad
\operatorname{div}_x(S_KF)=S_K\operatorname{div}_KF.
\tag{1.C2}
\]

If \(\mathbb P_x\) denotes the ordinary Leray projector in \(x\), applied
at each fixed \(\vartheta\), the exact extended projector is therefore

\[
\boxed{\mathbb P_K=S_K^{-1}\mathbb P_xS_K.}
\tag{1.C3}
\]

At \(\vartheta=0\), (1.C3) is precisely the physical Leray projection of
\(F(x,K\phi(x))\).  The conjugation preserves each integer angle charge.
What it does **not** preserve is a \(K\)-uniform slow-envelope norm:
\(S_K\) inserts \(x\)-oscillations of frequency \(K\).  The expansion
below is the tame semiclassical form of this exact but non-tame identity.

All coefficients in (1.1)--(1.3) are independent of \(\vartheta\).
Therefore angle averaging commutes with every operator.  The charge-zero
pressure is the ordinary slow elliptic problem and must be included in the
retained low field.  This note applies to

\[
QF:=F-\frac1{2\pi}\int_{\mathbb T}F\,d\vartheta.
\]

On the range of \(Q\),

\[
A_2^{-1}=a^{-1}\partial_\vartheta^{-2}
\tag{1.4}
\]

is bounded in every fixed angle Sobolev norm.  Put \(h=K^{-1}\).  Then

\[
h^2\Delta_K=A_2+hA_1+h^2A_0
\tag{1.5}
\]

and

\[
h^2\operatorname{div}_K F=hB_1F+h^2B_0F.
\tag{1.6}
\]

Seek

\[
p^{(N)}=h\sum_{n=0}^{N}h^np_n.
\tag{1.7}
\]

The coefficients are determined explicitly by

\[
\boxed{
\begin{aligned}
A_2p_0&=B_1F,\\
A_2p_1&=B_0F-A_1p_0,\\
A_2p_n&=-A_1p_{n-1}-A_0p_{n-2},
\qquad n\ge2.
\end{aligned}}
\tag{1.8}
\]

In particular,

\[
p_0
=a^{-1}\partial_\vartheta^{-1}(\xi\cdot F).
\tag{1.9}
\]

The corresponding pressure-gradient expansion is

\[
\begin{aligned}
{\cal D}p^{(N)}
={}&\xi\partial_\vartheta p_0\\
&+\sum_{n=1}^{N}h^n
\left(\xi\partial_\vartheta p_n+\nabla p_{n-1}\right)
+h^{N+1}\nabla p_N.
\end{aligned}
\tag{1.10}
\]

Its principal term is exactly

\[
\xi\partial_\vartheta p_0
=\frac{\xi\otimes\xi}{|\xi|^2}F.
\tag{1.11}
\]

Thus the leading high-charge Leray operator is the pointwise projection
onto \(\xi^\perp\), as geometric optics predicts.  It is order zero, not
order \(K\).

Substitution gives the exact finite-order residual

\[
\begin{aligned}
&(A_2+hA_1+h^2A_0)p^{(N)}
 -(hB_1+h^2B_0)F\\
&\quad =
h^{N+2}(A_1p_N+A_0p_{N-1})
+h^{N+3}A_0p_N,
\end{aligned}
\tag{1.12}
\]

with \(p_{-1}=0\).  Equivalently, in the unscaled pressure equation the
residual begins at \(h^N\).

---

## 2. The derivative budget is \(n\), not \(2n\)

It is useful to write (1.8) as

\[
\begin{aligned}
p_0&=a^{-1}\partial_\vartheta^{-1}(\xi\cdot F),\\
p_1&=a^{-1}\partial_\vartheta^{-2}\operatorname{div}F+T_1p_0,\\
p_n&=T_1p_{n-1}+T_2p_{n-2},
\end{aligned}
\tag{2.1}
\]

where

\[
\begin{aligned}
T_1u
&=-a^{-1}\partial_\vartheta^{-1}
\left(2\xi\cdot\nabla+\operatorname{div}\xi\right)u,\\
T_2u
&=-a^{-1}\partial_\vartheta^{-2}\Delta u.
\end{aligned}
\tag{2.2}
\]

\(T_1\) costs one slow derivative and advances the \(h\)-index by one.
\(T_2\) costs two slow derivatives and advances the \(h\)-index by two.
Every path contributing to \(p_n\) therefore costs at most \(n\) slow
derivatives, apart from a harmless fixed offset in the initial
coefficients.  Counting two derivatives at every recursion step would
incorrectly assign order \(2n\); the \(A_0\) term carries \(h^2\), not
\(h\).

Here is a coefficient estimate sufficient for the cascade ledger.  Suppose
\(s\ge1\), \(a\ge c_0^2\), and, on a fixed normalized stage,

\[
\|\partial^\alpha \xi\|_\infty
+\|\partial^\alpha F\|_{H^m_\vartheta}
\le C_0L^{|\alpha|}(|\alpha|!)^s
\tag{2.3}
\]

for every slow multi-index \(\alpha\).  Gevrey classes are closed under
products, and \(a^{-1}\) has the same Gevrey order because \(a\) is bounded
away from zero.  Induction on \(n+|\alpha|\) in (2.1) then gives

\[
\boxed{
\|\partial^\alpha p_n\|_{H^m_\vartheta}
\le C_1C_2^{n+|\alpha|}
\big((n+|\alpha|+1)!\big)^s.}
\tag{2.4}
\]

The proof uses only the elementary factorial inequality

\[
r!\,q!\le(r+q)!
\tag{2.5}
\]

and the fact that the number of product and recursion partitions is at
most exponential.  The inverse angle derivatives in (2.2) are bounded on
nonzero integer charges, so there is no small-divisor loss.

The \(h^n\) coefficient of (1.10) obeys the same Gevrey order.  The first
term uses \(p_n\) without another slow derivative; the second uses one
derivative of \(p_{n-1}\).  Similarly, (1.12) and (2.4) imply, in every
fixed slow seminorm,

\[
\|\operatorname{div}_K(F-{\cal D}p^{(N)})\|
\le C(C_2h)^N(N!)^s.
\tag{2.6}
\]

Fixed physical derivatives add only fixed powers of \(K\), which do not
alter the optimal exponential scale.

This order is also generally sharp.  The leading part of \(T_1^n\)
contains \(n\) repeated directional derivatives of a Gevrey-\(s\)
coefficient.  Without stronger regularity, those derivatives can grow
like \(C^n(n!)^s\).

### No hidden doubling with the nonlinear hierarchy

Suppose the order-\(r\) nonlinear coefficient already obeys a
Gevrey-\(\sigma\) majorant and the slow envelopes have order \(s\).
Let

\[
\rho=\max\{\sigma,s\}.
\]

Applying the pressure coefficient of order \(n\) to a hierarchy
coefficient of order \(r\) gives, after absorbing exponential constants,

\[
(n!)^\rho(r!)^\rho\le((n+r)!)^\rho.
\tag{2.7}
\]

Summing over \(n+r=N\) adds only an exponential factor.  Thus the combined
hierarchy remains Gevrey-\(\rho\); the two indices do not multiply the
Gevrey exponent.  In the preferred construction,

\[
\sigma=s=2
\]

is therefore still Gevrey two after pressure is included.

---

## 3. From the formal profile to the exact physical pressure

The finite expansion need not be promoted to a global inverse of
\(\Delta_K\).  Evaluate it on the physical phase:

\[
p_{\rm app}^K(x)
=p^{(N)}(x,K\phi(x)).
\]

Let

\[
r^K
=\operatorname{div}\left(
F(x,K\phi(x))-\nabla p_{\rm app}^K(x)
\right).
\tag{3.1}
\]

This physical residual has zero spatial mean exactly, because it is a
divergence.  Solve the ordinary physical problem

\[
\Delta q^K=r^K,\qquad \int q^K=0.
\tag{3.2}
\]

On the torus, standard elliptic estimates give

\[
\|\nabla q^K\|_{H^m}
\le C_m\|r^K\|_{H^{m-1}}.
\tag{3.3}
\]

Consequently, if (2.6) is \(e^{-cK^{1/s}}\) (or the weaker
\(e^{-cj^2}\) required by the stage ledger), the difference between the
formal material-phase Leray correction and the exact physical Leray
correction has the same exponential size, up to fixed polynomial factors.
The exact correction need not have a pure charged-profile representation;
it belongs in the flat residual/wake bookkeeping.

On \(\mathbb R^3\), the same argument is local plus a Calderón--Zygmund
estimate.  Compactly supported high-charge data can produce a noncompact
pressure tail, but its failure to be represented by the local WKB series
is exponentially small under the nonstationary Gevrey estimate below.
The tail is not literally compact and must not be silently discarded.

### Analytic Gaussian core followed by a flat cutoff

There is a useful alternative to imposing a compact Gevrey-two envelope on
the active transition equations.  Suppose the normalized active packet,
material phase, and low flow are analytic on a strip of radius bounded
below independently of the stage, and the packet has Gaussian spatial
decay

\[
|F(y,\vartheta)|\le Ce^{-c|y|^2}.
\tag{3.4}
\]

Then (2.4) has \(s=1\).  Optimal pressure truncation and the
nonstationary-phase argument in Section 4 give

\[
\text{pressure error}\le Ce^{-cK}.
\tag{3.5}
\]

One should apply this analytic parametrix to the **uncut** Gaussian core.
Only afterwards choose a physical cutoff radius

\[
R_j\asymp j\ell_j,
\qquad\text{equivalently}\qquad |y|\asymp j.
\tag{3.6}
\]

Every fixed derivative of the field changed by that cutoff is bounded by
a polynomial in \(j,\ell_j^{-1}\), times

\[
e^{-cj^2}.
\tag{3.7}
\]

The divergence correction in the cutoff annulus and its ordinary Leray
pressure correction have the same size in fixed scaled Sobolev norms.
They may have a noncompact pressure tail, but the order-zero physical
Leray projector does not destroy the \(e^{-cj^2}\) amplitude.  Thus there
is no need to run the all-order analytic parametrix through a nonanalytic
compact cutoff.

For \(K_j=j^A\), (3.5) is smaller than (3.7) for every \(A>2\); the
existing \(A>4\) choice has ample margin.  This route is compatible with
an analytic local transition theory on an arbitrarily short normalized
stage: choosing the outer ratio \(r>1\) close to one shortens the required
one-step deformation.  What remains to be proved is that the analytic
radius of the coupled low material flow stays uniformly positive on that
stage.  If the radius collapses with \(j\), the constant \(C_j\) gate in
Section 5 reappears.

There is also a spatial-domain caveat.  A cutoff at \(|y|\asymp j\) is
flat, but it requires the uncut analytic transition equations and their
background to be valid on a normalized ball whose radius grows like
\(j\).  A Gaussian tail at a fixed normalized radius \(C\) is only
\(e^{-cC^2}\), not \(e^{-cj^2}\).  Therefore a construction available only
on one fixed \(O(1)\) cell cannot discard everything outside that cell.
It must either extend the analytic transition uniformly to the growing
ball or include the intermediate Gaussian tail in the global wake.  Making
the normalized time interval short does not by itself solve this spatial
compatibility problem.

---

## 4. The exact gauge kernel and physical-frequency aliases

Mean zero in \(\vartheta\) does not make \(\Delta_K\) globally invertible.
For every mean-zero periodic \(H\),

\[
g_K(x,\vartheta)=H(\vartheta-K\phi(x))
\tag{4.1}
\]

satisfies

\[
{\cal D}g_K=0,\qquad \Delta_Kg_K=0.
\tag{4.2}
\]

This is not a contradiction to the expansion above.  The slow derivatives
of \(g_K\) are of size \(K^{|\alpha|}\); it lies outside every
\(K\)-uniform slow-envelope class.  Upon physical evaluation,

\[
g_K(x,K\phi(x))=H(0),
\]

so the ambiguity is a physically constant pressure gauge.

For an affine integer phase \(\phi=w\cdot x\), a profile Fourier mode
\((n,\ell)\) has physical frequency

\[
n+K\ell w.
\tag{4.3}
\]

The gauge alias is \(n=-K\ell w\).  A slow band \(|n|\le M\) is exactly
separated from it when

\[
M<cK.
\tag{4.4}
\]

The polynomial construction has
\(M_j\asymp j^2/\log j\) and \(K_j=j^A\), so (4.4) has enormous margin.

For a variable material phase there is no exact slow Fourier band after
composition, but the replacement is a Gevrey nonstationary-phase bound.
For a charge \(\ell\ne0\), amplitude \(A_\ell\), and physical Fourier
frequency \(n\), the phase is

\[
\Psi(x)=K\ell\phi(x)-n\cdot x.
\]

If \(|n|\le c_0K|\ell|/2\), then

\[
|\nabla\Psi|\ge c_0K|\ell|/2.
\]

Repeated integration by parts with

\[
\frac{\nabla\Psi}{i|\nabla\Psi|^2}\cdot\nabla
\]

and Gevrey-\(s\) derivative bounds gives

\[
\left|
\widehat{A_\ell e^{iK\ell\phi}}(n)
\right|
\le
C^{N+1}(N!)^s(K|\ell|)^{-N}.
\tag{4.5}
\]

Optimizing \(N\) yields

\[
\left|
\widehat{A_\ell e^{iK\ell\phi}}(n)
\right|
\le Ce^{-c(K|\ell|)^{1/s}}.
\tag{4.6}
\]

Thus physical low-frequency leakage is exponentially small, not zero.
The honest theorem is a microlocal parametrix separated from
\(\eta=-K\ell\xi\), followed by the exact correction (3.2).  Any claim of
an exact inverse on all mean-zero profiles would be false.

---

## 5. Material coordinates do not change the derivative count

Let \(x=X(t,y)\) be the low flow, \(g^{\alpha\beta}\) the inverse material
metric, and \(J=\det DX\).  With the fixed material direction \(w\), the
lifted Laplace--Beltrami operator is

\[
L_Ku
=J^{-1}(\partial_\alpha+Kw_\alpha\partial_\vartheta)
\left[
Jg^{\alpha\beta}
(\partial_\beta+Kw_\beta\partial_\vartheta)u
\right].
\tag{5.1}
\]

It again has the form

\[
L_K=K^2\widetilde A_2+K\widetilde A_1+\widetilde A_0,
\tag{5.2}
\]

with

\[
\begin{aligned}
\widetilde A_2
&=(g^{\alpha\beta}w_\alpha w_\beta)
\partial_\vartheta^2,\\
\widetilde A_1
&=\left[
2g^{\alpha\beta}w_\beta\partial_\alpha
+J^{-1}\partial_\alpha(Jg^{\alpha\beta}w_\beta)
\right]\partial_\vartheta,\\
\widetilde A_0&=\Delta_g.
\end{aligned}
\tag{5.3}
\]

Uniform ellipticity of \(DX\) gives

\[
g^{\alpha\beta}w_\alpha w_\beta\ge c|w|^2.
\]

The recurrence and its weighted derivative budget are therefore unchanged:
\(\widetilde A_1\) costs one derivative per one power of \(h\), while
\(\widetilde A_0\) costs two per two powers.  Incompressibility gives
\(J=1\), but this simplification is not needed.

The actual gate is uniform control of the material coefficients.  If their
Gevrey constant on stage \(j\) is \(C_j\), the pressure remainder behaves
schematically as

\[
\left(\frac{C_j}{K_j}\right)^{M_j}(M_j!)^s.
\tag{5.4}
\]

Uniform \(C_j\) closes.  Polynomial growth \(C_j\lesssim j^B\) also closes
after increasing the carrier exponent to

\[
A>B+2s.
\tag{5.5}
\]

Exponential growth \(C_j\gtrsim e^{cj}\) is different: with
\(M_j\asymp j^2/\log j\), its contribution is of order
\(e^{c j^3/\log j}\), which no fixed polynomial \(K_j=j^A\) can absorb.
The low-flow theorem must therefore give uniform or at worst polynomial
normalized Gevrey constants.  Pressure itself introduces no such growth,
but it cannot repair it.

---

## 6. Polynomial-carrier verdict

With a stage-uniform Gevrey-\(s\) estimate, (2.6) gives

\[
C^{M_j}(M_j!)^sK_j^{-M_j}.
\]

For

\[
M_j=\left\lfloor\eta\frac{j^2}{\log(e+j)}\right\rfloor,
\qquad
K_j=j^A,
\]

Stirling's formula yields

\[
\log\!\left(
C^{M_j}(M_j!)^sK_j^{-M_j}
\right)
=-\eta(A-2s)j^2+o(j^2).
\tag{6.1}
\]

Hence:

* compact Gevrey-\(s\) envelopes are allowed for every \(s>1\);
* Gevrey two gives the required \(e^{-cj^2}\) residual for \(A>4\);
* pressure does not force \(A>8\);
* increasing \(A\) absorbs polynomial, but not exponential, deterioration
  of the material Gevrey radius.

The variable-phase pressure/Leray operation is therefore not a fatal
obstruction to the single-carrier route.  The remaining theorem is now
more sharply stated:

1. construct the coupled low-flow/high-charge hierarchy with uniform or
   polynomially tame normalized Gevrey-two bounds;
2. keep its slow spectrum/microlocal support separated from the gauge set
   through order \(M_j\);
3. prove that the exponentially small exact pressure correction can be
   included in the global wake and terminal-force ledger;
4. retain the charge-zero pressure in the low material flow rather than
   trying to invert it with the high-charge parametrix.

These are substantial analytic tasks, but the feared repeated pressure
derivative loss does not occur.
