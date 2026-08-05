# Three material phases: lattice gap and all-order charge majorant

**Date:** 2026-08-03

**Status:** exact material-phase and combinatorial lemmas; local positive
stress chart by openness; zero-charge/global-wake equation still open.

**Scope:** the all-order multiphase endpoint inverse left by C53.  This note
does not solve the zero mode, the viscous wake, or the complete transition.

## 1. Outcome

Using three phases instead of the pointwise-minimal two removes two
structural problems at once.

1. Three coordinate phases are exactly curl-free and exactly transported.
2. Every nonzero integer cross-charge has a uniform material lattice gap;
   there are no small divisors at any order.

Moreover, the number of charges through nonlinear order \(r\) is only
\(O(r^3)\), and the exact Gevrey-2 convolution estimate absorbs this
polynomial without changing the factorial class.

Thus the all-order phase algebra splits sharply:

- every \(k\ne0\) charge is an elliptic high-frequency pressure/Hodge
  channel with a uniform divisor; and
- \(k=0\) is the sole non-elliptic channel and must carry the positive stress,
  energy transfer, and global viscous pressure wake.

The remaining theorem is difficult, but its difficulty is no longer phase
combinatorics or Diophantine resonance.

## 2. Exact material coordinate phases

Let \(X(t,a)\) be the flow map of a smooth incompressible transition and
\[
                         F(t,a)=D_aX(t,a).
\]
On the periodic box, take the three circle-valued initial phases
\[
                         \Phi_i(0,x)=x_i\pmod{2\pi}.
\]
Transport them by
\[
                         \Phi_i(t,X(t,a))=a_i.                  \tag{2.1}
\]
Then
\[
 (\partial_t+u\cdot\nabla)\Phi_i=0,\qquad
 \nabla_x\Phi_i=F^{-T}e_i.                                     \tag{2.2}
\]

For an integer charge \(k=(k_1,k_2,k_3)\), put
\[
                         \Phi_k=k_1\Phi_1+k_2\Phi_2+k_3\Phi_3.
\]
Equation (2.2) gives the exact identity
\[
                         \nabla\Phi_k=F^{-T}k.                  \tag{2.3}
\]

If
\[
                         \|F\|+\|F^{-1}\|\le K,                 \tag{2.4}
\]
then
\[
             K^{-1}|k|\le|\nabla\Phi_k|\le K|k|.               \tag{2.5}
\]
In particular,
\[
                         k\ne0\quad\Longrightarrow\quad
                         |\nabla\Phi_k|\ge K^{-1}.              \tag{2.6}
\]

The leading pressure denominator on charge \(k\) is
\(|\nabla\Phi_k|^2\), so (2.6) gives a uniform inverse bound
\(O(K^2)\).  If \(K_j\) is polynomial in the cascade index, that loss is
absorbed by the polynomial carrier ledger.

## 3. Positive stress chart with transported kernels

At \(F=I\), define
\[
\begin{aligned}
 R_1^0&={q\over2}(e_2\otimes e_2+e_3\otimes e_3),\\
 R_2^0&={q\over2}(e_1\otimes e_1+e_3\otimes e_3),\\
 R_3^0&={q\over2}(e_1\otimes e_1+e_2\otimes e_2).
\end{aligned}                                                  \tag{3.1}
\]
Then
\[
 R_i^0e_i=0,\qquad
 R_i^0>0\ \hbox{on }e_i^\perp,\qquad
 R_1^0+R_2^0+R_3^0=qI.                                        \tag{3.2}
\]

For fixed kernels \(e_i\), the sum map from the three symmetric transverse
blocks to \(\operatorname{Sym}_3\) is onto.  Indeed:

- \(R_1\) controls the \(23\) entry;
- \(R_2\) controls the \(13\) entry;
- \(R_3\) controls the \(12\) entry; and
- their six diagonal variables span all three diagonal entries.

Surjectivity and transverse positivity are open.  Therefore there are
\(\delta,c>0\) such that, whenever
\[
                         \|F-I\|<\delta,\qquad
                         \|Q-qI\|<\delta q,                     \tag{3.3}
\]
one can solve
\[
 Q=R_1+R_2+R_3,\qquad
 R_iF^{-T}e_i=0,\qquad
 R_i\ge cq\ \hbox{on }(F^{-T}e_i)^\perp.                       \tag{3.4}
\]

Each positive transverse block has a two-polarization factorization
\[
                         R_i=a_i\otimes a_i+b_i\otimes b_i,
 \qquad a_i,b_i\perp F^{-T}e_i.                                \tag{3.5}
\]
Thus it is the averaged covariance of one real oscillatory colour with two
transverse polarizations; no negative stress coefficient is introduced.

Thus the three fixed coordinate phases give an actual curl-free,
materially transported positive-stress chart on every sufficiently short
transition interval.  The pointwise-minimal two-direction chart remains
interesting, but it is not needed to solve phase integrability.

## 4. Exact charge count

Starting from the six fundamental real charges
\(\{\pm e_1,\pm e_2,\pm e_3\}\), a product of at most \(r\) factors has
integer charge
\[
                         |k|_1\le r.                            \tag{4.1}
\]

The number of integer lattice points in the three-dimensional
\(\ell^1\) ball is the cross-polytope Ehrhart polynomial
\[
\begin{aligned}
 N_3(r)
 &=\#\{k\in\mathbb Z^3:|k|_1\le r\}\\
 &=1+6r+12{r\choose2}+8{r\choose3}.                            \tag{4.2}
\end{aligned}
\]
Equivalently, the shell \(|k|_1=r\), \(r\ge1\), has
\[
                         N_3(r)-N_3(r-1)=4r^2+2.               \tag{4.3}
\]

This growth is cubic, not exponential.  In particular,
\[
                         N_3(r)\le(2r+1)^3.                    \tag{4.4}
\]

## 5. Gevrey-2 convolution majorant

The exact factorial convolution identity is
\[
\begin{aligned}
 \sum_{a=0}^r {r\choose a}(a!)^2((r-a)!)^2
 &=(r!)^2\sum_{a=0}^r{1\over{r\choose a}}\\
 &\le3(r!)^2.                                                   \tag{5.1}
\end{aligned}
\]
The bound holds for every \(r\ge1\): the two endpoint terms give \(2\),
the two next terms give \(2/r\), and the remaining binomial reciprocals
fit in the unused margin.  More explicitly, for \(r\ge4\) each remaining
term is at most \({r\choose2}^{-1}\), so their sum is at most
\(2(r-3)/(r(r-1))\); the cases \(r\le3\) are immediate.

Every polynomial charge multiplicity can be absorbed into the exponential
constant of a Gevrey majorant.  For example,
\[
                         (r+1)^3\le8^r,\qquad r\ge1.            \tag{5.2}
\]
Consequently a recursion whose order-\(r\) coefficient is bounded by a
fixed polynomial number of bilinear convolutions preserves
\[
                         M_r=C^r(r!)^2                         \tag{5.3}
\]
after increasing \(C\).  The three-phase charge lattice does not consume
an additional factorial.

## 6. The exact remaining split

For \(k\ne0\), (2.6) makes the principal material pressure equation
elliptic.  A WKB pressure/Hodge recursion can be solved charge by charge,
with losses depending polynomially on \(K\) and on the differentiation
order.  Equations (4.2) and (5.1) show that summing all such charges remains
inside the existing Gevrey-2 budget.

For \(k=0\), the principal phase derivative vanishes.  This is not a small
divisor; it is a genuinely different equation.  It contains:

1. the positive covariance \(Q\);
2. the work-carrying mean interaction with the parent;
3. the incompressible low-frequency pressure; and
4. the second- and higher-jet global viscous wake.

The zero charge must therefore be retained as a coupled macroscopic
unknown.  Treating it as an elliptic oscillatory corrector would repeat the
one-carrier and local-pressure errors already found by audit.

## 7. Remaining theorem

The constructive target can now be stated without hidden phase qualifiers:

> Build the zero-charge material stress/wake evolution and the nonzero-charge
> elliptic WKB recursion simultaneously on a short three-coordinate
> material interval, with \(C^r(r!)^2\) bounds and a rank-five terminal
> endpoint derivative.

This note closes curl-free phase existence, all-order nonzero-charge
separation, small divisors, and charge-count factorial loss.  It does not
close the zero-charge evolution or its global pressure wake.
