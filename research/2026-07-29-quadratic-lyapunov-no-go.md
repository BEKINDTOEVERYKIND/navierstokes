# No positive Fourier-diagonal quadratic Lyapunov function beyond energy

Date: 2026-07-29

## Result

Consider mean-zero, real, divergence-free fields on \(\mathbb T^3\).  Write the
helical Fourier decomposition

\[
 \widehat u(\mathbf k)
   =u_+(\mathbf k)h_+(\mathbf k)+u_-(\mathbf k)h_-(\mathbf k),\qquad
 i\mathbf k\times h_s(\mathbf k)=s|\mathbf k|h_s(\mathbf k).
\]

Let

\[
 Q_W(u)=\sum_{\mathbf k\ne0}\sum_{s=\pm}
 W_{\mathbf k,s}|u_s(\mathbf k)|^2.
\]

Assume the weights are real and reality-compatible,
\(W_{-\mathbf k,s}=W_{\mathbf k,s}\).  On real fields only this even part
is observable in \(Q_W\).

The following theorem closes the helical-diagonal quadratic-Lyapunov route.

**Theorem.** Suppose \(Q_W\) is conserved by every smooth 3D Euler solution
whose initial data have finite Fourier support.  Then there are constants
\(A,B\), independent of \(\mathbf k\), such that

\[
 W_{\mathbf k,s}=A+B\,s|\mathbf k|.
\]

Thus the only such quantities are linear combinations of energy and helicity.
If \(Q_W\) is nonnegative on both helical polarizations at arbitrarily large
wave number, then \(B=0\) and \(Q_W\) is a nonnegative multiple of energy.
In particular, a scalar multiplier

\[
 Q_w(u)=\sum_{\mathbf k\ne0}w(\mathbf k)|\widehat u(\mathbf k)|^2
\]

is conserved for all Euler flows only when \(w\) is constant.

The same conclusion holds if “conserved” is replaced by “the Euler
contribution to \(\dot Q_W\) has one sign for every field.”  That contribution
is cubic and changes sign under \(u\mapsto-u\), so one-sidedness forces it to
vanish identically.

## Exact triad proof

Fix a noncollinear Fourier triad

\[
 \mathbf k+\mathbf p+\mathbf q=0,\qquad
 k=|\mathbf k|,\quad p=|\mathbf p|,\quad q=|\mathbf q|.
\]

For a helicity triple \((\sigma,\tau,\rho)\in\{\pm1\}^3\), the three modal
energy rates in the corresponding Waleffe triad have a common real amplitude
factor \(\Xi\):

\[
\begin{aligned}
 \frac d{dt}|u_\sigma(\mathbf k)|^2&=(\tau p-\rho q)\Xi,\\
 \frac d{dt}|u_\tau(\mathbf p)|^2&=(\rho q-\sigma k)\Xi,\\
 \frac d{dt}|u_\rho(\mathbf q)|^2&=(\sigma k-\tau p)\Xi.
\end{aligned}
\]

The geometric coupling is nonzero for every sign triple when the wavevector
triangle is nondegenerate.  By choosing the phase of one amplitude, \(\Xi\)
can take either sign.  Therefore conservation for every amplitude and phase
requires

\[
\begin{split}
0={}&W_{\mathbf k,\sigma}(\tau p-\rho q)
  +W_{\mathbf p,\tau}(\rho q-\sigma k)\\
 &+W_{\mathbf q,\rho}(\sigma k-\tau p)
\end{split}
\tag{1}
\]

for all eight helicity triples.

Set

\[
 W_{\mathbf k,\sigma}=A_k+\sigma B_k
\]

and similarly at \(\mathbf p,\mathbf q\).  Expanding (1) in the independent
Walsh characters on \(\{\pm1\}^3\) gives

\[
\begin{aligned}
0={}&\sigma k(A_q-A_p)+\tau p(A_k-A_q)+\rho q(A_p-A_k)\\
&+\sigma\tau(pB_k-kB_p)
 +\sigma\rho(kB_q-qB_k)
 +\tau\rho(qB_p-pB_q).
\end{aligned}
\]

Every coefficient must vanish.  Hence, on this triad,

\[
 A_k=A_p=A_q,\qquad
 \frac{B_k}{k}=\frac{B_p}{p}=\frac{B_q}{q}.
\]

Any two nonzero lattice wavevectors can be connected by two noncollinear
triads sharing a third wavevector: choose an integer vector \(\mathbf r\)
parallel to neither, and use

\[
(\mathbf k,\mathbf r,-\mathbf k-\mathbf r),\qquad
(\boldsymbol\ell,\mathbf r,-\boldsymbol\ell-\mathbf r).
\]

The two weights at the shared wavevector, \(W_{\mathbf r,+}\) and
\(W_{\mathbf r,-}\), fix both constants.  This proves the global formula
\(W_{\mathbf k,s}=A+Bs|\mathbf k|\).

Finally, if \(B\ne0\), choosing the helicity \(s=-\operatorname{sgn}B\) makes
\(A+Bs|\mathbf k|<0\) for all sufficiently large \(|\mathbf k|\).  Hence a
positive full-Euler invariant has \(B=0\).

## Explicit failure of the critical Sobolev candidate

Take

\[
\mathbf k=(1,0,0),\quad
\mathbf p=(0,1,0),\quad
\mathbf q=(-1,-1,0)
\]

and helicities \((\sigma,\tau,\rho)=(+,-,+)\).  For the homogeneous Sobolev
weight \(w(\mathbf k)=|\mathbf k|^{2\alpha}\), the coefficient multiplying
\(\Xi\) is

\[
 2(2^\alpha-1).
\]

It is nonzero for every \(\alpha\ne0\).  In particular, the critical
\(\dot H^{1/2}\) weight gives \(2(\sqrt2-1)\).  Changing one modal phase
reverses the sign, so the Euler nonlinearity can instantaneously increase or
decrease \(\dot H^{1/2}\).

## Stronger Navier--Stokes consequence

For viscosity \(\nu>0\),

\[
 \frac d{dt}Q_w(u)
 =T_w(u)-2\nu\sum_{\mathbf k}w(\mathbf k)|\mathbf k|^2
 |\widehat u(\mathbf k)|^2,
\]

where \(T_w\) is the cubic Euler transfer.  If \(w\) is nonconstant, the
theorem supplies a finite triad state \(u\) with \(T_w(u)>0\).  Replacing
\(u\) by \(\lambda u\) makes the positive term scale as \(\lambda^3\) and
the viscous term as \(\lambda^2\).  For sufficiently large \(\lambda\),
\(\dot Q_w(0)>0\).  Thus no nonconstant positive scalar Fourier weight is
even a universal quadratic Lyapunov function for the viscous equation.

## What remains outside the theorem

The theorem completely rules out scalar Fourier multipliers and, more
generally, multipliers diagonal in curl polarization.  It does not by itself
classify a translation-invariant multiplier that mixes the two helicities at
each wavevector, nor a nonlinear or solution-dependent functional.

There are nevertheless two hard restrictions on those possible escapes.

1. For any quadratic functional, a universal Navier--Stokes Lyapunov
   inequality implies exact Euler conservation: amplitude scaling makes the
   cubic Euler term dominate viscosity, and \(u\mapsto-u\) reverses that cubic
   term.
2. The positive \(H^{1/2}\)-like helicity mechanism really does work after
   projecting away one helical sector, but that is a different, decimated
   equation.  Reintroducing both helicities restores the obstruction above.

Accordingly, searching over radial, anisotropic, optimized, or computer-found
scalar Fourier weights cannot solve the full problem.  A viable
global-regularity attack must use a genuinely nonquadratic, non-diagonal, or
conditional geometric mechanism rather than another weighted energy.

## Primary references

- F. Waleffe, *The nature of triad interactions in homogeneous turbulence*,
  Physics of Fluids A 4 (1992), 350--363:
  https://doi.org/10.1063/1.858309
- L. Biferale and E. S. Titi, *On the Global Regularity of a
  Helical-decimated Version of the 3D Navier--Stokes Equations*:
  https://arxiv.org/abs/1303.1215
- L. Biferale, S. Musacchio and F. Toschi, *Inverse Energy Cascade in
  Three-Dimensional Isotropic Turbulence*:
  https://arxiv.org/abs/1111.1412
