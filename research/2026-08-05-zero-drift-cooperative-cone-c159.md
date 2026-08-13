# C159: a cooperative-cone certificate for the zero-drift Kelvin orbit

**Date:** 2026-08-05

**Status:** exact reduction plus outward-rounded computer-assisted certificate

**Checker:**
[checks/zero_drift_cooperative_cone_c159.py](../checks/zero_drift_cooperative_cone_c159.py)

## 0. Claim boundary

This note stays on the \(\delta=4/5\) zero-drift orbit of C152.  It replaces
the floating trace experiment in C153 by a much smaller certificate.  No
amplitude column is integrated.  Instead, the transverse Kelvin equation
is put in a periodic frame in which it is cooperative and has a constant
strict subsolution.  The checker validates only the base orbit and one
scalar covector quadrature.

The conclusion is

\[
  \boxed{\rho(M)>1,\qquad \det M=1,\qquad \operatorname {tr}M>2.}
  \tag{0.1}
\]

Here \(M\) is the real two-dimensional Kelvin-amplitude monodromy on the
off-plane periodic covector specified in C152.  This is a principal-cocycle
result.  It does **not** prove a finite-frequency localized packet, a
relative leakage bound, or the nonlinear one-cell stage map.

## 1. Exact scalar reconstruction of the periodic covector

Use

\[
 n={N\over\sqrt3},\qquad g=\nabla f,\qquad
 U=N\times g,\qquad h=|g|^2,qquad c=|U(X_0)|^2={378\over25}.
 \tag{1.1}
\]

On \(f=0\), \(|U|^2=3h\).  Along the orbit, both \(k\cdot U\) and
\(k\cdot n\) are
constant.  Write

\[
 m=k\cdot n=\sqrt3\,\beta,
 \qquad
 k={c\over3h}U+\beta N+\gamma g.                  \tag{1.2}
\]

The exact velocity gradient is

\[
 A=C_NH-\sqrt2\,Ng^T,\qquad H=\nabla^2f.           \tag{1.3}
\]

Since \(AN=0\), \(A^TN=-3\sqrt2\,g\), and
\(g'=-A^Tg\), substitution of (1.2) into \(k'=-A^Tk\) gives the scalar
equation.  The simplification here is exact: on \(f=0\), one has
\(HN=0\) and \(\operatorname{tr}H=\Delta f=-2f=0\).  In the orthogonal
basis \((g,U,N)\), this implies
\(U^THU=-3g^THg\), and hence
\(g^T(A+A^T)U=-2U^THU\).  Therefore

\[
 \boxed{
 \gamma'=3\sqrt2\,\beta
       +{2c\over3}{U^THU\over h^2},
 \qquad \gamma(0)=0.}                              \tag{1.4}
\]

The C152 periodicity condition is exactly \(\gamma(T)=0\).  Thus the
three-component covector validation reduces to the phase orbit and one
scalar quadrature.

For completeness, the period parameters used by the checker have an AGM
form which removes all numerical quadrature.  Put

\[
 a=y_-,\quad c_*=y_+-y_-,\quad d_*=y_+-1,
 \quad \mathfrak m=1-{a d_*\over c_*}.              \tag{1.5}
\]

Indeed, \(t=\tan\theta\) turns the C152 period integral into

\[
 \int_0^\infty {dt\over\sqrt{(t^2+a)(d_*t^2+c_*)}},
\]

and \(t=\sqrt a\tan\phi\) turns this exactly into
\(K(1-ad_*/c_*)/\sqrt{c_*}\).  Thus

\[
 T={4\over3\delta}{K(\mathfrak m)\over\sqrt{c_*}},
 \tag{1.6}
\]

and, with primes denoting energy derivatives at zero,

\[
 {T'\over T}=
 \left\{{E/K\over2\mathfrak m(1-\mathfrak m)}
             -{1\over2\mathfrak m}\right\}\mathfrak m'
 -{c_*'\over2c_*}.                                  \tag{1.7}
\]

The standard AGM sequences give

\[
 K={\pi\over2\operatorname {AGM}(1,\sqrt{1-\mathfrak m})},
 \qquad
 {E\over K}=1-\sum_{j\ge0}2^{j-1}c_j^2,             \tag{1.8}
\]

where \(c_0=\sqrt{\mathfrak m}\) and
\(c_{j+1}=(a_j-b_j)/2\).  Finally

\[
 \beta=-{T'\over T}{c\over3\sqrt2}.                 \tag{1.9}
\]

The checker evaluates (1.5)--(1.9) with rational outward rounding; the AGM
tail is bounded after quadratic contraction, and \(\pi\) is enclosed by the
Machin identity.  The complete parameter intervals--not rounded midpoint
values--are then inserted into the residual certificate below.

## 2. Exact cooperative transverse frame

Let \(p=k-mn\) and \(D=|p|^2\), so \(q=|k|^2=D+m^2\).  The two vectors

\[
 E_1=P_k n=n-{m\over q}k,
 \qquad E_2=k\times n                                  \tag{2.1}
\]

are orthogonal, span \(k^\perp\), and return periodically with \(k\).
Writing a Kelvin amplitude as \(v=E_1z_1+E_2z_2\) gives \(z'=Bz\).

In coordinates with third axis \(n\), the Beltrami and axial-invariance
identities put

\[
 A=\begin{pmatrix}S&0\\ l^T&0\end{pmatrix},
 \qquad S^T=S,\quad \operatorname {tr}S=0,
 \quad l\cdot(k\times n)=-\sqrt2\,(k\cdot U)=-\sqrt2c.
 \tag{2.2}
\]

A direct differentiation of (2.1), using \(p'=-Sp-ml\), gives

\[
\begin{aligned}
 B_{11}&={m\,l\cdot p\over D},\\
 B_{22}&={2p\cdot Sp+m\,l\cdot p\over D},\\
 B_{21}&={m^2\sqrt2c\over qD},\\
 B_{12}&={2m\,p\cdot S(k\times n)+\sqrt2c(D-m^2)\over D}.
\end{aligned}                                           \tag{2.3}
\]

These are physical-time coefficients.  The checker works on unit time
\(s=t/T\) and therefore multiplies (2.3) by \(T\).

## 3. The finite outward certificate

The phase system is augmented by
\(\cos a,\sin a,\cos b,\sin b\).  On each of
32 equal unit-time panels the checker constructs a degree-40 decimal-rational
Taylor polynomial.  It treats every generated decimal coefficient as an
exact rational number and recomputes the polynomial residual with directed
outward rounding.  The exact algebraic initial sine is separately enclosed,
so reference generation supplies no premise of the proof.  The scalar
\(\gamma\) polynomial is checked against (1.4), with \(h^{-2}\) kept as a
rational expression rather than trusted through a floating inverse variable.

The certified budget is

\[
 \int_0^1\!|R_{\rm phase}|_\infty ds<7\times10^{-25},
 \qquad
 \int_0^1\!|R_\gamma|ds<4\times10^{-17}.             \tag{3.1}
\]

The directed panel jumps and the algebraic initial enclosure are charged in
these same componentwise budgets; the reported bounds are for the total
defects used by Gronwall comparison, not just the smooth panel interiors.

On the two-micro-unit phase tube, panelwise interval automatic
differentiation gives

\[
 \mu_\infty(D F_{\rm phase})<42,
 \qquad |D_{\rm phase}F_\gamma|_1<300.               \tag{3.2}
\]

Together with \(e^{42}<2\times10^{18}\), the residuals, algebraic initial-data
enclosures, and AGM parameter enclosures imply the self-consistent tube

\[
 \|\Delta(\cos a,\sin a,\cos b,\sin b)\|_\infty<2\times10^{-6},
 \qquad |\Delta\gamma|<8\times10^{-4}.                \tag{3.3}
\]

The checker then splits the same path into 1024 evaluation cells and
inserts (3.3) directly into the rational formulas (1.2) and (2.3).  Its
outward lower bounds are deliberately rounded down to

\[
\begin{aligned}
 B_{12}&>32, & B_{21}&>{9\over10},\\
 B_{11}+{3\over20}B_{12}&>{7\over10},
 &{20\over3}B_{21}+B_{22}&>{1\over5}.
\end{aligned}                                           \tag{3.4}
\]

The raw unrounded certificate margins are larger; (3.4) is the theorem
premise.

## 4. Cone conclusion

Let \(w=(1,3/20)^T\).  The first line of (3.4) says that \(B\) is Metzler.  The
second says componentwise

\[
 Bw>{1\over5}w.                                        \tag{4.1}
\]

Comparison for cooperative systems therefore yields

\[
 M w\ge e^{1/5}w>w.                                    \tag{4.2}
\]

The positive cone is invariant, so the period map has a positive real
eigenvalue \(\rho>1\).  C153's exact Liouville identity gives \(\det M=1\).
Hence the other transverse multiplier is \(\rho^{-1}>0\), and

\[
 \operatorname {tr}M=\rho+\rho^{-1}>2.                 \tag{4.3}
\]

This proves (0.1) without using the former floating matrix entries.  The
large numerical trace remains useful orientation only; it is not a premise
of the certificate.

## 5. What remains open

C159 closes the principal zero-drift Floquet sign.  It does not close the
one-cell theorem.  The remaining load-bearing tasks are unchanged:

1. construct a localized finite-frequency packet that follows this
   periodic ray for \(O(\log q)\) returns while respecting C154's fiber shear;
2. prove a relative complement/leakage estimate at the C147 microseed
   scale; and
3. integrate the retained packet and its wake into the full unforced
   Navier--Stokes stage map.
