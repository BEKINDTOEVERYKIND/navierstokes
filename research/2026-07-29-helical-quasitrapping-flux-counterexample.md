# Counterexample to the advertised helical quasi-trapping flux estimate

Date: 2026-07-29

## Scope and verdict

A July 2026 preprint record by Luca Eliseo Pavesi, *Global Regularity for
the Three-Dimensional Incompressible Navier--Stokes Equations via
Geometric Frustration and Helical Quasi-Trapping*, advertises the
universal estimate

\[
 |\Pi(K)|
 \leq
 C\,\frac{E_{>K}^{1/2}E^{1/2}}{K}
\tag{0.1}
\]

with an absolute constant \(C\), where \(\Pi(K)\) is nonlinear spectral
energy flux, \(E_{>K}\) is energy above the cutoff, and \(E\) is total
energy.

As stated, (0.1) is false.  A single real divergence-free Fourier triad
violates it by an arbitrarily large factor.  The obstruction is already
visible from homogeneity: flux is cubic in velocity, while the right-hand
side of (0.1) is quadratic.  The explicit family below also shows that
the factor \(K^{-1}\) has the wrong frequency scaling.

This note audits the central estimate advertised in the public preprint
record.  It is not a page-by-page audit of the PDF, which was unavailable
from the host during this check.  A differently normalized or
solution-dependent estimate would be a different statement, but it
would not be (0.1) with an absolute constant.

The calculation is reproduced by
[`checks/quasitrapping_flux_counterexample.py`](../checks/quasitrapping_flux_counterexample.py).

---

## 1. An exact real triad

Work on \(\mathbb T^3\).  For an integer \(N\geq1\), set

\[
 p=(N,0,0),\qquad
 q=(0,N,0),\qquad
 k=p+q=(N,N,0).
\tag{1.1}
\]

For \(A>0\), prescribe the nonzero Fourier coefficients

\[
\begin{aligned}
 \widehat u(p)&=A e_2,&
 \widehat u(-p)&=A e_2,\\
 \widehat u(q)&=A e_3,&
 \widehat u(-q)&=A e_3,\\
 \widehat u(k)&=-iA e_3,&
 \widehat u(-k)&= iA e_3.
\end{aligned}
\tag{1.2}
\]

They obey the reality condition
\(\widehat u(-n)=\overline{\widehat u(n)}\), and every coefficient is
orthogonal to its wavevector.  In physical variables, up to the Fourier
normalization,

\[
 u(x)=
 2A e_2\cos(Nx_1)
 2A e_3\cos(Nx_2)
 2A e_3\sin(N(x_1+x_2)).
\tag{1.3}
\]

It is smooth, real, periodic, and divergence free.

---

## 2. Exact nonlinear flux through \(K=N\)

Use the Leray-projected Euler nonlinearity

\[
 \widehat{\mathcal N(u)}(n)
 =
 iP_n
 \sum_{r+s=n}
 \big(\widehat u(r)\cdot s\big)\widehat u(s).
\tag{2.1}
\]

At \(n=k\), the only nonzero pair in the sum is \(p+q=k\), together
with its reversed ordering.  The reversed scalar product vanishes, and

\[
\begin{aligned}
 \widehat u(p)\cdot q&=AN,\\
 \widehat u(q)\cdot p&=0.
\end{aligned}
\tag{2.2}
\]

Since \(e_3\perp k\),

\[
 \widehat{\mathcal N(u)}(k)
 =iNA^2e_3.
\tag{2.3}
\]

The Euler contribution to the time derivative is
\(-\mathcal N(u)\).  Therefore

\[
\begin{aligned}
 2\operatorname{Re}
 \left(
 \overline{\widehat u(k)}
 \cdot[-\widehat{\mathcal N(u)}(k)]
 \right)
 &=
 2\operatorname{Re}
 \left(
 iAe_3\cdot[-iNA^2e_3]
 \right)\\
 &=2NA^3.
\end{aligned}
\tag{2.4}
\]

The conjugate mode \(-k\) contributes the same amount.  At that instant,
all other modes above \(N\) have zero amplitude and hence zero
instantaneous energy derivative.  Thus the nonlinear flux into
\(\{|n|>N\}\) has magnitude

\[
 \boxed{|\Pi(N)|=4NA^3}
\tag{2.5}
\]

under the convention \(E=\sum_n|\widehat u(n)|^2\).  Conventions with a
factor \(1/2\) change only the fixed factor \(4\).

The energies are

\[
 E=6A^2,\qquad
 E_{>N}=2A^2.
\tag{2.6}
\]

Consequently the right-hand side of (0.1) is

\[
 C\,\frac{\sqrt{12}\,A^2}{N}.
\tag{2.7}
\]

The claimed estimate would force

\[
 4NA^3
 \leq C\frac{\sqrt{12}A^2}{N},
\qquad\text{hence}\qquad
 N^2A\leq C\frac{\sqrt3}{2}.
\tag{2.8}
\]

This fails either by sending \(A\to\infty\) at fixed \(N\), or by sending
\(N\to\infty\) at fixed \(A\).

---

## 3. Why helicity does not repair the universal claim

Helical coordinates are a unitary change of basis on each transverse
Fourier plane.  They reorganize the coefficient in (2.1); they do not
change the physical flux (2.5).  The triad can be decomposed into its
helical signs, and at least one resulting helical interaction must carry
the nonzero total transfer.

Special homochiral or maximally aligned data may have very small or zero
flux.  That does not imply a bound for arbitrary divergence-free data.
Numerical cancellation in one aligned test likewise cannot establish
(0.1).

Any valid universal estimate must respect both homogeneities:

* cubic scaling in the velocity amplitude; and
* one spatial derivative in the nonlinear transfer.

For example, standard estimates contain a gradient or a positive power
of the active frequency.  An absolute \(K^{-1}\) gain using only two
energy factors cannot hold.

---

## Reference

* Luca Eliseo Pavesi,
  [*Global Regularity for the Three-Dimensional Incompressible
  Navier--Stokes Equations via Geometric Frustration and Helical
  Quasi-Trapping*](https://zenodo.org/records/21194906), version 3,
  July 2026.  The record states (0.1) in its abstract.

