# A full-mode analytic tail for the retained-wake return derivative

**Date:** 2026-09-06  
**Claim:** C199  
**Status:** positive PDE estimate with explicit constants; conditional finite-section inverse test. No nonzero profile or finite-section singular-value certificate is supplied.

This executes the analytic part of Route D's pre-registered endpoint
derivative test. It distinguishes two derivatives that have different
high-frequency behavior. Prescribing a complete terminal value requires
inverting the forward derivative \(S\); a same-space return equation has
derivative \(CS-I\). Backward heat obstructs the first inverse. It does
not obstruct the second merely because \(S\) contains forward heat.
No existing prescribed-endpoint obstruction is withdrawn.

## 1. Full linearized Navier--Stokes estimate

Use the normalized \(2\pi\)-periodic torus and the real Hilbert space
\(H=L^2_{\sigma,0}\) of mean-zero divergence-free vector fields. Let
\(U(t,x)\) be a smooth divergence-free reference solution on \([0,T]\),
with fixed forcing if a force is present, and let \(\mu,T>0\). Define

\[
 M_0=\sup_{0\le t\le T}\|U(t)\|_\infty,\qquad
 M_1=\sup_{0\le t\le T}\|\operatorname{sym}\nabla U(t)\|_{
 L^\infty(\mathrm{op})}.
 \tag{1.1}
\]

The full endpoint linearization \(S\) is the solution operator of

\[
 \partial_t v-\mu\Delta v+
 \mathbb P\operatorname{div}(U\otimes v+v\otimes U)=0,
 \qquad v(0)=h.                                      \tag{1.2}
\]

On a suitable strong-solution domain this is the endpoint derivative;
here its bounded linear extension to \(H\) is the object being estimated.
No nonlinear Navier--Stokes endpoint map on an open \(L^2\) ball is
asserted. There is no charge, shell, or wake deletion in (1.2). Let \(P_K\) project
onto all Fourier modes \(0<|k|\le K\), and set \(Q_K=I-P_K\), for
\(K\ge1\). Then

\[
 \boxed{\|S\|_{H\to H}\le e^{M_1T},\qquad
 \|S Q_K\|_{H\to H}
 \le {e^{M_1T}\over K}
       \sqrt{{1\over\mu T}+{4M_0^2\over\mu^2}}.}
 \tag{1.3}
\]

The estimate applies to the full linearized PDE at every fixed positive
viscosity. It is polynomial in \(\mu^{-1}\) at fixed \(M_0,M_1,T\),
and is not a viscosity-uniform statement.

### Proof

Write \(f(t)=\|v(t)\|_2^2\) and
\(y(t)=\|(-\Delta)^{-1/2}v(t)\|_2^2\). Transport is skew in the
\(L^2\) pairing, so

\[
 f'(t)\le 2M_1 f(t),\qquad
 f(s)\ge e^{-2M_1(T-s)}f(T)\quad(0\le s\le T).          \tag{1.4}
\]

Pair (1.2) with \((-\Delta)^{-1}v\). The Leray projection is
self-adjoint and fixes this test field. Since
\(\|U\otimes v+v\otimes U\|_2\le2M_0\|v\|_2\) and
\(\|\nabla(-\Delta)^{-1}v\|_2=y^{1/2}\),

\[
 {1\over2}y'+\mu f\le 2M_0 f^{1/2}y^{1/2}
 \le {\mu\over2}f+{2M_0^2\over\mu}y.
 \tag{1.5}
\]

Thus, with \(a=4M_0^2/\mu\),

\[
 y'+\mu f\le ay,\qquad
 \mu\int_0^T e^{-as}f(s)\,ds\le y(0).                 \tag{1.6}
\]

Keep the weighted dissipation integral; discarding it and bounding
\(y(T)\) instead would introduce an unnecessary \(e^{aT}\). Combining
(1.4) and (1.6) gives the sharper inequality

\[
 f(T)\le
 {e^{2M_1T}y(0)\over
  \mu\int_0^T e^{-(a-2M_1)s}\,ds}.
 \tag{1.7}
\]

Because \(M_1\ge0\), the denominator integral is at least
\(\int_0^T e^{-as}ds\ge T/(1+aT)\). The last inequality follows
from \(e^x\ge1+x\) for \(x\ge0\), including \(a=0\) by continuity.
For \(h=Q_Kh\), Parseval gives \(y(0)\le K^{-2}\|h\|_2^2\).
Substitution proves (1.3). The calculation first holds for smooth data
and extends to \(H\) by the energy estimate.

## 2. A certified finite block suffices for the whole inverse

Let \(C:H\to H\) be the actual bounded linear exit chart, or the
derivative of an exit chart at the endpoint, and suppose a rigorous bound
\(\|C\|\le c\) is available. Put

\[
 L=CS,\qquad M=c e^{M_1T},\qquad
 \delta_K={M\over K}
     \sqrt{{1\over\mu T}+{4M_0^2\over\mu^2}}.
 \tag{2.1}
\]

Suppose a validated finite-dimensional enclosure establishes

\[
 \sigma_{\min}\big(P_K(I-L)P_K|_{P_KH}\big)\ge s>0.
 \tag{2.2}
\]

This must be the compression of the full evolution, including leakage
that returns to the retained block. Integrating a truncated equation
without a certified truncation error does not prove (2.2).

Define \(\Lambda=1+M/s\). If

\[
 \boxed{\delta_K\Lambda<1,}
 \tag{2.3}
\]

then the full infinite-dimensional linearized return operator is invertible and

\[
 \boxed{\|(I-CS)^{-1}\|_{H\to H}
 \le {\Lambda\over1-\delta_K\Lambda}.}
 \tag{2.4}
\]

Indeed let \(A_K=P_K(I-L)P_K|_{P_KH}\) and \(B=I-LP_K\).
The exact identity

\[
 B^{-1}=I+LP_K A_K^{-1}P_K
 \tag{2.5}
\]

gives \(\|B^{-1}\|\le\Lambda\). Since
\(I-L=B-LQ_K=B(I-B^{-1}LQ_K)\), (2.3) yields a convergent Neumann
series and (2.4). No common flag, self-adjointness, normality, or
commutation of \(C\) with \(P_K\) is required. In particular, all
off-block transport is charged in \(M\) and \(\delta_K\).

For fixed \(M_0,M_1,T,c,s\), a sufficient explicit cutoff is

\[
 K\ge 2M\left(1+{M\over s}\right)
       \sqrt{{1\over\mu T}+{4M_0^2\over\mu^2}},
 \qquad \|(I-CS)^{-1}\|\le 2\left(1+{M\over s}\right).
 \tag{2.6}
\]

The number \(s\) still has to be certified at the chosen cutoff. Raising
\(K\) does not guarantee that (2.2) remains true.

## 3. Exact zero-background destruction test

At \(U=0\), \(S=e^{\mu T\Delta}\), so (2.1) improves to

\[
 M=c e^{-\mu T},\qquad
 \delta_K=c e^{-\mu T K^2}.                            \tag{3.1}
\]

By contrast, prescribing every terminal Fourier coefficient requires
\(S^{-1}\), whose multiplier is \(e^{\mu T|k|^2}\). This is the
backward-heat factor in the older prescribed-endpoint tests.

For a complete PDE example take the deliberately simple same-space
chart \(C=3I\) and \(\mu T=1\). Then

\[
 \boxed{\|(I-3e^\Delta)^{-1}\|_{H^r_{\sigma,0}\to
 H^r_{\sigma,0}}\le11\quad\text{for every real }r.}
 \tag{3.2}
\]

This operator is not a contraction perturbation on its lowest shell:
\(3/e>1\). Nevertheless its inverse is bounded on every mode. Indeed
\(8/3<e<11/4\), whence

\[
 \left|1-{3\over e}\right|={3\over e}-1>{1\over11},
 \qquad
 1-3e^{-|k|^2}\ge1-{3\over(8/3)^2}={37\over64}
 \quad(|k|^2\ge2).
 \tag{3.3}
\]

The constant 11 therefore follows by Parseval, on the full torus PDE
space rather than a Galerkin section. This is a test of the inverse
mechanism, not a nonzero return profile or the actual stage chart.

Viscosity-uniform invertibility is a separate question. With \(C=I\),

\[
 {1\over\mu T}\le
 \|(I-e^{\mu T\Delta})^{-1}\|
 ={1\over1-e^{-\mu T}}\le1+{1\over\mu T}.
 \tag{3.4}
\]

With \(C=cI\), \(c>1\), there are exact kernels at
\(\mu T=(\log c)/m^2\), \(m=1,2,\ldots\), using the transverse
mode \(k=(m,0,0)\). These resonances accumulate at zero viscosity.
For \(0\le c<1\), the full inverse is instead bounded by
\(1/(1-c)\), uniformly in positive viscosity. None of these scalar
examples determines the spectrum of a physical dilation chart.

## 4. What Route D must actually certify

The next finite-section experiment must specify whether its unknown is
a prescribed initial datum with an independent prescribed endpoint, a
same-space return state, or a nonautonomous stage graph. They have
different derivatives. For a graph with viscosity update \(\mu'\), the
linearization contains \(h(\mu')\), not an automatic same-fiber
\(-h(\mu)\). Equations (2.2)--(2.4) do not certify that graph operator.

Similarly the torus rescaling \(u(x)\mapsto u(x/q)\) is not a
well-defined global periodic chart for general \(q>1\). On
\(\mathbb R^3\), the physical chart
\(Cu(y)=g^{-1}Q^Tu(a+Qy/q)\) has exact \(L^2\) norm
\(q^{3/2}/g\), but a bounded Fourier ball there is infinite dimensional.
The PDE estimate (1.3) can be proved there on its homogeneous Sobolev
domain; it does not turn that low-frequency ball into a finite matrix.
Spatial localization and an exterior estimate would also be needed.

At the actual stage viscosities, (2.6) deteriorates as viscosity tends
to zero, and neither \(s\), \(M_0\), nor \(M_1\) is currently certified
uniformly. No nonlinear inverse-function theorem on \(L^2\),
growing-order Gevrey bound, terminally flat force, nonlinear trapping, or
singularity follows. The positive result is the explicit
full-mode tail estimate and the correct finite-block-to-PDE inverse
criterion; applying it to the actual retained-wake chart remains open.

**Checker:** `checks/retained_wake_endpoint_inverse_tail_c199.py` verifies
the exact Young-inequality square, the rational constants, the
finite-rank inverse identity with both off-block couplings present,
and the all-mode heat example. It does not validate a nonzero PDE
trajectory or supply the missing finite-section enclosure (2.2).
