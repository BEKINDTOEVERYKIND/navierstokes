# C198: full-wake return energy and circulation demands

Date: 2026-09-06. Status: proved identities and necessary conditions for
smooth, exact, unforced returns; explicitly conditional on their stated
norms for approximate returns and inviscid limits. This does not construct
a return or rule out Navier--Stokes blow-up.

This audits Route B of
`research/2026-08-30-counterfactual-success-portfolio.md` before its proposed
axisymmetric computation. A full finite-energy state cannot have a
nonzero regular Euler fixed point under the proposed chart. A viable graph
must carry nonvanishing viscous dissipation as its normalized viscosity
tends to zero, or fail the stated compactness/energy assumptions. The
circulation identity gives a second, independently measurable restriction
in the axisymmetric class. Three finite returns alone are not excluded.

## 1. Exact full-state identities

Let a smooth finite-energy solution on all of R^3 solve

\[
 \partial_t u+(u\cdot\nabla)u+\nabla p=\mu\Delta u,
 \qquad \nabla\cdot u=0,
\]

on a finite interval [0,T]. Assume enough decay/integrability to justify
the energy identity. Write E(v)=||v||_2^2, using unnormalized Lebesgue
measure, and

\[
 I=\int_0^T\|\nabla u(t)\|_2^2\,dt,
 \qquad ({\cal C}v)(y)=g^{-1}Q^Tv(a+q^{-1}Qy),
 \qquad Q\in SO(3).
\]

The Jacobian is q^3, hence

\[
 \boxed{E({\cal C}u(T))=k\{E(u(0))-D\},\quad
 k=\frac{q^3}{g^2},\quad D=2\mu I.}                 \tag{1.1}
\]

No discarded child band, wake cutoff, or localized energy can replace the
full E in this identity. The identity follows by pairing the equation
with 2u: incompressibility removes transport and pressure and integration
by parts gives E(T)+2mu I=E(0).

For an exact graph with mu_{n+1}=rho mu_n, rho=q/g, let X_n be its stage
states and u_n the corresponding full evolutions. Then

\[
 E_{n+1}=k(E_n-D_n),\qquad
 \boxed{\sum_{j=0}^{N-1}k^{-j}D_j
      =E_0-k^{-N}E_N.}                               \tag{1.2}
\]

In particular the actual three-return diagnostic is

\[
 D_0+k^{-1}D_1+k^{-2}D_2=E_0-k^{-3}E_3.              \tag{1.3}
\]

A zero-dissipation three-return orbit is compatible with E_3=k^3 E_0.
Thus neither (1.1) nor bounded finite-stage output by itself invalidates
the registered three-return experiment. If E_N<=B E_0, (1.2) gives the
explicit demand sum k^{-j}D_j >= (1-B k^{-N})E_0, useful when its right-hand
side is positive.

## 2. Approximate return: include the measured norm ratio

Let R_n=X_{n+1}-C u_n(T_n), e_n=||R_n||_2. The reverse triangle inequality
implies the exact interval

\[
 \boxed{
 E_n-k^{-1}(\sqrt{E_{n+1}}+e_n)^2
 \le D_n\le
 E_n-k^{-1}\max\{0,\sqrt{E_{n+1}}-e_n\}^2.}           \tag{2.1}
\]

Intersect it with [0,E_n]. This is a necessary consistency check, with
explicit constants, for an independently replayed unforced computation.
Equivalently, set

\[
 b_n=E_{n+1}-k(E_n-D_n),\quad
 h_n=2e_n\sqrt{E_{n+1}}+e_n^2.
\]

Then |b_n|<=h_n, and the three-return enclosure is

\[
 \boxed{\left|
 D_0+k^{-1}D_1+k^{-2}D_2-E_0+k^{-3}E_3
 \right|\le k^{-1}h_0+k^{-2}h_1+k^{-3}h_2.}          \tag{2.2}
\]

The portfolio's X_mu norm is inhomogeneous and dominates L2. A reported
relative X_mu residual epsilon does not automatically mean e_n/sqrt(E_n)
<=epsilon: multiply by the reported normalization factor
||X_reference||_{X_mu}/sqrt(E_n), or report the relative L2 residual itself.
The weighted circulation residual below requires its own measurement.

An equal-L2 illustration, not an assumption on the whole graph, is
particularly transparent. For

\[
 q=\frac65,\qquad \frac76<\gamma=\frac{\log g}{\log q}
          <\frac{37}{25},
 \qquad E_{n+1}=E_n>0,
\]

we have

\[
 k^{-1}=q^{2\gamma-3}
 <(5/6)^{1/25}<\frac{993}{1000},
 \quad\text{since }(993/1000)^{25}>5/6.              \tag{2.3}
\]

If in addition e_n/sqrt(E_n)<=1/100000, then (2.1) yields

\[
 \boxed{\frac{D_n}{E_n}
 >1-\frac{993}{1000}\left(\frac{100001}{100000}\right)^2
 =\frac{69801399007}{10^{13}}
 >\frac{69}{10000}.}                                 \tag{2.4}
\]

For exact equal-L2 returns the stronger D_n/E_n>7/1000 holds. If
T_n<=T_max and ||grad u_n(t)||_2<=M, the approximate equal-L2 assumptions
therefore require

\[
 \boxed{\mu_n>\frac{69E_n}{20000T_{\max}M^2}.}         \tag{2.5}
\]

Thus a fixed resolved-gradient cap cannot certify these equally normalized
returns to arbitrarily small viscosity. For a general graph, use (2.1)
with its actual E_{n+1}/E_n; (2.4) must not be applied by silently resetting
that ratio to one.

## 3. What a bounded inviscid face must carry

Assume q<g<q^(3/2), exact returns for all n, mu_n=rho^n mu_0, and
E_n->E_*>0. Formula (1.1), without any PDE limiting argument, proves

\[
 \boxed{D_n\longrightarrow (1-k^{-1})E_*,\qquad
 \mu_n I_n\longrightarrow\tfrac12(1-g^2/q^3)E_*>0.}   \tag{3.1}
\]

At the registered parameter range the limiting fractional loss exceeds
7/1000. In particular a uniformly bounded T_n and uniformly bounded
stage-integrated squared gradient cannot coexist with this energy limit.
Nonvanishing dissipation is a requirement, not a result of the existing
construction. Endpoint-only X_mu bounds do not control intermediate
evolutions. A bound throughout the stage has the stronger consequence
proved next.

More generally, if 0<E_min<=E_n<=E_max and D_n->0, choose n_0 with
D_n<[(1-k^{-1})/2]E_min for n>=n_0. Then

\[
 E_{n+1}>\frac{k+1}{2}E_n,
 \qquad
 E_{n_0+N}>\left(\frac{k+1}{2}\right)^N E_{\min},     \tag{3.2}
\]

contradicting E_max for large N. Thus uniform positive and finite energies
also exclude vanishing dissipation, even without convergence of E_n.

The registered norm supplies a more direct quantitative test without
assuming a regular Euler limit. For every frequency magnitude s>=0,

\[
 2\sqrt\mu\,s^2\le s+\mu s^3,
 \qquad s(\sqrt\mu s-1)^2\ge0.
\]

Plancherel and the inhomogeneous Sobolev norm convention therefore give

\[
 \boxed{D_n\le\sqrt{\mu_n}
     \int_0^{T_n}\|u_n(t)\|_{{\cal X}_{\mu_n}}^2\,dt
 \le\sqrt{\mu_n}T_n K_n^2,\quad
 K_n=\sup_{0\le t\le T_n}\|u_n(t)\|_{{\cal X}_{\mu_n}}.} \tag{3.3}
\]

Here the nonnegative mu^2 H^(5/2) term is simply available in X_mu;
the first two terms already prove the bound, with constant one. For the
equal-L2 approximate return of (2.4), this yields the explicit requirement

\[
 \boxed{K_n>
   \left(\frac{69E_n}{10000T_n}\right)^{1/2}
       \mu_n^{-1/4}.}                                  \tag{3.4}
\]

More generally use its measured lower bound from (2.1) in place of
69E_n/10000. An exact infinite orbit with 0<E_min<=E_n<=E_max, bounded T_n,
and uniformly bounded full-trajectory X_mu norm is impossible by (3.2)
and (3.3). Bounded endpoint X_mu norms alone do not establish the required
full-trajectory bound. The search must resolve the transient peaks or
increasing stage duration required by (3.3).

Finally suppose the exact graph/evolutions converge strongly enough to
produce a smooth finite-energy Euler return X_*=C S_T^0 X_* and to retain
its full energy. Euler energy conservation gives

\[
 \boxed{\|X_*\|_2^2=k\|X_*\|_2^2,\quad k>1
        \quad\Longrightarrow\quad X_*=0.}             \tag{3.5}
\]

This obstruction is independent of axisymmetry. It kills a nonzero regular
finite-energy Euler fixed face. It does not kill a singular-viscosity
graph, an energy defect in a weak limit, energy escaping every fixed
similarity ball, or a local inner profile of infinite similarity-space
energy. Each such alternative requires its corresponding full-state or
outer matching analysis. None is already a certified construction.

## 4. The axisymmetric circulation restriction

Restrict to the portfolio's chart a=0, Q e_z=e_z, and smooth axisymmetric
fields. Put G=r u^theta. The theta equation, multiplied by r, is

\[
 (\partial_t+u^r\partial_r+u^z\partial_z)G
 =\mu(\partial_r^2-r^{-1}\partial_r+\partial_z^2)G.   \tag{4.1}
\]

Assume bounded G and boundary/decay conditions permitting the maximum
principle. Smoothness gives G=0 at the axis. At an interior extremum with
r>0 the first derivatives vanish and the second-order term has the correct
sign, yielding ||G(T)||_infinity<=||G(0)||_infinity. This standard principle
is also explicitly stated by Lei--Zhang, *Criticality of the Axially
Symmetric Navier--Stokes Equations*,
[arXiv:1505.02628v2](https://arxiv.org/abs/1505.02628v2). No global regularity
claim for arbitrary swirl is imported from that paper.

The exact chart satisfies G_new(r,z)=rho G_T(r/q,z/q). Therefore

\[
 \boxed{A_{n+1}\le\rho A_n,\quad
 A_n=\|rX_n^\theta\|_\infty,\quad
 A_n\le\rho^n A_0=\frac{\mu_n}{\mu_0}A_0.}           \tag{4.2}
\]

For approximate returns let a_n=||r R_n^theta||_infinity. Then

\[
 A_{n+1}\le\rho A_n+a_n,
 \quad A_3\le\rho^3 A_0+\rho^2a_0+\rho a_1+a_2.     \tag{4.3}
\]

The parameter range gives

\[
 \rho=(5/6)^{\gamma-1}<(5/6)^{1/6}<\frac{971}{1000},
 \qquad(971/1000)^6>5/6.
\]

Hence the explicit three-return diagnostic is

\[
 \boxed{A_3\le\frac{915498611}{10^9}A_0+
 (971/1000)^2a_0+(971/1000)a_1+a_2.}                  \tag{4.4}
\]

A proposed equally normalized circulation A_{n+1}=A_n>0 necessarily has
weighted residual a_n/A_n>29/1000. Nonzero swirl at each of only three
stages remains compatible with (4.4). Nonzero normalized swirl at every
stage also remains compatible when it decreases to zero; a uniform
positive lower bound on A_n is impossible. The ratio A_n/mu_n need not
tend to zero, so (4.2) is not a small-swirl Navier--Stokes regularity theorem.

If an exact regular Euler fixed face exists with bounded circulation,
(4.2) at the fixed point gives G_*=0. If also
eta_*=omega_*^theta/r is bounded, its swirl-free Euler evolution is passive
transport and its chart multiplier is 1/(gq^2). The fixed point forces
||eta_*||_infinity<=[1/(gq^2)]||eta_*||_infinity, hence eta_*=0. A smooth
finite-energy divergence-free and curl-free field on R^3 is zero. This
independently checks (3.5) in the axisymmetric class, under the additional
stated circulation and vorticity norms.

## 5. Change to the actual search

Keep three full states and measure E_n, I_n, E_{n+1}, ||R_n||_2, A_n,
and ||rR_n^theta||_infinity on every replay; certify (2.1), (2.2), and
(4.4) before interpreting a small optimizer objective. Include full tail
bounds in these quantities. Do not fix the circulation to one at each
stage. Do not require a smooth nonzero finite-energy Euler endpoint at
mu=0. A proposed graph with E_n converging to E_*>0 must exhibit the dissipation
in (3.1); a graph with bounded nonconvergent positive energies must satisfy
(1.2) and (3.2). Both must resolve the full-trajectory norm cost in (3.3).
Energy accumulation
for three steps is allowed and must be distinguished from a bounded
infinite cascade. These are equations for the already proposed object,
not a successor mechanism.

Checker: `checks/full_wake_energy_circulation_c198.py` verifies all rational
constant comparisons, exact finite-stage identities, scaling powers, and
residual bounds on rational algebraic instances. The PDE integrations and
maximum-principle argument above are analytic proofs, not simulated by
the checker.
