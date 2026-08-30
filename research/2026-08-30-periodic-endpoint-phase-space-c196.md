# C196: exact periodic solenoidal endpoint profiles and the retained-band phase-space ceiling

**Date:** 2026-08-30

**Status:** exact real periodic curl construction with explicit constants;
exact Fourier-support ceilings; exact but conditional C193/C194 exponent
arithmetic; no common full-aperture Floquet/WKB witness, uniform
multi-beam parametrix, band retention, viscosity, nonlinear return, UVSR,
or singularity theorem

**Checker:**
[checks/periodic_endpoint_phase_space_c196.py](../checks/periodic_endpoint_phase_space_c196.py)

## 0. Verdict and correction to the positive chain

C182's \(bq^{9/8}\) estimate is an **upper bound**, not a constructed
baseline. C180 supplies \(q^3\) source--daughter labels, but its principal
action is norm-preserving transport. C193's displayed one-profile
principal witness has absolute endpoint concentration only \(bq^{3/8}\).
Thus “C193 fills C182's \(q^{3/8}\) deficit” was deficit arithmetic, not a
positive endpoint theorem.

There is a precise positive kinematic result. On a fixed-aperture periodic
band one can place \(q^2\) coherent carrier modes under a
three-dimensional envelope of bandwidth \(q^{1/3}\). The resulting field
is exactly real, periodic, and divergence free, has fixed \(L^2\) size, and
has point concentration at least \(c_*q^{3/2}\) with the explicit
constant (1.18). Replacing \(q^{1/3}\) by \(q^{1/12}\) gives at least
\(c_*q^{9/8}\). The two certified lower-bound exponents differ by
\(3/8\); no exact quotient of the full concentrations is claimed.

That construction is not yet a dynamic witness. C194 controls one compactly
supported \(\mathbb R^3\) beam in \(|f|\le1/10\), while the construction
below is a global periodic \(q^2\)-carrier superposition. Naively summing
the normalized beam errors loses \(\sqrt{q^2}=q\), which is larger than
both available error margins. A uniform Fourier-integral
almost-orthogonality theorem and a retained-band theorem are therefore
load-bearing, not bookkeeping.

Independently, every endpoint component projected into C180's narrower
retained slab has concentration at most

\[
                 8\delta^{3/2}\frac{q^{3/2}}{J^2}       \tag{0.1}
\]

for one oriented box, or \(\sqrt2\) times this constant after a disjoint
reality completion. The wider C176 slab similarly pays \(J^{-1/2}\).
Formal carrier--envelope labels cannot evade a bound on distinct output
Fourier modes.

The current \(A_2\) route therefore has a sharp phase-space specification:
prove a fixed-aperture dynamic band, prove an honest log-decorated
\(J^{-2}\)-taxed child, or prove a nonlinear mechanism that creates the
missing band. This is not an architecture-wide trichotomy.

## 1. Exact real periodic curl construction

Use normalized Haar measure on \(\mathbb T^3\). Let \(S_+\) contain \(M\)
integer carriers with \(S_+\cap(-S_+)=\varnothing\). Suppose their
directions lie within angle \(\delta\) of a unit vector
\(\widehat p_0\). Fix a real unit vector \(e\perp\widehat p_0\) and set

\[
 a_p=\frac{P_{p^\perp}e}{|P_{p^\perp}e|},\qquad
 a_{-p}=a_p.                                             \tag{1.1}
\]

Then \(p\cdot a_p=0\), \(e\cdot a_p\ge\cos\delta\), and

\[
 W(x)=\frac1{\sqrt{2M}}\sum_{p\in S_+}
             a_p\big(e^{ip\cdot x}+e^{-ip\cdot x}\big)   \tag{1.2}
\]

is exactly real and divergence free, with

\[
 \|W\|_2=1,\qquad e\cdot W(0)\ge\cos\delta\sqrt{2M}.      \tag{1.3}
\]

Define the real vector potential by

\[
 B_p=i\frac{p\times W_p}{|p|^2}.                         \tag{1.4}
\]

Because \(p\cdot W_p=0\),

\[
 ip\times B_p
 =-\frac{p\times(p\times W_p)}{|p|^2}=W_p,              \tag{1.5}
\]

so \(\operatorname{curl}B=W\) exactly.

For an integer \(L\ge1\), take the symmetric normalized Dirichlet envelope

\[
 \eta_L(x)=(2L+1)^{-3/2}
   \prod_{j=1}^3\sum_{r=-L}^L e^{irx_j}.                 \tag{1.6}
\]

It obeys

\[
 \|\eta_L\|_2=1,\qquad
 \eta_L(0)=(2L+1)^{3/2},\qquad \nabla\eta_L(0)=0.        \tag{1.7}
\]

Set

\[
                         V_L=\operatorname{curl}(\eta_LB). \tag{1.8}
\]

The field is exactly real, periodic, and solenoidal, and

\[
                         V_L(0)=\eta_L(0)W(0).           \tag{1.9}
\]

Suppose all blocks

\[
                 \pm p+[-L,L]^3,\qquad p\in S_+,        \tag{1.10}
\]

are pairwise disjoint, and let \(p_{\min}=\min|p|\). In the block of
\(p\), a sideband \(s\) changes \(W_p\) to

\[
 W_p-s\times\frac{p\times W_p}{|p|^2}.                  \tag{1.11}
\]

The second term has norm at most \(\alpha\), where

\[
                              \alpha=\frac{\sqrt3L}{p_{\min}}. \tag{1.12}
\]

Disjointness and Parseval yield

\[
       1-\alpha\le\|V_L\|_2\le1+\alpha,                 \tag{1.13}
\]

and hence the completely explicit lower bound

\[
 \frac{|V_L(0)|}{\|V_L\|_2}
 \ge \frac{\cos\delta}{1+\alpha}
       \sqrt{2M}\,(2L+1)^{3/2}.                         \tag{1.14}
\]

The quantity \(\alpha\) is a relative sideband/polarization correction;
there is no solenoidality error because (1.8) is an exact curl.

### 1.1 An explicit full-aperture family

Take \(q=m^{24}\), \(m\ge2\),

\[
 L_+=m^8=q^{1/3},\qquad L_-=m^2=q^{1/12},\qquad
 d=2L_++1.                                               \tag{1.15}
\]

Let \(K=\lfloor q/(4d)\rfloor\), and use

\[
 S_+=\{(10q+dr_1,dr_2,dr_3):0\le r_1,r_2,r_3<K\}.       \tag{1.16}
\]

Since \(d\le3m^8\) and \(\lfloor x\rfloor\ge x/2\) for \(x\ge2\),

\[
 M=K^3\ge\frac{q^2}{13824},\qquad p_{\min}\ge10q.        \tag{1.17}
\]

The translated blocks are disjoint. Their directions obey
\(\tan\delta\le\sqrt2/40\), and choosing \(e=e_2\) gives
\(e\cdot a_p\ge39/40\). Moreover

\[
 \alpha_+\le\frac{\sqrt3}{10}q^{-2/3},\qquad
 \alpha_-\le\frac{\sqrt3}{10}q^{-11/12}.
\]

Since both are at most \(\sqrt3/10\), (1.14) gives

\[
 {\cal C}(V_{L_-})\ge c_*q^{9/8},\qquad
 {\cal C}(V_{L_+})\ge c_*q^{3/2},\qquad
 c_*=\frac{39}{10\sqrt{13824}(1+\sqrt3/10)}.             \tag{1.18}
\]

Periodicity, reality, solenoidality, and both endpoint exponents are thus
simultaneously compatible. What is not proved is that the \(A_2\) Floquet
dynamics transports these profiles in covariant stable/expanding bundles.

## 2. Distinct-mode ceilings

Let a trigonometric vector field have Fourier support \(\Sigma\). At every
point,

\[
 |v(x)|\le\sum_{k\in\Sigma}|\widehat v(k)|
       \le\sqrt{|\Sigma|}\,\|v\|_2.                     \tag{2.1}
\]

This counts distinct output modes, not algebraic labels.

For C180's one-sided oriented box

\[
 |e_1\cdot(k-q\bar k)|\le\delta q/J,\quad
 |e_2\cdot(k-q\bar k)|\le\delta q/J,\quad
 |e_3\cdot(k-q\bar k)|\le\delta q/J^2,                  \tag{2.2}
\]

attach the disjoint unit cube centered at each lattice point. In the
\((e_1,e_2,e_3)\) frame every cube has half-width at most \(\sqrt3/2\).
When the smallest box half-width is at least \(\sqrt3/2\),

\[
 \#(\mathbb Z^3\cap\mathcal B^\sharp_{q,J})
 \le8\prod_{r=1}^3(w_r+\sqrt3/2)
 \le64\delta^3\frac{q^3}{J^4}.                          \tag{2.3}
\]

Equations (2.1)--(2.3) prove (0.1). A reality-completed pair of disjoint
boxes changes the constant by at most \(\sqrt2\), not the exponent. The
original C176 box of volume proportional to \(q^3/J\) analogously pays
\(J^{-1/2}\).

The phrase “\(q^{-1/4}\) tube” must specify its radial coordinate:

* a three-coordinate relative covector cube has three physical half-widths
  \(cq^{3/4}\), hence at most
  \(8(cq^{3/4}+\sqrt3/2)^3\) modes; when
  \(cq^{3/4}\ge\sqrt3/2\), this is at most \(64c^3q^{9/4}\) and its
  concentration is at most \(8c^{3/2}q^{9/8}\);
* a genuinely projective angular \(cq^{-1/4}\) tube in the fixed annulus
  \(0<|k|\le bq\) has two transverse half-widths \(bcq^{3/4}\) and one
  radial half-width \(bq\), hence at most
  \(8(bq+\sqrt3/2)(bcq^{3/4}+\sqrt3/2)^2\) modes; when all three
  half-widths are at least \(\sqrt3/2\), this is at most
  \(64b^3c^2q^{5/2}\) and its concentration is at most
  \(8b^{3/2}c\,q^{5/4}\); and
* C193's actual Dirichlet side length is \(4q^{1/4}\), so it has exactly
  \(64q^{3/4}\) modes and the displayed \(8q^{3/8}\) concentration.

C195's normalized-direction tube does not by itself impose a radial
\(q^{-1/4}\) restriction, so only the projective \(q^{5/4}\) ceiling
attaches automatically. Both ceilings remain below raw \(q^{3/2}\).

If \(P\) formal carrier--envelope pairs land in \(N\) output modes with
multiplicities \(m_\xi\), then

\[
          \sum_\xi m_\xi=P,\qquad
          \sum_\xi m_\xi^2\ge\frac{P^2}{N}.             \tag{2.4}
\]

For coherent equal-sign pairs the left side is their Parseval energy.
This is another proof that relabeling collisions does not create new
fixed-energy concentration.

The \(J^{-2}\) ceiling is kinematically attainable under an explicit
matching support hypothesis. If a retained reality-complete set contains
\(N\) positive modes and their negatives, and
\(e\cdot P_{k^\perp}e/|P_{k^\perp}e|\ge\kappa>0\) on the positive set,
choosing
\[
 a_k=P_{k^\perp}e/|P_{k^\perp}e|
\]
with equal reality-paired coefficients gives the exact lower bound
\(\kappa\sqrt{2N}\). Establishing \(N\) comparable to \(q^3/J^4\) for the
exact retained evolution is open. The corresponding same-energy child
volume is \(J^4q^{-3}\), with volume-equivalent isotropic scale
\(J^{4/3}/q\).

A variable-coefficient linear evolution can populate new modes. Therefore
these ceilings apply to the endpoint component actually retained or
projected in the stated \(\Sigma\); the principal fiber multiplier supplies
no out-of-band lower bound.

## 3. Conditional one-beam bridge arithmetic

This section is **not** a composition theorem. It records the powers that
would survive only if a future \(M\)-uniform Fourier-integral
almost-orthogonality estimate, spatial-annulus localization, and endpoint
band/tail theorem upgraded C194 from one beam to the full block.

Let the filter gain be \(G\ge q^g\). Since C192 gives \(\rho>e^8\), take

\[
 R_g=\left\lceil\frac g8\log q\right\rceil+1.
\]

Using \(T<76/25\),

\[
 t<\frac{19g}{50}\log q+\frac{152}{25},\qquad
 e^{6t}<e^{912/25}q^{57g/25}.                           \tag{3.1}
\]

For \(g=3/8\), a stable width \(q^{-1/12}\) and an expanding width
\(q^{-1/3}\), seeded by \(G^{-1}\), give C194 first-term margins

\[
 \mu_-=\frac{37}{600},\qquad \mu_+=\frac{14}{75},        \tag{3.2}
\]

and second-term margins

\[
 \widetilde\mu_-=\frac{29}{200},\qquad
 \widetilde\mu_+=\frac{13}{25}.                          \tag{3.3}
\]

The first terms leave the conditional interval

\[
                 \frac38\le g<\frac{275}{684},          \tag{3.4}
\]

of width \(37/1368\). For the explicit grid,
\(M\ge q^2/13824\). The normalized triangle majorant contains
\(\sqrt M\), hence one full power of \(q\), and therefore destroys every
margin in (3.2).

Even under the unproved uniform block theorem, every use would still need
the explicit finite-\(q\) condition

\[
\begin{aligned}
 e^{912/25}\big[&
 4{,}199{,}040(1+t)^3C_{\nabla}q^{-\mu_\pm}\\
 &+2{,}898{,}006{,}000{,}000{,}000(1+t)^7C_0
 q^{-\widetilde\mu_\pm}\big]\le\frac1{100},             \tag{3.5}
\end{aligned}
\]

where \(C_\nabla,C_0\) must be uniform block-profile constants proved by
that theorem; they are not supplied by the Dirichlet construction.

There is also a precise sufficient-majorant collision with C182. In the
principal fixed-band model, a raw \(bq^{3/2}\) endpoint with gain \(q^g\)
has expanding entrance scale \(bq^{3/2-g}\). Keeping it below \(bq\)
requires \(g\ge1/2\). C194's displayed broad-stable sufficient majorant
certifies decay only for

\[
                         g<\frac{25}{57}<\frac12.        \tag{3.6}
\]

These two certified windows do not overlap. This is not a theorem that
every improved bridge fails for \(g\ge25/57\). At \(g=3/8\), the expanding
seed has point scale \(bq^{9/8}\), so C182's entrance hypothesis is simply
inapplicable.

## 4. Claim boundary

C196 proves:

* exact real periodic solenoidal fixed-aperture profiles with explicit
  \(q^{9/8}\) and \(q^{3/2}\) concentration lower bounds;
* exact support-cardinality ceilings for three-coordinate, projective, and
  C180/C176 retained bands; and
* the explicitly conditional one-beam exponent ledger (3.1)--(3.6).

It does not prove that any one solution has both C194 evolution and the
C196 endpoint, and it does not prove band retention, viscosity, nonlinear
closure, UVSR, or a Navier--Stokes singularity.

The dependency-free checker verifies representative exact curl/sideband
algebra and every exponent, rational margin, grid-count inequality,
phase-space power, and clock constant used above. The general Fourier and
lattice proofs are the displayed arguments, not interval numerics hidden
inside the checker.
