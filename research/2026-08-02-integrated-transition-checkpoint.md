# Integrated transition checkpoint after branch reconciliation

**Date:** 2026-08-02

**Status:** theorem-oriented research checkpoint; exact algebra and scale
claims are separated from analytic hypotheses.  This is not a solution of
the Navier--Stokes Millennium problem.

## 1. Outcome

The continuation branch and the later quasimode/wake work now fit into one
coherent architecture.  Several former blockers are closed or corrected:

1. an explicit Batchelor vortex column is matched to the full BAS edge;
2. the AO edge deficit is compatible with a \(j^2\) gain window at the
   endpoint \(A=4\), conditional on a uniform PDE propagator and a
   first-order curved residual; the separate Gevrey-2 ledger requires the
   strict integrated choice \(A>4\);
3. the pointwise positive-stress algebra needs two variable directions, not
   three, while gradient integrability and material transport remain open;
4. the affine viscous capture ODE is robust on a fixed incoming section;
5. the false envelope-pressure identity is replaced by an exact remainder
   formula and a finite-energy projected-tube candidate; and
6. the second thin-torus exterior jet is aspect-uniform, while the third jet
   displays the first exact \(\varepsilon^{-2}\) moment; and
7. a compact pressure modulation in the actual locked-pitch Gavrilov family
   now has a unique full BAS edge, repairing the C67 base mismatch at
   principal-symbol level.

The unresolved theorem is no longer an unspecified finite-frequency
rectifier.  It is an all-order localized Leray/heat/pressure parametrix,
coupled to a finite-curvature edge-matched AO carrier.

The certified Batchelor column is not the straight-tube limit of the
standard Gavrilov bubble.  Section 8 records the subsequent repair: redesign
the common Gavrilov pressure modulation and certify its edge directly.  The
remaining pre-curvature question is the independently unaudited q-free
transfer of the AO eigenmode proof.

## 2. Corrections forced by audit

| Former statement | Exact correction |
|---|---|
| Three phases are minimal near isotropy. | False as pointwise covariance algebra when kernel directions vary.  The two-direction sum is a submersion at an explicit positive decomposition of \(qI\).  The selected direction fields need not be curl-free or transported; three fixed directions remain an explicit sufficient chart. |
| The scalar envelope-dark equation makes the common high sum pressure. | False.  A transverse product envelope has a nonzero Leray charge even when that scalar equation holds. |
| The four primary modes prove every simultaneous finite-colour cell impossible. | Too broad.  They exclude only that support; richer shared-output correctors are not ruled out. |
| The finite real BAS doublet is an exact rectifier. | False.  Reality forces an infinite outward sideband ladder; every nonzero finite extreme leaks. |
| The first explicit AO profile is at the full BAS edge. | False by a strict \(>1.0003\) transverse-ratio counterexample.  A nearby joint optimizer repairs it exactly. |
| The affine endpoint has uniform upper bound \(0.1981\) for all \(r\ge10\). | The uniform upper bound is \(0.2000\); \(0.1981\) applies only near \(r=10\).  The terminal damping margin is unchanged. |
| The second thin-torus wake might lose an aspect power. | The extracted exterior field is aspect-uniform.  The raw local pressure is not and must be recombined with the local nonlinear jet. |

## 3. Spectral carrier module

For a Batchelor column

\[
 V(r)={Q_*\over r}(1-e^{-r^2}),\qquad W(r)=e^{-r^2},
\]

let \(x_*\) be the unique positive zero of

\[
 (2x+1)e^{2x}-(7x+2)e^x+(5x+1)=0.
\]

Define

\[
 h(x)={e^x-1-x\over x^2},\qquad
 g(x)={2-e^x\over (x+1)e^x-(2x+1)},
\]

\[
 \beta_*=\sqrt{g(x_*)},\qquad Q_*={\beta_*\over h(x_*)}.
\]

Exact interval arithmetic gives

\[
\begin{aligned}
0.59671214&<x_*<0.59671216,\\
0.5101679&<\beta_*<0.5101682,\\
0.8278572&<Q_*<0.8278581,\\
0.4571745&<\sqrt{b_*}<0.4571753.
\end{aligned}
\]

The full column BAS admits the exact reduction

\[
 (K^2x)'=A_0y,\qquad y'=Cx.
\]

Every non-resonant Lyapunov exponent is zero.  At resonance,

\[
 \lambda^2=b_\beta(r){K_h^2\over \ell^2+K_h^2}.
\]

The constructed \((r_*,\beta_*)\) is the unique global maximum over all
resonant radii and ratios, so every BAS exponent is at most
\(\sqrt{b_*}\), with equality at the AO ring.

For a normalized quasimode with residual \(M^{-q}\), edge deficit
\(M^{-r}\), gain \(G_j\asymp j^g\), carrier \(M_j=j^A\), and propagator
prefactor \((1+t)^d\), the relative finite-action error is

\[
 E_j\lesssim
 j^{-Aq+g(d+1)}\exp(Cj^{g-Ar}).
\]

Hence \(Ar\ge g\) and \(Aq>g(d+1)\) suffice.  Equality in the first
condition leaves only a bounded edge exponential.  AO has \(r=1/2\);
with conditional curved residual \(q=1\), \(g=2\), and \(d=0\), this gives

\[
                              A\ge4.
\]

At \(A=4\), the relative error is \(O(j^{-2})\), hence summable.  The
independent Gevrey-2 carrier gate is strict and therefore retains \(A>4\)
for the full construction.

The original compatibility obstruction is still exact.  For the standard
thin Gavrilov seed,

\[
 V_G(s)=\kappa s g(s),\qquad W_G(s)=V_G(s)/\sqrt2.
\]

At the certified Batchelor AO radius \(r_*\), exact interval arithmetic gives
\(V_B'(r_*)>0>W_B'(r_*)\) and
\(-W_B'(r_*)/V_B'(r_*)\in(2.95,2.96)\).  Thus even axial reflection and a
Galilean shift cannot identify the first jets.  Section 8 replaces the
Batchelor carrier by a compatible compact locked-pitch profile with its own
unique full edge.  Its remaining obligations are the q-free AO proof audit,
a frequency-uniform PDE-sector propagator at that exact edge, a real-mode
lower bound, and the normalized finite-curvature residual.  At
principal-symbol level the prefactor issue is closed for both profiles: the
physical BAS propagator satisfies a uniform edge exponential with no
polynomial prefactor.  The remaining lift is genuinely a PDE/lower-order
problem.

## 4. Local transition module

The corrected pointwise two-direction stress chart begins at

\[
 A_0=\operatorname{diag}(0,q,a),\qquad
 B_0=\operatorname{diag}(q,0,q-a),\qquad 0<a<q.
\]

Five fixed-block variations give all symmetric coordinates except the
\(12\) entry; rotating the kernel of \(A_0\) toward \(e_2\) supplies that
sixth coordinate.  The implicit-function theorem therefore gives a smooth
positive pointwise covariance chart near \(qI\).  For a spatial stress field,
its selected kernel directions need not be curl-free gradients and are not
yet realized as transported material phases.

The affine carrier/parent/child ODE has a captured vanishing-seed orbit.
For \(\alpha=0.9\), \(r\ge10\),

\[
 0.1489<p_\infty<0.2000,\qquad
 {Z_\infty\over P_0}\asymp r^{-1},\qquad
 1-{P_\infty\over P_0}\asymp r^{-2}.
\]

Its robust form starts at a fixed incoming section and requires a weighted
affine-defect norm.  An exponentially small incoming seed cannot tolerate
arbitrary additive error; it needs multiplicative or backward-weighted
forcing control.

For one symmetric two-wave colour with exactly divergence-free real scalar
envelopes, put \(P=AB\).  After subtracting the nominal pressure, the exact
high-sum source is

\[
 R=2(W^2+\delta^2)P_e e+2\delta^2P_vv.
\]

Thus exact pressure requires \(P_e=P_v=0\); no nonzero compact lateral
product is exactly dark.

The leading finite-energy replacement is the compact projected tube

\[
 u_\pm=\mathbb P\!\left[
 g((e\cdot x)/L_e)
 \rho((t\cdot x)/R_t,(v\cdot x)/R_v)
 a_\pm e^{ik_\pm\cdot x}\right].
\]

Its exact per-envelope-mode solenoidal correction has size

\[
 {\delta|\xi|\over|k+\xi t|}
 =O\!\left({Q\over\Lambda^2R}\right).
\]

The first central symbol errors are relative

\[
 O((\Lambda R_t)^{-1})
 +O\!\left({Q\over\Lambda^2R_v}\right),
\]

and the remote axial endcaps cost \(O((QL_e)^{-1})\).  The missing theorem
is an all-order weighted pseudodifferential, pressure, and heat estimate
showing that the global Leray tails do not erase BAS gain or the preseeded
child.

## 5. Global wake module

For a thin torus of major radius \(R\), aspect \(\varepsilon\), actual
speed \(v\), and duration \(\tau=TR/v\), the second-jet source obeys

\[
 \left\|\nu(U\otimes\Delta U+\Delta U\otimes U)\right\|_1
 \lesssim \nu v^2R
\]

uniformly in \(\varepsilon\).  The extracted solenoidal exterior field
satisfies

\[
 \|\nabla^m(\tau^2Z_2)\|_2
 \lesssim_m \nu T^2R^{1/2-m}.
\]

For \(K^3\) packed microtori, the inter-stage absolute cost is \(K^2\);
same-stage nearest neighbors cost \(K^3\).  Three orthogonal orientations
cancel the zeroth trace-free moment.  Three centrally inverted pairs also
cancel first moments and reduce the clustered outer loss to \(K^0\).

At the third jet,

\[
 \int{\mathcal A}_2\,dx
 =4\int\Delta U\otimes\Delta U\,dx,
\]

so a generic \(\varepsilon^{-2}\) trace-free quadrupole is unavoidable.
After normalization, conditional on the global one-viscosity tubular
projection lemma for the first term,

\[
 \|\nabla^m(\tau^3Z_3)\|_2
 \lesssim_m
 \nu T^3R^{1/2-m}
 +\nu^2T^3v^{-1}R^{-1/2-m}\varepsilon^{-2}.
\]

Writing

\[
 \mu={\nu\over vR},\qquad
 \Theta={T\mu\over\varepsilon^2},
\]

the relative exterior size is

\[
 T^3\mu+T^3\mu^2\varepsilon^{-2}
 =T^2\varepsilon^2\Theta+T\varepsilon^2\Theta^2.
\]

Thus the correct all-order expansion parameter is the minor-scale
\(\Theta\), not the macro Reynolds number alone.  If
\(\varepsilon_j=\ell_j^\beta\), \(v_j=\ell_j^{-\gamma}K_j^\gamma\),
and \(K_j,T_j\) are polynomial, the compatible high-Reynolds gate is

\[
 1<\gamma<\frac32,\qquad
 0<\beta<{\gamma-1\over2}.
\]

## 6. Integrated theorem target

A prize-level construction along this route would follow from one linked
parametrix theorem with these uniform properties:

1. independently validate the q-free AO eigenmode proof transfer for the
   compact edge-matched Gavrilov profile, or construct its local Weber
   quasimode directly;
2. bend that compatible mode into a thin ring with normalized residual
   \(O(M^{-1})\);
3. lift the uniform BAS envelope to a frequency-uniform Fourier-block PDE
   propagator over the \(G_j\asymp j^2\) gain window;
4. evolve the compact projected tube, or a proven superior multiphase
   replacement, through the robust affine capture section;
5. retain every lateral collar, axial endcap, mixed sideband, periodic
   image, pressure wake, and same-stage neighbor field;
6. close a local/exterior Gevrey induction in powers of \(\Theta_j\)
   through \(M_j\asymp j^2/\log j\); and
7. obtain residual and seam errors \(e^{-cj^2}\), making the terminal force
   \(C^\infty\)-flat.

No item in the present repository proves this linked theorem.  The exact
results above reduce its degrees of freedom and prevent the next proof
attempt from relying on a finite sideband closure, a false pressure
cancellation, or an aspect-uniformity claim beyond the second jet.

## 7. Reproducibility

All dependency-free scripts in the checks directory pass.  Their output
certifies algebra, rational interval comparisons, and exponent ledgers only;
it does not certify the missing analytic parametrix.

## 8. Carrier addendum: C67 repaired inside the Gavrilov family

The base-profile gate identified in Sections 1 and 3 now has a
principal-symbol repair.  For the actual straight Gavrilov relation

\[
                 V=r\Omega,\qquad W=V/\sqrt2,
\]

write `h=r Omega'/Omega`.  Its exact full resonant BAS rate is

\[
 \lambda^2=\Omega^2
 \left[-2h{3h^2+6h+1\over3h^2+2h+1}\right]
 {K_h^2\over\ell^2+K_h^2}.
\]

The compactified hollow log-normal profile centered at

\[
 h_*={-5-\sqrt{22}\over3},\qquad
 \dot h_*={8h_*^2\over3h_*+1}
\]

has a unique unrestricted edge

\[
                 \lambda_*^2=-h_*\Omega_0^2,
\]

a strict local fixed-sector maximum, and an isolated helical-phase minimum.  It is
obtained by choosing the common pressure modulation in the standard
Gavrilov seed, so no Batchelor-to-Gavrilov profile matching is needed.

Literal AO Assumption A cannot hold for a compact locked-pitch bump because
`q_AO=-(rV)'/W'` has a pole at a nonzero interior extremum of `V`.  The pole
is artificial: the physical Rayleigh coefficients are exactly

\[
 a=r\partial_r\left[{\beta r^2W'-\Gamma'\over
 r(1+\beta^2r^2)}\right],\qquad
 b=-{2\beta V(W'+\beta\Gamma')\over1+\beta^2r^2}.
\]

They are smooth and vanish near the axis for the hollow profile.  A
line-by-line audit of AO v3 finds no use of the quotient after these
coefficients are introduced; the gluing proof uses their regularity, the
local trapping data, and isolation of the `Lambda` level.  The resulting
q-free proof transfer is a high-priority self-derived claim, not yet an
independently audited theorem.

The revised spectral chain is now:

1. independently audit the q-free AO proof transfer, or construct its local
   Weber packet as a quasimode;
2. prove a frequency-uniform Fourier-block semigroup bound for this compact
   locked-pitch profile;
3. bend the compatible real packet into the thin Gavrilov torus with a
   normalized first-order curvature residual; and
4. feed it into the unchanged localized transition and all-order wake
   modules.

The full derivation and source audit are in
`2026-08-02-locked-pitch-gavrilov-edge-profile.md`; the independent profile
flexibility/literature gate is in
`2026-08-02-compact-ring-profile-flexibility-audit.md`.
