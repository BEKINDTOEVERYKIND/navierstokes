# Return-cell no-go map and the surviving infinite-tail target

Date: 2026-07-29

## Status

No solution of the Navier--Stokes Millennium problem is claimed here.

This note records three kill-or-advance results for the forced periodic
breakdown program:

1. a fixed finite Fourier/helical transfer cell is impossible, even as an
   asymptotic terminal-stage template with vanishing normalized force and
   vanishing strong off-template tail;
2. known universal steady Euler flows program particle Poincaré maps, not
   Euler state-to-state return maps, and one Beltrami curl eigenspace is
   exactly dark to the projected Euler nonlinearity;
3. allowing the scale ratio to grow repairs the scalar smoothness and
   localization ledger only by demanding an unbounded-gain nonlinear return
   with a super-algebraically accurate shape reset.

The only version of the present route not eliminated by these tests has a
genuinely infinite or localized tail as part of its leading active profile
and a scale-renormalized nonlinear recurrence theorem.

## 1. Asymptotically finite active cells are impossible

For a lattice-compatible active ansatz that can be placed on one common
normalized torus, rescale carrier frequency \(N_j\), velocity amplitude
\(A_j\), and nonlinear clock \(A_jN_j\):

\[
 y=N_jx,\qquad
 \tau=A_jN_j(t-t_j),\qquad
 u(t,x)=A_jv_j(\tau,y).
\]

The forced Navier--Stokes equation becomes

\[
 \partial_\tau v_j+B(v_j,v_j)
 =\varepsilon_j\Delta v_j+g_j,\qquad
 \varepsilon_j=\frac{\nu N_j}{A_j},\qquad
 g_j(\tau,y)
 =\frac{\mathbb P f(t_j+\tau/(A_jN_j),\,y/N_j)}
        {A_j^2N_j},
 \tag{1.1}
\]

where \(B(v,v)=\mathbb P(v\cdot\nabla v)\).
An arbitrary periodic stage is not automatically put on a fixed torus by
\(y=N_jx\); applying the result requires this lattice compatibility (or
an analogous common-domain localization argument). Smoothness of the
physical Clay force also does not by itself imply \(g_j\to0\) in the
normalized norm. That is a separate active-stage hypothesis.

Let \(S\subset\mathbb Z^3\setminus\{0\}\) be finite and symmetric, fix
\(s>5/2\), and work on a compact normalized time interval \(I\). Suppose

\[
\begin{aligned}
 &\sup_j\|v_j\|_{C(I;H^s)}<\infty,\\
 &\|(I-P_S)v_j\|_{C(I;H^s)}\to0,\\
 &\|g_j\|_{C(I;H^{s-1})}\to0,\\
 &\varepsilon_j\to\varepsilon\in[0,\infty).
\end{aligned}
\tag{1.2}
\]

The finitely many coefficients in \(P_Sv_j\) are equicontinuous by
(1.1). After taking a subsequence they converge in \(C(I;H^s)\), and the
tail hypothesis upgrades this to \(v_j\to v\) in \(C(I;H^s)\). The
Sobolev algebra estimate

\[
\|B(v_j,v_j)-B(P_Sv_j,P_Sv_j)\|_{H^{s-1}}
\lesssim
(\|v_j\|_{H^s}+\|P_Sv_j\|_{H^s})
\|(I-P_S)v_j\|_{H^s}
\tag{1.3}
\]

makes the nonlinear term converge in \(C(I;H^{s-1})\), while the
Laplacian term converges in \(C(I;H^{s-2})\). Passing the integrated
equation in \(H^{s-2}\), including at frequencies outside \(S\), shows
that the limit \(v\) is a real finite-Fourier-support solution of

\[
\partial_\tau v+B(v,v)=\varepsilon\Delta v.
\tag{1.4}
\]

Kishimoto--Yoneda's classification, including its viscous extension,
then gives

\[
\boxed{B(v,v)\equiv0.}
\tag{1.5}
\]

Every limiting coefficient evolves only by heat. Since every subsequence
has a further subsequence with a nonlinear-dark limit, the full sequence
of nonlinear transfer functionals tends uniformly to zero. Hence no
fixed-fraction autonomous gate with an order-one normalized energy or
transfer denominator can satisfy (1.2). The nondegeneracy condition is
essential: a ratio alone can remain nonzero while numerator and
denominator both vanish.

If instead \(\varepsilon_j\to\infty\), the normalized energy identity
forces \(\nabla v_j\to0\) in spacetime \(L^2\). Because \(0\notin S\),
the tail hypothesis also gives \(\overline v_j\to0\); Poincaré applied to
\(v_j-\overline v_j\) then yields \(v_j\to0\) on every interior time
interval, apart from a possible vanishing initial boundary layer. That
does not yield a controlled nontrivial gate either.

The theorem is deliberately conditional on strong normalized compactness.
A prospective cascade can escape it only if at least one of the following
is essential:

- the number of leading active modes grows with the generation;
- an infinite/off-template tail stays order one on the nonlinear clock;
- the normalized active force stays order one;
- the normalized template degenerates rather than converging; or
- no uniform strong bound is available.

The third escape is incompatible with the intended role of the Clay force:
it would make the force implement the active cascade rather than protect a
super-algebraically small dormant seed.

Source:
[Kishimoto--Yoneda, *Characterization of three-dimensional Euler flows
supported on finitely many Fourier modes*](https://arxiv.org/abs/2110.08039).

## 2. A pairwise reality-sideband obstruction

Let \(u_p=a\), \(u_q=b\), \(u_{-p}=\overline a\), and
\(u_{-q}=\overline b\), with

\[
p\cdot a=q\cdot b=0,\qquad
\alpha=a\cdot q,\qquad \beta=b\cdot p.
\]

The unprojected symmetrized Euler outputs from this pair at \(p+q\) and
\(p-q\) are

\[
T_+=\alpha b+\beta a,\qquad
T_-=-\alpha\overline b+\overline\beta a.
\tag{2.1}
\]

They obey

\[
T_+\cdot(p-q)=0,\qquad T_-\cdot(p+q)=0.
\tag{2.2}
\]

If \(|p|\ne|q|\) and \(P_{p-q}T_-=0\), then \(T_-\) is parallel to
\(p-q\). Equation (2.2) and

\[
(p-q)\cdot(p+q)=|p|^2-|q|^2\ne0
\]

force \(T_-=0\). This implies
\(\alpha\overline b=\overline\beta a\). For nonzero coefficients, either
both scalars vanish, or
\(a=(\alpha/\overline\beta)\overline b\), which contradicts
\(a\cdot p=0\). Thus

\[
\boxed{
 |p|\ne|q|,\quad P_{p-q}T_-=0
 \ \Longrightarrow\ T_+=0.}
\tag{2.3}
\]

The sum/difference-exchanged statement is identical. This is pairwise:
other decompositions of the same output frequency can cancel the total
coefficient. Section 1 is what rules out a complete fixed finite polyad,
including one using such cancellations.

Equal-radius pairs are the narrow exception. They can delete one companion
while retaining an output normal to \(\operatorname{span}\{p,q\}\). This
is useful for programmable high--high-to--low stress, but it is not a
low--high recurrent amplifier.

## 3. Additive Fourier closure has only three outcomes

For

\[
\Lambda=\operatorname{span}_{\mathbb Z}\{k_1,\ldots,k_m\}
\subset\mathbb Z^3,
\]

the complete Euler and Navier--Stokes convolution preserves \(\Lambda\).
Its rank gives an exhaustive trichotomy.

- Rank one is a generalized shear and \(B(u,u)=0\).
- Rank two is an exact 2D3C system: the planar velocity is two-dimensional
  and the normal component is passive. Smooth evolution is global.
- Rank three is a finite-index sublattice of \(\mathbb Z^3\), hence the
  full three-dimensional PDE on a finite cover/quotient torus rather than
  a sparse shell closure.

For the two-chain candidate previously tested,

\[
p=(6,3,-6),\quad q_1=(2,3,6),\quad q_2=(-6,3,-2),
\]

\[
\det[p,q_1,q_2]=-384.
\]

Its exact nonlinear completion is therefore volume-filling. Adding all
sidebands does not close a reduced gate; it restores the unresolved 3D
equation.

This also changes the decision on the pending high-resolution
charge-lattice causal fork. Its `diagonal_source_off` branch lives in

\[
\mathbb Zp+\mathbb Z(q_1+q_2),
\]

which has rank two and is an invariant 2D3C system. Even a perfect
autonomous-growth pass cannot certify a self-reproducing 3D daughter: it
contains only one transverse charge. A prize-relevant handoff must
produce a new pump together with **two independent daughter charge
chains**. The existing single-diagonal \(N=384/512\) refinement therefore
has no prize-level decision value.

There is a useful positive boundary. A maintained unidirectional/Beltrami
pump and a type-I viscous Bloch chain can have a rigorous unstable
eigenvalue. Parabolic unstable-manifold theory then gives an exact local
forced-Navier--Stokes infinite-lattice amplifier. It is not recurrent:
the old pump is maintained by order-one active forcing, and no theorem
returns the nonlinear endpoint to a scaled daughter pump.

Source:
[Vasudevan, *Instability of unidirectional flows for the 2D
Navier--Stokes equations and related alpha-models*](https://arxiv.org/abs/2011.02244).

## 4. Particle universality is not state recurrence

Berger--Florio--Peralta-Salas construct stationary Beltrami fields whose
particle Poincaré maps approximate arbitrary area-preserving disk maps.
The velocity field itself remains stationary.

The distinction is exact. If

\[
\operatorname{curl}B=\lambda B,\qquad
\operatorname{curl}V=\lambda V,
\]

then the cross interaction is a gradient:

\[
\mathbb P\big((B\cdot\nabla)V+(V\cdot\nabla)B\big)=0.
\tag{4.1}
\]

One curl eigenspace is therefore a flat manifold of Euler equilibria,
despite potentially universal particle dynamics. For two different curl
eigenvalues,

\[
\mathbb P\big((B_\lambda\cdot\nabla)C_\mu+
(C_\mu\cdot\nabla)B_\lambda\big)
=(\lambda-\mu)\mathbb P(B_\lambda\times C_\mu),
\tag{4.2}
\]

but particle universality does not control this projected cross product
or its infinite off-shell output.

The usable positive interface is finite-stage only: inverse localization
can place a robust particle router into a high-curl torus Beltrami
eigenfield, whose scalar heat decay is an exact Navier--Stokes parent. A
prize route would still need a new cocycle-enhanced localization theorem
that controls the full pressure-projected child and its nonlinear return
to a non-Beltrami cone.

Sources:
[Berger--Florio--Peralta-Salas, *Steady Euler flows on
\(\mathbb R^3\) with wild and universal dynamics*](https://arxiv.org/abs/2202.02848)
and
[Enciso--Luque--Peralta-Salas, high-frequency Beltrami localization on
the three-torus](https://arxiv.org/abs/1909.07448).

## 5. Growing ratios repair bookkeeping, not recurrence

Set \(A_j=\log N_j\). One explicit scalar schedule is

\[
\varepsilon=\frac18,\quad \delta=\frac14,\quad
r_j=N_j^{1/32},\quad h_j=A_j^2,\quad
Q_j=r_j^{1/2},\quad K_j\asymp A_j.
\tag{5.1}
\]

With strain gain \(\mu_j=r_j^{2+\varepsilon}\), it gives

\[
\begin{aligned}
E_j&\asymp N_j^{-1/2},&
D_j&\asymp A_j^2N_j^{-5/8},\\
\Delta t_j&\asymp A_j^2N_j^{-17/8},&
\text{active heat error}&\asymp A_j^2N_j^{-1/16},\\
\|u_j\|_3&\asymp N_j^{5/24},&
\operatorname{Re}_j&\asymp N_j^{1/8}.
\end{aligned}
\tag{5.2}
\]

Dormant seeds and their protection force are
\(\exp[-A_j^2+O(A_j)]\), which is flat in every fixed derivative order.
Conditional on tame all-order homological estimates, choosing
\(K_j\asymp A_j\) also makes localization/WKB residuals
\(\exp(-cA_j^2)\).

The unavoidable chain is

\[
h_j\gg A_{j+1},\qquad L_j\gtrsim h_j,\qquad
L_j\ll Q_j\ll r_j.
\tag{5.3}
\]

It forces \(r_j\to\infty\), hence

\[
\mu_j=r_j^{2+\varepsilon}\to\infty.
\tag{5.4}
\]

Thus growing ratios remove the fixed-ratio localization obstruction only
by replacing a finite clone with an arbitrary-gain nonlinear return whose
shape defect is \(O(N^{-\infty})\). A linear spectral instability supplies
the seed gain but not the order-one nonlinear endpoint. Subdividing the
large jump restores the original fixed-ratio clone problem.

## 6. The surviving theorem target

The surviving object cannot be another finite triad or a fixed finite
helical polyad. It must have all of the following properties.

1. **Leading infinite/localized tail.** The sideband tail is part of the
   normalized profile, not a perturbation that vanishes in strong topology.
2. **Spectral seed gain without material frequency inflation.** A dormant
   child grows on a genuine eigenline/Floquet line while active frequencies
   remain polynomial in the carrier frequency.
3. **Nonlinear scale return.** At an order-one endpoint, the complete
   infinite-tail state, after scaling/rotation/translation/phase, lies in
   the same parent-plus-unstable-mode class.
4. **Vanishing normalized active force.** The force protects dormant tails
   but does not implement the nonlinear handoff.
5. **All-order terminal accuracy.** Drain, wakes, pressure tails, and shape
   defects are exact or super-algebraically small through the Zeno cascade.

This is stronger and more specific than “find a turbulent cascade.” It is
a renormalized heteroclinic/unstable-manifold theorem in a weighted
infinite-dimensional profile space.

The high--high-to--low half is no longer wholly speculative:
Cheskidov--Dai--Palasek construct a complete inverse cascade in a classical
Navier--Stokes flow, although their solution injects energy
instantaneously from infinite frequency and therefore does not give a Clay
finite-energy forward breakdown. The unresolved prize-level half is a
forward, scale-changing, low--high nonlinear return with the five
properties above.

Sources:
[Cheskidov--Dai--Palasek, *Instantaneous Type I blow-up and
non-uniqueness of smooth solutions of the Navier--Stokes
equations*](https://arxiv.org/abs/2511.09556) and
[Palasek, *Finite-time blow-up in an elementary model of the 3D
Navier--Stokes equations*](https://arxiv.org/abs/2605.13827).

## 7. Computational decision

A larger generic DNS is not a prize-relevant next step. A GPU experiment
becomes decisive only after a candidate leading infinite-tail profile and
an output section are specified. Its endpoint must contain a new pump and
two independent transverse daughter charges; a single diagonal child is
rank-two and cannot reproduce the 3D handoff. The run must test:

- daughter unstable projection;
- rescaled endpoint distance in a weighted tail norm;
- old-parent drain;
- off-profile energy;
- maximum active frequency inflation;
- normalized active-force residual; and
- convergence under both resolution and domain refinement.

A pass would identify a concrete theorem target. A fail would kill the
candidate profile. Until such a profile is specified, more gain
optimization measures a local amplifier already known to exist rather
than the missing recurrence.
