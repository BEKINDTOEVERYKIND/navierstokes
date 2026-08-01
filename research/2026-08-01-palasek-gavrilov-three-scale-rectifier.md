# A Palasek–Gavrilov three-scale rectifier

Date: 2026-08-01

## Result and claim boundary

This note records a promising route and two exact obstructions. It does **not** construct a Navier–Stokes singularity.

The positive result is that thin pressure-modulated Gavrilov rings reproduce Palasek's intermittency volume, the `L2`/`L∞` conversion, the amplifier coefficient, and the viscous parameter window exactly. In the straight-tube limit, the outer flank of a standard compact pressure modulation also satisfies a strict centrifugal-instability sign condition.

The negative result is that the exponentially growing field cannot simultaneously be a real invariant eigenline and a steady/self-dark clone of the next Gavrilov parent. Euler energy polarization forbids positive growth in that case. The recurrent object must therefore have at least three roles:

1. a steady low-frequency parent `P`;
2. a two- or multi-colour unstable carrier bundle `A`;
3. a distinct steady child `Z`, written by the carrier's low beat/mean stress.

At amplitude level the correct conservative target is

\[
 P'=-dA^2,\qquad
 A'=(dP-cZ-\kappa)A,\qquad
 Z'=cA^2.
\tag{0.1}
\]

This is the earlier parent--active--child gate, but with an important interpretation: `A` is an auxiliary short-wave bundle and is **not** the next parent. The load-bearing PDE problem is a multi-phase Floquet-to-stress rectifier lemma.

Everything below is `VERIFIED-SELF` unless explicitly identified as a literature theorem. It needs independent cross-audit before entering `CLAIMS.md` as established.

## 1. Exact geometric match to Palasek intermittency

Let `N` denote the inverse minor radius of a thin toroidal Gavrilov bubble. For

\[
 2<\alpha<\frac52,
\]

choose minor and major radii

\[
 r_N=N^{-1},\qquad R_N=N^{4-2\alpha}.
\tag{1.1}
\]

Then

\[
 \delta_N:=\frac{r_N}{R_N}=N^{2\alpha-5}\longrightarrow0
\tag{1.2}
\]

and the torus volume is

\[
 |\operatorname{supp}u_N|
 \asymp R_Nr_N^2
 =N^{-2(\alpha-1)}.
\tag{1.3}
\]

This is exactly the fractal-support volume postulated in Palasek's Obukhov model. If

\[
 Y_N\asymp\|u_N\|_\infty,
 \qquad X_N\asymp\|u_N\|_2,
\]

then (1.3) gives

\[
 X_N\asymp N^{-(\alpha-1)}Y_N,
 \qquad
 Y_N\asymp N^{\alpha-1}X_N.
\tag{1.4}
\]

The dominant gradient is across the minor radius, hence

\[
 \|\nabla u_N\|_\infty
 \asymp \frac{Y_N}{r_N}
 =NY_N
 =N^\alpha X_N.
\tag{1.5}
\]

Thus a parent at frequency `N_{k-1}` supplies the precise low--high growth rate

\[
 \sigma_k\asymp
 N_{k-1}Y_{k-1}
 =N_{k-1}^{\alpha}X_{k-1}.
\tag{1.6}
\]

For Palasek target amplitudes

\[
 N_k=N_{k-1}^{b},\qquad
 X_k=N_k^{\beta-\alpha},
 \qquad 1<b,\quad 2b<\beta<\alpha<\frac52,
\tag{1.7}
\]

we have

\[
 \sigma_k\asymp N_{k-1}^{\beta},
 \qquad
 \nu N_k^2\ll\sigma_k
 \quad\Longleftrightarrow\quad
 2b<\beta
\tag{1.8}
\]

(up to subpower carrier factors). The Reynolds number of the parent is

\[
 \operatorname{Re}_{k-1}
 \asymp\nu^{-1}N_{k-1}^{\beta-2}\to\infty.
\tag{1.9}
\]

So the same inequality which clears viscosity in the shell model clears it for the thin Gavrilov geometry.

## 2. Strict centrifugal sign in the thin-tube limit

Gavrilov's analytic pressure coordinate has expansion near the central circle

\[
 a_G=2\left[
 \left(\frac{\rho-R}{R}\right)^2+
 \left(\frac zR\right)^2
 \right]+O(\delta_N^3).
\tag{2.1}
\]

After multiplying the local flow by a smooth pressure cutoff `g` and passing to straight thin-tube coordinates, the leading helical column has profiles

\[
 V(s)=\kappa_0 s g(s),
 \qquad
 W(s)=\frac{\kappa_0 s g(s)}{\sqrt2},
 \qquad
 \kappa_0\asymp\frac{Y_N}{r_N}.
\tag{2.2}
\]

Here `V` is azimuthal velocity and `W` is axial velocity. Put

\[
 \Omega=\frac Vs=\kappa_0 g,
 \qquad
 \Gamma=sV=\kappa_0s^2g,
 \qquad
 h=\frac{sg'}g.
\tag{2.3}
\]

A direct substitution into the Leibovich--Stewartson discriminant gives

\[
 V\Omega'\left(\Omega'\Gamma'+(W')^2\right)
 =\frac{\kappa_0^4g^4}{2}
 h(3h^2+6h+1).
\tag{2.4}
\]

Every ordinary compact bump has an outer descending flank on which `h<-2`. There

\[
 3h^2+6h+1>0,
 \qquad
 h(3h^2+6h+1)<0.
\tag{2.5}
\]

Moreover

\[
 \Gamma'=\kappa_0sg(2+h)<0.
\tag{2.6}
\]

Thus the limiting column satisfies a strict centrifugal-instability sign condition, and the projected poloidal circulation decreases outward as in Bayly's criterion. The expected short-wave growth rate is

\[
 \sigma_{\rm WKB}\asymp c_g\kappa_0
 \asymp c_gNY_N
 =c_gN^\alpha X_N,
\tag{2.7}
\]

with `c_g>0` depending only on the chosen rescaled modulation.

What is proved by the literature and what is not:

- Bayly proves a short-wave instability for two-dimensional flows with convex closed streamlines and outward-decreasing circulation.
- Lifschitz--Hameiri reduce stability of axisymmetric vortex rings with swirl to a bicharacteristic/Floquet system.
- Albritton--Ożański rigorously construct the related high-wavenumber ring eigenmodes for vortex columns under their global column hypotheses.
- It remains to prove that the exact finite-curvature, compactly modulated Gavrilov ring retains a Floquet multiplier greater than one. Strictness of (2.5) and `delta_N -> 0` makes this a credible perturbative lemma, not an established theorem here.

## 3. Exact no-go for “unstable child = steady next parent”

Let

\[
 B(u,v):=\mathbb P((u\cdot\nabla)v),
 \qquad Q(u):=B(u,u),
\tag{3.1}
\]

where `P` is the Leray projection. Linearization about a steady parent `P` is

\[
 L_PA=B(P,A)+B(A,P).
\tag{3.2}
\]

For real divergence-free fields, integration by parts gives the exact polarization identity

\[
 \langle A,L_PA\rangle
 =-\langle P,Q(A)\rangle.
\tag{3.3}
\]

If `A` spans a real growing invariant line,

\[
 -L_PA=\lambda A,
 \qquad \lambda>0,
\tag{3.4}
\]

then

\[
 \lambda\|A\|_2^2
 =\langle P,Q(A)\rangle.
\tag{3.5}
\]

If `A` is already a steady Euler child, `Q(A)=0`, and (3.5) forces `lambda=0`. Hence a real scale-uniform unstable daughter cannot already be a steady/self-dark clone.

For a complex Floquet mode `A=C+iD`, the corresponding phase-averaged identity involves

\[
 Q(C)+Q(D).
\]

Positive growth is excluded if the entire real phase circle `C cos theta-D sin theta` is self-dark. One specially chosen steady snapshot does not by itself exclude a rotating two-dimensional bundle. This is why the conclusion is “multi-phase/nonstationary bundle required,” not “all complex modes impossible.”

A related conservation filter applies to a pure scale return of one material Gavrilov tube. If a distinguished circulation loop is returned to its homologous loop, Kelvin circulation gives `a ell = constant`, while energy gives `a^2 ell^3 = constant`. Without a wake these force no change of scale. With a wake, the same material tube is still locked to the critical scaling `a proportional ell^{-1}`. The cascade must activate distinct preloaded descendants and retain or dissipate an ancestral wake; it cannot simply shrink one bubble into its clone.

## 4. The exact conservative rectifier gate

Consider

\[
 P'=-dA^2,
 \qquad
 A'=(dP-cZ-\kappa)A,
 \qquad
 Z'=cA^2,
\tag{4.1}
\]

with `d,c>0`, carrier damping `kappa>=0`, and incoming state

\[
 (P,A,Z)\longrightarrow(P_0,0,0)
 \quad\text{as }t\to-\infty.
\tag{4.2}
\]

The exact energy identity is

\[
 \frac d{dt}(P^2+A^2+Z^2)=-2\kappa A^2.
\tag{4.3}
\]

There is also the linear invariant

\[
 Z+\frac cdP=\frac cdP_0.
\tag{4.4}
\]

Let

\[
 s=dP_0-\kappa>0,
 \qquad
 P_c=\frac{c^2P_0+d\kappa}{d^2+c^2},
 \qquad
 q_0=\frac{ds}{d^2+c^2}.
\tag{4.5}
\]

Then the exact heteroclinic is

\[
 \begin{aligned}
 P(t)&=P_c-q_0\tanh(s(t-t_0)),\\
 A(t)&=\frac{s}{\sqrt{d^2+c^2}}\operatorname{sech}(s(t-t_0)),\\
 Z(t)&=\frac cd(P_0-P(t)).
 \end{aligned}
\tag{4.6}
\]

Its outgoing state is

\[
 \begin{aligned}
 P_f&=\frac{(c^2-d^2)P_0+2d\kappa}{d^2+c^2},\\
 Z_f&=\frac{2c(dP_0-\kappa)}{d^2+c^2},\\
 A_f&=0.
 \end{aligned}
\tag{4.7}
\]

When `kappa << dP_0` and `c >> d`,

\[
 \frac{Z_f}{P_0}=\frac{2d}{c}+O\left(\frac{d^3}{c^3}+\frac\kappa{cP_0}\right),
 \qquad
 \frac{P_0-P_f}{P_0}=2\frac{d^2}{c^2}+o\left(\frac{d^2}{c^2}\right).
\tag{4.8}
\]

Thus a small child can be created while depleting only a still smaller fraction of the parent.

## 5. Mapping the gate to the superlacunary cascade

Set

\[
 P_0=X_{k-1},\qquad Z_f=X_k,
 \qquad
 r_k:=\frac{X_k}{X_{k-1}}
 =N_{k-1}^{-(b-1)(\alpha-\beta)}.
\tag{5.1}
\]

The parent-to-carrier coefficient should obey

\[
 d_k\asymp N_{k-1}^{\alpha},
 \qquad
 d_kP_0\asymp N_{k-1}^{\beta}.
\tag{5.2}
\]

Equation (4.8) asks for

\[
 \frac{c_k}{d_k}\asymp\frac{2}{r_k}
 \asymp 2N_{k-1}^{(b-1)(\alpha-\beta)}.
\tag{5.3}
\]

Equivalently,

\[
 c_k\asymp
 N_{k-1}^{\alpha+(b-1)(\alpha-\beta)}
 =N_{k-1}^{b\alpha-(b-1)\beta}.
\tag{5.4}
\]

This is below the largest child-scale coefficient `N_k^alpha=N_{k-1}^{b alpha}` by the factor

\[
 N_{k-1}^{-(b-1)\beta}.
\tag{5.5}
\]

So the requested child coupling is large relative to `d_k` but not larger than the dimensional interaction available at the child scale. It must be selected by overlap, phase, and stress geometry.

Use a two- or multi-colour carrier with frequencies

\[
 \Lambda_k\asymp K_kN_k,
\tag{5.6}
\]

where `K_k` is subpower (polynomial in the stage index is enough). Two nearby carrier phases can have difference frequency of order `N_k`, so their beat stress can write the child. The carrier damping is

\[
 \kappa_k\asymp\nu\Lambda_k^2.
\tag{5.7}
\]

The gate is amplifying precisely when

\[
 \nu K_k^2N_k^2\ll N_{k-1}^{\beta},
\tag{5.8}
\]

which again reduces to `2b<beta` for subpower `K_k`.

Truncate (4.6) at `s|t-t_0|=G_k`. The incoming carrier is `O(e^{-G_k})`, and the gate duration is

\[
 \tau_k\asymp\frac{G_k}{N_{k-1}^{\beta}}.
\tag{5.9}
\]

For a `C∞` incoming seed it is enough to take

\[
 \frac{G_k}{\log\Lambda_k}\longrightarrow\infty,
\tag{5.10}
\]

for example `G_k=(log Lambda_k)^2`. Then

\[
 \nu N_{k-1}^2\tau_k
 \asymp G_kN_{k-1}^{2-\beta}\longrightarrow0,
\tag{5.11}
\]

and the gate times are summable under the superlacunary schedule. This is substantially less demanding than forcing one linear packet to survive all the way to Palasek's terminal profile without a capture mechanism.

## 6. The missing Floquet-to-stress rectifier lemma

The remaining PDE target can now be stated concretely.

> **Multi-colour Gavrilov rectifier target.** For every sufficiently large scale ratio, construct inside a thin pressure-modulated Gavrilov parent a localized carrier bundle with at least two transported phases such that:
>
> 1. all colours lie in a common exponentially unstable Floquet cone with rate `d_k P_0`;
> 2. their phase-averaged parent work and low stress realize the reciprocal coefficient `d_k`;
> 3. their transported beat covariance spans a prescribed compact child-writing stress with effective coefficient `c_k` from (5.3);
> 4. nonlinear feedback rotates/captures the bundle onto the heteroclinic (4.6), leaving a steady Gavrilov child and an admissible wake rather than a growing carrier;
> 5. at least two independent phase gradients realize positive-definite covariance; incompressibility makes one phase rank at most two;
> 6. viscosity, transport, off-shell products, compact pressure tails, and local momentum constraints are smaller in a weighted norm compatible with the infinite cascade.

A common unstable cone is not enough: the quadratic covariance of that cone must meet the child stress-synthesis cone. Conversely, an arbitrary high--high-to-low stress construction is not enough: its phases must be reachable as the output of the unstable Floquet propagator. The decisive condition is an endpoint submersion from incoming unstable phase/polarization data to outgoing child covariance.

The direct one-doublet Fourier realization already failed because constant-coefficient triad geometry cannot make the same active line parent-unstable and child-stable. The present target escapes only by using a time-dependent Floquet bundle and at least two colours. That escape must be proved, not assumed.

## 7. What is worth computing

A GPU search is not yet the bottleneck. The next finite computation should be a narrow certificate, after the exact Gavrilov modulation is fixed:

1. integrate the exact Lifschitz--Hameiri bicharacteristic/Floquet system around the finite-curvature ring;
2. certify a multiplier `>1` uniformly for small but nonzero `delta`;
3. propagate two independent unstable phases;
4. compute the rank and condition number of the map from their quadratic beat covariances to the finite child-stress dictionary.

A positive result would validate items 1--3 of the rectifier target. A negative result would kill this particular modulation without spending effort on the all-order corrector.

## Primary sources

- S. Palasek, *Finite-time blow-up in an elementary model of the 3D Navier--Stokes equations*: https://arxiv.org/abs/2605.13827
- A. V. Gavrilov, *A steady Euler flow with compact support*: https://arxiv.org/abs/1810.08020
- B. J. Bayly, *Three-dimensional centrifugal-type instabilities in inviscid two-dimensional flows*: https://doi.org/10.1063/1.867002
- A. Lifschitz and E. Hameiri, *Localized instabilities of vortex rings with swirl*: https://doi.org/10.1002/cpa.3160461005
- D. Albritton and W. S. Ożański, *Linear and nonlinear instability of vortex columns*: https://arxiv.org/abs/2310.20674
- A. Cheskidov, M. Dai, and S. Palasek, *Instantaneous Type I blow-up and non-uniqueness of smooth solutions of the Navier--Stokes equations*: https://arxiv.org/abs/2511.09556
