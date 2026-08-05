# Long-window Kelvin conditioning: passive no-go and terminal relaunch gate

**Date:** 2026-08-03

**Status:** exact passive-transport no-go, exact diagonal control ledger, and
an exact pointwise finite-dimensional terminal-pulse parameterization.  The
latter is not a Gevrey endpoint theorem or an unforced Navier--Stokes launch
theorem: a genuinely new terminal bath requires a source, reservoir, or
nonlinear transfer.

**Scope:** the \(e^{8G}\) conditioning qualification in C96.  This note does
not edit the claim register and does not construct the localized source which
would activate a new oscillatory bath.

## 1. Outcome

There are two different operations which should not be called the same
``reparameterization.''

1. Expressing a bath which has been passively transported through the whole
   gain window in terminal coordinates removes the small determinant from the
   displayed Jacobian, but transfers the same exponential loss into the
   initial phase and covariance data.  It is only a coordinate change.
2. Specifying a fresh bath at a time \(\tau\) within a bounded interval of the
   target time \(T\), and transporting it only on \([\tau,T]\), gives a
   uniformly conditioned pointwise Kelvin chart.  Fixed torus cohomology is
   preserved exactly.  However, the homogeneous Kelvin equation cannot create
   this bath from zero.  The relaunch is useful only after a nonlinear
   source/reservoir theorem is supplied.  Uniform \(C^m\) or Gevrey control
   additionally requires corresponding bounds for the full flow path and its
   derivatives on a phase neighborhood.

The distinction is quantitative.  If the phase wave number gains a factor

\[
 s={|k(T)|\over |k(0)|},
\]

then every passively transported positive rank-two covariance with a terminal
transverse gap at least \(q>0\) has initial covariance energy at least

\[
                    \operatorname{tr}_{k(0)^\perp}Q(0)
                    \ge 2qs.                                \tag{1.1}
\]

Thus \(s=e^G\) already forces an \(e^G\) physical cost, independently of
the choice of control coordinates.  The explicit C89 diagonal strain has
still sharper directional losses: a current \(e_3\)-kernel angle pulls back
by \(e^{(2\alpha+\beta)G}\), and a uniform current \(e_3\)-covariance gap
pulls back by \(e^{2(\alpha+\beta)G}\).

A bounded terminal pulse avoids these factors because its remaining
deformation has bounded condition number.  It does **not** supply the pulse's
oscillatory energy or its activation source.

## 2. An invariant passive-energy no-go

Let \(F\in SL(3,\mathbb R)\), let

\[
 k=F^{-T}k_0,
\]

and let \({\cal L}_F(k_0):k_0^\perp\to k^\perp\) be the exact Kelvin velocity
map from C96.  Its transverse determinant is

\[
 \left|\det_\perp {\cal L}_F(k_0)\right|
 ={ |k_0|\over |k|}=s^{-1}.                                 \tag{2.1}
\]

For a positive covariance \(Q_0\) on \(k_0^\perp\), passive inviscid
transport gives

\[
 Q_T={\cal L}_F Q_0 {\cal L}_F^T.                           \tag{2.2}
\]

Taking two-dimensional determinants and using (2.1),

\[
             \det_\perp Q_T=s^{-2}\det_\perp Q_0.          \tag{2.3}
\]

If \(Q_T\ge qI\) on \(k^\perp\), then
\(\det_\perp Q_T\ge q^2\).  For a positive two-by-two matrix,
\(\operatorname{tr}Q_0\ge2\sqrt{\det Q_0}\).  Equations
(2.3) therefore give (1.1).

Equivalently, if the available launch covariance obeys
\(\operatorname{tr}Q_0\le E\), then

\[
          \lambda_{\min}(Q_T|_{k^\perp})\le {E\over2s}.     \tag{2.4}
\]

The local pressure-gauge chart can remain algebraically full rank for every
finite \(G\), as C96 proves, while its uniformly positive target
neighborhood collapses at least like \(s^{-1}\) under a bounded launch-energy
budget.  Positive viscosity only decreases the two labelled amplitudes and
cannot improve (1.1).

This no-go is independent of whether \(F\) is diagonal, normal, or generated
by a time-dependent strain.  A bounded terminal deformation changes \(s\)
by only a bounded factor, so it cannot undo a wave-number gain accumulated
before the terminal interval.

## 3. Exact diagonal control ledger

For the C89 family write

\[
 S=\operatorname{diag}(-\alpha,-\beta,\alpha+\beta),
 \qquad \alpha,\beta>0,
\]

and evolve for time \(G\).  Then

\[
 F_G=\operatorname{diag}
       (e^{-\alpha G},e^{-\beta G},e^{(\alpha+\beta)G}).     \tag{3.1}
\]

At \(k_0=e_1\), the Kelvin map on the transverse velocity plane is

\[
 {\cal L}_{F_G}|_{e_1^\perp}
 =\operatorname{diag}
      (e^{\beta G},e^{-(\alpha+\beta)G})                   \tag{3.2}
\]

in the \((e_2,e_3)\) basis.  Its determinant is
\(e^{-\alpha G}\), exactly the reciprocal wave-number gain.

Let the launch phase direction be represented by

\[
 k_0=e_1+x e_2+y e_3.
\]

The terminal projective kernel coordinates have first derivatives

\[
 \theta_2=e^{(\beta-\alpha)G}x,
 \qquad
 \theta_3=e^{-(2\alpha+\beta)G}y.                           \tag{3.3}
\]

For launch covariance block
\(B_0=\left(\begin{smallmatrix}A&C\\C&B\end{smallmatrix}\right)\),
the terminal block is

\[
 B_T=
 \begin{pmatrix}
 e^{2\beta G}A&e^{-\alpha G}C\\
 e^{-\alpha G}C&e^{-2(\alpha+\beta)G}B
 \end{pmatrix}.                                             \tag{3.4}
\]

At a diagonal launch base \(\operatorname{diag}(a,b)\), the six stress
columns, apart from signs and the pressure column, have exponential factors

\[
\begin{array}{c|ccccc}
\text{control}&A&B&C&x&y\\ \hline
\text{factor}
 &e^{2\beta G}
 &e^{-2(\alpha+\beta)G}
 &e^{-\alpha G}
 &a e^{(3\beta-\alpha)G}
 &b e^{-(4\alpha+3\beta)G}.
\end{array}                                                  \tag{3.5}
\]

Their product is \(ab e^{-8\alpha G}\), reproducing C96.

Now use terminal variables \((A_T,B_T,C_T,\theta_2,\theta_3)\).  The stress
Jacobian in those variables is the static boundary-chart Jacobian, but the
physical launch data are

\[
\begin{aligned}
 A&=e^{-2\beta G}A_T,&
 B&=e^{2(\alpha+\beta)G}B_T,&
 C&=e^{\alpha G}C_T,\\
 x&=e^{(\alpha-\beta)G}\theta_2,&
 y&=e^{(2\alpha+\beta)G}\theta_3.                           \tag{3.6}
\end{aligned}
\]

Therefore current-time coordinates do not make the forward control map
tame.  For the fixed-cohomology core family
\(\Phi=x_1+y\psi_3\), any fixed launch neighborhood \(|y|\le y_*\) reaches
only

\[
 |\theta_3|\le y_*e^{-(2\alpha+\beta)G}.                    \tag{3.7}
\]

Likewise a fixed positive terminal \(B_T\) needs launch covariance
\(B=e^{2(\alpha+\beta)G}B_T\).  Re-labeling (3.6) as a new control norm
hides, rather than removes, these physical costs.

For the numerical C89 choice \(\alpha=1\), \(\beta=5/4\), the relevant
inverse exponents are

\[
 2\alpha+\beta={13\over4},\qquad
 2(\alpha+\beta)={9\over2},\qquad
 4\alpha+3\beta={31\over4}.                                \tag{3.8}
\]

## 4. Exact bounded-terminal parameterization

There is nevertheless a clean finite-dimensional terminal-pulse theorem.
Let \(X=X(T,\tau):\mathbb T^3\to\mathbb T^3\) be the volume-preserving flow
map on \([\tau,T]\), put \(H(y)=D X(y)\) in the active core, and assume

\[
 \|H\|+\|H^{-1}\|\le C_H.                                  \tag{4.1}
\]

Choose a desired terminal circle-valued phase \(\Phi_T\) with

\[
 [d\Phi_T/(2\pi)]=[dx_1/(2\pi)].                            \tag{4.2}
\]

Define its launch value by pullback,

\[
 \Phi_\tau=\Phi_T\circ X(T,\tau).                           \tag{4.3}
\]

Transporting (4.3) forward gives exactly \(\Phi_T\) at time \(T\).  A flow
map is isotopic to the identity, so

\[
 [d\Phi_\tau/(2\pi)]=[d\Phi_T/(2\pi)]                      \tag{4.4}
\]

and fixed integral cohomology is preserved.  Bound (4.1) makes the pointwise
phase gradients, and hence the finite-dimensional \(C^1\) phase coordinates,
comparable.  It does **not** by itself compare \(C^m\) or Gevrey norms.  Those
require uniform derivative bounds for \(X,X^{-1}\) through the requested
order (or one Gevrey scale), as well as a neighborhood on which \(d\Phi\)
does not vanish.

Given a positive rank-two terminal covariance field \(Q_T\) with kernel
\(d\Phi_T\), its exact pointwise pullback is

\[
 Q_\tau(y)={\cal T}_{H(y)}^{-1}
             \bigl(Q_T(X(y))\bigr).                          \tag{4.5}
\]

If the two harmonic damping factors also obey

\[
                         d_1,d_2\ge d_*>0,                  \tag{4.6}
\]

then the viscous fibre map in C96 has a smooth inverse with constants
depending pointwise on \(C_H,d_*\), and the terminal covariance gaps.  Hence
the phase, covariance, and pressure-gauge chart over this interval has a
uniform finite-dimensional pointwise inverse.  A uniform function-space
inverse additionally needs bounds along the deformation path and on the
derivatives, with respect to phase direction and space, of the damping
integral and a local transverse frame.  Those estimates are not supplied by
(4.1) and (4.6) alone.

For a strain bound \(\|S\|\le M\) and a terminal duration
\(T-\tau\le\delta_0\), (4.1) follows with
\(C_H\le2e^{M\delta_0}\).  In the inviscid affine determinant coordinates
of C96,

\[
 |\det J_{\tau,T}|
 ={ab\over |H^{-T}n|^8}
 \ge ab e^{-8M\delta_0}.                                    \tag{4.7}
\]

With the labelled viscous harmonics, C96 multiplies this determinant by
\((\delta_1\delta_2)^2c_D\), where \(\delta_j=d_j^2\) and
\(c_D\) is their positive weighted mean.  Under (4.6) the corresponding
crude lower bound is therefore

\[
 |\det J_{\tau,T,\nu}|
 \ge ab e^{-8M\delta_0}d_*^{10}.                             \tag{4.8}
\]

These pointwise constants are independent of the earlier gain length \(G\).
The heat condition (4.6), however, still requires the existing finite-band
small-heat ledger; bounded affine time alone does not control
\(\nu K^2(T-\tau)\).

## 5. Why relaunch is not a free operation

For each labelled harmonic the Kelvin amplitude equation, with or without
viscosity, is homogeneous.  If its amplitude is zero at \(\tau\), it stays
zero.  Multiplying a homogeneous Kelvin solution \(a(t)\) by a nonconstant
activation factor \(\chi(t)\) produces the exact extra term

\[
                         \chi'(t)a(t).                       \tag{5.1}
\]

Thus a bath cannot be turned on by a cutoff inside an exact unforced
solution.  A terminal relaunch requires at least one of:

1. a nonzero carrier/reservoir already present in the state and a proved
   nonlinear transfer into the two labelled harmonics;
2. a source which is later absorbed by the complete Reynolds-stress/wake
   construction; or
3. a separate local controllability theorem for the exact nonlinear
   transition.

In particular, a bounded terminal pulse is a valid way to localize the C96
chart **after** such a launch theorem exists.  It is not itself that theorem.
If the same nonzero bath must remain passively material through the whole
gain window, (1.1)--(3.7) apply and the exponential cost is unavoidable.

The wave-number area identity also shows why bounded flow feedback is not
enough.  A terminal deformation of bounded condition number can rotate and
redistribute two existing polarizations, but it changes
\(|\det_\perp{\cal L}|\) by only a bounded factor.  It cannot create a
uniform rank-two covariance from a bounded-energy bath after an unbounded
wave-number gain.

## 6. What control-dependent flow can and cannot repair

The difficult terminal phase angle can be generated by a bounded terminal
shear if the coarse flow itself is admitted as a control.  Linearize the
covector equation about the diagonal strain and divide by the base
\(k_1(T)=e^{\alpha T}\).  For perturbations \(\delta S_{12}\) and
\(\delta S_{13}\),

\[
\begin{aligned}
 \theta_2(T)
 &=e^{(\beta-\alpha)T}x_0
   -\int_0^T e^{(\beta-\alpha)(T-s)}\delta S_{12}(s)\,ds,\\
 \theta_3(T)
 &=e^{-(2\alpha+\beta)T}y_0
   -\int_0^T e^{-(2\alpha+\beta)(T-s)}\delta S_{13}(s)\,ds.
                                                               \tag{6.1}
\end{aligned}
\]

A bounded \(\delta S_{13}\) supported in \([T-\delta,T]\) therefore makes
an \(O(1)\) change in \(\theta_3(T)\), independently of \(T\).  Feedback
ending \(\delta\) time units earlier is suppressed by
\(e^{-(2\alpha+\beta)\delta}\).  This identifies a genuine terminal-control
channel for the phase direction.

It does not repair the passive covariance-area deficit (2.3).  Zero
harmonic amplitude remains zero under every homogeneous control-dependent
Kelvin propagator, and an accumulated gain \(s\) still imposes (1.1) for
each realized flow.  To obtain a uniformly positive second polarization,
the feedback must be accompanied by amplitude injection/nonlinear transfer,
not just by a change of \(F\).

Moreover, the affine shear \(\delta Sx\) is not a periodic finite-energy
velocity.  A localized realization creates the same divergence, pressure,
and collar wake obligations already isolated elsewhere in the repository.
Equation (6.1) is a control-kernel calculation, not a localized
Navier--Stokes endpoint map.

## 7. Endpoint and cascade consequences

The terminal pulse is compatible with the algebraic endpoint ledger only
in its factored form.  If physical amplitudes must have zero derivatives
through order \(M\) at a seam, take

\[
                         \chi(s)=s^{M+1}.
\]

Then the covariance carries the factor

\[
                         \chi(s)^2=s^{2M+2}.                 \tag{7.1}
\]

The normalized covariance chart may stay uniformly conditioned on the
active side, but the absolute chart necessarily degenerates at the seam.
Consequently the target residual/stress must have the doubled factor
(7.1), exactly as in C91.  A terminal pulse does not permit generic
independent stress data at a clean endpoint.

For the long cascade window \(G_j\asymp j^2\), passive transport with
\(s_j=e^{cG_j}\) requires launch covariance energy at least
\(2q e^{cG_j}\).  This grows faster than every fixed power of the repository's
geometric-in-\(j\) length/amplitude scales and polynomial carrier scales.
The passive long-window chart is therefore incompatible with a uniform
stage-normalized energy budget.  The logical escape routes are precise:

* let the target covariance gap shrink like \(s_j^{-1}\) and redo the full
  endpoint/stress scale ledger;
* retain a reservoir whose energy pays (1.1);
* prove a nonlinear terminal injection and use the bounded-pulse theorem;
  or
* separate the long AO gain carrier from the short one-phase stress bath.

The last option may introduce a second physical packet during overlap.  Its
charge interactions must be included rather than treated as a free reset.

## 8. Claim boundary

This note resolves the interpretation of the \(e^{8G}\) factor but not the
transition theorem.

* **Proved:** terminal coordinates alone cannot remove the physical
  long-window loss; (1.1) is a coordinate-invariant passive no-go.
* **Proved:** a bath genuinely launched within a bounded terminal interval
  has a uniformly conditioned pointwise finite-dimensional
  phase/covariance chart and preserves fixed cohomology exactly, subject to
  a positive finite-band damping floor.
* **Proved:** bounded terminal flow feedback can rotate the material phase,
  but passive flow feedback cannot repair the rank-two covariance-energy
  deficit.
* **Open:** prove uniform \(C^m\)/Gevrey path and damping estimates, and
  construct the localized nonlinear source/reservoir which launches the
  terminal bath, satisfies the doubled endpoint factor, and carries all
  charged and zero-charge wake terms.

The companion checker records the diagonal exponents, the exact
wave-number/covariance determinant identity, a rational current-to-launch
pullback, and the terminal flow-feedback kernels.
