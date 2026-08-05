# Endpoint-jet interpolation: what the three-phase chart does and does not solve

**Date:** 2026-08-03

**Status:** exact finite-jet compatibility and linearization; constructive
periodic Navier--Stokes--Reynolds splice; endpoint-rank no-go for vanishing
quadratic amplitudes.  Exact oscillatory realization remains open.

**Scope:** the endpoint splice left by C85--C88.  This note distinguishes a
Gevrey-2 Reynolds-stress interpolation theorem from the much stronger exact
three-phase transition theorem.  It does not construct a Navier--Stokes
singularity or an unforced oscillatory return map.

## 1. Outcome

There is no remaining finite-jet obstruction at the **coarse periodic
Navier--Stokes--Reynolds level**.  A prescribed outgoing collar jet through
order \(M\) can be joined to an incoming compatible jet by a divergence-free
Gevrey-2 path.  Its projected defect has zero spatial mean, the periodic
symmetric anti-divergence of C86 routes it exactly, and an isotropic pressure
gauge makes the resulting covariance positive.  On a short material
interval, the three transported coordinate phases factor that positive
covariance smoothly into two polarizations per phase.

The exact linear compatibility conditions are:

1. the two endpoint values have the same spatial mean and every positive
   time derivative has zero mean. This is the complete condition when the
   Reynolds stress is free; and
2. if that stress is required to vanish through order \(M-1\), the velocity
   jet must obey the Navier--Stokes recurrence through order \(M-1\).

The exact collar of C87 satisfies both conditions automatically.

There is one additional nonlinear endpoint price if the *physical*
oscillatory amplitudes, rather than only their stress, must vanish through
order \(M\).  A quadratic covariance of amplitudes which are
\(O(s^{M+1})\) is \(O(s^{2M+2})\).  Positivity has no signed cancellation
which can avoid this.  The static splice must therefore use the canonical
collar recurrence through order \(2M+2\).  C88's scalar gate has the same
exponent and changes only by the explicit ratio
\((2M+2)^2/M^2\), which tends to four and is at most \(16\) for \(M\ge1\):
\[
                         C\mu(2M+2)^2<1.                      \tag{1.1}
\]

This positive result does **not** supply the missing endpoint submersion.
The covariance map is quadratic in the physical amplitudes.  Where all fast
amplitudes vanish, its derivative is identically zero.  More quantitatively,
the norm of the interior amplitude right inverse grows like the reciprocal
of the amplitude and therefore blows up as a flat endpoint is approached.
Thus Borel interpolation, the pointwise three-phase chart, and periodic
anti-divergence alone cannot prove a robust exact terminal map.  That still
requires a dynamical controllability/right-inverse theorem for the complete
packet-plus-wake evolution.

## 2. Exact finite Navier--Stokes jet compatibility

Work on \(\mathbb T^3\) with the projected operator

\[
 \mathcal N_\nu(U)
 :=\partial_tU+\mathbb P\operatorname{div}(U\otimes U)-\nu\Delta U.
 \tag{2.1}
\]

Write \(U_n=\partial_t^nU(t_*)\).  Direct differentiation gives

\[
 \partial_t^n\mathcal N_\nu(U)(t_*)
 =C_n(U_0,\ldots,U_{n+1}),                                  \tag{2.2}
\]

where

\[
 C_n
 =U_{n+1}-\nu\Delta U_n
  +\mathbb P\sum_{a=0}^n{n\choose a}
       U_a\cdot\nabla U_{n-a}.                              \tag{2.3}
\]

Consequently, within the sufficiently regular divergence-free jet space, a
jet \(J_M=(U_0,\ldots,U_M)\) is the order-\(M\) jet of an unforced
Navier--Stokes germ exactly when

\[
                         C_n(J_M)=0,
 \qquad 0\le n\le M-1.                                      \tag{2.4}
\]

The recurrence in C87 is precisely (2.4).  In particular, the full jet is
not a collection of independent controls: \(U_0\) determines every
\(U_n\) recursively.

This also reduces the apparent endpoint burden.  If a transition is
unforced on an outgoing neighborhood, matching its terminal value to the
collar value matches the complete collar jet automatically.  More
generally, if its Reynolds stress is flat to all orders at the seam, the
left jet obeys the same recurrence as the right collar and equality of the
value again forces equality of every time derivative.  The \(M\) collar
coefficients are therefore consistency checks, not \(M\) independent
terminal controls.

The linearization at a compatible base jet is

\[
\begin{aligned}
 DC_n[h]
 &=h_{n+1}-\nu\Delta h_n\\
 &\quad+\mathbb P\sum_{a=0}^n{n\choose a}
   \big(h_a\cdot\nabla U_{n-a}
       +U_a\cdot\nabla h_{n-a}\big).                         \tag{2.5}
\end{aligned}
\]

As a map from \((h_0,\ldots,h_M)\) to
\((DC_0[h],\ldots,DC_{M-1}[h])\), (2.5) is block lower triangular with the
identity multiplying \(h_{n+1}\).  It is therefore onto the product of
the target function spaces.  Its kernel is exactly the graph of the
linearized Navier--Stokes recurrence: \(h_0\) is free and determines all
\(h_1,\ldots,h_M\).

This is the first important distinction:

* the **residual-jet map** is onto because each new time derivative enters
  with coefficient one;
* the **unforced solution-jet map** is not onto the independent jet space;
  its range is the recurrence graph (2.4).

## 3. The only periodic linear compatibility is spatial mean

Consider the projected Navier--Stokes--Reynolds equation

\[
                  \mathcal N_\nu(U)+
                  \mathbb P\operatorname{div}Q=0.            \tag{3.1}
\]

The spatial integral of (3.1) gives

\[
                         {d\over dt}\int_{\mathbb T^3}U=0.   \tag{3.2}
\]

Hence two endpoint values must have the same mean, and an admissible jet
has

\[
 \int U_0^-=\int U_0^+,
 \qquad \int U_n^\pm=0\quad(n\ge1).                          \tag{3.3}
\]

Conversely, let \(U(t)\) be any divergence-free interpolation with constant
mean.  Then \(F=\mathcal N_\nu(U)\) is divergence free and has zero mean.
The C86 operator gives

\[
                         S=\mathcal RF,
 \qquad \operatorname{div}S=F.                              \tag{3.4}
\]

Thus

\[
                         Q=\rho I-S                           \tag{3.5}
\]

solves (3.1), because
\(\mathbb P\operatorname{div}(\rho I)=0\).  No angular-momentum or other
periodic cokernel remains.

At the level of endpoint jets, one can set

\[
                         S_n=\mathcal R C_n(J).                \tag{3.6}
\]

The routed trace-free stress \(S\) is flat through order \(M-1\) if and only
if (2.4) holds.  A
target advertised as an exact collar jet but violating (2.4) cannot be
repaired by a stress which is required to disappear at that endpoint.

The work condition is not an additional scalar control.  Multiplying
(3.1) by \(U\) yields identically

\[
 {1\over2}{d\over dt}\|U\|_2^2+\nu\|\nabla U\|_2^2
 =\int_{\mathbb T^3}Q:\operatorname{sym}\nabla U.             \tag{3.7}
\]

The isotropic part drops out.  Equation (3.7) is forced by (3.1), so the
anti-divergence does not introduce a separate energy parameter.

## 4. Finite Gevrey-2 interpolation

For finite \(M\), no infinite Borel theorem is actually needed.  In any
linear spatial Gevrey-2 space, the two endpoint jet map

\[
 \mathcal H_M:P_{2M+1}\longrightarrow X^{2M+2},
 \qquad
 P\longmapsto
 \big(P^{(0)}(0),\ldots,P^{(M)}(0),
       P^{(0)}(1),\ldots,P^{(M)}(1)\big)                      \tag{4.1}
\]

is an isomorphism.  Indeed, an element of its kernel is divisible by
\(t^{M+1}(t-1)^{M+1}\), which has degree \(2M+2\), while \(P\) has degree at
most \(2M+1\).  Applying the inverse componentwise gives a time-analytic,
spatially Gevrey-2 interpolation.

Because divergence and spatial mean commute with (4.1), data satisfying
(3.3) give a divergence-free constant-mean path.  If the endpoint jets
obey the C88 bound, the standard non-quasianalytic Gevrey-2 Borel operator
can instead be used after a fixed enlargement of the Roumieu constant.
For the finite construction, (4.1) is the completely elementary version;
its inverse norm by itself is not asserted to be uniform in \(M\).

For a one-sided endpoint construction, let \(s\ge0\) be distance to the
seam and require the physical amplitudes to have zero derivatives through
order \(M\).  Put
\[
                         r=M+1,\qquad N=2r=2M+2.              \tag{4.2}
\]
Interpolate the coarse velocity using the canonical collar jet through
order \(N\), not merely \(M\).  Equations (2.2)--(2.4) and the Hadamard
division lemma give
\[
             F=\mathcal N_\nu(U)=s^N\widetilde F,\qquad
             S=\mathcal RF=s^N\widetilde S                    \tag{4.3}
\]
near the seam, with Gevrey-2 quotients when the interpolation is Gevrey-2.

Choose a nonnegative Gevrey-2 scalar \(\alpha\) such that
\[
                         \alpha(s)=s^r
 \quad\hbox{near the seam},                                  \tag{4.4}
\]
and \(\alpha>0\) in the active interior.  Then \(S/\alpha^2\) extends
smoothly across the seam.  Select \(q\) larger than the global operator norm
of this quotient (and later by the fixed three-phase chart constant), and
set
\[
\begin{aligned}
 \widetilde Q&=qI-\alpha^{-2}S,\\
 Q&=\alpha^2\widetilde Q=q\alpha^2I-S.                        \tag{4.5}
\end{aligned}
\]
Here \(\alpha^{-2}S\) denotes its smooth extension at \(\alpha=0\).
Then \(Q\ge c\alpha^2I\), it solves (3.1), and its square-root amplitudes
have the form
\[
                         A_i=\alpha\widetilde A_i.             \tag{4.6}
\]
All derivatives of \(A_i\) through order \(M=r-1\) vanish at the seam.
For the explicit three-phase factorization below, the enlarged \(q\)
ensures \(\|\widetilde Q-qI\|_{\mathrm{op}}\ll q\).

The doubling in (4.2) is unavoidable for a static positive covariance.
If every smooth amplitude is \(O(s^{M+1})\), then the trace of their
covariance is a sum of squares and is \(O(s^{2M+2})\).  Hence the projected
defect must have the same order.  There is no cancellation between
positive squares which can produce an earlier signed stress coefficient.

The square root does not leave the Gevrey-2 factorial class.  If the first
possibly nonzero stress derivative obeys
\[
 \|\partial_s^{2r}S(0)\|
 \le C^{2r}((2r)!)^2,                                       \tag{4.7}
\]
then the coefficient of \(s^{2r}\) is at most
\(C^{2r}(2r)!\).  An amplitude coefficient \(D_rs^r\) dominating its
square may be chosen with
\[
 D_r\le C^r\sqrt{(2r)!}\le(2C)^rr!,                         \tag{4.8}
\]
because \(\binom{2r}{r}\le4^r\).  Therefore
\[
                         |\partial_s^r(D_rs^r)(0)|
                         \le(2C)^r(r!)^2.                    \tag{4.9}
\]
The order doubling is exactly compatible with Gevrey 2.  This endpoint
coefficient ledger does not by itself prove a uniform global Borel
extension or a tame oscillatory corrector, but it rules out an additional
factorial-class loss.

If genuine two-sided exact collar germs in the chosen Gevrey class happen
to be available, \(F\) can instead be made zero on endpoint neighborhoods
and \(\alpha\) can be chosen flat.  The one-sided construction above does
not assume backward parabolic solvability and is the relevant general
statement.

The C88 estimate applied at \(N=2M+2\) requires (1.1).  Since the cascade
has \(\mu M^2\to0\), this constant-factor enlargement creates no new scalar
exponent restriction.  It still gives no cascade-uniform bound for the
exact oscillatory corrector: \(q\) controls a coarse pressure gauge, not
the viscous loss of the unresolved fast field.

## 5. Exact three-phase stress linearization

At undeformed coordinate phases, put

\[
 E_1=(e_2,e_3),\qquad E_2=(e_1,e_3),\qquad E_3=(e_1,e_2),     \tag{5.1}
\]

and let \(H_i\in\operatorname{Sym}_2\).  The transverse stress map is

\[
 \Sigma(H_1,H_2,H_3)=\sum_{i=1}^3E_iH_iE_i^T.                \tag{5.2}
\]

It has the following explicit right inverse.  For
\(G=(g_{ij})\in\operatorname{Sym}_3\), define

\[
\begin{aligned}
 T_1(G)&=
 \begin{pmatrix}0&0&0\\0&g_{22}/2&g_{23}\\0&g_{23}&g_{33}/2\end{pmatrix},\\
 T_2(G)&=
 \begin{pmatrix}g_{11}/2&0&g_{13}\\0&0&0\\g_{13}&0&g_{33}/2\end{pmatrix},\\
 T_3(G)&=
 \begin{pmatrix}g_{11}/2&g_{12}&0\\g_{12}&g_{22}/2&0\\0&0&0\end{pmatrix}.
                                                                    \tag{5.3}
\end{aligned}
\]

Then \(T_i(G)e_i=0\) and

\[
                         T_1(G)+T_2(G)+T_3(G)=G.              \tag{5.4}
\]

For transported directions \(n_i=F^{-T}e_i\), choose smooth transverse
frames \(E_i(n_i)\) and define \(\Sigma_n\) as in (5.2).  Surjectivity is
open.  A canonical smooth right inverse near the coordinate frame is

\[
                         \Sigma_n^*(\Sigma_n\Sigma_n^*)^{-1}.\tag{5.5}
\]

If \(\widetilde Q\) in (4.5) is sufficiently close to \(qI\), this right
inverse gives positive \(H_i\) close to \((q/2)I_2\).  Their positive square
roots supply two transverse polarizations per phase.  Short-time material
transport keeps \(n_i\) in this open chart.

The amplitude-level rank is equally explicit.  Set

\[
                         A_i^0=cE_i,
 \qquad c=\sqrt{q/2},                                        \tag{5.6}
\]

where the columns are the two polarizations.  For a desired stress
variation \(G\), let

\[
                         \delta A_i={1\over2c}T_i(G)E_i.      \tag{5.7}
\]

Then

\[
 \sum_i\left(\delta A_i(A_i^0)^T+A_i^0\delta A_i^T\right)=G.\tag{5.8}
\]

Thus the twelve real transverse amplitude variables have rank six onto
\(\operatorname{Sym}_3\) in the positive interior.  The right-inverse norm
is \(O(c^{-1})\).

## 6. The endpoint degeneracy is unavoidable

Let the physical covariance map be

\[
                         \Gamma(A_1,A_2,A_3)
                         =\sum_iA_iA_i^T.                     \tag{6.1}
\]

At a clean endpoint the fast oscillations must vanish, so \(A_i=0\).
But

\[
                         D\Gamma_{(0,0,0)}[\delta A]=0.       \tag{6.2}
\]

Therefore the pointwise amplitude map has rank zero at that endpoint.
Along the positive splice (4.5), \(c(t)\asymp\alpha(t)\), and (5.7) gives

\[
                         \|D\Gamma_{A(t)}^\dagger\|
                         \gtrsim\alpha(t)^{-1}.               \tag{6.3}
\]

The inverse is not uniform as the amplitudes turn off.  Equivalently, the
tangent cone of positive semidefinite stresses at zero is not all of
\(\operatorname{Sym}_3\), and a covariance built from smooth amplitudes
which vanish at the endpoint has no linear stress variation there.

Retaining a nonzero isotropic \(qI\) at the endpoint avoids (6.2) only in
the formal Reynolds tensor.  If it is actually realized by oscillatory
velocities, those velocities have nonzero terminal energy and the endpoint
is no longer the prescribed bare collar.  Treating \(qI\) solely as a
pressure gauge does not supply physical amplitude rank.

This is not a no-go for **interior dynamical control**.  Variations of
nonzero amplitudes earlier in the interval can propagate to the terminal
low-frequency state.  It is a no-go for deducing that control from static
Borel interpolation or the pointwise covariance chart.

## 7. The controlled Reynolds endpoint derivative is onto

If the three-phase amplitudes are treated as freely prescribable controls,
the preceding ingredients do give an exact finite-order endpoint right
inverse at the Reynolds level.  This is stronger than pointwise covariance
rank and weaker than the required oscillatory Navier--Stokes theorem.

Write the controlled equation as
\[
 \mathcal E(U,A)
 =\mathcal N_\nu(U)
  +\mathbb P\operatorname{div}\Gamma(A)=0.                   \tag{7.1}
\]
At a base splice \((U,A)\), its derivative is
\[
\begin{aligned}
 D\mathcal E_{(U,A)}[h,\dot A]
 &=\mathcal L_Uh+
   \mathbb P\operatorname{div}\big(D\Gamma_A[\dot A]\big),\\
 \mathcal L_Uh
 &=\partial_th+\mathbb P\operatorname{div}
       (U\otimes h+h\otimes U)-\nu\Delta h.                  \tag{7.2}
\end{aligned}
\]

Let \(g_0\) be any divergence-free, zero-mean terminal perturbation in the
chosen Gevrey space.  Generate its terminal jet
\((g_0,\ldots,g_N)\), \(N=2M+2\), by the linearized collar recurrence
(2.5).  Use the Hermite/Borel interpolation to choose \(h\) with zero
incoming jet and this outgoing jet.  Then
\[
                         f=\mathcal L_Uh
                         =s^N\widetilde f                   \tag{7.3}
\]
at the clean endpoint.  Set
\[
                         \dot Q=-\mathcal Rf.                \tag{7.4}
\]
It has zero mean compatibility automatically and satisfies
\[
                         \mathcal L_Uh+
                         \mathbb P\operatorname{div}\dot Q=0.\tag{7.5}
\]

Since \(A_i=\alpha\widetilde A_i\), with
\(\alpha=s^{M+1}\), and \(\dot Q=\alpha^2\widetilde{\dot Q}\),
the scaled version of (5.7) gives
\[
                         \dot A_i
                         =\alpha\,\mathfrak T_i
                           (\widetilde{\dot Q}),              \tag{7.6}
\]
where \(\mathfrak T_i\) is the smooth interior amplitude right inverse.
Thus \(\dot A_i\) is smooth, has zero derivatives through order \(M\), and
\[
                         D\Gamma_A[\dot A]=\dot Q.            \tag{7.7}
\]

Equations (7.2)--(7.7) show that the terminal derivative
\[
 \{\hbox{three-phase amplitude paths flat through \(M\)}\}
 \longrightarrow
 \{\hbox{zero-mean compatible terminal velocity jets through \(M\)}\}
                                                                    \tag{7.8}
\]
is onto for the controlled periodic Reynolds system.  Restricting the
target to the five affine-capture coordinates remains onto.

This construction does not contradict (6.2): the pointwise endpoint
amplitude map has rank zero, while (7.8) uses the complete time-dependent
control path and the equation.  Nor does it close the original problem.
In the actual WKB construction the amplitudes are not arbitrary controls;
transport, incompressibility, every charged corrector, viscosity, and the
zero-charge wake constrain them simultaneously.  A cascade-uniform bound
for the composite right inverse in (7.6) is also not proved here.

## 8. Exact conclusion

The endpoint problem now has a clean two-level answer.

> **Finite Reynolds-splice theorem.**  Equal spatial mean is the complete
> linear compatibility condition for a finite periodic Reynolds splice
> with free stress; the Navier--Stokes recurrence is the additional exact
> condition at every order where the endpoint stress is required to vanish.
> To make its physical
> three-phase amplitudes vanish through order \(M\) while preserving
> positivity, use the canonical collar jet through order \(2M+2\); then
> (4.5)--(4.6) give a positive Gevrey-2 principal covariance with the
> required endpoint flatness.  The doubling changes only the constant in
> the C88 heat gate.  With amplitudes regarded as controls, the resulting
> finite-order terminal derivative is onto every zero-mean compatible
> velocity jet.

> **Missing exact theorem.**  The preceding construction is not an exact
> oscillatory Navier--Stokes transition and its endpoint amplitude
> linearization has rank zero.  One must still prove that interior
> three-phase controls, after every nonzero charge, zero-charge wake,
> pressure correction, and viscous term is evolved, have a uniformly tame
> terminal derivative onto the affine-capture variables.

Consequently C85--C88 remove the finite-jet interpolation, periodic
compatibility, positivity, and controlled principal endpoint-rank
obstacles.  They do not show that the actual phase-resolved
Navier--Stokes/WKB dynamics realizes those controls with a tame norm.
