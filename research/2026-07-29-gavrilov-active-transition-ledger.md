# A localized Euler--Reynolds transition between Gavrilov bubbles

Date: 2026-07-29

## Result and claim boundary

This note resolves the finite-dimensional compatibility part of the
proposed active transition between disjoint scaled Gavrilov bubbles.
It does **not** construct a smooth Navier--Stokes singularity.

There are four concrete conclusions.

1. Gavrilov's pressure modulation can be chosen so that one nonzero
   compact steady bubble has exactly zero linear and angular momentum.
   Linear momentum is automatic for every compact divergence-free field;
   zero angular momentum is obtained by one scalar condition on the
   pressure modulation.
2. The same scalar modulation cannot cancel the helicity of a sufficiently
   thin Gavrilov shell: the helicity is a positive quadratic functional of
   the modulation.  One orientation-reversing reflected copy cancels it
   exactly.  Thus a two-element symmetry orbit is enough for simultaneous
   momentum, angular-momentum, and helicity closure.
3. Because each member of the reflected pair can already be made
   zero-angular-momentum, the compact symmetric anti-divergence can be
   constructed separately inside each small bubble ball.  No
   moment-forced stress bridge between the reflected copies, or between
   adjacent scales, is necessary.
4. For fixed disjoint bubbles and time-dependent amplitudes, the
   Euler--Reynolds defect is exactly the time derivative.  At turnover
   time its compact symmetric stress has the correct inertial size
   \(a^2\), and the viscous correction is smaller by the inverse Reynolds
   number.  The remaining prize-level problem is not a finite moment
   condition.  It is the all-order oscillatory realization and viscous
   endpoint-jet handoff.

The key new design choice is the zero-angular-momentum modulation in
Section 3.  A reflection pair alone cancels total angular momentum, but
then a symmetric anti-divergence of either individual component is
impossible when its angular momentum is nonzero.  The stress must connect
the two components.  Removing the angular momentum of each component
separately eliminates this unnecessary nonlocality.

## 1. The exact moment conditions

Let \(F\in C^\infty_c(\mathbb R^3;\mathbb R^3)\).  If a compactly
supported symmetric tensor \(R=R^T\) satisfies

\[
 \operatorname{div}R=F,
\tag{1.1}
\]

then integration by parts gives

\[
 \int_{\mathbb R^3}F\,dx=0,
 \qquad
 \int_{\mathbb R^3}x\times F\,dx=0.
\tag{1.2}
\]

The second condition follows because

\[
\begin{aligned}
 \left(\int x\times\operatorname{div}R\,dx\right)_k
 &=\epsilon_{k i j}\int x_i\partial_\ell R_{j\ell}\,dx\\
 &=-\epsilon_{k\ell j}\int R_{j\ell}\,dx=0.
\end{aligned}
\tag{1.3}
\]

These are also sufficient.  The compact symmetric-divergence theorem of
Isett and Oh gives a linear operator

\[
 {\cal R}_B:
 \left\{
 F\in C^\infty_c(B;\mathbb R^3):
 \int F=0,\ \int x\times F=0
 \right\}
 \longrightarrow C^\infty_c(B;\operatorname{Sym}_3)
\tag{1.4}
\]

such that

\[
 \operatorname{div}{\cal R}_B F=F.
\tag{1.5}
\]

Here \(B\) may be any ball containing the support with a positive support
margin.  Under a dilation by \(\ell\), the operator gains exactly one
factor of \(\ell\):

\[
 \|\nabla^m{\cal R}_B F\|_\infty
 \lesssim_m
 \ell^{1-m}\|F\|_\infty
\tag{1.6}
\]

when the derivatives of \(F\) occur on scale \(\ell\).  The explicit
averaged-ray formula in Isett--Oh also preserves Gevrey-\(\sigma\) bounds
if both the datum and averaging bump are Gevrey-\(\sigma\), \(\sigma>1\).
Indeed, their formulas are finite sums of parameter integrals of products
of derivatives of \(F\), a fixed averaging bump, and polynomials of
degree at most two.  If

\[
 \|\nabla^mF\|_\infty+
 \|\nabla^m\zeta\|_\infty
 \le C_0C_1^m(m!)^\sigma,
\tag{1.6a}
\]

Leibniz' rule is controlled by

\[
 \binom{m}{k}(k!)^\sigma((m-k)!)^\sigma
 =(m!)^\sigma\binom{m}{k}^{1-\sigma}
 \le(m!)^\sigma.
\tag{1.6b}
\]

The parameter integrals are exactly those bounded in the proof of the
Isett--Oh \(C^m\) estimates.  Summing the at most \(2^m\) derivative
splittings only changes \(C_1^m\).  After spatial rescaling, this gives

\[
 \|\nabla^m{\cal R}_BF\|_\infty
 \le C_0'C_1'{}^m(m!)^\sigma\ell^{1-m}.
\tag{1.6c}
\]

For a compactly supported divergence-free velocity \(V\), its total
linear momentum vanishes without any design:

\[
 \int V_i\,dx
 =\int\partial_j(x_iV_j)\,dx=0.
\tag{1.7}
\]

The symmetric first moment vanishes as well,

\[
 \int(x_iV_j+x_jV_i)\,dx
 =\int\partial_k(x_ix_jV_k)\,dx=0,
\tag{1.8}
\]

but the antisymmetric first moment, equivalently

\[
 J(V):=\int x\times V\,dx,
\tag{1.9}
\]

need not vanish.

## 2. Scaling and rigid-motion ledger

Let \(W\) be one compact steady Euler bubble with pressure \(P\).  For
\(a,\ell>0\), \(c\in\mathbb R^3\), and \(Q\in O(3)\), put

\[
\begin{aligned}
 W_{a,\ell,c,Q}(x)
 &=aQW\!\left(Q^T\frac{x-c}{\ell}\right),\\
 P_{a,\ell,c,Q}(x)
 &=a^2P\!\left(Q^T\frac{x-c}{\ell}\right).
\end{aligned}
\tag{2.1}
\]

This is again a steady Euler solution.  Since \(\int W=0\), translation
does not change angular momentum.  Direct change of variables gives

\[
\begin{aligned}
 E(W_{a,\ell,c,Q})
 &=a^2\ell^3E(W),\\
 J(W_{a,\ell,c,Q})
 &=a\ell^4\det(Q)\,QJ(W),\\
 H(W_{a,\ell,c,Q})
 &=a^2\ell^2\det(Q)\,H(W),
\end{aligned}
\tag{2.2}
\]

where

\[
 E(V)=\frac12\int|V|^2,\qquad
 H(V)=\int V\cdot\operatorname{curl}V.
\tag{2.3}
\]

The determinant in the last two lines is essential: angular momentum is
an axial vector and helicity is a pseudoscalar.  In particular, central
reflection \(Q=-I\) flips helicity but leaves \(J\) unchanged.  A generic
orientation-reversing reflection cancels total helicity, but it does not
by itself guarantee that the two individual components have zero angular
momentum.

There is a precise support consequence.  Suppose \(F=F_1+F_2\), with
\(F_i\) supported in disjoint open sets \(D_i\), and suppose a symmetric
solution of \(\operatorname{div}R=F\) has support in
\(D_1\cup D_2\).  Then \(R=R_1+R_2\) with
\(\operatorname{supp}R_i\subset D_i\), so (1.2) applies to each \(F_i\)
separately.  Therefore

\[
 \int x\times F_1\ne0
 \quad\Longrightarrow\quad
 \operatorname{supp}R\ \hbox{cannot stay in }D_1\cup D_2.
\tag{2.4}
\]

For the time derivative of a reflection pair with nonzero individual
angular momentum, the anti-divergence must cross the gap between the
copies.  Global cancellation alone is not a localized construction.

## 3. A zero-angular-momentum Gavrilov modulation

Gavrilov's local axisymmetric solution has cylindrical form

\[
 u=\frac1{\rho}
 \left(
 p_z e_\rho-p_\rho e_z+b(p)e_\varphi
 \right),
\tag{3.1}
\]

on a thin neighborhood of a circle.  It obeys

\[
 u\cdot\nabla p=0,\qquad |u|^2=3p,
\tag{3.2}
\]

and consequently, for every smooth scalar function \(\chi\),

\[
 W=\chi(p)u
\tag{3.3}
\]

is another steady Euler flow, with modified pressure determined by

\[
 dP_\chi=\chi(p)^2\,dp.
\tag{3.4}
\]

Taking \(\chi\) compactly supported in a sufficiently small positive
pressure interval makes \(W\) smooth and compactly supported.

Axisymmetry makes \(J(W)\) parallel to the symmetry axis.  Its axial
component is

\[
 J_z(W)
 =\int \rho W_\varphi\,dx
 =\int \chi(p)b(p)\,dx.
\tag{3.5}
\]

On a sufficiently small positive pressure interval \(I\), \(b(p)>0\).
Choose two nonzero, nonnegative Gevrey-\(\sigma\) bumps
\(\chi_1,\chi_2\) with disjoint supports contained in \(I\), and define

\[
 L_i=\int\chi_i(p)b(p)\,dx>0,
 \qquad
 \chi=\chi_1-\frac{L_1}{L_2}\chi_2.
\tag{3.6}
\]

Then \(\chi\not\equiv0\), and (3.5) gives

\[
 J(W)=0.
\tag{3.7}
\]

Thus:

> **Zero-angular-momentum seed lemma.**
> For every \(\sigma>1\), Gavrilov's construction contains a nonzero
> compact Gevrey-\(\sigma\) steady Euler velocity \(W\) satisfying
> \[
> \int W\,dx=0,\qquad \int x\times W\,dx=0.
> \tag{3.8}
> \]

The sign change in \(\chi\) is harmless: the Euler pressure sees
\(\chi^2\), while angular momentum sees \(\chi\).

## 4. Why the same modulation does not cancel helicity

The derivative term in

\[
 \operatorname{curl}(\chi(p)u)
 =\chi(p)\operatorname{curl}u+
 \chi'(p)\nabla p\times u
\tag{4.1}
\]

does not contribute to helicity, since
\(u\cdot(\nabla p\times u)=0\).  Therefore

\[
 H(W)
 =\int\chi(p)^2\,u\cdot\operatorname{curl}u\,dx.
\tag{4.2}
\]

For (3.1), a direct cylindrical calculation gives

\[
 u\cdot\operatorname{curl}u
 =
 \frac{b\,\Delta_*p-b'|\nabla p|^2}{\rho^2},
 \qquad
 \Delta_*p=p_{\rho\rho}+p_{zz}-\rho^{-1}p_\rho.
\tag{4.3}
\]

Gavrilov's identities

\[
 \Delta_*p+bb'=\frac52\rho^2,
 \qquad
 |\nabla p|^2+b^2=3p\rho^2
\tag{4.4}
\]

reduce this to the scalar function

\[
 h(p):=u\cdot\operatorname{curl}u
 =\frac52b(p)-3p\,b'(p).
\tag{4.5}
\]

Near the seed circle,

\[
 b(p)=C\sqrt p\,(1+O(p)),
 \qquad
 h(p)=C\sqrt p\,(1+O(p))>0.
\tag{4.6}
\]

Hence every nonzero modulation supported in a sufficiently thin positive
pressure shell has

\[
 H(\chi(p)u)=\int\chi(p)^2h(p)\,dx>0.
\tag{4.7}
\]

This proves a small but exact obstruction: changing the sign of the
pressure modulation can kill angular momentum but cannot kill helicity.

Let \(S\in O(3)\) be orientation reversing, and place a rigidly moved
copy of \(W\) in a small ball disjoint from its image under \(S\).  Define

\[
 W^S(x)=S W(S^Tx),\qquad P^S(x)=P(S^Tx).
\tag{4.8}
\]

The disjoint sum

\[
 G=W+W^S
\tag{4.9}
\]

is again a compact steady Euler field.  Equations (2.2), (3.7), and
\(\det S=-1\) give

\[
 \int G=0,\qquad J(G)=0,\qquad H(G)=0.
\tag{4.10}
\]

Thus the two-element orbit \(\{I,S\}\) is enough.  Choosing the small
seed ball generically away from the dilation centre also makes the
resulting scaled lattice nonaxisymmetric about that centre.

## 5. An exact disjoint-bubble interpolation

First consider an arbitrary smooth compactly supported divergence-free
path \(v(t)\), choose a compactly supported pressure \(\pi(t)\) (in
particular, \(\pi=0\) is allowed), and let

\[
 F=\partial_tv+\operatorname{div}(v\otimes v)+\nabla\pi.
\tag{5.1a}
\]

The nonlinear tensor \(v\otimes v\) is symmetric, as is \(\pi I\).
Therefore

\[
 \int F\,dx=\frac d{dt}\int v\,dx=0,
 \qquad
 \int x\times F\,dx=\frac d{dt}J(v).
\tag{5.1b}
\]

Thus a compact symmetric anti-divergence exists along a general path
precisely when its angular momentum is constant.  This is the complete
finite-dimensional compatibility condition; helicity is not in the
cokernel of the symmetric-divergence operator.

Let \((B_m,P_m)\), \(m=0,\ldots,M\), be finitely many steady Euler
bubbles whose **joint velocity-pressure supports** are pairwise
disjoint.  They may represent one parent, one child, and finitely many
outgoing wake bubbles.
Assume each \(B_m\) is a rigidly moved and scaled copy of the seed from
Section 3, so

\[
 \int B_m=0,\qquad J(B_m)=0.
\tag{5.1}
\]

Let \(\lambda_m\in C^\infty([0,1])\), set
\(s=(t-t_0)/\tau\), and define

\[
\begin{aligned}
 v(t,x)&=\sum_{m=0}^M\lambda_m(s)B_m(x),\\
 \pi(t,x)&=\sum_{m=0}^M\lambda_m(s)^2P_m(x).
\end{aligned}
\tag{5.2}
\]

Because the supports are disjoint, all quadratic cross terms vanish and
each amplitude-scaled bubble remains steady Euler:

\[
 \operatorname{div}(v\otimes v)+\nabla\pi=0.
\tag{5.3}
\]

The exact Euler defect is therefore

\[
 F_E
 :=\partial_tv+\operatorname{div}(v\otimes v)+\nabla\pi
 =\frac1\tau\sum_{m=0}^M\lambda_m'(s)B_m.
\tag{5.4}
\]

Every summand in (5.4) separately satisfies both conditions (1.2).
Choose mutually disjoint balls \(D_m\) containing the joint
velocity-pressure support of each bubble, and put

\[
 R_E
 =\frac1\tau\sum_{m=0}^M
 \lambda_m'(s){\cal R}_{D_m}B_m.
\tag{5.5}
\]

Then

\[
 \operatorname{div}R_E=F_E,
 \qquad
 \operatorname{supp}R_E
 \subset\bigcup_mD_m.
\tag{5.6}
\]

This is the promised removal of the finite moment obstruction.  The
stress does not need to cross from the parent annulus to the child
annulus, nor from one reflected copy to the other.

If \(B_m\) has velocity size \(a_m\) and spatial scale \(\ell_m\), then

\[
 \|\nabla^nR_{E,m}\|_\infty
 \lesssim_n
 \frac{a_m\ell_m^{1-n}}{\tau}
 \max_{s}|\lambda_m'(s)|.
\tag{5.7}
\]

At turnover time

\[
 \tau\asymp\frac{\ell}{a},
\tag{5.8}
\]

and for child/wake scales comparable by a fixed scale ratio, (5.7) is

\[
 \|R_E\|_\infty\asymp a^2.
\tag{5.9}
\]

Thus the Reynolds stress has exactly the inertial size required for a
principal oscillatory realization; there is no extra small denominator.

## 6. Positivity and the exact work identity

Use the Euler--Reynolds convention

\[
 \partial_tv+\operatorname{div}(v\otimes v)+\nabla\pi
 =\operatorname{div}R_E.
\tag{6.1}
\]

Choose a nonnegative Gevrey bump \(\rho\), supported in the same disjoint
balls and strictly positive on \(\operatorname{supp}R_E\), large enough
that

\[
 Q:=\rho I-R_E>0
\tag{6.2}
\]

on the active region.  Then

\[
 \partial_tv+\operatorname{div}(v\otimes v)
 +\nabla(\pi-\rho)
 =-\operatorname{div}Q.
\tag{6.3}
\]

The positive tensor \(Q\) is the covariance to be supplied by fast
oscillations.  Equation (6.2) is only an algebraic positive-covariance
reduction.  It does not prove that a compact smooth divergence-free
oscillation with that covariance, transport law, and endpoint phase
exists.

The isotropic gauge changes neither the projected defect nor its work.
Indeed,

\[
\begin{aligned}
 E_v'(t)
 &=\int v\cdot\partial_tv\,dx\\
 &=\int v\cdot\operatorname{div}R_E\,dx\\
 &=-\int R_E:\operatorname{sym}\nabla v\,dx\\
 &=\int Q:\operatorname{sym}\nabla v\,dx.
\end{aligned}
\tag{6.4}
\]

For the disjoint path (5.2),

\[
 E_v(t)=\sum_m\lambda_m(s)^2E(B_m),
\qquad
 E_v'(t)=\frac2\tau\sum_m
 \lambda_m(s)\lambda_m'(s)E(B_m).
\tag{6.5}
\]

Thus work compatibility imposes no additional hidden scalar condition:
it is exactly the derivative of the chosen coarse energy.  What is
required for a Navier--Stokes stage is that parent energy minus
child-plus-wake energy equal the actual viscous loss of the full
low-plus-high field.  A wake amplitude supplies a scalar parameter for
that balance.  Constant total energy is not allowed once viscosity is
present.

## 7. Viscosity and the bare-endpoint obstruction

Write the Navier--Stokes--Reynolds equation as

\[
 \partial_tv+\operatorname{div}(v\otimes v)+\nabla\pi
 -\nu\Delta v
 =\operatorname{div}R_\nu.
\tag{7.1}
\]

Since \(v\) is divergence free,

\[
 \operatorname{div}
 \big(\nabla v+(\nabla v)^T\big)=\Delta v.
\tag{7.2}
\]

Therefore the exact stress is

\[
 R_\nu
 =R_E-\nu\big(\nabla v+(\nabla v)^T\big).
\tag{7.3}
\]

The viscous term satisfies the same moment conditions automatically.
At bubble scale \((a,\ell)\),

\[
 \|R_{\rm visc}\|_\infty
 \asymp\frac{\nu a}{\ell},
\qquad
 \frac{\|R_{\rm visc}\|_\infty}{a^2}
 \asymp\frac{\nu}{a\ell}
 =\operatorname{Re}^{-1}.
\tag{7.4}
\]

Thus it is perturbative in the cascade window where
\(\operatorname{Re}\to\infty\).

Multiplying (7.1) by \(v\) gives

\[
 E_v'+\nu\|\nabla v\|_2^2
 =-\int R_\nu:\operatorname{sym}\nabla v.
\tag{7.5}
\]

Using

\[
 2\int|\operatorname{sym}\nabla v|^2
 =\int|\nabla v|^2
\tag{7.6}
\]

for compact divergence-free \(v\), (7.3) reduces (7.5) exactly to
(6.4).  There is no missing factor in the work ledger.

There is, however, an unavoidable seam obstruction.  If a nonzero steady
bubble is held constant in an endpoint collar, then \(R_E=0\) but

\[
 R_\nu=-\nu\big(\nabla v+(\nabla v)^T\big)\ne0.
\tag{7.7}
\]

No isotropic pressure gauge can turn this trace-free anisotropic tensor
into zero.  Equivalently, the bare endpoint residual is
\(-\nu\Delta v\), which is not flat along the shrinking cascade.

Consequently an endpoint cannot consist only of a steady Gavrilov
bubble.  It must carry one of the following equivalent pieces of data:

* the full forward Navier--Stokes jet;
* a nonzero oscillatory bath whose stress begins with (7.7); or
* a wake/centre correction that absorbs the same jet.

This is an infinite-dimensional endpoint condition.  It is not repaired
by the finite symmetry orbit.

## 8. Reflection symmetry controls helicity through the transition

Choose all parent, child, and wake bubbles in reflected pairs and use the
same amplitude on the two members of each pair.  Then

\[
 v(t,Sx)=Sv(t,x)
\tag{8.1}
\]

for an orientation-reversing \(S\).  The Isett--Oh stresses can be paired
by

\[
 R^S(t,x)=S R(t,S^Tx)S^T,
\tag{8.2}
\]

so the full Euler--Reynolds or Navier--Stokes--Reynolds system preserves
the same symmetry.

Curl is an axial vector under \(S\):

\[
 \operatorname{curl}v(t,Sx)
 =-\;S\operatorname{curl}v(t,x).
\tag{8.3}
\]

Hence the helicity density is odd and

\[
 H(v(t))=0
\tag{8.4}
\]

at every transition time, not merely at the endpoints.  The same
symmetry cancels every pseudoscalar helicity-production integral.  There
is therefore no separate scalar helicity parameter left to tune.

The local angular-momentum condition can also be built into the fast
carrier rather than repaired afterward.  If

\[
 w=\operatorname{curl}\operatorname{curl}\Psi,
 \qquad \Psi\in C^\infty_c(D;\mathbb R^3),
\tag{8.5}
\]

then \(w\) is compactly supported and divergence free.  Writing
\(A=\operatorname{curl}\Psi\), an integration by parts gives

\[
 \int x\times\operatorname{curl}A\,dx=2\int A\,dx.
\tag{8.6}
\]

Since \(\int\operatorname{curl}\Psi=0\),

\[
 \int w\,dx=0,\qquad \int x\times w\,dx=0.
\tag{8.7}
\]

Every localized polarized plane wave can be produced to principal order
in the form (8.5), with the envelope terms entering lower WKB orders.
Thus a double-curl carrier enforces the two compact anti-divergence
moments separately in every bubble ball.  The higher-order construction
must retain this representation, but it does not spend a new global wake
parameter.

## 9. Concrete transition lemma and the remaining theorem

The calculations prove the following finite-dimensional reduction.

> **Localized Gavrilov transition lemma.**
> Fix \(\sigma>1\).  Let a parent, child, and finitely many wake states be
> finite sums of mutually disjoint, scaled and rigidly moved
> Gavrilov bubbles built from the zero-angular-momentum modulation in
> Section 3 and paired by one improper reflection.  Let their scalar
> amplitudes follow any Gevrey-\(\sigma\) path.
>
> Then there is a compactly supported, reflection-equivariant symmetric
> stress \(R_E\) realizing the exact Euler defect, separately localized
> inside the disjoint bubble balls.  At turnover time it has inertial size
> \(a^2\) and Gevrey-\(\sigma\) derivative bounds.  After adding the
> explicit stress \(-2\nu\operatorname{sym}\nabla v\), the same statement
> holds for the Navier--Stokes defect, with relative viscous size
> \(\operatorname{Re}^{-1}\).  An isotropic gauge makes the associated
> covariance algebraically positive without changing work.  Linear
> momentum, angular momentum, and total helicity vanish identically.

This is a finite-dimensional anti-divergence lemma, not a realization
lemma.  In particular, it does not realize the positive covariance by an
exact smooth velocity.  The remaining theorem is now sharply isolated:

> **One-carrier Gevrey realization target.**
> Realize \(Q=\rho I-R_\nu\) by a phase-resolved, localized oscillatory
> field at the polynomial carrier \(K_j=j^A\), solve its transport,
> pressure, resonance, and viscous corrector equations through
> \(M_j\asymp j^2/\log j\), and hand the full oscillatory and pressure jet
> to the next stage with residual \(e^{-c j^2}\).

The leading and carrier-level finite moment and helicity equations no
longer obstruct this target.  Every subsequent WKB correction still has
to be projected onto the same local moment-free class.  The hard parts
are:

1. a Gevrey-tame inverse for the multiwave cell equations without
   frequency inflation;
2. phase and circulation provenance for the child carrier;
3. the viscous endpoint invariant graph forced by (7.7); and
4. control of the pressure/centre tails created by localization.

These are analytic, not GPU-search, bottlenecks.

## Primary sources

* A. V. Gavrilov, *A steady Euler flow with compact support*:
  https://arxiv.org/abs/1810.08020
* P. Isett and S.-J. Oh, *On nonperiodic Euler flows with Hölder
  regularity*, especially Section 11 on the compact symmetric divergence
  equation:
  https://arxiv.org/abs/1402.2305
* C. De Lellis and L. Székelyhidi Jr., *Dissipative continuous Euler
  flows*:
  https://arxiv.org/abs/1202.1751
