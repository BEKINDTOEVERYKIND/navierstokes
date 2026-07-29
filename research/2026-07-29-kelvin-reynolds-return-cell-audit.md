# Kelvin--Reynolds localized return cell: exact algebra, two obstructions, and the corrected theorem target

Date: 2026-07-29

## Claim boundary

This note does **not** construct a Navier--Stokes singularity or even a
localized Euler return cell.  It audits the most concrete surviving
candidate:

1. an approximately affine strain amplifies a Kelvin packet from a parent
   scale to a smaller active scale;
2. the amplified packet's Reynolds stress drains the old strain and creates
   the next strain packet;
3. the complement is retained as a finite-energy wake;
4. stage viscosity is removed to increasing order.

Three exact conclusions result.

* There is no energy-sign obstruction.  The affine strain has an explicit
  positive-semidefinite, rank-two strict-drain stress made only from its two
  compressive polarizations.
* Kelvin circulation and zero helicity can be reconciled at the level of the
  amplifier.  Central-reflection symmetry makes total helicity exactly zero,
  while the apparently growing wavelength-scale circulation is exactly the
  circulation of a wider material preimage.
* Two shortcuts fail.  A single exactly polarized shear is nonlinearly dark
  and cannot reset the parent.  Also, pressure-quadrupole cancellation
  imposed throughout the affine amplification makes the affine energy
  transfer identically zero.

There is a further functional-analytic correction: a linearized Euler
operator with both complete endpoints fixed has an infinite-dimensional
cokernel.  The viable theorem cannot ask for that operator to have a right
inverse after quotienting only finitely many neutral modes.  It must leave
the outgoing wake free and impose only a projected child endpoint.

The remaining step is still prize-level: prove a localized nonlinear Euler
endpoint submersion, in a packet-plus-wake space, near the algebraic cell
described below.

## 1. Exact affine Kelvin amplifier

Let

\[
 S=\operatorname{diag}(-\alpha,-\beta,\alpha+\beta),
 \qquad \alpha>0,\qquad
 \frac{\beta}{\alpha}=\gamma\in(1,3/2).
\]

The affine field \(U=Sx\) is an exact steady Euler solution with pressure
\(-x\cdot S^2x/2\).  The field

\[
 w_2(x,t)=A(t)e_2\sin(k(t)x_1)
\]

gives another exact Euler solution \(U+w_2\) when

\[
 k'=\alpha k,\qquad A'=\beta A.
\]

Consequently, during

\[
 h=\frac{\log r}{\alpha},
\]

one has

\[
 k(h)=rk(0),\qquad A(h)=r^\gamma A(0).
\]

The second compressive polarization is

\[
 w_1(x,t)=B(t)e_1\sin(m(t)x_2),
\qquad
 m'=\beta m,\qquad B'=\alpha B.
\]

Each field \(U+w_i\) is exact separately.  Their uncoloured sum is not in
general exact, because the two shears have a cross interaction once
\(k(t)\ne m(t)\).  Thus the pair is an algebraic/WKB dictionary, not an
exact multiwave solution.

For a general Kelvin mode

\[
 w=\operatorname{Re}\{a(t)e^{i\xi(t)\cdot x}\},
\qquad \xi\cdot a=0,
\]

the affine equations are

\[
 \xi'=-S^T\xi,\qquad
 a'=-Sa+2\xi\frac{\xi\cdot Sa}{|\xi|^2}.
\]

The two axis-aligned fields above are the two growing elementary solutions.
This is the diagonal special case of the Craik--Criminale construction.

With viscosity, the first mode remains exact after replacing

\[
 A'=\left(\beta-\nu k^2\right)A.
\]

Its damping over the interval \(h\) is

\[
 \exp\!\left(
   -\frac{\nu k(0)^2(r^2-1)}{2\alpha}
 \right).
\]

At physical stage scale \(\ell\), with strain
\(\alpha\asymp a/\ell\) and carrier \(k(0)\asymp\Lambda/\ell\), the exponent
is \(O(\Lambda^2\nu/(a\ell))\).  Thus a fixed internal carrier ratio
\(\Lambda\) is compatible with the small stage viscosity
\(\varepsilon=\nu/(a\ell)\to0\).

More usefully, the internal ratio may grow with the cascade.  If

\[
 \ell_j=r^{-j},\qquad
 \Lambda_j=\ell_j^{-\kappa},
\]

then carrier viscosity is still perturbative provided

\[
 \varepsilon_j\Lambda_j^2
 \asymp\ell_j^{\gamma-1-2\kappa}\longrightarrow0,
 \qquad
 0<\kappa<\frac{\gamma-1}{2}.
\]

If a WKB/localization hierarchy can be solved to order \(M_j\asymp j\),
then

\[
 \Lambda_j^{-M_j}\le e^{-c j^2}.
\]

The balanced choice
\(\kappa=(\gamma-1)/3\) makes the leading localization and carrier-heat
exponents equal.  This removes a purely scale-theoretic objection to an
all-order localized Kelvin construction; it does not construct the
stage-dependent higher-carrier seed.

## 2. An exact strict-drain cone

The low flow in an Euler--Reynolds decomposition obeys

\[
 \partial_tV+\mathbb P\operatorname{div}(V\otimes V+R)=0.
\]

In a region where \(\nabla V=S\), its energy exchange with \(R\) is

\[
 \frac12\frac{d}{dt}\|V\|_2^2
 =\int R:S.
\]

Define, for \(\kappa>0\),

\[
 \begin{aligned}
 R_{\rm drain}
 &:=\kappa\big((\alpha+\beta)I-S\big)\\
 &=\kappa\operatorname{diag}
       (2\alpha+\beta,\ \alpha+2\beta,\ 0).
 \end{aligned}
\]

Then

\[
 R_{\rm drain}\succeq0,\qquad
 R_{\rm drain}:S=-\kappa|S|^2<0.
\]

Moreover,

\[
 R_{\rm drain}
 =
 \kappa(2\alpha+\beta)e_1\otimes e_1
 +\kappa(\alpha+2\beta)e_2\otimes e_2.
\]

These are precisely the covariance rays of the two growing compressive
Kelvin polarizations in Section 1.  The isotropic part in the definition of
\(R_{\rm drain}\) is pressure gauge:

\[
 \mathbb P\operatorname{div}
 \left(\rho(\alpha+\beta)I\right)=0
\]

for every scalar envelope \(\rho\).  Hence
\(\rho R_{\rm drain}\) is dynamically equivalent to
\(-\kappa\rho S\), while retaining a positive-semidefinite rank-one
decomposition.

This disposes of two possible no-go arguments:

1. positivity of a realizable Reynolds covariance does not forbid draining
   this strain;
2. an \(e_3\)-polarized growing wave is not needed to represent the diagonal
   drain modulo pressure.

It does not prove that the stress follows the required time history or
creates the desired child field.

### 2.1 The same tensor has the correct affine jet

There is also no local linear-jet obstruction.  On a region where

\[
 \rho(x)=\rho_0-\frac q2|x|^2,
 \qquad q>0,
\]

one has

\[
 \operatorname{div}(\rho R_{\rm drain})
 =-qR_{\rm drain}x.
\]

Since

\[
 \operatorname{tr}R_{\rm drain}
 =3\kappa(\alpha+\beta),\qquad
 R_{\rm drain}
 -\frac{\operatorname{tr}R_{\rm drain}}3I
 =-\kappa S,
\]

the polynomial Leray projection is

\[
 \mathbb P\operatorname{div}(\rho R_{\rm drain})
 =q\kappa Sx.
\]

Thus, in the affine Euler--Reynolds equation, this choice gives

\[
 S'+q\kappa S=0:
\]

the same positive stress both has strict negative work and drains the affine
jet.  Reversing the curvature reverses the jet direction.

This is only an exact polynomial germ.  A compact cutoff changes the
canonical pressure by a harmonic interior field determined by the boundary
stress.  Controlling or incorporating that harmonic Hessian is part of the
localized pressure/wake problem; the calculation is not a compactly
supported Euler--Reynolds solution.

## 3. Exact localization germ and the single-shear no-go

Let \(\eta,\theta\in\mathbb S^2\) satisfy
\(\eta\cdot\theta=0\), put \(\zeta=\eta\times\theta\), and let \(a\) be a
smooth compactly supported envelope.  The vector-potential block

\[
 W_{\Lambda}[a,\theta,\eta]
 :=
 \Lambda^{-1}
 \nabla\times
 \left(a(x)\zeta\cos(\Lambda\eta\cdot x)\right)
\]

is exactly divergence-free and has expansion

\[
 W_{\Lambda}
 =
 a\theta\sin(\Lambda\eta\cdot x)
 +\Lambda^{-1}(\nabla a\times\zeta)
      \cos(\Lambda\eta\cdot x).
\]

After projection below the carrier, its self-covariance has leading term

\[
 \frac12a^2\theta\otimes\theta.
\]

Taking \(a_i\) to be constants times a common cutoff \(\chi\), the two
choices

\[
 (\theta,\eta)=(e_1,e_2),\qquad
 (\theta,\eta)=(e_2,e_1)
\]

therefore realize \(\chi^2R_{\rm drain}\) at leading low frequency.
Different carrier radii or physical colours are needed to prevent their
cross interaction from entering the child band.

There is an exact obstruction behind the qualifier "leading."  If

\[
 W=\theta f,\qquad \theta\ \hbox{constant},\qquad
 \nabla\cdot W=0,
\]

then \(\theta\cdot\nabla f=0\), and hence

\[
 (W\cdot\nabla)W
 =\theta f\,(\theta\cdot\nabla f)=0.
\]

Thus every exactly divergence-free, single-polarization shear is
nonlinearly dark.  It can be amplified by the strain but cannot create the
reset by itself.  A reset requires at least one of:

1. envelope variation along the polarization together with its
   divergence correction;
2. cross-polarization interactions;
3. a nonlinear background/wake component.

The first option is exactly where the low rank-one stress appears, but the
finite-\(\Lambda\) correction and its pressure are part of the nonlinear
Euler cell.  They cannot be discarded as a terminally flat force.

## 4. Helicity is removable by an exact invariant symmetry

Let \(J=-I\).  Restrict the full velocity to the polar-equivariant subspace

\[
 u(-x,t)=-u(x,t).
\]

Euler and Navier--Stokes preserve this subspace.  Since the curl of an odd
velocity is even,

\[
 \omega(-x,t)=\omega(x,t),
\]

and the helicity density is odd:

\[
 (u\cdot\omega)(-x,t)=-(u\cdot\omega)(x,t).
\]

Therefore

\[
 \int u\cdot\omega\,dx=0
\]

exactly at every time.  The same parity makes every coefficient of a formal
viscous correction helicity-neutral when the correction is constructed
inside the odd subspace.

The affine strain, the two sine shears, and
\(W_\Lambda[a,\theta,\eta]\) with even \(a\) are all odd.  Consequently no
48-copy isotropization is required merely to cancel helicity.  A single
improper symmetry suffices and leaves the anisotropic drain tensor
available.

This proves only total helicity neutrality.  It does not settle vortex-line
topology or the endpoint return.

## 5. Exact Kelvin-circulation provenance

The material flow map of the **full** solution \(Sx+w_2\), not merely of its
affine part, can be integrated explicitly.  If \(X=(X_1,X_2,X_3)\) is the
initial label, then at the amplification time

\[
\begin{aligned}
 x_1(h)&=r^{-1}X_1,\\
 x_2(h)&=r^{-\gamma}X_2+
 \frac{A(0)}{2\beta}
 (r^\gamma-r^{-\gamma})\sin(k(0)X_1),\\
 x_3(h)&=r^{1+\gamma}X_3.
\end{aligned}
\]

The extra \(x_2\) displacement is a phase-dependent shear, but the
derivative of the map in the \(X_2\) direction is exactly
\(r^{-\gamma}\).

Consider the first shear.  At the output, take a rectangle in the
\(x_1x_2\)-plane whose two \(x_2\)-directed sides lie at phases
\(\pi/2\) and \(3\pi/2\), and whose \(x_2\)-length is \(L_0/r\).
The affine field \(Sx\) is a gradient and contributes no circulation.
The shear circulation is

\[
 \Gamma_{\rm out}
 =2A(h)\frac{L_0}{r}
 =2A(0)L_0r^{\gamma-1}.
\]

The material preimage of each \(x_2\)-directed side has \(X_2\)-length

\[
 r^\gamma\frac{L_0}{r}
 =L_0r^{\gamma-1},
\]

and the same two wave phases.  The two other sides of the preimage are
sheared curves; their shear-wave circulation contributions cancel because
they traverse the same phase interval in opposite directions.  Its initial
circulation is therefore

\[
 \Gamma_{\rm in}
 =2A(0)L_0r^{\gamma-1}
 =\Gamma_{\rm out}.
\]

So the factor \(r^{\gamma-1}\) relative to a geometrically scaled
parent-size loop is not circulation creation.  The output loop is the image
of a wider input loop.  In a localized design, this yields a precise
footprint requirement: the incoming mediator must contain
\(r^{\gamma-1}\) times the transverse circulation width, or an equivalent
number of circulation strands.

The mechanism is therefore compatible with Kelvin's theorem, but it is not
free.  If the incoming packet does not contain that wider strip/strand
family, the proposed handoff fails.  The calculation also does not prove
that the vorticity supporting the new low strain is on the correct
coadjoint orbit; that remains an endpoint condition.

## 6. A pressure-darkness/amplification incompatibility

Let \(W\) be a sufficiently decaying perturbation of the constant affine
strain and define its integrated covariance

\[
 C(t)=\int W(x,t)\otimes W(x,t)\,dx.
\]

The linearized Euler energy identity is

\[
 \frac12\frac{d}{dt}\|W(t)\|_2^2=-S:C(t).
\]

For a compact packet, the trace-free part \(C^\circ\) is also the coefficient
of the leading canonical pressure quadrupole.  In particular, the usual
finite-symmetry way to kill that quadrupole is to impose

\[
 C(t)=c(t)I.
\]

But then

\[
 S:C(t)=c(t)\operatorname{tr}S=0.
\]

Hence:

> A perturbation whose integrated stress is pressure-quadrupole-dark at
> every instant cannot gain energy from a globally constant trace-free
> affine strain.

This rules out symmetrizing the amplifier by the full cubic reflection group
throughout the active interval: that symmetry isotropizes \(C(t)\) and
removes the very affine work needed for growth.

It does **not** rule out:

1. the single-reflection helicity cancellation in Section 4;
2. terminal-only pressure multipole cancellation after the energy transfer;
3. a non-affine parent whose strain varies among the symmetry copies;
4. retaining the anisotropic pressure quadrupole as part of a closed
   packet-plus-wake state.

The fourth option appears unavoidable for the simplest affine cell.
Ordinary support separation is not enough: a newly shed wake is only
order-one stage distances from the child, so its quadrupolar pressure jet is
also order one in child-normalized variables.

## 7. What the 2025 rank-one technology does and does not give

Cheskidov--Dai--Palasek use directions
\(\theta_j\) whose tensors \(\theta_j\otimes\theta_j\) span the symmetric
matrices and the smooth geometric identity

\[
 R=\sum_j\Gamma_j(R)^2\theta_j\otimes\theta_j
\]

near a positive matrix.  Their amplitudes are chosen so that the
high-frequency self-stress prescribes a lower-frequency symmetric tensor
modulo an isotropic pressure.  Velocity potentials provide exact
incompressibility and control nonlocal projection errors.

This imports two useful ingredients here:

1. the local stress map is smoothly right-invertible around a strictly
   positive covariance;
2. vector-potential localization is compatible with a prescribed
   rank-one low stress.

The affine drain matrix in Section 2 lies on the boundary of the positive
cone but has an explicit two-ray factorization.  If a full tensor
neighborhood is needed, one may add an isotropic pressure gauge and use the
interior rank-one lemma.

What does not import is the dynamical orientation.  Their proved hierarchy
is a forward, heat-assisted high-to-low cascade beginning with an
infinite-frequency reservoir.  It does not show that the localized
Kelvin-amplify/reset ansatz is one exact smooth Euler trajectory, nor that
its endpoint map is a return.

In particular, an \(O(\Lambda^{-1})\) WKB error at a fixed internal carrier
ratio is an order-one normalized error at every cascade stage.  It is not a
flat external force.  There are two logically viable remedies:

1. make the Euler cell exact and then expand only in stage viscosity;
2. use the growing carrier
   \(\Lambda_j=\ell_j^{-\kappa}\),
   \(0<\kappa<(\gamma-1)/2\), and solve both the WKB and viscous hierarchies
   to order \(M_j\asymp j\).

The second option gives
\(\Lambda_j^{-M_j}+\left(\varepsilon_j\Lambda_j^2\right)^{M_j}
\le e^{-cj^2}\) after decreasing the order constant.  It replaces one
fixed Euler return profile by a stage-dependent, increasingly oscillatory
annular state.  The missing endpoint and wake estimates are not weakened.

## 8. Two rigidity gates for the outer wake and steady parent

### 8.1 A strictly homogeneous axisymmetric wake is excluded

In child variables, an old annular wake at normalized radius
\(|y|\asymp\ell_i/\ell_j\) has the natural size

\[
 \frac{a_i}{a_j}
 =
 \left(\frac{\ell_i}{\ell_j}\right)^{-\gamma}
 \asymp |y|^{-\gamma}.
\]

Thus a fixed renormalized outer state would have the homogeneous tail

\[
 V(\lambda y)=\lambda^{-\gamma}V(y).
\]

For the generalized Euler similarity exponents

\[
 a_*=\frac{\gamma}{1+\gamma},\qquad
 b_*=\frac1{1+\gamma},
\]

the similarity linear terms cancel on this tail:

\[
 a_*V+b_*y\cdot\nabla V
 =(a_*-b_*\gamma)V=0.
\]

A strictly homogeneous stationary outer tail must therefore solve the
stationary Euler equations on \(\mathbb R^3\setminus\{0\}\).

Shvydkoy proves that there are no \(C^1\)-smooth **axisymmetric**
homogeneous stationary Euler fields of degree \(-\gamma\) for
\(0<\gamma<2\).  This covers the entire cascade window here.  The same paper
conjectures much broader nonexistence for smooth homogeneous fields, apart
from irrotational fields at special integer homogeneities, but that broader
statement is not a theorem.

Consequently an axisymmetric strictly homogeneous wake ansatz is closed.
A survivor must use at least one of:

1. a genuinely three-dimensional, non-axisymmetric homogeneous tail not
   covered by the theorem;
2. a log-periodic/DSS tail, for which a logarithmic radial derivative
   survives in the similarity equation;
3. a nonhomogeneous discrete annular wake whose shells never converge to
   one stationary spherical profile.

The tail is not in unweighted \(L^2_y\):

\[
 \int_{1<|y|<R}|V(y)|^2\,dy
 \asymp R^{3-2\gamma},
 \qquad \gamma<3/2.
\]

This does not contradict finite physical energy.  Every physical cascade
is truncated at each preterminal time, and its shell energies
\(\ell_i^{3-2\gamma}\) form a summable series toward small scales.  Moreover
the local velocity-gradient and pressure-jet contributions of distant
annuli can converge even while the normalized \(L^2_y\) norm diverges.
The correct fixed-point space is therefore a weighted annular space, not
finite-energy \(L^2_y\).

### 8.2 Gavrilov localization does not give a generic steady parent

The 2026 theorem of Peralta-Salas and Slobodeanu shows, under its stated
analyticity, regular-boundary, and localizability hypotheses, that a
localizable steady Euler flow in a bounded domain is axisymmetric and the
domain is rotationally symmetric toroidal (or a toroidal annulus) with
convex transverse boundary.

This theorem does not exclude every localized steady Euler field, nor does
it prove that a biaxial affine jet cannot occur at one point of a toroidal
flow.  It does show that Gavrilov's localizable compact-support technology
cannot be treated as a free device for realizing an arbitrary localized
three-dimensional strain/reset geometry.  The most plausible survivor is
therefore a time-dependent or non-localizable parent, unless one proves that
the required Kelvin drain and endpoint submersion occur inside the permitted
axisymmetric toroidal class.

## 9. The full two-endpoint right inverse is impossible

This point is independent of the fluid-specific algebra.

Let \(X\) be an infinite-dimensional phase space, and let
\(\mathcal U(t,s)\) be the evolution family of the linearized Euler equation
about a smooth finite-time trajectory.  For

\[
 \mathcal Lh=\partial_t h+A(t)h=g,\qquad h(0)=0,
\]

Duhamel's formula gives

\[
 h(S)=\int_0^S\mathcal U(S,s)g(s)\,ds.
\]

If one also fixes the complete outgoing endpoint \(h(S)=0\), the range of
\(\mathcal L\) consists only of forcings satisfying the \(X\)-valued
compatibility condition

\[
 \int_0^S\mathcal U(S,s)g(s)\,ds=0.
\]

This is an infinite-dimensional cokernel.  Removing translations,
rotations, phase, energy, helicity, and finitely many other neutral
directions does not make \(\mathcal L\) onto.

Therefore a hypothesis that fixes the complete incoming packet-plus-wake
and the complete outgoing packet-plus-wake, while asking for a right inverse
on generic spacetime residuals, cannot be correct as stated.

## 10. Correct endpoint-submersion target

Let \(Z\subset X\) be the allowed incoming seed/modulation space, and let
\(\Pi:X\to Y_{\rm child}\) retain only the outgoing child constraints.  The
outgoing complement is declared to be the updated wake.  The relevant
linear map is

\[
 \Gamma
 :=
 \Pi\mathcal U(S,0)\big|_Z.
\]

The correct hypothesis is:

\[
 \Gamma:Z\longrightarrow Y_{\rm child}
\quad\hbox{is onto and has a scale-uniform tame right inverse.}
\]

Indeed, for

\[
 \mathcal Lh=g,\qquad h(0)=z,
\]

one has

\[
 \Pi h(S)
 =
 \Gamma z+
 \Pi\int_0^S\mathcal U(S,s)g(s)\,ds.
\]

A prescribed child correction can be attained by choosing \(z\) through a
right inverse of \(\Gamma\); the unconstrained part
\((I-\Pi)h(S)\) becomes wake.

This is the right formulation for the all-order viscous equations

\[
 \mathcal LV_n
 =
 \Delta V_{n-1}
 -\sum_{p+q=n}
   \mathbb P\operatorname{div}(V_p\otimes V_q).
\]

If the endpoint right inverse and the nonlinear estimates have Gevrey-tame
bounds, stage order \(M_j\asymp j\) converts
\(\varepsilon_j^{M_j}\) into \(e^{-cj^2}\), which is sufficient for a
terminally flat force.  Odd symmetry keeps every order helicity-neutral.

There remains a global compatibility: the chosen incoming correction
\(z_j\) must itself be the child/wake output of stage \(j-1\).  Thus even the
projected endpoint submersion is only a local theorem; closing the infinite
renormalized wake chain is an additional fixed-point problem.

## 11. Sharpened candidate theorem

A useful breakthrough theorem would now have the following form.

There exist \(r>1\), \(\gamma\in(1,3/2)\), an odd Euler stage trajectory
\(V_0(\sigma)\), and a weighted packet-plus-annular-wake state space, with
finite energy for every physical truncated cascade, such that:

1. in an active core, \(V_0\) contains the affine strain \(S\) and the two
   vector-potential Kelvin families of Section 3;
2. their time-integrated low stress has a strict component in the
   \(R_{\rm drain}\) direction;
3. the outgoing child is a smaller strain packet plus its next Kelvin seed;
4. the outgoing complement contains the old strain, anisotropic pressure
   tail, mediator remnants, and vorticity connectors as an admissible wake;
5. the material preimage of the selected child circulation satisfies the
   transverse-width/strand condition of Section 5;
6. the projected endpoint map \(\Gamma\), not the full two-endpoint
   operator, has a uniform Gevrey-tame right inverse;
7. the renormalized infinite wake has a fixed point in weighted annular
   seminorms, while its physical shell energies are summable; no false
   \(L^2_y\) requirement is imposed on the homogeneous outer tail;
8. either an exact Euler endpoint cell exists, or the growing-carrier WKB
   hierarchy and the viscous hierarchy are jointly solvable to order
   \(M_j\asymp j\).

The exact calculations here make items 1--2 and the helicity/circulation
part of item 5 plausible.  They do not prove items 3, 6, or 7.  No current
GPU experiment can decide those functional-analytic and Lagrangian gates.

## Primary sources

* A. D. D. Craik and W. O. Criminale, *Evolution of wavelike disturbances
  in shear flows: a class of exact solutions of the Navier--Stokes
  equations*, Proc. R. Soc. Lond. A 406 (1986), 13--26:
  https://doi.org/10.1098/rspa.1986.0061
* A. Cheskidov, M. Dai, and S. Palasek, *Instantaneous Type I blow-up and
  non-uniqueness of smooth solutions of the Navier--Stokes equations*:
  https://arxiv.org/abs/2511.09556
* A. V. Gavrilov, *A steady Euler flow with compact support*:
  https://arxiv.org/abs/1810.08020
* R. Shvydkoy, *Homogeneous solutions to the 3D Euler system*:
  https://arxiv.org/abs/1510.03378
* D. Peralta-Salas and R. Slobodeanu, *A symmetry theorem for localizable
  steady solutions of the 3D Euler equations*:
  https://arxiv.org/abs/2606.13462
