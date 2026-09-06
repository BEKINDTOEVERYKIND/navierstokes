# Navier--Stokes frontier

**Branch:** `agent/aug2-integrated-transition`
**Registry frontier:** C204 (2026-09-06 checkpoint; C189 is the ingested
auditor cross-audit)
**Boot rule:** every research session starts from this file, `CLAIMS.md`, and
the referenced artifacts on the branch. Chat history is not a premise.
Every checkpoint commit must update this file in the same commit.

## Current architecture target

Define the exact unforced nonautonomous stage-renormalization skew-product

\[
 (X,\mu,n)\longmapsto
 \left(\mathcal R_{n,\mu}X,{q_n\over g_n}\mu,n+1\right),
 \qquad
 \mathcal R_{n,\mu}=\mathcal C_{{\rm exit},n}\circ
 S^{NS,\mu}_{T_n},
\]

on one complete structured state space containing the active packet, carried
wake, phase/modulation variables, the normalized-viscosity coordinate, and
the exit chart.  Construct a nonautonomous invariant trapping tube with one
expanding direction.  An autonomous inviscid-face fixed point is a further
conditional target only after proving a normalized limit
\(\mathcal R_{n,0}\to\mathcal R_{\infty,0}\).  C188 proves that a fixed
supercritical scaling \(g>q\) cannot have an autonomous fixed point at
positive normalized viscosity because \(\mu'=(q/g)\mu<\mu\); the actual
factorial map has no finite-\(n\) fixed point.
Leak, pressure, and wake channels are components of the single residual
\(\mathcal R_{n,\mu}(X_n)-X_{n+1}\); they are not independent modules to be
closed one at a time.

The unmatched leading object is the **unforced viscous physical-velocity
scale return (UVSR)**:

> a nonzero smooth boundaryless three-dimensional unforced Navier--Stokes
> orbit/profile whose renormalized endpoint is a smaller-copy structured
> state and whose physical velocity normalization realizes a same-energy
> focus \(F\) and net endpoint gain \(g=bF\), with the admitted wake
> included in the state.  In the equal-shape ledger
> \(F=q^{3/2}\) and \(g=q^\gamma\); C188 records the explicit correction
> for unequal normalized \(L^2\) shapes.

UVSR is a demand-side property, not a unique implementation.  The phrase
“one realistic path left” is withdrawn.  The active program now keeps four
separately falsifiable mechanisms in parallel: the \(A_2\) same-witness
cascade, a direct full-wake Navier--Stokes invariant graph, a
nonaxisymmetric Euler-dominant Type-II relative profile, and an all-order
terminally-flat forced construction.  Their exact first discriminators are
recorded below; failure of one does not imply exhaustiveness of the others.

C188 replaces the hard-coded C127 net exponent \(\gamma=5/4\).  If the
bounded-profile C176 worst-case upper envelope
\(C(1+\Lambda)q^{-1}J^{7/2}/b^3\to0\) remains in the specification, its
exact polynomial window is \(7/6<\gamma<3/2\), with a sufficient
slowly-diverging logarithmic correction at the lower boundary.  This is
not a physical necessity theorem because C176 supplies no matching lower
bound.  Here and below \(C_{\rm col}\) absorbs the declared uniform bound
on \(1+\Lambda_j\).  These intervals describe normalized pure powers.  The
exact equal-shape all-sequence demand specifications are
\[
 g=q\rho,\quad 1<\rho=bq^{1/2}<q^{1/2},\quad
 {\log\rho\over\log q}\longrightarrow0
 \qquad\hbox{(direct class)},
\]
and, with \(K_*=C_{\rm col}C_J^{7/2}\),
\[
 g=K_*^{1/3}q^{7/6}(\log q)^{7/6}\omega,\quad
 \omega\longrightarrow\infty,\quad
 {\log\omega\over\log q}\longrightarrow0,\quad b<1
 \qquad\hbox{(declared C176 envelope)}.
\]
Both sequence specifications must also satisfy the exact global energy
condition (2.0a) in C188.  The explicit boundary witnesses there do.
In particular \(q=n^4\), \(g=2q\), \(b=2n^{-2}\) doubles Reynolds at
every return while attaining direct polynomial order one; it is a scalar
witness and does not reuse C180's \(q=n^8\) shell theorem.
One convenient lower-demand scalar schedule which keeps C180's
proved \(q=n^8\) shell is
\[
 b=n^{-5/2},\qquad g=n^{19/2}=q^{19/16}.
\]
It requires C161's explicit respecification
\(H=n^{51/2}\), \(\varepsilon=n^{-28}\), and
\(J_{\rm split}=\lceil n^{5/2}\rceil\), as recorded and checked in C188.
If a direct full-residual UVSR certificate cancels or absorbs the separate
C176 collar term, the broader normalized pure-power window is
\(1<\gamma<3/2\), with the direct all-sequence boundary stated above.
Changing the equal-shape \(q^{3/2}\) power would require nonuniform profile
constants or an explicit energy-disposal/wake ledger.

The fixed-point theorem itself is not the novel component. Once a UVSR
approximate profile with a certified full residual exists, the intended
imports are ABC's full-operator spectral splitting and Duhamel contraction,
Elgindi's symmetry-adapted coercivity and compactness, and Chen--Hou's
coercive-bulk/validated-finite-rank trapping method. The detailed primary-
source map is recorded in
`research/2026-08-23-fixed-point-literature-map-c185.md`.

## C188 pre-registration and verdict

The following outcome rule was fixed before the C188 checker was run.

1. Optimize the net demand over a declared schedule class. If no schedule
   meets all imported scalar constraints, issue an architecture-level
   no-go and stop UVSR. If the class is nonempty, its exact exponent window
   becomes the UVSR specification.
2. Test that window against CKN/Lin epsilon regularity, the local energy
   inequality, and Tao's quantitative critical-norm rate. If those theorems
   make the corridor empty under their actual hypotheses, issue the
   architecture-level no-go. Otherwise retain the exact surviving window.

**Verdict:** nonempty at the landed scalar-ledger level.  The
bounded-profile legacy-collar worst-case-envelope window is
\((7/6,3/2)\); the direct-full-residual window is \((1,3/2)\).  These are
the normalized pure-power windows; the exact decorated boundary
specifications are displayed above.  No further
numerical power restriction follows from CKN/Tao without the open
singular-center, time-occupancy, pressure, and wake hypotheses.  The
proposed inference that every smooth intermediate stage center must remain
above one epsilon threshold is false: C188 gives explicit \(r^2,r^3,r^4\)
upper bounds there.  This does not control backward cylinders centered at
an eventual singular point.  No criterion was weakened after seeing the
result.

## C190 pre-registered PPRG verdict

The C186 rotating-gradient orbit has now completed outcome (b) of the
unchanged witness test.  Along its fixed origin, any two consecutive maps
between three distinct endpoint times can both have determinant one only
when the conserved vertical charge is \(m=0\); this includes the two
quarter episodes.  A genuine \(2\pi\) coefficient return also forces
\(m=0\).  For the two equal quarter episodes, in the canonical common
co-rotating orthonormal fiber frame,

\[
 \Phi_1=\Phi_2=I+{\pi\over2}\beta E_{12},\quad |\beta|\le1,\qquad
 (\Phi_i-I)^2=0,\qquad
 \operatorname {tr}(\Phi_2\Phi_1)=2.
\]

Thus the ordered \(1/100\) boxes about \(I+E_{12}\) and \(I+E_{21}\)
cannot both hold, and the direct \(2+\delta\) branch fails for every
\(\delta>0\).  The genuine full-return block
\(S=I+2\pi\beta E_{12}\) obeys the explicit bound

\[
                         \|S^N\|_2\le1+{44\over7}N.
\]

This is a theorem-grade obstruction for the chosen orbit, not a failed
float search and not a no-go for the complete admissible Kelvin class.
Accordingly it feeds but does not fire the architecture trigger.  No
criterion was weakened and no successor gate was created.

## C191 C185-deficit verdict and class-scope correction

C185's allowed operator-norm floor has enough abstract scalar exponent for
the raw C182 power.  On C188's \(q=n^8\) schedule,

\[
 R_\Delta=\left\lceil{15\over8}\log q\right\rceil
 \quad\Longrightarrow\quad
 \|G_{R_\Delta T}\|_{2\to2}\ge q^{3/8},
\]

and the required inertial pump action is enclosed explicitly by

\[
 {45\over8}\log q<TR_\Delta
 <{57\over10}\log q+{76\over25}.
\]

If the same floor is also assigned to C188's formal preparation factor,
\[
R_*=\left\lceil{285\over16}\log q\right\rceil
\quad\Longrightarrow\quad
\|G_{R_*T}\|_{2\to2}\ge Hq^{3/8}.
\]
Conditionally on C176/C188's declared collar constants, this longer clock
has nonviscous normalized collar at most
\(C_{\rm col}856^{7/2}n^{-1/4}\) and viscous-collar factor at most
\(C_{\nu{\rm col}}\nu856^{3/2}(j!)^{-3/2}n^{437/28}\).  The explicit
C191 heat condition retains at least \(99/100\) of the pump and charges at
most a \(100/99\) factor in normalized time; on C188's physical scale,
\[
(j!)^{-35/2}TR_*\le t_*\le{100\over99}(j!)^{-35/2}TR_*.
\]

This does not close the physical endpoint.  The corpus supplies no lower
stage-action theorem guaranteeing those returns.  C182 is an
\(L^\infty\) upper bound on one packet class, while C185 is an unrestricted
supremal \(L^2\) operator lower bound with no common finite-\(q\) retained
band.  Treating \(q^{3/8}\) as a common scalar multiplier overruns child
energy by \(q^{3/4}=n^6\); dividing the entrance seed by \(q^{3/8}\) restores
the final energy and cancels the point multiplier.  Under uniform
rescaling, C147's coherent-\(q^3\)-packet writer diagnostic remains
divergent, \((57/2)n^8\log n\), but its transfer to C188/C161's displayed
\(q^2\)-source coefficient is not proved.  The exact reservoir
specification is therefore coherent fixed-final-energy
\(L^\infty/L^2\) concentration by \(q^{3/8}\) on one retained-band witness,
with window, viscosity, active-retention, depletion, and wake losses
charged explicitly.  Bare \(L^2\) growth receives zero focus credit.

C191 also corrects the proposed class dichotomy.  C152's background is
itself passive 2D3C:
\[
U=N\times\nabla f-\sqrt2fN=v_h+\Theta n,\qquad
n=N/\sqrt3,\quad\Theta=-\sqrt6f,\quad v_h\cdot\nabla\Theta=0.
\]
C159/C183/C184/C185 therefore already supplies an \(m\ne0\) exact return
and growth inside the broad C179/C183 class.  Universal secular no-return
is false on that class.  Calling this member the active pump rather than an
auxiliary reservoir is not a mathematical class exclusion, and no narrower
class is invented here.  The two requested forms are accepted program
verdict forms, not an exhaustive theorem about the landed class.  Form (i)
is unavailable; no third form is introduced.

## C192 strong-growth verdict

The outward C159 tube contains much more quantitative gain than its former
constant-subsolution floor recorded.  A 2048-cell directed Metzler
comparison, with no unstable amplitude-column integration, proves

\[
 M(1,3/20)^T>3000(1,3/20)^T,
 \qquad \rho(M)>3000>e^8.
\]

The C192 adversarial pass also repaired C159's inherited Decimal
sign/absolute-value context defect and reran both complete certificates;
the displayed lower bound is post-repair.

Using only C189's approved abstract spectral-inclusion step, the robust
infinite-dimensional conclusion is now

\[
             \|G_{rT}\|_{L^2\to L^2}\ge3000^r>e^{8r}.
\]

No essential-radius form is used.  On \(q=n^8\), the raw deficit clock and
action are

\[
 R_\Delta^{(192)}=\left\lceil{3\over64}\log q\right\rceil,
 \qquad
 {9\over64}\log q<TR_\Delta^{(192)}
 <{57\over400}\log q+{76\over25}.
\]

If the same floor is assigned to the formal \(Hq^{3/8}\) factor, then

\[
 R_*^{(192)}=\left\lceil{57\over128}\log q\right\rceil,
 \qquad
 TR_*^{(192)}<{1083\over800}\log q+{76\over25}.
\]

Both logarithmic action coefficients are exactly one fortieth of C191's
allowed-floor coefficients. Under the same conditional C176/C188 scalar
ledger, the scalar collar bound improves from \(856^{7/2}\) to
\(23^{7/2}\), with the corresponding \(23^{3/2}\) viscous-collar
constant.

This removes the former **power-level** objection to a first-order
finite-frequency bridge.  On the raw clock, a remainder of entrance-unit
size \(C_0(1+t)^Pq^{-1/2}e^{\Gamma t}\) is smaller than the
\(q^{3/8}\) signal whenever

\[
                         \Gamma<{350\over57};
\]

an order-\(-1\) remainder has the weaker threshold
\(\Gamma<550/57\).  These were conditional arithmetic thresholds in C192.
C194 now lands the first local pressure-resolved remainder with exponent
\(6\), but not its periodic/off-ray/same-witness completion.  A narrow
C159 Gaussian at this scale has spatial
width \(q^{-1/4}\), not the \(q^{-1}\) child width; squeezing it to the
child scale destroys the narrow-cone premise.  Bare scalar growth also
still cancels under fixed-energy normalization.  Thus C192 sharply moves
the frontier but does not overturn C191's endpoint verdict.

## C193 fixed-energy filter verdict

C193 certifies the contracting line of the same C159 monodromy rather than
treating its \(>3000\) multiplier as a scalar.  In the returning coefficient
frame the two eigenlines have slopes \(\pm x\), with
\[
 \frac{13}{100}<x<\frac15.
\]
After conversion from the orthogonal C159 frame to a physical orthonormal
frame, their slopes obey \(1/2<y=x|k_0|<2\).  Hence the physical spectral
projectors have norm below \(5/4\) and the eigenbasis condition number is
below \(2\).  The same claim proves the sharp global \(A_2\) jets
\[
 \|DU\|_{\rm op}\le6,\qquad
 \|D^2U\|_{\rm mult}\le3\sqrt6,\qquad
 \|D^3U\|_{\rm mult}\le9.
\]

In the exact complex principal two-fiber model, with
\[
 R_{\rm filt}=\left\lceil\frac38\log n\right\rceil+1,
 \qquad q=n^8,\qquad G=\rho^{R_{\rm filt}}>3000n^3,
\]
C193 gives an explicit localized expanding profile and broad contracting
profile whose entrance and endpoint \(L^2\) norms both equal \(b\), with no
overshoot at any discrete return, and
\[
 \frac{\mathcal C(v_{R_{\rm filt}})}{\mathcal C(v_0)}
 >\frac{1750}{251}n^3>q^{3/8},
 \qquad \mathcal C(v)=\frac{\|v\|_\infty}{\|v\|_2}.
\]
This positively resolves C191's finite-dimensional scalar-normalization
collision.  It does not control continuous within-period energy and is not
yet a real, divergence-free, off-ray, periodic finite-frequency witness.

## C194 pressure-resolved local bridge verdict

C194 constructs the exact curl parametrix
\[
 v_{\rm app}=\frac{\hbar}{i}\operatorname {curl}(e^{i\phi/\hbar}c),
 \quad c=-\frac{k\times b}{|k|^2},\quad
 p_{\rm app}=\hbar e^{i\phi/\hbar}\,2i\frac{k^TAb}{|k|^2},
\]
on the \(\mathbb R^3\) \(A_2\) invariant annulus \(|f|\le1/10\).  It is
exactly divergence free, resolves pressure explicitly, and has the exact
residual identity
\[
 (D_t+DU)v_{\rm app}+\nabla p_{\rm app}
 =\hbar e^{i\phi/\hbar}(D_td+DU\,d+\nabla\pi).
\]
For the exact linearized-Euler solution with the same initial data,
\[
\begin{aligned}
\|v-v_{\rm app}\|_2
\le\frac{\hbar e^{6t}|b_*|}{|k_0|}\big[&
4{,}199{,}040\,\varepsilon^{-1}(1+t)^3\|\nabla\chi\|_2\\
&+2{,}898{,}006{,}000{,}000{,}000(1+t)^7\|\chi\|_2\big].
\end{aligned}
\]
Thus the full local symbol/pressure error has \(\Gamma=6<350/57\).  With
\(\hbar=q^{-1}\) and the C193 concentration width
\(\varepsilon=q^{-1/4}\), its first term has relative power margin
\(q^{-27/100}\) against a hypothetical same-packet \(q^{3/8}\) signal;
the narrower C192 stress test \(\varepsilon=q^{-1/2}\) retains margin
\(q^{-1/50}\).  C193's extra return changes only the fixed clock factor
from \(e^{456/25}\) to \(e^{912/25}\).

C194 supplies the upper-error half only.  It does not prove off-ray
expansion/stable transport, an integer periodic phase, reality completion,
initial or endpoint band retention, the C193 concentration endpoint for
the same exact solution, viscosity, or nonlinear closure.

## C195 finite-horizon dominated-cone verdict

C195 closes the finite-horizon \(C^0\) off-ray dominated-cone problem
requested after C194.  The repaired off-level Kelvin block retains the skew term
\[
              \frac{m(t\cdot Sp-p\cdot St)}{QD}
              =\frac{2\sqrt3\,mf}{Q},
\]
and carries \(|f|\le10^{-11}\).  On the fat 2048-cell tube, both the forward
and reflected-inverse coefficient cones gain by more than \(3000\).  Their
complete slope images lie respectively in
\[
 [0.1405737,0.1898012],\qquad[0.1402257,0.1903237],
\]
any forward-cone line and sign-reflected backward-cone line have physical
angle sine at least \(520/569\), their two-line oblique projector norms are
at most \(569/520\), and each block's projective slope contraction is below
\(1/4{,}500{,}000\).  Continuous one-block physical propagation is below
\(3^{55}\).

On a common Euclidean lift of the initial torus neighborhood, set
\(r=\max\{\|\delta X_0\|_2,\|\delta k_0\|_2\}\). The honest finite-horizon
closeness condition is
\[
 r(1+TR_{\rm filt})^3\le1/(8.7\cdot10^{13}).
\]
For \(q=n^8\), it certifies \(r=q^{-1/4}\) at \(n\ge10^{10}\),
\(r=q^{-1/12}\) at \(n\ge10^{39}\), and \(r=q^{-1/3}\) at
\(n\ge10^8\); hence the combined \((q^{-1/12},q^{-1/3})\) two-width box is
certified at \(n\ge10^{39}\).  Cone-field and selected-line first derivatives were removed
after audit and remain open.  No invariant or canonical stable/unstable
bundle, closed phase, finite-frequency packet, viscosity, or nonlinear
return follows.

## C196 endpoint phase-space verdict

C196 gives an explicit real periodic exact-curl construction on
\(q=m^{24}\).  Its carrier grid has
\[
 M\ge q^2/13824,\qquad p_{\min}\ge10q,
\]
and the two envelope widths \(q^{1/12}\), \(q^{1/3}\) give fixed-energy
concentrations at least
\[
 c_*q^{9/8},\qquad c_*q^{3/2},\qquad
 c_*=\frac{39}{10\sqrt{13824}(1+\sqrt3/10)}.
\]
This is exact periodic/reality/solenoidality kinematics, not a common
Floquet/WKB solution.

The support count is now dimensionally explicit.  A three-coordinate
relative \(q^{-1/4}\) cube has ceiling \(q^{9/8}\); a genuinely projective
angular tube with full annular radial width has ceiling \(q^{5/4}\);
C193's one \(q^{-1/4}\) spatial envelope has concentration \(q^{3/8}\).
C180's retained one-sided box obeys
\[
 \#\Sigma\le64\delta^3q^3/J^4,\qquad
 {\cal C}(v)\le8\delta^{3/2}q^{3/2}/J^2,
\]
with at most a \(\sqrt2\) constant for a disjoint reality completion.

C194 controls one compactly supported beam, whereas C196 uses
\(M\ge q^2/13824\) global periodic carriers.  Its normalized triangle
majorant contains \(\sqrt M\), hence a full power of \(q\), destroying the
one-beam margins.  Thus a uniform
multi-beam Fourier-integral almost-orthogonality theorem plus spatial
localization and endpoint band/tail retention is load-bearing.  C196
records the single-beam powers only conditionally and makes no composition
claim.

C197 below now resolves the compact common-lattice upper-error sum. The
global periodic endpoint and its dynamical lower bounds remain open.

## C197 uniform compact multi-beam upper-error verdict

C197 removes the carrier-count loss for C194-type compact beams with a
common envelope and carriers in one separated lattice. The ambient Kelvin
generator is \((-I+2kk^T/|k|^2)DU\), so its operator norm is at most \(6\).
Four covector derivatives and the required mixed spatial derivatives cost
explicit polynomial factors in time and retain the exponent \(6\).
Parameter Fourier synthesis costs \(45\); the compact-envelope Bessel bound
costs \(27\), independently of the number of carriers.

With a real unit-ball-supported envelope \(\chi\), \(0<\varepsilon\le1/40\),
an initial zero-level center, \(d\varepsilon\ge1\), and normalized covectors
in the stated fixed cube, the exact \(\mathbb R^3\) linearized-Euler
solution with the same curl initial datum satisfies
\[
 \|u-u_{\rm app}\|_2\le\hbar e^{6t}\|a\|_{\ell^2}
 \left[10^{94}(1+t)^{37}\|\chi\|_\infty+
 10^{80}\varepsilon^{-1}(1+t)^{31}\|\nabla\chi\|_\infty\right].
\]
The constants are conservative and explicit; no \(\sqrt M\) is present.
Reality completion costs at most \(2\) relative to the positive-box
coefficient norm. This proves the compact multi-beam upper-error estimate.
On the stated C193 clock and the two C196 envelope/scaling choices, the
real error divided by the declared profile/coefficient scale is bounded by
\[
 4\cdot10^{112}(\log q)^{37}q^{-37/600}<1/100
 \quad\text{for }q\ge10^{10000}.
\]
This is a sufficient finite-scale budget, not physical-energy or signal
normalization.
It does not prove a common expanding aperture, compare coefficient norm
to the physical entrance norm from below, realize C196's global periodic
endpoint, retain the output band, or control viscosity.

## C198 full-wake energy and circulation verdict

For exact full-state returns under the already proposed free-space chart
and smooth unforced stages,
write \(E_n=\|X_n\|_2^2\), \(k=q^3/g^2>1\), and
\(D_n=2\mu_n\int_0^{T_n}\|\nabla u_n\|_2^2dt\). C198 proves
\[
 E_{n+1}=k(E_n-D_n),\qquad
 D_n\le\sqrt{\mu_n}\int_0^{T_n}
                  \|u_n(t)\|_{\mathcal X_{\mu_n}}^2dt.
\]
Full energies bounded above and away from zero cannot coexist with vanishing dissipation.
In particular, bounded stage durations and uniformly bounded full-trajectory
\(\mathcal X_\mu\) norms cannot produce this trapping graph as \(\mu\to0\).
Endpoint-only bounds and singular graphs remain open. A nonzero regular
finite-energy Euler fixed face is also excluded by energy conservation.

The three-return search has not been run or ruled out: finite-stage
normalized energy growth is permitted. C198 provides exact measured-energy
and residual enclosures for its replay. Only for equally normalized \(L^2\)
stages with relative \(L^2\) residual at most \(10^{-5}\), each stage must
dissipate more than \(69/10000\) of its energy. The axisymmetric circulation
contracts by at most \(971/1000\) per exact return, hence by at most
\(915498611/10^9\) over three returns; weighted residual corrections are
explicit.

## C199 full linearized-PDE inverse-tail verdict

For the full linearized Navier--Stokes evolution \(S\) on the mean-zero
periodic \(L^2\) space, with viscosity \(\mu>0\), velocity bound \(M_0\),
and symmetric-gradient bound \(M_1\),
\[
 \|SQ_K\|\le\frac{e^{M_1T}}K
       \sqrt{\frac1{\mu T}+\frac{4M_0^2}{\mu^2}}.
\]
This positive PDE tail estimate follows by retaining weighted
\(H^{-1}\) dissipation. For a bounded same-space chart \(C\), put
\[
 M=\|C\|e^{M_1T},\quad
 \delta_K=\frac M K\sqrt{\frac1{\mu T}+\frac{4M_0^2}{\mu^2}},
 \quad \Lambda=1+M/s.
\]
A certified \(s>0\) for the smallest singular value of the full compression
\(P_K(I-CS)P_K\), with \(\delta_K\Lambda<1\), gives
\[
 \|(I-CS)^{-1}\|\le\frac{\Lambda}{1-\delta_K\Lambda}.
\]
The full-mode example \(\|(I-3e^\Delta)^{-1}\|_{H^r\to H^r}\le11\)
shows why backward heat alone does not reject this return derivative.
It is an example with an artificial chart, not a nonzero stage candidate.
No actual stage finite block, viscosity-uniform inverse, nonlinear \(L^2\)
endpoint map, or nonautonomous graph inverse is certified.

## C200--C204 same-solution viscous endpoint verdict

C204 completes the finite-frequency linear concentration discriminator on
the periodic A2 pump. For
\[
 n\ge\max(10^{12500},\nu),\quad q=n^8,\quad
 \mu=\nu((n-1)!)^{-3/2},\quad
 R=\lceil(3/64)\log q\rceil+1,\quad L=1+RT,
\]
there is one smooth real mean-zero divergence-free solution of the actual
linearized Navier--Stokes equation around the unforced heat-decaying A2
pump, at its actual action-return time \(S\), with
\[
 \|v(0)\|_2=\|v(S)\|_2=b,\qquad
 \|P_{\mathcal B_{\mathbb R}}v(S)\|_\infty
 \ge10^{-157}b\,q^{3/2}/L^2,
\]
\[
 \frac{\|v(S)\|_\infty/\|v(S)\|_2}
      {\|v(0)\|_\infty/\|v(0)\|_2}
 >\frac{25}{24}q^{3/8}.
\]
The same absolute endpoint bound holds without projection. The retained
band has explicit half-widths \(8\cdot10^{-100}q/L\),
\(8\cdot10^{-100}q/L\), \(8\cdot10^{-100}q/L^2\), centered at a nearest
integer to the C159 reference covector times \(q\). Its compatibility
with C180's particular center, frame, splitter, and exit chart is open.

The proof uses C200's smooth finite-horizon growing/contracting selectors,
C201's rank-one-shear-adapted integer lattice and physical Riesz bounds,
C202's periodic general-symbol viscous bridge, and C203's finite Fourier
test. The full parameter/support tube is verified together; all curl,
viscous, reality, stable-ballast and physical-normalization losses are
charged on the same solution. C185's supremal operator norm is not
multiplied into a packet bound.

The linear perturbation's budgets are also explicit:
\[
 \|v(0)\|_\infty\le20000bq^{9/8},\quad
 \sup_s\|v(s)\|_2<10^{16}b,\quad
 \int_0^S\|v(s)\|_2^2ds<4\cdot10^{33}b^2,\quad
 \mu\int_0^S\|\nabla v(s)\|_2^2ds<2\cdot10^{34}b^2.
\]
The constants and sufficient scale are intentionally conservative and
have not been numerically instantiated.

**Current A2 target:** evaluate and control the full nonlinear defect
\(P(v\cdot\nabla v)\), the pump response and retained wake, and the actual
exit-chart return for this constructed family. This is the existing
nonlinear UVSR obligation, not another PPRG gate. C204 proves a linear
concentration endpoint and a genuine same-energy gain; it does not prove
a nonlinear stage, a repeated return, or C182's stronger entrance
\(\|v(0)\|_\infty\le bq\). No independent passive reservoir is used in
the proved linear concentration step.

## Active research portfolio

No item below is described as the globally “last path.”

1. **\(A_2\) same-witness cascade.** C204 has now proved the real periodic
   viscous linear endpoint, physical energy equality and \(L^{-2}\)-taxed
   retained concentration for one solution, with its actual clock and
   full linear trajectory charged. Attack its full nonlinear defect and
   pump/wake response at the stage amplitude, then test the actual exit
   chart. The constructed frequency center must be matched explicitly;
   no automatic C180 seed/splitter compatibility is credited.
2. **Direct full-wake Navier--Stokes graph.**  Search
   \[
   {\cal C}S_{T(\mu)}^\mu X(\mu)-X((q/g)\mu)=0
   \]
   in axisymmetric-with-swirl free-space variables with the entire endpoint
   retained.  First run: \(q=6/5\), \(7/6<\log g/\log q<1.48\), three
   consecutive returns, \(128\times256\) then \(192\times384\), and
   12--20 restarts.  Candidate thresholds are relative full residual
   \(<10^{-5}\), independent replay \(<10^{-4}\), and top-third/outer-collar
   fractions each \(<10^{-7}\).  Boundary-hitting and one-component
   collapse are failed searches. C198 requires measured full-energy and
   circulation residual ledgers on each replay. A bounded-time trapping
   tube with uniformly bounded full-trajectory \(\mathcal X_\mu\) norm and
   \(L^2\) energy bounded above and away from zero is excluded as \(\mu\to0\). The finite
   three-return test remains meaningful, but an infinite graph must resolve
   the required transient norm growth or increasing stage durations.
3. **Euler-dominant Type-II relative profile.**  At
   \((\alpha,\beta)=(11/20,9/20)\), search a genuinely nonaxisymmetric,
   non-outgoing rotating inner profile with a \(|y|^{-11/9}\)
   logarithmically twisted large-similarity-radius law.  Weighted residual
   \(<10^{-10}\), centered-ball flux ratio \(1/8\) at three radii, a stable
   tail law, and a resolution-stable finite bad spectrum qualify only as an
   inner-profile discovery candidate.  Promotion requires a certified
   time-dependent finite-energy cutoff/outer matching, the complete induced
   matching defect, and explicit pressure/Biot--Savart coupling to the outer
   flow.  For the unforced route that defect must be cancelled dynamically,
   not declared as forcing.  CIV/Chae--Wolf guardrails are source-scoped;
   Pineau--Vicol is methodological context only.
4. **All-order forced route.**  Test finite sections of the full
   retained-wake endpoint derivative before any profile search.
   C199 proves an explicit analytic high-frequency tail for the full
   linearized PDE and a finite-block inverse criterion for \(I-CS\).
   Backward-heat growth of \(S^{-1}\) alone does not decide that return
   derivative. Any singular-value rejection must use the actual derivative
   of the prescribed-endpoint, same-space return, or nonautonomous graph
   equation being tested. A uniform
   \(C^M(M!)^2\)-type inverse growth bound warrants a Nash--Moser/Borel
   continuation.  The final defect must be Clay-admissible and terminally
   flat.

The detailed equations and failure criteria are in
research/2026-08-30-counterfactual-success-portfolio.md.

## Proved frontier

- C200--C204: smooth finite-horizon selectors and compact periodic
  synthesis now compose into an actual linearized-Navier--Stokes
  retained endpoint with equal physical endpoint norms, absolute
  \(10^{-157}bq^{3/2}/L^2\) concentration and concentration gain
  \(>(25/24)q^{3/8}\), under the displayed finite-scale threshold.
  Explicit linear trajectory/action/dissipation budgets are proved;
  the nonlinear defect and return remain open.

- C121: the homochiral \(A_2\) pump is an exact unforced heat-decaying
  Navier--Stokes background.
- C152/C159/C192: one exact zero-drift periodic orbit has a rigorously
  certified returning \(m\ne0\) Kelvin cone, and the sharpened directed
  comparison gives one-period multiplier strictly larger than \(3000\).
- C161/C176/C179/C180: the required support cardinality, correlated tube,
  exact passive 2D3C background, and factorial-shell tight-star arithmetic
  exist, but no physical velocity endpoint follows from them.
- C181: a fixed-cone static vertical shear has an explicit uniform bound.
- C182: under its stated entrance regularity hypotheses, fixed-time endpoint
  growth is at most \(bq^{9/8}\), short of \(bq^{3/2}\).
- C183: the exact passive 2D3C Lagrangian gauge gives the return-resonance and
  linear/quadratic covector-drift formulas; a common-frame
  \(O(q^{-1/2})\) perturbation gives no extra \(\sqrt q\).
- C184: degree-\(q\) stationary polynomial palettes pay an exponential
  off-line cost and cannot meet a polynomial reservoir budget.
- C185/C189/C192: the certified C159 multiplier plus the abstract-level
  Shvydkoy spectral inclusion gives
  the robust infinite-dimensional operator-norm estimate
  \(\|G_{nT}\|_{L^2\to L^2}\ge3000^n>e^{8n}\).  The
  essential-spectral-radius form is citation-held until the paper-body
  Theorem 4.1 check lands.  C189 independently confirms C186 and C187 in
  full and C185 with precisely this split.
- C186: conditionally, candidate unipotent polarization blocks do not
  algebraically force a common flag. The exact pair
  \(I+E_{12},I+E_{21}\) has per-episode exponent \(>12/25\), robustly
  \(>9/20\) in its entrywise \(1/100\) boxes. For separately established
  square-zero generators, the exact decision scalar is
  \(\tau=\operatorname{tr}(N_1N_2)\). A bounded exact passive 2D3C flow can
  rotate its scalar gradient through transverse lines, but this does not
  determine the Kelvin polarization flags.
- C187: ordinary finite-stage linearized Navier--Stokes Duhamel continuity
  holds in \(H^3\) with the explicit constant
  \(44\sqrt{2T/\nu}\exp(7744V^2T/\nu)\). It is not scale-uniform.
- C188: equal-normalized-shape same-energy transfer into a \(q^{-3}\) child
  has ledger focus \(q^{3/2}\); unequal shapes contribute the explicit
  square root of their \(L^2\)-constant ratio.  The net gain has sharp
  polynomial infimum \(7/6\) within the declared bounded-profile C176
  worst-case envelope, and infimum \(1\) for a direct full-residual UVSR.
  The direct floor has the exact scalar witness
  \(q=n^4,g=2q,b=2n^{-2}\), with total energy at most \(16/3\) in its
  displayed normalization.
  The exact fixed-\(q=n^8\) scalar schedule \(b=n^{-5/2}\) lowers \(5/4\)
  to \(19/16\), with C161 seed/split normalization respecified explicitly.
  At a smooth stage center, the standard scaled energy/cubic/enstrophy
  quantities are bounded by explicit constants times \(r^2,r^3,r^4\); no
  numerical CKN/Tao power follows from the scalar ledger alone.  Also
  \(\mu'=(q/g)\mu\), and the factorial schedule is an augmented
  nonautonomous map, not an autonomous fixed point.
- C190: on the explicit C186 rotating-gradient orbit, the determinant-one
  gate for any two consecutive episodes forces \(m=0\).  In the canonical
  common co-rotating frame both maps equal
  \(I+(\pi/2)\beta E_{12}\), \(|\beta|\le1\); their product trace is
  exactly two, their
  square-zero generators share a flag, and the full-return powers satisfy
  \(\|S^N\|_2\le1+(44/7)N\).  This is the pre-registered orbit-specific
  obstruction, not the complete-class PPRG verdict.
- C191: the C185 floor supplies the raw scalar \(q^{3/8}\) after
  \(R_\Delta=\left\lceil(15/8)\log q\right\rceil\) returns and is compatible with the
  declared scalar collar/heat powers after an explicit logarithmic-window
  allocation.  No actual lower stage coverage is proved.  The landed
  unrestricted \(L^2\) operator norm cannot be multiplied by C182's
  \(L^\infty\) upper bound; common scalar reuse overruns energy by
  \(q^{3/4}\), fixed-energy normalization cancels the gain, and C185
  supplies no C125 relative-return cancellation.  Also, the
  C159/C185 returning growing covector lies in the broad passive-2D3C
  class, falsifying universal \(m\ne0\) secular no-return on that class.
- C192: a 2048-cell outward Metzler product proves
  \(Mw>3000w\), upgrades the allowed operator-norm floor to \(3000^r\),
  and divides both C191 logarithmic action coefficients by forty.  The raw
  action is below \((57/400)\log q+76/25\); consequently a first-order
  \(q^{-1/2}\) finite-frequency remainder has the explicit nonempty power
  window \(\Gamma<350/57\).  No such remainder bound or physical endpoint
  is claimed.
- C193: the exact \(A_2\) jets satisfy the sharp bounds
  \(6,3\sqrt6,9\); the central C159 stable line has physical projector norm
  below \(5/4\); and an exact principal two-polarization filter preserves
  endpoint \(L^2\), has no discrete-return overshoot, and improves
  concentration by more than \(q^{3/8}\).  It is complex, central-fiber,
  and finite-dimensional.
- C194: one exactly solenoidal pressure-resolved local \(A_2\) WKB beam
  has the explicit linearized-Euler error displayed above, with
  exponential rate \(6\) and margins \(q^{-27/100}\) at width \(q^{-1/4}\)
  and \(q^{-1/50}\) at width \(q^{-1/2}\), relative to a hypothetical
  same-packet signal.  Periodicity, off-ray growth, band, and endpoint
  concentration are not included.
- C195: an explicit fat off-ray finite-horizon tube carries strictly mapped
  forward and sign-reflected inverse coefficient cones.  Both cone gains
  exceed \(3000\); every forward/backward cone-line pair has two-line
  oblique projector norms at most \(569/520\); each block's slope contraction
  is below \(1/4{,}500{,}000\); and the exact initial-radius criterion is
  \(r(1+TR_{\rm filt})^3\le1/(8.7\cdot10^{13})\).  This is a finite-horizon
  \(C^0\) dominated-cone field, not an invariant or canonical bundle;
  cone-field derivatives and finite-frequency realization remain open.
- C196: exact real periodic divergence-free fixed-aperture profiles attain
  explicit fixed-energy concentration lower bounds \(c_*q^{9/8}\) and
  \(c_*q^{3/2}\).  Endpoint support in C180's retained box imposes the exact
  upper tax \(8\delta^{3/2}q^{3/2}/J^2\) one-sided.  A projective angular
  \(q^{-1/4}\) tube has ceiling \(q^{5/4}\), while a three-coordinate tube
  has ceiling \(q^{9/8}\).  C194 is one-beam only; no uniform multi-beam
  composition or dynamic endpoint follows.
- C197: the compact common-lattice multi-beam curl/WKB upper error has
  explicit constants \(10^{94},10^{80}\), time powers \(37,31\), and
  exponential rate \(6\), independently of carrier count. It is measured
  against coefficient \(\ell^2\), with envelope and separation hypotheses;
  no lower-growth, periodic endpoint, or physical normalization is proved.
- C198: full-state energy and circulation identities yield exact
  approximate-return enclosures. The constant-one estimate
  \(D_n\le\sqrt{\mu_n}\int\|u_n\|_{\mathcal X_{\mu_n}}^2dt\)
  excludes uniformly bounded full trajectories with bounded times and
  energies bounded above and away from zero. Endpoint-only singular graphs remain open.
- C199: the full linearized Navier--Stokes operator has an explicit
  polynomial high-frequency tail at fixed positive viscosity. A validated
  finite compression can certify its same-space return inverse by the
  displayed criterion; the full-mode example has inverse norm at most
  \(11\). The actual stage compression is not certified.

## Open

- UVSR: no exact or rigorously trapped unforced viscous scale-return profile
  with physical velocity focus is known.
- C198 excludes the bounded-time, uniformly full-trajectory
  \(\mathcal X_\mu\)-bounded full-energy graph with energy bounded above
  and away from zero. Endpoint-only graphs require resolved transient
  norms or longer stages; neither has been constructed.
- C199 supplies a fixed-positive-viscosity tail test, but the actual
  full-stage finite compression and its singular-value enclosure remain
  open. Its same-space inverse theorem is not a theorem about the
  viscosity-shifting graph.
- PPRG realization: the broad landed passive-2D3C class contains both
  C190's secular orbit and C159/C185's exact returning growing orbit, so a
  universal secular-lock theorem is false there.  The only unresolved
  accepted witness form **within the PPRG line** is one genuinely
  incommensurate non-fixed orbit
  with an \(m\ne0\) return or transported-frame Lyapunov growth.  It must
  carry the same-witness retained-band, viscosity, fixed-final-energy
  \(q^{3/8}\) concentration, C125, depletion, and wake estimates.
- C204 completes the C193/C194 linear same-witness bridge with finite
  carrier scale, reality, periodicity, retained band, viscosity and
  physical energy normalization. Its quadratic nonlinear defect,
  pump depletion, retained wake and repeating exit-chart return are
  open. C200's smooth finite-horizon selected lines do not constitute
  an invariant or canonical infinite-time splitting.
- A complete state space and normalization for \(\mathcal R\), including the
  retained wake and normalized-viscosity coordinate, have not been fixed.
- No normalized convergence \(\mathcal R_{n,0}\to\mathcal R_{\infty,0}\)
  has been proved, so there is not yet an autonomous inviscid-face fixed-
  point problem for the actual factorial schedule.
- No coercive bulk estimate or spectral gap has yet been proved for the
  proposed renormalized Navier--Stokes operator.
- **UVSR standing terminal obligation:** terminal singular-center tracking
  has not been proved.  C188's CKN/Lin/Tao clearance is confined to the
  scalar ledger and smooth intermediate stage centers.  Any terminal UVSR
  object must still control backward cylinders at the proposed singular
  center, including time occupancy, pressure/wake contributions, the local
  energy inequality, and the relevant critical-norm comparison.  No
  residual-stable forced epsilon criterion or torus transfer of Tao's
  quantitative theorem has been proved.
- C125's relative return, the old BAFL split, and RIGM remain unproved, but
  are to be absorbed into one full residual/trapping inequality rather than
  promoted into further architecture gates.

## Pre-registered kill criteria

These verdicts are fixed before the next computation.

1. **PPRG architecture-trigger scope.** The proposed universal
   \(m\ne0\) secular-lock branch is unavailable on the broad C179/C183
   class because C159/C185 is a landed counterexample.  Architectural role
   labels do not remove it from that mathematical class, and no narrower
   reservoir class is introduced to rescue the statement.  Failure to find
   the remaining incommensurate witness is not a theorem and does not fire
   the architecture trigger.  No third verdict form or successor gate is
   permitted.
2. **Witness boundary.** An abstract matrix pair or unrestricted
   operator-norm exponent does not validate PPRG.  Promotion requires one
   exact non-fixed incommensurate passive 2D3C orbit, its physically
   transported return/growth, and on the same witness a finite-frequency
   PDE estimate proving fixed-final-energy concentration by \(q^{3/8}\)
   after every explicit loss.
   The periodic C159/C192 pump asset is subject to the same physical
   endpoint rule even though it is not the incommensurate PPRG reservoir.
3. **Fixed-time regular branch.** Any candidate satisfying C182's entrance
   hypotheses on only \(O(1)\) normalized time is rejected as a source of
   \(bq^{3/2}\) focus.
4. **Stationary palette branch.** A candidate whose outer selected
   coefficient forces super-polynomial collateral norm, as in C184, is
   rejected unless the collateral modes are explicitly part of the useful
   endpoint state.
5. **No sixth PPRG gate.** If the remaining accepted PPRG witness search
   fails, record a failed search.  Do not rename that residue, narrow its
   class ad hoc, or create another PPRG gate; this says nothing exhaustive
   about the three independent portfolio routes above.
6. **Anti-formalism checkpoint rule.** A checkpoint must deliver an
   explicit-constant estimate, a strict corridor narrowing, or a completed
   dichotomy branch. State-space/operator definitions alone are not a
   checkpoint. Numerical residual minima count only with a full interval
   certificate.
7. **Phase-space honesty.** Count distinct endpoint Fourier modes, not
   carrier--envelope labels. C204 meets the logarithmically taxed linear
   endpoint specification using C201's physical norm comparison and
   C203's actual Fourier projection. Its explicit band must not be
   identified with C180's prescribed center/chart without a proof.
   Neither the compact envelopes nor the full exact solution are
   assumed bandlimited.
8. **Direct full-wake candidate.**  The thresholds in Active research
   portfolio item 2 are fixed.  Selected-shell recurrence, \(T\to0\),
   cutoff/collar pile-up, or vanishing swirl/poloidal energy is a failed
   search, not a weakened residual target. The measured norm ratio must
   appear in C198's energy enclosure; the \(69/10000\) loss floor applies
   only to its stated equal-\(L^2\) corollary. Circulation cannot be reset to
   one at each stage.
9. **Type-II candidate.**  The inner-discovery thresholds in item 3 are
   fixed.  Axisymmetric or locally outgoing profiles, missing
   \(|y|^{-11/9}\) inner tail/twist, failure of the \(1/8\) flux identity, or
   proliferating unstable eigenvalues are failed searches.  Passing those
   tests does not advance a Clay claim without certified finite-energy outer
   matching, matching-annulus defect control, pressure/nonlocal coupling,
   and the standing terminal singular-center obligations.
10. **Forced tame inverse.** Test the derivative of the actual residual.
    Exponential inverse growth or a stable adjoint compatibility for that
    derivative kills the tested formulation. C199's same-space
    \(I-CS\) criterion cannot be substituted for a viscosity-shifting graph
    operator. Only a validated finite compression of the full evolution
    plus its analytic tail can promote the same-space inverse.

## Audit obligations

- **Resolved by C187:** C136's (5.1) is folded into the single full
  residual/trapping hypothesis. The explicit \(H^3\) substitute is a
  finite-stage continuity theorem only; no uniform structured-state
  constant is claimed.
- **Withdrawn by C187:** the session-only infinite-ladder spectral
  enclosures \(0.66855\ldots\) and \(2.63707\ldots\) have no landed
  operator, normalization, interval certificate, tail resolvent bound, or
  checker. They carry zero evidentiary weight. C120's certified finite
  \(6\times6\) enclosure is the only relevant landed eigenvalue interval.
- Every load-bearing inequality added from this checkpoint onward has an
  explicit numerical constant. Bare \(O(\cdot)\) and \(o(\cdot)\) may
  describe non-load-bearing orientation only.
- No obstruction may be the sole advance of two consecutive checkpoints.
- **Resolved by C188:** the demand-side scalar schedule optimization is
  nonempty and the CKN/Tao “every smooth stage center” premise is
  withdrawn.  The bounded-profile C176 worst-case-envelope infimum is
  \(7/6\), not the current \(5/4\); the equal-shape ledger multiplier is
  \(q^{3/2}\), with an explicit physical shape-factor correction.
- **Standing UVSR obligation after C188:** the CKN/Lin/Tao conclusion is
  scalar-ledger-scoped.  It does not clear backward-cylinder regularity at
  a terminal singular center.  Terminal center tracking, time occupancy,
  pressure/wake terms, the local energy inequality, and critical-norm
  comparison remain attached to UVSR and may not be dropped at profile
  certification.
- **C185/C192 citation restriction from C189:** until the Shvydkoy Theorem 4.1
  paper-body check lands, every downstream use cites only
  operator norm.  C192's allowed strengthened form is
  \(\|G_{nT}\|_{L^2\to L^2}\ge3000^n>e^{8n}\); the
  \(r_{\rm ess}\) form remains source-held and is not a premise.
- **Completed by C190 for the chosen orbit:** the PPRG two-box/trace test
  landed as pre-registered outcome (b), with exact square-zero and trace
  calculations.  C191's class-scope correction does not alter that
  orbit-specific theorem.  A partial float search is not a checkpoint.
- **Resolved by C191:** the C185 floor has sufficient abstract scalar
  exponent in an explicitly long logarithmic action allocation, but the
  direct C182 reconciliation is negative on the landed statements.  There
  is no lower stage-coverage theorem, no common retained-band concentration
  witness or viscous multiplier, fixed-energy scalar reuse cancels, and
  C125 remains open.  Do not credit
  C185's \(L^2\) exponent toward the raw \(q^{3/8}\) concentration demand
  before a same-witness theorem.
- **Class correction from C191:** universal \(m\ne0\) secular no-return is
  false on the broad C179/C183 passive-2D3C class because C159/C185 is a
  returning growing member.  The incommensurate transported-frame witness
  is the only unresolved accepted form.  Failed search is not a class
  theorem, and no third formulation is permitted.
- **C192 boundary:** the short clock removes the generic power mismatch but
  does not by itself prove a remainder, retained band, viscosity,
  child-scale localization, concentration, or C125.  C194 proves only the
  local upper-error half.  Do not
  multiply its operator norm into C182's upper bound or treat a
  \(q^{-1/4}\)-width Gaussian as a \(q^{-1}\)-width child.
- **C193/C194 boundary:** the fixed-energy filter is an exact complex
  principal two-fiber lemma and the WKB estimate is an \(\mathbb R^3\)
  local-beam upper bound.  Neither may be quoted as the real periodic
  off-ray same-witness endpoint.  Composite uses pay the extra-return fixed
  factor \(e^{912/25}\), not C194's raw-clock \(e^{456/25}\).
- **C195 boundary:** only a \(C^0\) finite-horizon forward/reflected-inverse
  dominated-cone field is certified.  The separation and oblique-projector
  bounds hold for arbitrary lines chosen from the two cones; they do not
  select persistent \(E^u,E^s\).  The draft cone/projector-derivative claim
  was removed after audit.  Do not promote C195 to an invariant or canonical
  bundle, a closed phase, a uniform FIO propagator, or a viscous packet.
- **C196 boundary:** the exact periodic profiles are kinematic endpoints,
  not C194 solutions.  A projective angular tube has two narrowed
  directions, not three.  Every phase-space ceiling is conditional on
  actual endpoint support, and C194's block powers remain conditional on a
  uniform multi-beam theorem missing at that checkpoint. C197 subsequently
  proves the compact common-lattice upper estimate with distinct envelope
  hypotheses; it does not complete the periodic dynamic composition.
- **Type-II promotion boundary:** the \(|y|^{-11/9}\) law belongs only to the
  inner similarity-profile discovery problem and gives a global
  infinite-energy profile.  A weighted residual, flux identity, and finite
  bad spectrum do not by themselves produce finite-energy data or a Clay
  solution.  Promotion requires a certified time-dependent cutoff/outer
  match, its full induced defect, explicit pressure/Biot--Savart coupling,
  and dynamic cancellation of the defect in the unforced equation.
- **Portfolio engine portability:** engines/return_map_opt.py now retains
  the seed-dependent time in autograd, uses vector-speed diagnostics, and
  snapshots score/state/shift atomically.  PyTorch is absent in the current
  container, so this patch is syntax-checked only and supplies no numerical
  candidate.  The historical independent QP audit remains unchanged for
  provenance.
- **C198 boundary:** use complete energies and actual norm ratios in the
  residual ledger. Endpoint \(\mathcal X_\mu\) bounds alone do not bound
  intermediate norms. The \(69/10000\) floor is an equal-\(L^2\) corollary;
  circulation needs its own weighted supremum residual. No general
  axisymmetric or Navier--Stokes singularity exclusion is claimed.
- **C197 boundary:** the upper-error theorem uses compact
  \(\mathbb R^3\) beams, a common lattice with \(d\varepsilon\ge1\), and
  coefficient \(\ell^2\) normalization. It does not combine C195's
  shrinking cone with C196's fixed aperture, nor turn a bandlimited global
  envelope into an annulus-supported one. The retained-band absolute
  endpoint and all viscous losses remain open.
- **C199 boundary:** \(S\) is the bounded \(L^2\) extension of a linearized
  PDE operator, not a proved nonlinear endpoint map on an open \(L^2\)
  ball. A finite block must compress the full evolution with validated
  tail error. Neither a global torus dilation nor a same-fiber replacement
  of \(h(\mu')\) is permitted.
- **C204 boundary:** equal perturbation endpoint energies, a retained
  linear concentration gain and bounded linear trajectory costs do not
  control \(P(v\cdot\nabla v)\), the nonlinear pressure/pump/wake
  response, or a repeated scale return. The sufficient amplitude
  bound is \(\|v(0)\|_\infty\le20000bq^{9/8}\), not C182's \(bq\)
  premise. C180 center/frame/seed/splitter compatibility and the standing
  terminal singular-center obligations remain open. No external
  cross-audit or numerical nonlinear residual certificate is asserted.
