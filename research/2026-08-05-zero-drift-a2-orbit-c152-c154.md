# C152--C154: a zero-drift periodic orbit on the existing \(A_2\) pump

**Date:** 2026-08-05

**Status:** C152 and C154 exact; C153 is a **NUMERICAL CANDIDATE only**

**Checker:**
[checks/zero_drift_a2_orbit_c152_c154.py](../checks/zero_drift_a2_orbit_c152_c154.py)

## 0. Claim boundary

This note stays on the C121/C149 three-root \(A_2\) Beltrami pump.  It does
not introduce a new geometry.  Its purpose is to remove the axial-drift
defect of the vertical C149 ray and then to state exactly what is and is not
known about the resulting Kelvin cocycle.

* **C152 (EXACT):** for \(\delta=4/5\), the regular level \(f_\delta=0\)
  is one closed streamline in the \((a,b)\) phase torus (its physical lifts
  are closed translates).  Its period has an explicit complete integral
  representation.  The full three-dimensional *linearized* return at any
  chosen lift is an exact rank-one shear, and its periodic covectors form an
  explicit two-plane containing off-plane covectors.
* **C153 (NUMERICAL CANDIDATE, not a theorem):** ordinary high-accuracy
  integration gives transverse Kelvin trace
  \(16716.8837799\ldots\) and a positive exponent
  \(3.20281068\ldots\); Liouville's formula separately gives the exact
  determinant identity \(\det M=1\).  The checker reproduces the numerical
  values with a dependency-free RK4 convergence diagnostic.  It is not an
  outward-rounded interval or Taylor certificate, so C153 must not be
  registered as an
  established instability.
* **C154 (EXACT):** the same rank-one return shears generic
  three-dimensional Kelvin-covector bandwidth in the fiber over this orbit,
  linearly in the number of returns.  Consequently an \(O(q)\)
  normal-charge width becomes \(O(q\log q)\) during \(O(\log q)\) returns.
  An \(O(q)\) final fiber band requires either a correlated two-dimensional
  packet or an \(O(q/\log q)\) initial width in the shearing direction.
  Turning this fiber statement into a localized Fourier-packet theorem is
  not part of C154.

No finite-frequency localized Navier--Stokes amplifier, leakage estimate,
or nonlinear stage map is proved here.

## 1. C152: the exact zero-drift orbit

Keep

\[
 N=(1,1,1),\qquad r_1=(1,-1,0),\qquad r_2=(0,1,-1),
 \qquad r_3=-r_1-r_2=(-1,0,1),                         \tag{1.1}
\]

and put

\[
 f(a,b)=\cos a+\cos b+\delta\cos(a+b),\qquad
 U=N\times\nabla_xf-\sqrt2 fN,\qquad \delta={4\over5},       \tag{1.2}
\]

where \(a=r_1\cdot x\), \(b=r_2\cdot x\).  As in C121,
\(\nabla\times U=\sqrt2 U\).  The phase equations are

\[
 \dot a=-3f_b=3\{\sin b+\delta\sin(a+b)\},\qquad
 \dot b= 3f_a=-3\{\sin a+\delta\sin(a+b)\}.             \tag{1.3}
\]

They preserve \(f\).  On the level \(f=0\), the axial term in (1.2)
vanishes pointwise, not merely on average.  Thus every trajectory on this
level has exactly zero axial drift.

### 1.1 Regularity and topology of the level

The complete critical-point list on the phase torus is

\[
\begin{array}{c|c}
 (a,b)&f(a,b)\\ \hline
 (0,0)&14/5\\
 (\pi,\pi)&-6/5\\
 (0,\pi),(\pi,0)&-4/5\\
 a=b,\ \cos a=-5/8&-57/40.
\end{array}                                                   \tag{1.4}
\]

Hence zero is a regular value.  There is only one critical point above
zero, the strict nondegenerate maximum at \((0,0)\): its phase Hessian has
eigenvalues \(-1\) and \(-(1+2\delta)\).  A superlevel immediately below
that maximum is a disk, and there is no critical value between it and zero.
The regular-level isotopy theorem (equivalently, the elementary Morse
superlevel argument) therefore gives that \(\{f\ge0\}\) is a disk and
\(\{f=0\}\) is its single smooth boundary circle in the phase torus.  Since
this circle is contractible, each lift relevant to a fixed physical orbit
also closes; “single” here refers to the phase-torus level, not to all its
translated physical copies.

One algebraic point on that circle is

\[
 (a_0,b_0)=(A,-A),\qquad \cos A=-{2\over5},\qquad
 \sin A={\sqrt{21}\over5}.                                  \tag{1.5}
\]

At the corresponding point \(X_0\),

\[
 g_0:=\nabla f(X_0)={\sqrt{21}\over5}(-1,2,-1),\qquad
 U_0:=U(X_0)={3\sqrt{21}\over5}r_3.                         \tag{1.6}
\]

### 1.2 Exact period representation

Set \(p=(a+b)/2\), \(q=(a-b)/2\).  On the nearby level \(f=E\),

\[
 \dot p^{\,2}={9\over4}\left[4\cos^2p-
                 \{E-\delta\cos(2p)\}^2\right].             \tag{1.7}
\]

Define

\[
 y_\pm(E)=
 {1+\delta(\delta+E)\pm\sqrt{1+2\delta(\delta+E)}
  \over2\delta^2}.                                          \tag{1.8}
\]

At \(E=0\),

\[
 y_\pm(0)={41\pm5\sqrt{57}\over32},\qquad
 p_+=\arccos\sqrt{y_-(0)}
     =\arccos {\sqrt{57}-5\over8}.                          \tag{1.9}
\]

The endpoint-regularized complete period is

\[
 \boxed{
 T(E)={4\over3\delta}\int_0^{\pi/2}
 {d\theta\over\sqrt{A_E(\theta)B_E(\theta)}}},             \tag{1.10}
\]

where

\[
 A_E=y_-(E)+\{1-y_-(E)\}\sin^2\theta,
 \quad
 B_E=y_+(E)-1+\{1-y_-(E)\}\cos^2\theta.                   \tag{1.11}
\]

This follows from

\[
 4\cos^2p-\{E-\delta\cos2p\}^2
 =4\delta^2\{\cos^2p-y_-(E)\}\{y_+(E)-\cos^2p\}          \tag{1.12}
\]

and \(\cos^2p=y_-+(1-y_-)\sin^2\theta\).  Moreover,

\[
 y_\pm'(E)={1\pm[1+2\delta(\delta+E)]^{-1/2}\over2\delta},\tag{1.13}
\]

so \(A_E'\ge0\), \(B_E'>0\) near zero and

\[
 T'(E)=-{2\over3\delta}\int_0^{\pi/2}
 {A_E'/A_E+B_E'/B_E\over\sqrt{A_EB_E}}\,d\theta<0.          \tag{1.14}
\]

For orientation only (not as an interval certificate),

\[
 T(0)=3.03613777009396\ldots,\qquad
 T'(0)=-1.79832641310817\ldots.                              \tag{1.15}
\]

### 1.3 Exact rank-one return and the periodic-covector plane

Let \(F=D\Phi_{T(0)}(X_0)\) be the derivative of the full
three-dimensional flow after one period.  At first variation, write
\(E=g_0\cdot h\) for the level displacement.  The planar phase slip is
\(-T'(0)U_0E\), while differentiating the exact axial displacement
\(-\sqrt2 E T(0)N\) gives \(-\sqrt2T(0)NE\).  Therefore the following is an
exact identity for the derivative \(F\):

\[
 \boxed{F=I+u\otimes g_0},\qquad
 u=-T'(0)U_0-\sqrt2T(0)N.                                   \tag{1.16}
\]

Since \(g_0\cdot U_0=g_0\cdot N=0\), one has
\((u\otimes g_0)^2=0\), and hence

\[
 F^{-T}=I-g_0\otimes u.                                     \tag{1.17}
\]

The periodic covectors are exactly the two-plane

\[
 \boxed{\mathcal P=u^\perp}.                                \tag{1.18}
\]

It contains the explicit off-plane vector

\[
 k_0=U_0+\beta N,\qquad
 \beta=-{T'(0)|U_0|^2\over\sqrt2T(0)|N|^2}>0.               \tag{1.19}
\]

Numerically, \(\beta=2.11087534389\ldots\) and
\(k_0/|k_0|=(-0.1196603532,0.3954907233,0.9106417998)\).

## 2. C153: a large Kelvin multiplier -- numerical candidate only

Along the orbit let \(g=\nabla f\), \(H=\nabla^2f\), and let \(C_Nv=N\times
v\).  The exact velocity gradient is

\[
 A(t)=\nabla U=C_NH-\sqrt2Ng^T,                              \tag{2.1}
\]

with

\[
 H=-\cos a\,r_1r_1^T-\cos b\,r_2r_2^T
   -\delta\cos(a+b)(r_1+r_2)(r_1+r_2)^T.                   \tag{2.2}
\]

For the periodic covector (1.19), the Kelvin system is

\[
 \dot k=-A^Tk,\qquad
 \dot v=-Av+2{k\over|k|^2}(k^TAv),\qquad k\cdot v=0.       \tag{2.3}
\]

Using the deterministic orthonormal frame obtained by projecting
\((1,0,0)\) to \(k_0^\perp\), direct integration gives

\[
 M_{\rm num}=\begin{pmatrix}
 2028.59222878&6784.66646127\\
 4391.74913973&14688.29155113
 \end{pmatrix},                                             \tag{2.4}
\]

and therefore

\[
 \operatorname{tr}M_{\rm num}=16716.88377991\ldots,
 \quad \rho_{\rm num}=16716.88372010\ldots,
 \quad {\log\rho_{\rm num}\over T}=3.20281068\ldots.     \tag{2.5}
\]

For an exactly periodic covector, Liouville's formula does give the exact
identity

\[
 \det M={|k(0)|\over|k(T)|}=1,                              \tag{2.6}
\]

because the transverse trace of the generator in (2.3) is
\(-d(\log|k|)/dt\).  Equation (2.6) does **not** turn the floating-point
trace in (2.5) into a proof that \(\operatorname{tr}M>2\).

### 2.1 Failed certification attempt and the exact missing certificate

The checker integrates (1.3) and (2.3) with 4096, 8192 and 16384 fixed RK4
steps.  The trace stabilizes under these refinements at the scale quoted in
(2.5), and the computed covector return residual is tiny.  This is a strong
target for validation, but agreement under mesh refinement is not an a
posteriori truncation bound and IEEE roundoff is not directed interval
arithmetic.  In particular, none of the floating-point inequalities in the
checker is usable as a premise in C152, C154, or a later theorem.

A dependency-free proof can be made finite as follows: use the smooth
\(\theta\)-representation (1.10) on four orbit quarters; enclose
\(T,T',\beta\) with rational outward rounding; propagate the polynomial
sin--cos augmentation of (1.3), together with (2.3), by a validated Taylor
method; and bound the terminal frame projection.  A successful certificate
only needs

\[
 \inf\operatorname{tr}M>2                                  \tag{2.7}
\]

(the numerical margin suggests the stronger bound \(>10^4\)).  Such an
outward-rounded remainder certificate is not present here.  Componentwise
interval RK without a Taylor-model/QR wrapping control was also rejected:
the unstable amplitude column makes its box width swamp the final trace.
Thus C153 remains deliberately unregistered as a theorem.

## 3. C154: exact fiberwise Kelvin-bandwidth shear

Write \(K=F^{-T}\).  From (1.17) and \(u\cdot g_0=0\),

\[
 \boxed{K^\ell=I-\ell g_0\otimes u}\qquad(\ell\in\mathbb Z).\tag{3.1}
\]

Consequently a covector displacement in the fiber over the chosen periodic
orbit obeys

\[
 \Delta k_\ell=\Delta k_0-\ell g_0(u\cdot\Delta k_0).      \tag{3.2}
\]

In particular,

\[
 u\cdot N=-\sqrt2T(0)|N|^2\ne0.                            \tag{3.3}
\]

If an initially \(O(q)\)-wide fiber support contains a displacement with
\(|u\cdot\Delta k_0|\ge c q\), for a fixed \(c>0\), (3.2) gives a
\(\Theta(q\ell)\) component in the \(g_0\) direction after \(\ell\)
returns (up to the initial \(O(q)\) term).  A fixed per-return multiplier
would require \(\ell\asymp\log q\) for polynomial gain, so at the exact
linear-cocycle level

\[
 \boxed{\text{nondegenerate generic final bandwidth}
        =\Theta(q\log q),\quad\text{not }O(q).}              \tag{3.4}
\]

Two immediate linear repairs at this level are:

1. correlate the packet support so every displacement lies in
   \(u^\perp\), a two-dimensional band fixed by the shear; or
2. retain a three-dimensional packet but narrow its initial width in the
   shearing direction to \(O(q/\log q)\).

If the numerical multiplier in C153 is eventually certified, the C147
gain \(H=n^{26}\), \(q=n^8\) would use

\[
 \ell_H={26\log n\over\log\rho}
        =0.334\ldots\log q.                                 \tag{3.5}
\]

The coefficient in (3.5) is numerical/conditional on C153.  The exact
\(q\log q\) conclusion in (3.4) is a fiberwise conclusion conditional only
on taking \(\ell\asymp\log q\); C154 does not prove that a spatially
localized packet follows this monodromy for that long.  The extra logarithm
is factorially harmless for viscosity, but it is not harmless for the
promised \(O(q)\) Fourier support and packet-coherence ledger.  That support
repair, the rigorous Kelvin certificate, and the localized nonlinear
embedding are the surviving obligations.
