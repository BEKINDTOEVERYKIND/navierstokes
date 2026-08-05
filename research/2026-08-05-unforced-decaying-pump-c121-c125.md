# C121--C125: an exact unforced Beltrami pump and the decaying-pump gain gate

**Date:** 2026-08-05  
**Status:** exact Navier--Stokes background; exact linear and nonlinear
finite-dimensional normal-form gain ledgers; exact full-perturbation
retained-amplitude identity; no one-cell embedding theorem  
**Scope:** the C118--C120 A2 geometry.  No new geometry and no blow-up claim.

## 1. Verdict

The A2 pump itself needs no force.  A real field supported on the six A2
roots with one helicity is Beltrami, so viscosity only multiplies it by an
exponential.  Along any unstable singular direction of the finite hexagon
normal form, the perturbation is therefore exactly solvable.  If

\[
 \dot z=(\sigma P_0e^{-\delta t}-d)z,
\tag{1.1}
\]

then its unique gain maximum occurs at

\[
 t_*={1\over\delta}\log {\sigma P_0\over d}
\quad\hbox{when}\quad \sigma P_0>d,
\tag{1.2}
\]

and

\[
 \boxed{
 \log {z(t_*)\over z(0)}
 ={d\over\delta}\bigl(R-1-\log R\bigr),
 \qquad R={\sigma P_0\over d}>1.}
\tag{1.3}
\]

This is a genuine unforced decaying-pump gain formula for the ideal retained
block.  Its nonlinear finite normal form retains almost all of the gain for
a sufficiently small seed, with an explicit seed condition below.

The result is not yet the requested full Navier--Stokes stage map.  The
six-leaf support admits centre and second-sideband outputs at the same
bilinear order as the desired hexagon edge.  C125 gives an exact relative
Duhamel gate for retained-mode gain, including the quadratic perturbation
term.  It neither proves that gate nor supplies the terminal chart and
pump-to-pump conversion required by a stage map.

The formulas and rational tests are checked by
[`checks/unforced_decaying_pump_c121_c125.py`](../checks/unforced_decaying_pump_c121_c125.py).

## 2. C121: the A2 Beltrami pump is an exact unforced NS solution

Let \(n=(1,1,1)\), let \(r\) be any root from C118, and put

\[
 t_r=n\times r,\qquad h_r=t_r+i\sqrt2\,n.
\tag{2.1}
\]

The A2 identities give

\[
 r\cdot h_r=0,\qquad i r\times h_r=\sqrt2\,h_r.
\tag{2.2}
\]

At integer pump frequency \(K>0\), choose arbitrary complex coefficients
\(c_r\) on the three displayed roots and set

\[
 U_*(x)=\sum_{r\in\{r_1,r_2,r_3\}}
 \left(c_rh_r e^{iKr\cdot x}
       +\overline{c_rh_r}\,e^{-iKr\cdot x}\right).
\tag{2.3}
\]

The resulting real periodic field satisfies

\[
 \nabla\cdot U_*=0,\qquad
 \operatorname{curl}U_*=\kappa U_*,\qquad
 \Delta U_*=-\kappa^2U_*,\qquad \kappa=\sqrt2K.
\tag{2.4}
\]

Since

\[
 (U_*\cdot\nabla)U_*
 =\nabla\frac{|U_*|^2}{2}-U_*\times\operatorname{curl}U_*,
\tag{2.5}
\]

its Leray-projected Euler nonlinearity vanishes.  Therefore

\[
 \boxed{U(t)=P_0e^{-\delta t}U_*,\qquad
 \delta=\nu\kappa^2=2\nu K^2}
\tag{2.6}
\]

is an exact unforced Navier--Stokes solution on the torus.

This exact statement concerns the pump only.  It does not make a chosen
finite perturbation ladder invariant.

## 3. C122: exact linear gain on a retained ladder eigenmode

Let \(L\) be the frozen finite ladder matrix and let
\(Lv=\sigma v\), \(\sigma>0\).  Suppose all retained leaf wavevectors
have squared length \(k_\ell^2\), as in the A2 shell identity.  In the
ideal retained (truncated) block, linearization about (2.6) is exactly
(1.1), with

\[
 d=\nu k_\ell^2.
\tag{3.1}
\]

Consequently

\[
 z(t)=z_0\exp G(t),\qquad
 G(t)={\sigma P_0\over\delta}(1-e^{-\delta t})-dt.
\tag{3.2}
\]

For the geometric leaf
\(q_m\pm Kr_i\),

\[
 k_\ell^2=3m^2+2K^2,\qquad
 {d\over\delta}={3m^2+2K^2\over2K^2}.
\tag{3.3}
\]

The symbol \(\sigma\) in (3.2) is the **physical** inviscid eigenvalue per
unit pump amplitude.  C120 encloses \(\sigma_*\) for its scalar weighted
normal form.  Equating \(\sigma\) with \(\sigma_*\) in the PDE requires
the still-open polarization/coefficient identification.

## 4. C123: the threshold and optimum are sharp

Differentiating (3.2) gives

\[
 G'(t)=\sigma P_0e^{-\delta t}-d,\qquad
 G''(t)=-\delta\sigma P_0e^{-\delta t}<0.
\tag{4.1}
\]

If \(R=\sigma P_0/d\le1\), there is no positive-time gain.  If \(R>1\),
the unique maximum is (1.2), and direct substitution gives (1.3).  Since

\[
 R-1-\log R>0\qquad(R>1),
\tag{4.2}
\]

the maximum gain is strict.  Near threshold it is quadratic,

\[
 R-1-\log R=\frac12(R-1)^2+O((R-1)^3),
\tag{4.3}
\]

while for large \(R\) it is \(R+O(\log R)\).  Thus pump heat decay does
not kill a sufficiently supercritical one-stage gain; it fixes the
available integrated exponent.

## 5. C124: nonlinear finite-normal-form persistence with an explicit seed gate

In the normalized singular-direction normal form, include pump feedback
and leaf viscosity:

\[
\begin{aligned}
 \dot p&=-\delta p-2\sigma R_\ell^2,\\
 \dot R_\ell&=(\sigma p-d)R_\ell,
 \qquad p(0)=P_0>0,\quad R_\ell(0)=\varepsilon.
\end{aligned}
\tag{5.1}
\]

For three equal copies, \(R_\ell^2\) is simply the sum of the three
squared singular amplitudes; their ratios remain fixed.  The exact energy
identity is

\[
 {d\over dt}(p^2+2R_\ell^2)
 =-2\delta p^2-4dR_\ell^2.
\tag{5.2}
\]

Let \(\bar p=P_0e^{-\delta t}\) and
\(\bar R=\varepsilon e^{G(t)}\) denote the undepleted solution.  As long
as \(p\ge0\), comparison gives

\[
 0\le q:=\bar p-p,\qquad R_\ell\le\bar R,\qquad
 q(t)\le2\sigma t M(t)^2,
\tag{5.3}
\]

where, for \(0\le t\le t_*\),

\[
 M(t)=\varepsilon\exp\bigl(\max_{0\le s\le t}G(s)\bigr).
\tag{5.4}
\]

Moreover

\[
 \log{\bar R(t)\over R_\ell(t)}
 =\sigma\int_0^tq(s)\,ds
 \le \sigma^2M(t)^2t^2.
\tag{5.5}
\]

Since \(G\) is increasing on \([0,t_*]\), put
\(M_*:=M(t_*)=\varepsilon e^{G(t_*)}\).  The two explicit small-seed
conditions

\[
 2\sigma t_*M_*^2<\frac{d}{2\sigma},
 \qquad
 \sigma^2M_*^2t_*^2\le\eta
\tag{5.6}
\]

close the continuation assumption in (5.3): on \([0,t_*]\),
\(\bar p(t)\ge d/\sigma\) and
\(q(t)<d/(2\sigma)\), so \(p(t)>0\).  They also give

\[
 \boxed{R_\ell(t_*)\ge
 \varepsilon\exp\bigl(G(t_*)-\eta\bigr).}
\tag{5.7}
\]

Hence arbitrary finite ideal gain survives this nonlinear finite normal
form by choosing the incoming seed sufficiently small.  Identifying (5.1)
with a physical Galerkin block still requires the coefficient and
polarization realization left open in C118--C120.  This is not a uniform
seed theorem for an iterated PDE cascade: exponentially small seeds make
the relative leakage requirement in the next section harder, not easier.

## 6. C125: an exact full-perturbation retained-amplitude gate

Write the exact perturbation equation about the decaying Beltrami pump as

\[
 v'=A(t)v+\mathcal N(v,v),
 \qquad \mathcal N(v,v)=-\mathbb P(v\cdot\nabla v),
\tag{6.1}
\]

and let \(P_H\) project onto the retained six-leaf state.  Put
\(A_H(t)=P_HA(t)P_H\) on that state.  Thus

\[
\begin{aligned}
 v_H'&=A_H(t)v_H+P_HA(t)v_\perp+P_H\mathcal N(v,v),\\
 v_\perp'&=P_\perp A(t)(v_H+v_\perp)
             +P_\perp\mathcal N(v,v).
\end{aligned}
\tag{6.2}
\]

Let \(A_H^0(t)\) be the ideal retained normal-form block and let \(\phi\)
be a fixed left eigenfunctional satisfying

\[
 \phi A_H^0(t)=(\sigma P_0e^{-\delta t}-d)\phi,
 \qquad E_H(t)=A_H(t)-A_H^0(t).
\tag{6.3}
\]

Normalize the corresponding right eigenvector \(v_+\) by
\(\phi(v_+)=1\), and set
\(z=\phi(v_H)\).  Variation of constants gives the exact scalar identity

\[
 e^{-G(T)}z(T)-z(0)
 =\int_0^T e^{-G(s)}
 \phi\!\left(P_HA(s)v_\perp(s)+E_H(s)v_H(s)
              +P_H\mathcal N(v,v)(s)\right)ds
\tag{6.4}
\]

for every time interval on which the smooth perturbation exists.  Thus a
sufficient gate for the retained growing amplitude is

\[
 \boxed{
 {1\over|z(0)|}
 \int_0^{t_*}e^{-G(s)}
 \left|
 \phi\!\left(P_HA(s)v_\perp(s)+E_H(s)v_H(s)
              +P_H\mathcal N(v,v)(s)\right)
 \right|ds\le\theta<1.}
\tag{6.5}
\]

Under (6.5),

\[
 |z(t_*)|\ge(1-\theta)|z(0)|e^{G(t_*)}.
\tag{6.6}
\]

The linearized version is obtained by deleting \(\mathcal N\).  The A2
support calculation shows why even that restricted form of (6.5) is
load-bearing.  For each first leaf,

\[
 q_m+Kr_i\longrightarrow q_m,\quad q_m+2Kr_i
\tag{6.7}
\]

after one interaction with the pump pair \(\mp Kr_i\).  These modes are
allowed at the same bilinear order in pump amplitude and time as the
retained hexagon neighbours.  Support arithmetic alone neither proves that
their coefficients are nonzero nor provides a \(K^{-1}\), viscosity, or
small-seed advantage in the *relative* norm (6.5).

Therefore the exact named obstruction after C125 is:

> **same-order off-ladder slaving:** prove (6.5), uniformly in the stage
> schedule, for the centre, second sidebands, their descendants, and the
> nonlinear pump-shape deformation; or exhibit an exact cancellation which
> removes their return into the growing leaf.

The decaying Beltrami pump and all finite-dimensional gain/depletion
normal-form formulas are exact.  The gate controls only one retained
amplitude; it does not establish the active focusing or terminal chart
needed to identify that amplitude with the next pump.  The same-order
off-ladder/nonlinear estimate is not proved, so there is no full
Navier--Stokes one-cell closure and no singularity claim.
