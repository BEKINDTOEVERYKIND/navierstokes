# C160: exact unequal-radius detuning of the limiting DACR mean

**Date:** 2026-08-05
**Status:** exact limiting periodic-wake normal form; the causal finite-time
transient, finite \(\epsilon\), additive quartets, and localization remain
open
**Checker:**
[checks/two_radius_detuned_dacr_c160.py](../checks/two_radius_detuned_dacr_c160.py)

## 0. Claim boundary

This note stays on the C149 limiting elliptic ray of the C121 \(A_2\) pump.
It computes the persistent quadratic wake and directed cubic return for two
distinct radii on the resonant cone.

The exact asymptotic result is favorable: at every fixed nonzero radial
gap and every interior nonzero angular separation, the periodic sum-wake
and difference-wake secular means cancel.  The zero-initial causal
solution has zero directed **Cesàro** long-time DACR mean.

The limit is nonuniform.  The difference-wake inertial frequency is linear
in the relative radial gap.  Adjacent normalized charge layers have gap
\(O(q^{-1})\), so their detuning time is \(O(q)\), much longer than the
scheduled \(O(\log q)\) gain interval.  The zero-initial causal wake is a
long transient on that interval; the long-time zero cannot be inserted
into the stage ledger.

## 1. Two unequal resonant radii

Use

\[
 \eta_\phi=(\cos\phi,\sin\phi,1/\sqrt3),
 \qquad p=A\eta_0,qquad q=B\eta_\phi,
 \qquad A,B>0,\quad A\ne B,
 \tag{1.1}
\]

Choose the principal separation \(-\pi<\phi<\pi\), assume
\(\phi\ne0\), and put

\[
 t=\tan(\phi/2).
 \tag{1.2}
\]

The limiting selected rotating-frame parents are the C155 modes

\[
 b_0(s)=R(-s)a(s),\qquad b_\phi(s)=R(-s)a(s+\phi).
 \tag{1.3}
\]

Every quadratic forcing contains only the harmonics \(0,\pm2\).  For an
output \(k\), the rotating-frame Kelvin operator is

\[
 C_k=-2J+2{k k^T\over|k|^2}J.
 \tag{1.4}
\]

On \(k^\perp\),

\[
 C_k^2=-\omega_k^2I,
 \qquad
 \omega_k^2={4k_z^2\over|k|^2}.
 \tag{1.5}
\]

For \(k_+=p+q\) and \(k_-=p-q\), exact half-angle simplification gives

\[
 \boxed{
 \omega_+^2={(A+B)^2(1+t^2)\over D_+},\qquad
 \omega_-^2={(A-B)^2(1+t^2)\over D_-}}
 \tag{1.6}
\]

where, in two useful equivalent forms,

\[
 D_+=(1+t^2)(A^2+B^2)+AB(2-t^2)
     =(A+B)^2+t^2(A^2+B^2-AB),
 \tag{1.7}
\]

\[
 D_-=(1+t^2)(A^2+B^2)+AB(t^2-2)
     =(A-B)^2+t^2(A^2+B^2+AB).
 \tag{1.8}
\]

The strict nonresonance needed below is exact:

\[
 \omega_+^2-1={3ABt^2\over D_+},\qquad
 1-\omega_-^2={3ABt^2\over D_-},
 \tag{1.9}
\]

and

\[
 4-\omega_+^2
 ={3\{(A+B)^2+t^2(A-B)^2\}\over D_+}.
 \tag{1.10}
\]

Thus

\[
 0<\omega_-<1<\omega_+<2.
 \tag{1.11}
\]

At \(t=0\), both frequencies equal one and the periodic homogeneous
problem is resonant; at \(A=B\), \(\omega_-=0\).  The antipodal endpoint
\(|t|=\infty\) requires a separate limiting parametrization.  None of
these endpoints is part of C160.

## 2. Exact periodic wake solver

Write either projected quadratic forcing as

\[
 g(s)=g_0+g_c\cos2s+g_s\sin2s.
 \tag{2.1}
\]

Because of (1.11), the unique \(2\pi\)-periodic solution of

\[
 w'=C_kw+g
 \tag{2.2}
\]

is

\[
 w=w_0+w_c\cos2s+w_s\sin2s,
 \tag{2.3}
\]

with

\[
 \boxed{
 w_0={C_kg_0\over\omega_k^2},\qquad
 w_c={C_kg_c+2g_s\over\omega_k^2-4},\qquad
 w_s={C_kw_c+g_c\over2}.}
 \tag{2.4}
\]

Substitution uses only \(C_k^2=-\omega_k^2I\) on the transverse plane.
The checker reconstructs both sum and difference wakes from (2.4) in
\(\mathbb Q(\sqrt2,\sqrt3)\) and verifies (2.2) coefficient by coefficient.

For \(A\ne B\) and \(t\ne0\), both periodic wake solutions are unique: the
only possible integer frequency in the open interval \((0,2)\) is one,
which (1.9) excludes.  A common phase shift of the parent pair can
therefore be absorbed by the corresponding time shift and planar
rotation, so taking the target phase to be zero in (1.1) loses no
generality.  This argument fails on the excluded sets: at \(A=B\) the
difference frequency is zero, while at \(t=0\) both homogeneous
frequencies are one.

## 3. Exact cancellation of the directed periodic mean

Let \(\mathcal S\) be the symmetric projected Euler symbol from C155.  Feed
the periodic sum and difference wakes from (2.4) back to the target:

\[
 R_+=-\mathcal S(p+q,w_+;-q,b_\phi),
 \qquad
 R_-=-\mathcal S(p-q,w_-;q,b_\phi).
 \tag{3.1}
\]

Exact trigonometric averaging gives

\[
 \boxed{
 \left\langle b_0,R_+\right\rangle_{\rm per}
 =-{AB\,t^3(t^2+2)(t^2-t+2)(t^2+t+2)
       \over2(1+t^2)^5}}
 \tag{3.2}
\]

and

\[
 \boxed{
 \left\langle b_0,R_-\right\rangle_{\rm per}
 =+{AB\,t^3(t^2+2)(t^2-t+2)(t^2+t+2)
       \over2(1+t^2)^5}}
 \tag{3.3}
\]

Therefore

\[
 \boxed{
 \left\langle b_0,R_++R_-\right\rangle_{\rm per}=0
 \qquad(A\ne B,\ t\ne0)}
 \tag{3.4}
\]

This is stronger than pairwise energy antisymmetry: each physical directed
coordinate has zero periodic steady DACR mean across two unequal radii.
It concerns the limiting frozen parents and the persistent periodic wake,
not a short zero-initial causal stage.

## 4. The noncommuting gap/time limits

For a zero-initial quadratic wake, the exact solution is

\[
 w_{\rm causal}(s)=w_{\rm per}(s)-e^{sC_k}w_{\rm per}(0).
 \tag{4.1}
\]

Thus (3.4) removes the periodic steady contribution, but the homogeneous
transient remains.  The return rows contain only the harmonics
\(0,\pm2\).  Since (1.11) excludes resonance with those harmonics, each
homogeneous transient has zero Cesàro mean.  Consequently the complete
zero-initial directed return satisfies

\[
 \lim_{T\to\infty}{1\over T}\int_0^T
 \langle b_0,R_{+,\rm causal}+R_{-,\rm causal}\rangle\,ds=0.
 \tag{4.2}
\]

This is a time-average statement; the inviscid homogeneous wake does not
converge to the periodic wake.

The difference-channel frequency is \(\omega_-\).  On parameter sets with
\(t\) bounded away from zero and infinity and with \(A/B\) bounded above
and below, (1.6) gives the scale-invariant estimate

\[
 \omega_-\asymp {|A-B|\over A+B}.
 \tag{4.3}
\]

If the common radius is additionally normalized to order one, the right
side is comparable to \(|A-B|\).  Direct substitution in (2.4) shows on
the same compact sets that \(w_{-,\rm per}(0)\) stays bounded as
\(A\to B\): its constant forcing is \(O(|A-B|)\), while \(C_-\) is a
skew-adjoint rotation of norm \(\omega_-\) on \((p-q)^\perp\), with

\[
 e^{sC_-}=\cos(\omega_-s)I
       +{\sin(\omega_-s)\over\omega_-}C_-.
\]

Integrating the difference transient against its bounded \(0,\pm2\)
return row then gives only the upper estimate

\[
 \left|\int_0^T h_-(s)e^{sC_-}w_{-,\rm per}(0)\,ds\right|
 \le C_{\mathcal K}\min\{T,\omega_-^{-1}\}.
 \tag{4.4}
\]

Here \(\mathcal K\) is the chosen compact normalized parameter set; the
nonresonant \(2\pm\omega_-\) contributions are uniformly bounded.
Equation (4.4) is only an upper scale statement.  It is neither a lower
bound nor a proof that a coefficient of a useful sign remains coherent
until \(T\simeq\omega_-^{-1}\).

For adjacent normalized layers, take \(A=1+q^{-1}\), \(B=1\).  At quarter
separation,

\[
 \omega_-^2={2(A-B)^2\over2A^2-AB+2B^2}\asymp q^{-2}.
 \tag{4.5}
\]

Consequently

\[
 T_{\rm detune}\asymp\omega_-^{-1}\asymp q,
 \qquad
 \omega_-T_{\rm stage}\asymp{\log q\over q}\longrightarrow0.
 \tag{4.6}
\]

The two limits do not commute:

- at fixed \(A\ne B\), sending the number of periods to infinity yields
  the zero mean (3.4);
- on the scheduled \(O(\log q)\) interval with
  \(|A-B|\asymp q^{-1}\), the solution remains in the unresolved causal
  transient and approaches the equal-radius dynamics rather than its
  long-time detuned average.

Thus radial thickening supplies a genuine asymptotic cancellation but not
the needed short-stage estimate.

## 5. What remains load-bearing

C160 closes the unequal-radius **periodic steady diagonal** calculation.
It does not classify:

1. the causal transient with a signed, backward-weighted C125 estimate;
2. shared outputs with more than one additive quartet;
3. finite-\(\epsilon\) slowly growing parent lines, where the rigid-rotation
   periodic solver acquires errors over \(O(\log q)\);
4. a positive balanced distribution across \(q\) adjacent charge layers;
5. lattice discretization, the neutral-wake invariant graph, or
   localization/export.

The exact next question is whether the fixed-ring positive balance can be
kept during the long gain and radial/charge thickness introduced only in a
just-in-time terminal splitter.  No one-cell Navier--Stokes stage or
Millennium conclusion is claimed.
