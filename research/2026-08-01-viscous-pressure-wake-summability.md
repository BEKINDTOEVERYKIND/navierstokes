# Scaling and summability of the unavoidable viscous pressure wake

**Date:** 2026-08-01
**Status:** self-derived continuation of the second-jet obstruction; not
cross-audited
**Scope:** whole-space and local periodic estimates for the first nonlocal
velocity jet.  No all-order endpoint theorem is claimed.

## 1. Result

The second-jet obstruction in
[`2026-07-29-gavrilov-viscous-endpoint-jet-obstruction.md`](2026-07-29-gavrilov-viscous-endpoint-jet-obstruction.md)
forces a noncompact velocity tail.  The obstruction is real, but its scale
does not by itself destroy the high-Reynolds cascade ledger.

For one bubble of amplitude \(a\) and diameter \(\ell\), the nonlocal part
of the second velocity jet obeys, a distance \(d\gtrsim\ell\) away,

\[
        |\nabla^m u_2^{\rm nl}(x)|
        \lesssim_m \nu a^2\ell\,d^{-4-m}.                 \tag{1.1}
\]

Over the inertial time \(\tau=\ell/a\), its Taylor coefficient is

\[
        |\nabla^m(\tau^2u_2^{\rm nl})(x)|
        \lesssim_m \nu\ell^3d^{-4-m}.                    \tag{1.2}
\]

The amplitude cancels exactly.  On the next comparable separated scale,
the wake-to-main ratio is \(O({\rm Re}^{-1})\).

For the packed carrier with \(N\asymp K^3\) bubbles of diameter
\(\delta=\ell/K\), but the same macro-stage time \(\ell/a\), the crude
no-cancellation sum is instead

\[
        |\nabla^m Z^{(2)}(x)|
        \lesssim_m \nu\ell^3K^2d^{-4-m}.                 \tag{1.3}
\]

At a comparable separated scale its zeroth-order ratio to the next main
carrier is

\[
       {\nu K^2/\ell\over a}
       ={K\over {\rm Re}_{\rm carrier}},
       \qquad
       {\rm Re}_{\rm carrier}={a\ell\over\nu K}.          \tag{1.4}
\]

For \(\ell_j=r^{-j}\), \(K_j=(j+j_0)^A\), and
\(a_j=\ell_j^{-\gamma}K_j^\gamma\), \(1<\gamma<3/2\), this ratio is

\[
       \nu\ell_j^{\gamma-1}K_j^{2-\gamma}\longrightarrow0.  \tag{1.5}
\]

Thus the mandatory second-jet wake is perturbative and summable in the
existing geometric/polynomial ledger, even in the crude packed-bubble
bound.  This does **not** prove that all later pressure jets, wake--core
cross stresses, or endpoint derivatives obey the Gevrey-2 majorant.

## 2. One-bubble kernel estimate

Let

\[
 W_{a,\ell,c,Q}(x)
 =aQW\!\left(Q^T{x-c\over\ell}\right),                   \tag{2.1}
\]

where \(W\in C_c^\infty\) is a fixed steady Euler bubble.  At a bare
steady endpoint, the first pressure derivative is

\[
 -\Delta p_1=\partial_i\partial_jS_{ij},\qquad
 S=\nu(W_{a,\ell}\otimes\Delta W_{a,\ell}
       +\Delta W_{a,\ell}\otimes W_{a,\ell}).             \tag{2.2}
\]

The source has the exact scaling

\[
                  \|S\|_1\lesssim \nu a^2\ell.           \tag{2.3}
\]

With \(\Gamma(x)=(4\pi|x|)^{-1}\),

\[
 p_1=\partial_i\partial_j\Gamma*S_{ij},
 \qquad u_2^{\rm nl}=-\nabla p_1                         \tag{2.4}
\]

outside the bubble.  If
\(d=|x-c|\ge C_0\ell\), the Newton-kernel derivative bound gives

\[
\begin{aligned}
 |\nabla^m u_2^{\rm nl}(x)|
 &\le C_m d^{-4-m}\|S\|_1\\
 &\le C_m\nu a^2\ell d^{-4-m},                           \tag{2.5}
\end{aligned}
\]

which proves (1.1).  Multiplication by
\(\tau^2=\ell^2/a^2\) proves (1.2).

The leading coefficient is generally nonzero.  Indeed,

\[
 \int S=-2\nu{\mathsf G}(W_{a,\ell,c,Q}),\qquad
 {\mathsf G}(W_{a,\ell,c,Q})
 =a^2\ell Q{\mathsf G}(W)Q^T,                            \tag{2.6}
\]

and a thin Gavrilov seed has nonisotropic \({\mathsf G}\).  Estimate
(2.5) is therefore not merely bounding an artifact that symmetry already
removes for each individual bubble.

There is also an exact global scaling for the nonlocal nonlinear part.
For a fixed profile \({\cal F}\in H^\infty_\sigma(\mathbb R^3)\),

\[
       \tau^2u_2^{\rm nl}(x)
       ={\nu\over\ell}Q{\cal F}
          \!\left(Q^T{x-c\over\ell}\right).              \tag{2.7}
\]

Consequently

\[
              \|\tau^2u_2^{\rm nl}\|_2
              \lesssim \nu\ell^{1/2}.                   \tag{2.8}
\]

A geometric family therefore has an absolutely summable sequence of
\(L^2\) norms.  Equation (2.8) concerns the pressure-generated nonlinear
part; the compact \(\nu^2\Delta^2W\) part of the second jet has a different
amplitude scaling and is not being hidden in this statement.

## 3. Separated geometric cascade

Assume the stage-\(j\) bubbles lie in separated annuli of radius
\(\ell_j=r^{-j}\).  There is a fixed \(c_*>0\) such that an active point in
annulus \(n\) has distance at least \(c_*\ell_j\) from every stored outer
annulus \(j<n\).  The use of one stage time in (1.2) is also consistent
with the remaining cascade time: for the ledger's geometrically decreasing
\(\tau_j\),

\[
                    \sum_{k\ge j}\tau_k\lesssim_r\tau_j.  \tag{3.0}
\]

Thus a second Taylor coefficient evolved all the way to the terminal time
still carries only a fixed multiple of \(\tau_j^2\).  This observation does
not control the higher Taylor coefficients.  Applying (1.2) gives

\[
 \sum_{j<n}|\nabla^m Z_j^{(2)}(x)|
 \lesssim_{m,r,c_*}
 \nu\sum_{j<n}\ell_j^{-1-m}
 \lesssim_{m,r,c_*}\nu\ell_n^{-1-m}.                    \tag{3.1}
\]

For a one-bubble active field whose natural derivative size is
\(a_n\ell_n^{-m}\), (3.1) gives

\[
 {\sum_{j<n}|\nabla^m Z_j^{(2)}(x)|
  \over a_n\ell_n^{-m}}
 \lesssim {\nu\over a_n\ell_n}
 ={1\over {\rm Re}_n}.                                   \tag{3.2}
\]

The constant deteriorates as the fixed annular gap closes, for example as
\(r\downarrow1\).  The construction only needs one fixed \(r>1\), but this
loss has to be included in any near-identity inverse theorem.

## 4. Packed-bubble carrier

At stage \(j\), take

\[
 N_j\lesssim K_j^3,qquad
 \delta_j\asymp{\ell_j\over K_j},qquad
 \tau_j={\ell_j\over a_j}.                               \tag{4.1}
\]

For one microbubble, (2.5) multiplied by the **macro** time \(\tau_j^2\)
has coefficient

\[
 \nu a_j^2\delta_j\tau_j^2
 =\nu\delta_j\ell_j^2.                                  \tag{4.2}
\]

Summing absolute values over all \(N_j\) bubbles proves

\[
 |\nabla^m Z_j^{(2)}(x)|
 \lesssim_m \nu\ell_j^3K_j^2d_j(x)^{-4-m}.              \tag{4.3}
\]

No orientational or multipole cancellation is used.  In the separated
active annulus \(n\), the polynomial factor does not spoil geometric
dominance:

\[
 \sum_{j<n}|\nabla^mZ_j^{(2)}(x)|
 \lesssim_{m,r,A,j_0}
 \nu K_n^2\ell_n^{-1-m}.                                 \tag{4.4}
\]

The active carrier derivative scale is
\(a_n(K_n/\ell_n)^m\).  Hence

\[
 {\sum_{j<n}|\nabla^mZ_j^{(2)}(x)|
  \over a_n(K_n/\ell_n)^m}
 \lesssim
 \nu\ell_n^{\gamma-1}K_n^{2-\gamma-m}.                  \tag{4.5}
\]

For every fixed \(m\), the right side tends to zero exponentially in
\(n\), up to a polynomial.  At \(m=0\), it is exactly
\(K_n/{\rm Re}_{n,{\rm carrier}}\), modulo fixed adjacent-scale constants.

A crude global \(L^2\) triangle bound also remains summable.  One
microbubble's macro-time coefficient has \(L^2\) size
\(O(\nu\ell_j^{1/2}K_j^{3/2})\); summing over \(K_j^3\) bubbles gives

\[
                  O(\nu\ell_j^{1/2}K_j^{9/2}).            \tag{4.6}
\]

The series over \(j\) converges for polynomial \(K_j\).  This deliberately
uses the triangle inequality rather than assuming orthogonality of
overlapping pressure tails.

## 5. Cross stress is small, not absent

Before the global velocity wake is introduced, disjoint bubble velocities
have no pointwise cross stress:

\[
        W_i\otimes W_j=0\qquad(i\ne j).                   \tag{5.1}
\]

After a nonlocal wake \(z\) is retained, the tensors
\(W\otimes z+z\otimes W\) are present.  If on a core of volume
\(O(\ell^3)\),

\[
             \|z\|_\infty\le\epsilon a,qquad
             \|W\|_\infty\lesssim a,                    \tag{5.2}
\]

then

\[
        \|W\otimes z+z\otimes W\|_1
        \lesssim\epsilon a^2\ell^3.                      \tag{5.3}
\]

At a distance \(d\gtrsim\ell\), its pressure-gradient derivatives are
therefore bounded by

\[
        C_m\epsilon a^2\ell^3d^{-4-m}.                   \tag{5.4}
\]

Thus cross pressure is perturbative if the incoming wake is perturbative
on the core.  Equations (3.2) and (4.5) provide that property at the
second-jet level.  They do not establish the inductive all-order statement;
that remains part of the missing Gevrey endpoint theorem.

## 6. Periodic-image caveat

On \(\mathbb T^3\), the periodic Green function in a coordinate ball can
be written

\[
                      \Gamma_{\mathbb T^3}=\Gamma+H,      \tag{6.1}
\]

where \(H\) is smooth.  The singular part obeys the same separated estimate
as (2.5).  The image part satisfies, for each fixed \(m\),

\[
 |\nabla^m(\tau^2u_{2,H}^{\rm nl})|
 \lesssim_m \nu\ell^3                                  \tag{6.2}
\]

for one macro bubble, and \(O(\nu\ell^3K^2)\) for a packed stage.  These
smooth image contributions are summable over geometrically shrinking
stages, but they do **not** decay as a power of the distance between two
sets.  They must be retained among the global low/centre modes of a torus
endpoint map.

Accordingly, the whole-space estimate cannot simply be quoted unchanged on
the torus.  The correct periodic statement is “Newton singular part plus a
summable smooth global image field.”

## 7. Claim boundary and next test

This note resolves only the scale objection raised by the forced
\(|x|^{-4}\) second-jet velocity tail:

* the tail is mandatory;
* it is finite-energy and summable over the proposed wake;
* it is perturbative at the next active core in the high-Reynolds window;
* packing loses \(K\) relative to carrier Reynolds number but remains
  asymptotically small for polynomial \(K_j\); and
* periodic images are smooth and summable but must be kept as global data.

The next load-bearing task is an inductive version of (4.5): prove that all
pressure/velocity jets through
\(M_j\asymp j^2/\log j\), including every wake--core cross stress, obey a
uniform \(C^M(M!)^2\) majorant after the multi-colour replacement.  Nothing
here proves that theorem.
