# C169: a fixed-plane terminal charge star has a passive-scalar reservoir bound

**Date:** 2026-08-06
**Status:** exact full nonlinear 2D3C Navier--Stokes reduction, maximum
principle, and conditional signed-\(q\)-mode gate-budget no-go; off-plane
three-dimensional conversion remains open
**Checker:**
[checks/fixed_plane_passive_star_no_go_c169.py](../checks/fixed_plane_passive_star_no_go_c169.py)

## 0. Claim boundary

C164 proves a pointwise-modulus obstruction on the fixed tangential
**linearized** fiber of a pure-normal gate.  The full nonlinear statement is
stronger and simpler whenever the source packet and all of its normal-charge
daughters remain in one fixed Fourier plane.  The complete velocity is then
2D3C: its in-plane part evolves autonomously by two-dimensional
Navier--Stokes, while its transverse part is one passive scalar.  Consequently
the transverse component cannot acquire a larger physical-space supremum
than was already stored in that scalar at the beginning of the terminal
pulse.

For a real gate with \(q_{\rm s}\) **signed** pure-normal Fourier
coefficients, let \(G_t\) be its transverse coefficient-\(\ell^2\) norm and
\(G\geq G_t\) its full coefficient norm.  Then

\[
 \|u(t)\mathbin\cdot e_t\|_{L^\infty}
 \leq \sqrt {q_{\rm s}}\,G_t
 \leq \sqrt {q_{\rm s}}\,G.                        \tag{0.1}
\]

The first coefficient bound is sharp, already at the initial time.

The estimate is independent of the number and amplitude of the **in-plane**
source coefficients.  Hence a C161-type endpoint lower bound

\[
 \|u(t_*)\mathbin\cdot e_t\|_{L^\infty}
 \geq c_*bq^{3/2}                                      \tag{0.2}
\]

In the literal C161 Section 2 convention, \(q_{\rm s}=q\), so this would
require \(G_t\geq c_*bq\), hence \(G\geq c_*bq\).  Under the additional
normalized-gate hypothesis
\(G\leq C_Gb\), (0.1)--(0.2) are incompatible once
\(q>C_G/c_*\).  On the C161 schedule, the upper scale is
\(b\sqrt q=n^2\), whereas the requested coherent point scale is
\(bq^{3/2}=n^{10}\); the missing factor is exactly \(q=n^8\).

Here \(G_t\) and \(G\) are norms, not energies: the corresponding
coefficient energies are \(G_t^2\) and \(G^2\).  If instead \(q\) denotes
the number of independent
half-lattice gate frequencies, reality gives \(q_{\rm s}=2q\).  With the
full signed transverse norm \(G_t\), its sharp middle bound reads
\(\sqrt{2q}\,G_t\); with the half-lattice transverse norm \(G_{t,+}\),
for which \(G_t=\sqrt2G_{t,+}\), it reads
\(2\sqrt q\,G_{t,+}\).  These are exact convention factors and do not change
the powers below.

This is an exact implication, not an unconditional no-go for every C161
realization.  C161 explicitly does not identify its abstract star angle with
a physical gate coefficient.  A fixed-plane realization could satisfy
(0.2) only by prestoring transverse \(L^\infty\) of that size, and a
\(q\)-mode gate could do so only with transverse coefficient-\(\ell^2\)
reservoir norm \(G_t\gtrsim bq\), hence full norm \(G\gtrsim bq\).  Such a
reservoir is outside the \(G=O(b)\) hypothesis (and has transverse
coefficient energy \(G_t^2\gtrsim b^2q^2\)).
Localization may also add more than \(q\) gate frequencies, but it cannot
evade the underlying necessary condition that the terminal transverse
supremum was present initially.

The result applies only to a fixed-plane, pure-normal, 2D3C terminal
realization and only to its fixed transverse component.  It does not control
an in-plane point value.  An off-plane leaf, an actively rotated plane, or a
genuinely three-dimensional noncommuting conversion can leave this invariant
class and is not ruled out.

## 1. Exact full nonlinear triangular system

Let \((e_r,e_t,e_z)\) be an orthonormal frame and let every field be
independent of the transverse coordinate \(x_t=e_t\cdot x\).  Write

\[
 u(x,t)=v_r(x_r,x_z,t)e_r+v_z(x_r,x_z,t)e_z
          +\Theta(x_r,x_z,t)e_t
       =v+\Theta e_t,                               \tag{1.1}
\]

with

\[
 \partial_rv_r+\partial_zv_z=0.                    \tag{1.2}
\]

Because \(\partial_{x_t}=e_t\cdot\nabla=0\),

\[
 (u\cdot\nabla)u=(v\cdot\nabla_{rz})v
                  +(v\cdot\nabla_{rz}\Theta)e_t.  \tag{1.3}
\]

In particular there is no term containing \(\Theta\) in the in-plane
equation.  The pressure Poisson source also contains no \(\Theta\): every
summand with a transverse velocity index contains a transverse derivative
and vanishes.  The mean-zero periodic pressure is therefore independent of
\(x_t\).  Substitution into unforced Navier--Stokes gives exactly

\[
 \boxed{
 \begin{aligned}
  \partial_tv+v\cdot\nabla_{rz}v+\nabla_{rz}p
    &=\nu\Delta_{rz}v,\\
  \partial_t\Theta+v\cdot\nabla_{rz}\Theta
    &=\nu\Delta_{rz}\Theta.
 \end{aligned}}                                    \tag{1.4}
\]

Thus \(v\) is autonomous and \(\Theta\) is passive.  More explicitly,
divergence of the momentum equation gives

\[
 \Delta_{rz}p=-\sum_{i,j\in\{r,z\}}
    (\partial_i v_j)(\partial_jv_i),               \tag{1.5}
\]

and the unique mean-zero periodic solution is independent of \(x_t\).
This is the general fixed-plane reduction already underlying C140, now
applied to the C161/C164 terminal charge geometry.  Smooth periodic data
remain in this class for all time.  No smallness assumption is imposed on
\(v\).

For \(\nu=0\), if \(X(a,t)\) is the Lagrangian flow of \(v\), then

\[
 \Theta(X(a,t),t)=\Theta_0(a),
 \qquad \|\Theta(t)\|_{L^p}=\|\Theta_0\|_{L^p}
 \quad(1\leq p\leq\infty).                         \tag{1.6}
\]

For \(\nu>0\), the real maximum principle gives

\[
 \boxed{\|\Theta(t)\|_{L^\infty}
       \leq\|\Theta_0\|_{L^\infty}.}              \tag{1.7}
\]

For a complex Fourier representative the same conclusion follows from

\[
 (\partial_t+v\cdot\nabla_{rz}-\nu\Delta_{rz})|\Theta|^2
 =-2\nu|\nabla_{rz}\Theta|^2\leq0.                \tag{1.8}
\]

Here the advecting field \(v\) is real; only the passive scalar has been
complexified.  The real maximum principle follows directly at a spatial
maximum and minimum, while (1.8) applies the same principle to
\(|\Theta|^2\).

The \(L^2\) identity

\[
 {1\over2}{d\over dt}\|\Theta\|_2^2
 =-\nu\|\nabla_{rz}\Theta\|_2^2                  \tag{1.9}
\]

also shows that the plane flow can redistribute transverse Fourier
coefficients but cannot increase their squared \(\ell^2\) norm (the
transverse coefficient energy).

## 2. Application to a \(q_{\rm s}\)-shift gate and a plane source

Let \(S=-S\) be any finite signed set of wavevectors in
\(\operatorname{span}\{e_r,e_z\}\), and let the in-plane source be real:

\[
 V_S(x)=\sum_{s\in S}\widehat V_s e^{is\cdot x},
 \qquad \widehat V_s\in\operatorname{span}\{e_r,e_z\},
 \qquad s\cdot\widehat V_s=0,
 \qquad \widehat V_{-s}=\overline{\widehat V_s}.  \tag{2.1}
\]

The coefficients in (2.1) are arbitrary.  In particular, their number and
amplitude may depend on \(q_{\rm s}\) in any way.  Let the real pure-normal
gate be

\[
 W_G(x)=\sum_{a\in\mathcal A}
   (\alpha_a e_r+\gamma_a e_t)e^{iax_z},
 \qquad |\mathcal A|=q_{\rm s},                    \tag{2.2}
\]

where \(\mathcal A=-\mathcal A\),
\(\alpha_{-a}=\overline{\alpha_a}\), and
\(\gamma_{-a}=\overline{\gamma_a}\).  Thus
\(q_{\rm s}=|\mathcal A|\) is the number of signed nonzero modes and is
even.  This is the literal convention of C161 Section 2.

With initial datum \(u_0=V_S+W_G\), the initial transverse scalar is only

\[
 \Theta_0(x_z)=\sum_{a\in\mathcal A}\gamma_a e^{iax_z}.
                                                               \tag{2.3}
\]

Define the transverse and full-gate coefficient-\(\ell^2\) norms by

\[
 G_t^2=\sum_{a\in\mathcal A}|\gamma_a|^2,
 \qquad
 G^2=\sum_{a\in\mathcal A}
        (|\alpha_a|^2+|\gamma_a|^2).               \tag{2.4}
\]

Thus \(G_t^2\) and \(G^2\), rather than \(G_t\) and \(G\), are the
corresponding coefficient energies.

Parseval gives \(G_t=\|\Theta_0\|_2\) in the normalized-torus convention,
and Cauchy--Schwarz gives the sharp finite-mode Bernstein bound

\[
 \|\Theta_0\|_\infty
 \leq\sum_a|\gamma_a|
 \leq\sqrt {q_{\rm s}}\,G_t
 \leq\sqrt {q_{\rm s}}\,G.                        \tag{2.5}
\]

If \(\mathcal A_+\) contains one member of each reality pair and
\(h=|\mathcal A_+|=q_{\rm s}/2\), then, exactly,

\[
 G_t^2=2G_{t,+}^2,
 \qquad
 \|\Theta_0\|_\infty\leq2\sqrt h\,G_{t,+}
       =\sqrt{2h}\,G_t.                            \tag{2.6}
\]

Combining (1.7) and (2.5) proves (0.1).  Notice what does **not** occur in
this estimate: neither \(|S|\) nor any norm of \(V_S\) appears.  The
plane packet changes the advecting flow, potentially by an arbitrarily
large amount, but it cannot create or amplify the total transverse
supremum.

If the desired terminal signal is defined after subtracting a gate-only
reference solution, let \(\Theta^S\) and \(\Theta^0\) start from the same
\(\Theta_0\) but be advected by the source-plus-gate and gate-only plane
solutions, respectively.  The two maximum principles and the triangle
inequality give the unconditional finite-difference bound

\[
 \|\Theta^S(t)-\Theta^0(t)\|_\infty
 \leq2\|\Theta_0\|_\infty
 \leq2\sqrt {q_{\rm s}}\,G_t
 \leq2\sqrt {q_{\rm s}}\,G.                        \tag{2.7}
\]

This factor-two estimate is about two full nonlinear solutions.  It is not
a maximum principle for their difference.  The factor two cannot be
improved from the two individual maximum principles alone: in the inviscid
problem, two translations of \(\Theta_0(z)=\cos z\) differing by \(\pi\)
attain opposite values at the same point.

More generally, if the source already contains a transverse scalar
\(\Theta_{S,0}\), then the honest bound is

\[
 \|\Theta(t)\|_\infty
 \leq\|\Theta_{G,0}+\Theta_{S,0}\|_\infty
 \leq\|\Theta_{G,0}\|_\infty
     +\|\Theta_{S,0}\|_\infty.                    \tag{2.8}
\]

There is then no no-go unless the initial transverse reservoir is itself
budgeted.  In particular a coefficient-\(\ell^2\) norm bound implies an
\(L^\infty\) bound only after paying the finite-mode factor
\(\sqrt{q_{\rm s}}\) in (2.5).

## 3. The C161 exponent consequence

In the fixed-plane specialization under test, take C161's post-gain
coefficient ledger: \(q^2\) independent half-lattice in-plane source
coefficients have scale \(c_0=b/q\), and the abstract normalized star asks
for \(q^3\) independent transverse daughters of scale \(c_0/\sqrt q\).
Reality completion doubles each of those two cardinalities, a fixed factor
absorbed into \(c_*\) below.  If the daughters were physically coherent in
one fixed transverse component, their coefficient-\(\ell^1\) scale would be

\[
 q^3{c_0\over\sqrt q}=bq^{3/2}.                   \tag{3.1}
\]

Suppose a fixed-plane physical realization claims the lower bound

\[
 \|\Theta(t_*)\|_\infty\geq c_*bq^{3/2}.          \tag{3.2}
\]

Equations (1.7) and (2.5) force the convention-independent bound

\[
 \boxed{G_t\geq {c_*bq^{3/2}\over\sqrt{q_{\rm s}}},
        \qquad G\geq G_t.}                         \tag{3.3}
\]

With C161 Section 2's literal signed-shift convention
\(q_{\rm s}=q\), this is \(G_t\geq c_*bq\).  If \(q\) counts half-lattice
gate frequencies, then \(q_{\rm s}=2q\) and the full signed transverse norm
must obey \(G_t\geq c_*bq/\sqrt2\).  For a source-induced difference
satisfying the same target, (2.7) weakens each lower bound by exactly a
factor two.  Every version is a factor comparable to \(q\) above a gate
norm budget \(G=O(b)\).

The hypothesis \(G=O(b)\) is an additional physical gate normalization;
C161 does not prove it, identify the abstract star angle with a physical
gate coefficient, or supply the coherence lower bound (3.2).  Thus (3.3)
is a necessary-condition obstruction to that proposed realization, not a
contradiction internal to C161's abstract ledger.

On C161's schedule

\[
 q=n^8,\qquad b=n^{-2},                            \tag{3.4}
\]

the three relevant powers are

\[
 b\sqrt q=n^2,\qquad bq=n^6,\qquad
 bq^{3/2}=n^{10}.                                  \tag{3.5}
\]

Thus a \(q\)-signed-mode gate with coefficient-\(\ell^2\) norm \(O(b)\)
is short by the full norm factor \(q=n^8\).  Meeting the target requires
transverse coefficient norm \(G_t\gtrsim bq=n^6\), hence transverse
coefficient energy \(G_t^2\gtrsim b^2q^2=n^{12}\).  The hypothesized full
gate norm budget \(G=O(b)\) has energy \(O(b^2)=n^{-4}\), so the
squared-energy deficit is \(q^2=n^{16}\).

## 4. Why this does not contradict the finite-fiber pulse calculations

Linearization around a 2D3C base \((\bar v,\bar\Theta)\) gives

\[
 \partial_t\delta\Theta
 +\bar v\cdot\nabla\delta\Theta
 +\delta v\cdot\nabla\bar\Theta
 =\nu\Delta\delta\Theta.                          \tag{4.1}
\]

The forcing \(\delta v\cdot\nabla\bar\Theta\) means that
\(\delta\Theta\) need not obey a maximum principle.  Therefore a
linearized or Galerkin-compressed polarization gain, such as the finite
first-neighbour mechanism in C166, is not algebraically inconsistent with
C169.  It measures sensitivity of the rearrangement of a pre-existing
transverse gate field.  The **total** nonlinear scalar still satisfies
(1.7), and the finite difference between two rearrangements satisfies
(2.7).

This distinction is load-bearing.  A terminal star cannot count a large
linearized transverse response as newly created physical point amplitude
unless it also pays for the nonlinear base scalar whose rearrangement
supports that response.

## 5. Exact surviving target

C169 rules out the simplest full nonlinear completion of the C161/C164
terminal star:

> a \(q^2\)-mode in-plane source packet plus a \(q\)-signed-mode
> pure-normal gate,
> all confined to one fixed 2D3C plane, with transverse gate
> coefficient-\(\ell^2\) norm \(O(b)\), cannot generate a transverse endpoint
> point value comparable to \(bq^{3/2}\).

The same geometry can survive only by leaving at least one hypothesis:

1. prestore transverse \(L^\infty\gtrsim bq^{3/2}\) (and, for \(q\) gate
   signed modes, transverse coefficient norm \(G_t\gtrsim bq\), equivalently
   transverse energy \(G_t^2\gtrsim b^2q^2\));
2. target an in-plane component and prove a different focusing theorem;
3. introduce off-plane wavevectors or a genuinely rotating spatial plane,
   making the converter fully three-dimensional;
4. add an initially transverse source component, while budgeting its
   \(L^\infty\) reservoir honestly.

Options 1 and 4 merely pay the reservoir unless coupled to a different
in-plane or off-plane target.  Localization by itself is not a loophole: if
it preserves the fixed-plane 2D3C class, (1.7) still applies; if it introduces
transverse-coordinate dependence, it is part of option 3.  No one-cell stage
or Millennium conclusion is claimed.
