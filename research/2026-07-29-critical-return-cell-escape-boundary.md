# Critical parabolic return cells: the exact escape boundary

Date: 2026-07-29

## Verdict

The scaling-critical choice

\[
 \ell_j=r^{-j},\qquad a_j=\ell_j^{-1},\qquad
 \tau_j=\ell_j^2
\]

is a genuinely different possible route.  Viscosity, inertia, and the time
derivative all remain comparable at every stage, so a single exact
Navier--Stokes return cell would eliminate the vanishing-viscosity WKB
hierarchy used in the supercritical-amplitude program.

There is, however, a sharp import boundary.

> A compact, globally \(L^3\)-tight critical return cell which converges to
> a periodic orbit in logarithmic similarity time is ruled out by the
> backward discretely self-similar nonexistence theorem.

A smooth Clay force does not change that conclusion: after critical
rescaling its \(L^3\) norm tends to zero exactly like \(T-t\).

The only visible critical escape is a **non-tight persistent wake**.  Each
completed scale may retain finite energy of order \(\ell_j\), while
contributing an order-one amount to the cubed \(L^3\) norm.  The energy and
dissipation series converge, but the rescaled \(L^3\) mass escapes to
similarity infinity and grows linearly in the number of completed stages.
That regime lies outside the global-\(L^3\) profile hypothesis and is not a
fixed compact return cell.

This note proves neither existence of that wake transition nor a
Navier--Stokes singularity.  It says exactly which critical object remains
worth constructing and which numerical target is already excluded.

## 1. Exact forced similarity equation

Let \(s=T-t\), choose a prospective singular point \(x_*\), and set

\[
 y=\frac{x-x_*}{\sqrt{s}},\qquad
 \tau=-\log s,\qquad
 u(x,t)=s^{-1/2}U(y,\tau),\qquad
 p(x,t)=s^{-1}P(y,\tau).
\tag{1.1}
\]

For

\[
 \partial_tu-\nu\Delta u+(u\cdot\nabla)u+\nabla p=f,
 \qquad \operatorname{div}u=0,
\tag{1.2}
\]

direct differentiation gives

\[
\boxed{
 \partial_\tau U+\frac12U+\frac12y\cdot\nabla U
 -\nu\Delta U+(U\cdot\nabla)U+\nabla P
 =F(y,\tau),
}
\tag{1.3}
\]

where

\[
 F(y,\tau)
 =s^{3/2}f(x_*+\sqrt{s}\,y,T-s).
\tag{1.4}
\]

On the rescaled copy

\[
 \Omega_\tau=s^{-1/2}(\mathbb T^3-x_*),
\]

the change of variables \(x=x_*+\sqrt{s}\,y\) yields, for every
\(1\le p<\infty\),

\[
 \|F(\tau)\|_{L^p(\Omega_\tau)}
 =s^{\frac32-\frac{3}{2p}}\,
   \|f(T-s)\|_{L^p(\mathbb T^3)}.
\tag{1.5}
\]

In particular,

\[
\boxed{
 \|F(\tau)\|_{L^3(\Omega_\tau)}
 =s\,\|f(T-s)\|_{L^3(\mathbb T^3)}
 \longrightarrow0
}
\tag{1.6}
\]

for every force continuous at \(T\).  A \(C^\infty\) force gives the
corresponding local derivative convergence as well.  Thus a smooth
external force disappears from every critical blowup limit; it cannot
support a nonzero limiting periodic similarity orbit.

## 2. Tight critical recurrence is closed

The following is the precise consequence of the backward DSS
nonexistence theorem.

### Proposition 2.1 (tight critical-return obstruction)

Suppose a smooth forced solution on \(\mathbb T^3\times[0,T)\) has a
critical rescaling (1.1).  Assume that along shifts
\(\tau_n\to\infty\), the functions

\[
 U_n(y,\sigma)=U(y,\tau_n+\sigma)
\]

converge strongly enough to pass (1.3) to the limit, globally in
\(L^3_y\) and locally through two spatial derivatives and one time
derivative, to a time-periodic profile

\[
 U_\infty\in
 C^1(\mathbb R_\sigma;
 L^3(\mathbb R^3)\cap C^2(\mathbb R^3)).
\tag{2.1}
\]

Then \(U_\infty=0\).  In particular, it cannot be a nonzero compact
critical return cell.

#### Proof

By (1.6), the shifted forcing converges to zero in the critical
\(L^3\) space.  Passing to the limit in (1.3) gives an unforced periodic
solution of the backward Leray equation.  Chae's backward
discretely-self-similar nonexistence theorem applies under (2.1), and the
periodic profile is zero. \(\square\)

The role of global \(L^3\) convergence is essential.  Local convergence
alone does not prevent critical mass from moving to
\(|y|\to\infty\).

## 3. Why a fixed compact GPU cell is the wrong critical target

Suppose one searches for one normalized trajectory \(V(\sigma,y)\),
supported or rapidly decaying in a fixed \(y\)-region, whose outgoing
state is a scaled copy of its incoming state.  Repeating that trajectory
at \(\ell_j=r^{-j}\) makes the rescaled dynamics periodic in \(\tau\).
Compactness supplies global \(L^3\) tightness.  Subject to the regularity
needed to pass to the limit, Proposition 2.1 rules out a nonzero orbit.

Increasing spatial resolution or optimization time does not alter this
theorem.  A useful critical computation must therefore include a state
space whose similarity radius and retained history both grow with the
number of stages.  Optimizing only a fixed local box is aimed at the
excluded object.

This conclusion does not rule out:

1. a nonperiodic critical cascade;
2. a periodic local core with non-tight global wake;
3. a profile outside global \(L^3\); or
4. the \(\gamma>1\) inertial route, which is not critical DSS.

## 4. The non-tight wake has the correct scalar ledger

Let \(w_j\) be disjoint retained pieces with

\[
 \operatorname{diam}(\operatorname{supp}w_j)\asymp\ell_j,
 \qquad
 \|w_j\|_\infty\asymp\ell_j^{-1}.
\tag{4.1}
\]

Then

\[
\begin{aligned}
 E_j:=\|w_j\|_2^2
 &\asymp \ell_j,\\
 \|w_j\|_3^3
 &\asymp1,\\
 \tau_j
 &\asymp\ell_j^2,\\
 D_j:=\nu\tau_j\|\nabla w_j\|_2^2
 &\asymp\nu\ell_j.
\end{aligned}
\tag{4.2}
\]

Consequently,

\[
 \sum_jE_j<\infty,\qquad
 \sum_jD_j<\infty,\qquad
 \sum_j\tau_j<\infty,
\tag{4.3}
\]

while after \(n\) disjoint pieces have been retained,

\[
 \left\|\sum_{j\le n}w_j\right\|_3^3
 \asymp n.
\tag{4.4}
\]

The \(L^3\) growth is not a defect.  Boundedness of the critical norm would
exclude a first singularity.  Equation (4.4) supplies exactly the
non-tightness needed to evade Proposition 2.1.

To see it in similarity variables at \(s_n\asymp\ell_n^2\), an older piece
\(w_j\), \(j<n\), appears at similarity radius and diameter

\[
 \frac{\ell_j}{\ell_n}=r^{n-j}.
\tag{4.5}
\]

Critical scaling preserves its \(L^3\) contribution.  Thus every fixed
similarity ball sees only finitely many recent pieces, while an increasing
amount of \(L^3\) mass sits at radii tending to infinity.  Global
\(L^3\) compactness fails in the exact manner left open by the theorem.

## 5. Flat-force cost at critical scaling

Let a normalized approximate cell have residual of size \(\delta_j\).
At critical scaling, the undifferentiated physical equation has size
\(\ell_j^{-3}\), and \(m\) spatial derivatives give

\[
 \|\nabla_x^m f_j\|_\infty
 \lesssim
 \delta_j\ell_j^{-3-m}.
\tag{5.1}
\]

Since \(T-t_j\asymp\ell_j^2\), extension by zero with all time and space
derivatives flat requires, for every fixed \(m,N\),

\[
 \delta_j\ell_j^{-3-m}
 =O(\ell_j^{2N}).
\tag{5.2}
\]

Equivalently, \(\delta_j\) must beat every exponential in \(j\).
A bound

\[
 \delta_j\le e^{-c j^2}
\tag{5.3}
\]

is sufficient.  Therefore the critical route does not remove the
all-order accuracy requirement unless the transition is an exact
Navier--Stokes orbit.  It removes the small-viscosity hierarchy, but
replaces the compact return problem by a globally growing wake problem.

## 6. Revised critical theorem target

The smallest critical theorem not covered by Proposition 2.1 is:

> Construct a sequence of forward parabolic transitions on expanding
> similarity domains which transfers a fixed energy fraction to the next
> inner scale, retains an order-one cubed \(L^3\) wake at every completed
> scale, has summable energy and viscous loss, and matches the growing
> outer wake and pressure jets with residual \(e^{-c j^2}\).

This is substantially different from finding one recurrent compact cell.
It is also different from the current \(\gamma>1\) homochiral
Kelvin--Reynolds parametrix, where viscosity is asymptotically small and
the outer wake is built from exact steady Euler pieces.

No new GPU run is warranted until one of these two finite systems is
specified:

* the growing-domain critical wake return operator; or
* the finite charged homochiral transition operator in the inertial route.

## Primary sources

* D. Chae, *Remarks on the asymptotically discretely self-similar
  solutions of the Navier--Stokes and the Euler equations*,
  arXiv:1306.0305:
  https://arxiv.org/abs/1306.0305
* C. L. Fefferman, *Existence and smoothness of the Navier--Stokes
  equation*, official Clay problem description:
  https://www.claymath.org/wp-content/uploads/2022/06/navierstokes.pdf
