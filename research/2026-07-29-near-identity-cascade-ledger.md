# A uniformly near-identity polynomial-carrier cascade

Date: 2026-07-29

## Status

This note repairs the quantitative ledger behind the analytic
near-identity proposal.  It does **not** prove the missing endpoint
return theorem.

The important change is to shift the polynomial carrier:

\[
 \ell_j=r^{-j},\qquad K_j=(j+j_0)^A,\qquad
 k_j=\frac{K_j}{\ell_j}.
\tag{0.1}
\]

With \(r>1\) close to one and \(j_0\) large, the complete physical
carrier handoff is uniformly close to the identity from the very first
stage retained in the infinite tail.  One does not need to make the
normalized stage duration tend to zero, and the physical stage time is
the inertial time \(\ell_j/a_j\), not \(\ell_j^2\).

---

## 1. Exact one-step ratios

Write

\[
 \delta=\log r>0,\qquad n_j=j+j_0.
\tag{1.1}
\]

Then

\[
 q_j:=\frac{k_{j+1}}{k_j}
 =r\left(1+\frac1{n_j}\right)^A
\tag{1.2}
\]

and the logarithmic deformation required on stage \(j\) is exactly

\[
 \sigma_j:=\log q_j
 =\delta+A\log\left(1+\frac1{n_j}\right).
\tag{1.3}
\]

The elementary inequalities

\[
 \frac1{n+1}\leq\log\left(1+\frac1n\right)\leq\frac1n
\tag{1.4}
\]

give

\[
 \delta+\frac{A}{n_j+1}
 \leq \sigma_j
 \leq \delta+\frac{A}{n_j}.
\tag{1.5}
\]

Hence, for any prescribed \(\varepsilon_*>0\), choosing

\[
 0<\delta\leq\frac{\varepsilon_*}{2},
 \qquad
 j_0\geq\frac{2A}{\varepsilon_*}
\tag{1.6}
\]

ensures

\[
 0<\sigma_j\leq\varepsilon_*
\quad\text{for every }j\geq0.
\tag{1.7}
\]

This removes the finite collection of large polynomial-carrier
handoffs which was hidden by the unshifted choice \(K_j=j^A\).

For the Kelvin amplitude exponent \(1<\gamma<3/2\), set

\[
 a_j=k_j^\gamma
 =\ell_j^{-\gamma}K_j^\gamma.
\tag{1.8}
\]

The amplitude recurrence is then exact:

\[
 \frac{a_{j+1}}{a_j}=q_j^\gamma
 =e^{\gamma\sigma_j}.
\tag{1.9}
\]

Thus both the carrier and amplitude targets differ from the preceding
ones by \(O(\varepsilon_*)\).

---

## 2. Fixed normalized duration is preferable

Let \(T_0>0\) be a fixed, sufficiently short normalized interval.  Use
the affine Kelvin deformation already recorded in the repository,
rescaled so that its integrated logarithmic compression is
\(\sigma_j\).  Its normalized strain amplitude is then

\[
 \lambda_j\asymp \frac{\sigma_j}{T_0}.
\tag{2.1}
\]

For fixed \(T_0\), (1.7) makes this a uniformly small perturbation.
This is better suited to a local endpoint inverse theorem than taking
\(T_j\asymp\sigma_j\): in the latter parametrization, a low strain which
is itself generated during the stage can lose an additional power of
\(T_j\).

Suppose the normalized analytic evolution has the following two
properties on a packet class \({\cal A}_{\rho_*}\):

1. its endpoint map has a right inverse bounded uniformly in \(j\) on
   the five-dimensional strain target after division by the fixed
   factor \(T_0\);
2. its analytic-radius loss satisfies
   \[
   \rho(T_0)\geq \rho(0)-C_*T_0
   \tag{2.2}
   \]
   in a fixed neighbourhood of the parent packet.

Choose \(T_0<\rho_*/(4C_*)\), then choose \(\varepsilon_*\) below the
corresponding inverse-function radius.  The target (1.3) is then small
enough for a quantitative analytic inverse theorem at every stage.
The required controls are \(O(\sigma_j)\), and the analytic radius
remaining at the endpoint is bounded below independently of \(j\).

This is a conditional reduction, not a proof of either property.  Its
point is that no parameter in the desired endpoint theorem now
deteriorates merely because \(j\) increases.

---

## 3. Correct physical time and summability

The inertial time unit is

\[
 \tau_j=\frac{\ell_j}{a_j}
 =\ell_j^{1+\gamma}K_j^{-\gamma}.
\tag{3.1}
\]

For the fixed normalized interval in Section 2, the physical stage
duration is

\[
 \Delta t_j=T_0\tau_j,
\tag{3.2}
\]

not \(O(\delta\ell_j^2)\).  If one instead chooses a normalized
duration \(T_j\asymp\sigma_j\), the right formula is
\(\Delta t_j\asymp\sigma_j\tau_j\).

Because

\[
 \tau_j=r^{-(1+\gamma)j}(j+j_0)^{-A\gamma},
\tag{3.3}
\]

the total time is finite:

\[
 \sum_{j=0}^\infty\Delta t_j<\infty.
\tag{3.4}
\]

The energy scale remains

\[
 E_j\asymp a_j^2\ell_j^3
 =r^{-(3-2\gamma)j}(j+j_0)^{2A\gamma},
\tag{3.5}
\]

so

\[
 \sum_j E_j<\infty
\quad\Longleftrightarrow\quad
\gamma<\frac32.
\tag{3.6}
\]

The blowup diagnostic

\[
 a_j\ell_j
 =r^{(\gamma-1)j}(j+j_0)^{A\gamma}
 \longrightarrow\infty
\tag{3.7}
\]

still requires \(\gamma>1\).  Hence the sharp window

\[
 1<\gamma<\frac32
\tag{3.8}
\]

is unchanged.

The constants in (3.4)--(3.6) can grow as \(r\downarrow1\).  The Clay
problem asks for one finite-energy construction, not estimates uniform
as \(\delta\to0\), so this is harmless once a single sufficiently small
positive \(\delta\) is fixed.

---

## 4. Flat Gaussian cutoff with the corrected time scale

Assume the uncut normalized active packet has Gaussian decay and cut it
at normalized radius \(L(j+j_0)\).  The physical cutoff radius is

\[
 R_j=L(j+j_0)\ell_j.
\tag{4.1}
\]

The tail is bounded by

\[
 \exp\big(-cL^2(j+j_0)^2\big).
\tag{4.2}
\]

A fixed spatial derivative can cost a power of the physical carrier

\[
 k_j=\ell_j^{-1}K_j
 =r^j(j+j_0)^A,
\tag{4.3}
\]

while a fixed time derivative can cost a power of

\[
 \tau_j^{-1}
 =r^{(1+\gamma)j}(j+j_0)^{A\gamma}.
\tag{4.4}
\]

For fixed \(m,n\),

\[
\begin{aligned}
 \log\!\left(k_j^m\tau_j^{-n}\right)
 &=
 \big(m+n(1+\gamma)\big)j\log r\\
 &\quad+A(m+n\gamma)\log(j+j_0),
\end{aligned}
\tag{4.5}
\]

which is \(O_{m,n}(j+\log j)\).  Therefore

\[
 k_j^m\tau_j^{-n}
 e^{-cL^2(j+j_0)^2}
 \leq C_{m,n}e^{-c' L^2j^2}.
\tag{4.6}
\]

The same quadratic margin can absorb the \(e^{Cj^2}\) loss of an
order-\(M_j\) analytic/Gevrey construction after increasing \(L\).
Thus correcting the physical time scale does not damage terminal
flatness.

Finally,

\[
 \frac{R_{j+1}}{R_j}
 =r^{-1}\left(1+\frac1{j+j_0}\right)<1
\tag{4.7}
\]

whenever \(j+j_0>(r-1)^{-1}\).  This can be required already at
\(j=0\) by increasing \(j_0\).

---

## 5. Refined theorem gate

The remaining local theorem can now be posed with all scale parameters
in a compact set:

> **Uniform analytic near-identity return gate.**  For some fixed
> \(T_0,\rho_*,\varepsilon_*>0\), every target deformation
> \(0<\sigma\leq\varepsilon_*\) and every sufficiently large integer
> carrier \(K\) admit a one-carrier material-phase transition on
> \([0,T_0]\) which:
>
> 1. realizes the exact carrier and amplitude ratios
>    \(e^\sigma,e^{\gamma\sigma}\);
> 2. returns to a uniformly bounded analytic packet class with radius
>    at least \(\rho_*\);
> 3. retains the global pressure/wake state needed at the endpoint;
> 4. has an endpoint right inverse uniform in \(K\) and \(\sigma\);
> 5. admits the order-\(M\) pressure, viscosity, and wake expansion with
>    constants \(C^M(M!)^2\).

The exact rank-five instantaneous chart, the all-generation
\(K\)-cancellation, and the pressure coefficient count established in
the companion notes are the principal-symbol ingredients for this
theorem.  What remains is the nonlinear analytic endpoint construction
and its localization/wake coupling.

