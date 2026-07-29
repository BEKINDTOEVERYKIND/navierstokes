# Adversarial audit of the eight-sideband homochiral star

Date: 2026-07-29

## Verdict and claim boundary

This memo independently checks the proposed common-helicity,
near-common-radius repair in
`2026-07-29-forward-multiphase-parametrix.md`.  It does not edit or rely on
that note's checker.

The leading algebra is substantially correct, but the strongest natural
reading of the claim is not.

1. A real divergence-free homochiral field can be built from the listed
   integer modes.
2. Equal-radius homochiral interactions cancel exactly.  The proposed
   matched modes are therefore necessarily **not** on one exact shell.
   Only the central modes are on one shell; the partner modes occupy three
   additional exact shells in an \(O(1)\)-thick annulus.
3. For a near-opposite pair \(Kw,q-Kw\), the leading child is indeed the
   single real polarization \(\mathbb P_qw\).
4. The eight normalized matrices
   \((\mathbb P_{q_a}w_a)\otimes q_a\) have rational determinant
   \(-3/10\), and the stated six-column diagonal synthesis is correct.
   But \(-3/10\) is not the determinant of the actual Fourier symbol.
   The physical leading symbol has seven additional factors \(1/2\), so
   its raw determinant is \(-3/1280\), up to controllable phase and
   amplitude gauges.
5. The determinant is a matched-child **point-jet** calculation.  The
   eight children live in seven distinct Fourier sectors; their matrices
   do not synthesize a spatially constant affine flow.  With sine phases
   they synthesize the desired gradient only at the common zero \(x=0\).
6. Every first quadratic interaction is \(O(1)\), as claimed.  However,
   there are unmatched low sideband differences and unmatched high
   products at the same order.  At the proposed nonzero six-sideband pump,
   the true endpoint Jacobian is not the displayed eight-column matrix.
7. More seriously, an explicit two-interaction chain creates an
   \(O(K)\) coefficient at the second Picard step.  Thus the common-annulus
   cancellation removes the first \(O(K)\) source but does not remove
   carrier stiffness from the forward chain.

Accordingly:

> The eight-sideband star survives as a legitimate **principal jet
> parametrix**, not as an exact common-shell transition module.  Its
> advertised rank is real, but the unproved off-shell forward-chain
> theorem remains essential and cannot be inferred from first-step
> transparency.

---

## 1. Fourier convention, reality, and divergence

For a nonzero integer wave vector \(p\), choose a unit helical vector
\(h_s(p)\), \(s\in\{+1,-1\}\), satisfying

\[
 p\cdot h_s(p)=0,\qquad
 i p\times h_s(p)=s|p|h_s(p),\qquad
 |h_s(p)|=1.
\tag{1.1}
\]

The phase gauge can be chosen pairwise so that

\[
 h_s(-p)=\overline{h_s(p)}.
\tag{1.2}
\]

For any finite positive-frequency set \({\cal P}\) and complex
coefficients \(A_p\), the field

\[
 u(x)=
 \sum_{p\in{\cal P}}
 \left(
 A_ph_s(p)e^{ip\cdot x}
 +\overline{A_p}\,\overline{h_s(p)}e^{-ip\cdot x}
 \right)
\tag{1.3}
\]

is real and divergence-free.  Conjugation preserves the helicity sign:
the mode at \(-p\) is still an eigenvector with eigenvalue \(s|p|\).

Thus there is no reality or incompressibility defect in the proposed
homochiral star.

For two coefficients \(a\) at \(p\) and \(b\) at \(r\), define the
projected quadratic symbol at \(n=p+r\) by

\[
 {\cal B}_{n}(a,b)
 =
 \mathbb P_n
 \left(
 (a\cdot r)b+(b\cdot p)a
 \right).
\tag{1.4}
\]

The actual Euler Fourier coefficient has the additional common factor
\(-i\).

---

## 2. Exact homochiral identity

The vector identity

\[
 \nabla(u\cdot v)
 =
 (u\cdot\nabla)v+(v\cdot\nabla)u
 +u\times\operatorname{curl}v
 +v\times\operatorname{curl}u
\tag{2.1}
\]

gives, after Leray projection,

\[
 \boxed{
 i{\cal B}_{p+r}(h_s(p),h_s(r))
 =
 s(|p|-|r|)
 \mathbb P_{p+r}
 \left(h_s(p)\times h_s(r)\right).
 }
\tag{2.2}
\]

Two consequences are exact.

First,

\[
 |p|=|r|
 \quad\Longrightarrow\quad
 {\cal B}_{p+r}=0.
\tag{2.3}
\]

Any real superposition supported on one radius and one helicity is a
Beltrami curl eigenfield, so its entire projected Euler nonlinearity
vanishes.

Second,

\[
 |{\cal B}_{p+r}|
 \le \big||p|-|r|\big|.
\tag{2.4}
\]

Hence an \(O(1)\)-thick homochiral annulus suppresses the generic
\(O(K)\) first interaction to \(O(1)\).

---

## 3. Near-opposite child symbol from scratch

Take

\[
 p=Kw,\qquad r=q-Kw,
\qquad w,q\in\mathbb Z^3,\quad q\ne0.
\tag{3.1}
\]

Then

\[
 |p|-|r|
 =
 \frac{w\cdot q}{|w|}+O(K^{-1}).
\tag{3.2}
\]

Choose phases continuously so that, as \(K\to\infty\),

\[
 h_s(r)\longrightarrow
 \overline{h_s(p)}.
\tag{3.3}
\]

If \(e=w/|w|\), a helical frame gives

\[
 h_s(p)\times\overline{h_s(p)}
 =-is\,e.
\tag{3.4}
\]

Substituting (3.2)--(3.4) into (2.2) yields the actual raw symbol

\[
 \boxed{
 {\cal B}_{q}
 =
 -\frac{w\cdot q}{|w|^2}\,
 \mathbb P_qw
 +O(K^{-1}).
 }
\tag{3.5}
\]

Thus the proposed limiting polarization

\[
 c(w,q)=\mathbb P_qw
\tag{3.6}
\]

is correct, but it carries the nonzero scalar

\[
 g(w,q)=-\frac{w\cdot q}{|w|^2}.
\tag{3.7}
\]

The scalar can be absorbed into a sideband amplitude.  It cannot be
omitted when calling a determinant the determinant of the physical
symbol.

---

## 4. What is and is not on a common shell

The proposed data are

\[
 k_a=Kw_a,\qquad l_a=q_a-Kw_a,
\tag{4.1}
\]

where all \(|w_a|^2=2\).  For integer \(K\), every \(k_a,l_a\) is an
integer vector on \(\mathbb T^3\).  There is no lattice obstruction.

However,

\[
 |k_a|=|l_a|
 \quad\Longleftrightarrow\quad
 2K\,w_a\cdot q_a=|q_a|^2.
\tag{4.2}
\]

All listed \(w_a\cdot q_a\) are negative, whereas the right side is
positive.  Therefore no listed matched pair lies on one exact shell for
any positive \(K\).

More generally, even if (4.2) were arranged for some other data, (2.3)
would make the homochiral child exactly zero.  A nonzero homochiral beat
and an exact common shell are mutually exclusive.

For the listed modes the exact squared radii are:

\[
\begin{array}{c|c}
\text{modes}&|p|^2\\ \hline
k_1,k_2&2K^2\\
l_1,l_2&2K^2+2K+1\\
l_3,l_4,l_5,l_6,l_8&2K^2+2K+2\\
l_7&2K^2+4K+5.
\end{array}
\tag{4.3}
\]

Here \(k_1=\cdots=k_7=K(1,1,0)\) is one Fourier degree of freedom and
\(k_8=K(1,-1,0)\) is the other.  The seven repeated entries do not give
seven independent parent amplitudes; independence comes from the seven
distinct partners.

The accurate description is therefore a four-shell homochiral star in an
\(O(1)\)-thick annulus, with its two central parent modes on a common
shell.

---

## 5. The rank-eight calculation

Put

\[
 C_a=(\mathbb P_{q_a}w_a)\otimes q_a.
\tag{5.1}
\]

In the coordinates

\[
 (M_{11},M_{22},M_{12},M_{13},
   M_{21},M_{23},M_{31},M_{32})
\tag{5.2}
\]

on \(\mathfrak{sl}(3)\), direct rational elimination gives

\[
 \det[C_1\ \cdots\ C_8]=-\frac3{10}.
\tag{5.3}
\]

The displayed coefficients

\[
 \left(-2,-\frac52,1,1,\frac54,\frac54,0,0\right)
\tag{5.4}
\]

also satisfy

\[
 \sum_a\lambda_aC_a
 =
 \operatorname{diag}\left(-1,-\frac54,\frac94\right).
\tag{5.5}
\]

These two rational claims are correct.

For the physical leading symbols, (3.7) gives

\[
 (g_1,\ldots,g_8)
 =
 \left(
 \frac12,\frac12,\frac12,\frac12,
 \frac12,\frac12,1,\frac12
 \right).
\tag{5.6}
\]

Consequently the raw matched-symbol determinant tends to

\[
 \boxed{
 \det[g_1C_1\ \cdots\ g_8C_8]
 =
 -\frac3{1280},
 }
\tag{5.7}
\]

up to the freely chosen helical and complex-amplitude phase gauges.
It is still nonzero.  Thus this correction does not kill the principal
rank; it identifies what the exact \(-3/10\) actually measures:
amplitude-normalized limiting columns.

There are seven distinct child Fourier sectors, because

\[
 q_3=q_8=(-1,0,-1).
\tag{5.8}
\]

At that repeated sector the two limiting polarizations are independent.
Across distinct \(q_a\), however, Fourier orthogonality prevents the
matrices \(C_a\) from adding to one spatially constant matrix.  With sine
phases,

\[
 u_{\rm match}(x)
 =
 \sum_a\lambda_ac_a\sin(q_a\cdot x)
\tag{5.9}
\]

is real and divergence-free and obeys

\[
 u_{\rm match}(0)=0,\qquad
 \nabla u_{\rm match}(0)
 =
 \operatorname{diag}\left(-1,-\frac54,\frac94\right).
\tag{5.10}
\]

Away from \(x=0\), its gradient is not constant.  The calculation is a
full point-jet chart, not an exact affine child on a neighborhood.

---

## 6. Complete classification of first low products

Let the positive-frequency set consist of the two distinct central modes
and the eight partners, together with their reality conjugates.

The intended low products are

\[
 k_1+l_a=q_a,\quad 1\le a\le7,
\qquad
 k_8+l_8=q_8.
\tag{6.1}
\]

There are additional low products between partners with the same leading
direction:

\[
 l_a-l_b=q_a-q_b,
\qquad 1\le a,b\le7.
\tag{6.2}
\]

Their sizes follow immediately from (2.2) and the shell table (4.3).

* Within \(\{l_1,l_2\}\), they vanish exactly.
* Within \(\{l_3,l_4,l_5,l_6\}\), they vanish exactly.
* Between the first two and the next four, they are \(O(K^{-1})\), because
  the squared radii differ by only one.
* Between \(l_7\) and any of \(l_1,\ldots,l_6\), they are \(O(1)\).

Two of the last interactions feed an intended Fourier sector:

\[
 l_7-l_1=q_7-q_1=q_3,
\qquad
 l_7-l_3=q_7-q_3=q_1.
\tag{6.3}
\]

Their limiting symbols are nonzero.  Therefore the matched eight-column
matrix is not the full child Jacobian at a state where \(l_7\) and the
other partners are active.

At the particular diagonal synthesis (5.4), \(l_7\) and \(l_8\) have zero
coefficient.  The unmatched low output among the six active partners is
then only \(O(K^{-1})\).  This explains why the six-column pump is a
consistent **principal** state.  But the derivative in the seventh
control direction at that nonzero state includes its interactions with
all six active partners; those terms are absent from the determinant
(5.3).

All remaining nonzero pair sums have wave number \(O(K)\).  Formula (2.4)
and the annular bound show that their first symbols are \(O(1)\), while
same-exact-shell pairs vanish.  They are not absent: they form the first
off-shell wake.

---

## 7. An explicit second-step return of the \(O(K)\) stiffness

The loss of first-step \(O(K)\) terms is not inherited by generated
modes.  Here is a concrete chain using only the first parent/partner and
the second central direction.

Let

\[
 w_1=(1,1,0),\qquad
 w_2=(1,-1,0),\qquad
 q_1=(-1,0,0),
\tag{7.1}
\]

and start with the same positive helicity on

\[
 k_1=Kw_1,\qquad
 k_2=Kw_2,\qquad
 l_1=q_1-Kw_1.
\tag{7.2}
\]

The unmatched first product

\[
 n=k_2+l_1
 =K(w_2-w_1)+q_1
\tag{7.3}
\]

has an \(O(1)\), nonzero coefficient.  To see the leading value, use

\[
\begin{aligned}
 u&=(1,1,0)/\sqrt2,&
 v&=(1,-1,0)/\sqrt2,&
 z&=(0,0,1),\\
 h_u&=(-v+iz)/\sqrt2,&
 h_v&=(u+iz)/\sqrt2,&
 h_{-u}&=\overline{h_u}.
\end{aligned}
\tag{7.4}
\]

For \(s=+1\), direct substitution in (1.4) gives

\[
 {\cal B}_{n}(h_v,h_{-u})
 =
 \frac12 e_1+\frac{i}{2\sqrt2}e_3
 +O(K^{-1}).
\tag{7.5}
\]

The Euler-generated coefficient is

\[
 A_n'
 =
 -i{\cal B}_{n}
 =
 -\frac i2e_1+\frac1{2\sqrt2}e_3
 +O(K^{-1}).
\tag{7.6}
\]

Now interact this generated mode with the original \(k_1\) mode.  The
output wave vector is

\[
 m=n+k_1=Kw_2+q_1.
\tag{7.7}
\]

A second direct use of (1.4) yields

\[
 \boxed{
 {\cal B}_{m}(A_n',h_u)
 =
 \frac{iK}{4}(e_1+e_2)+O(1).
 }
\tag{7.8}
\]

Thus an \(O(1)\) first wake coefficient couples back to the annulus with
an \(O(K)\) coefficient at the next interaction.  In a Picard expansion,

\[
 A_n(t)=O(t),\qquad
 A_m(t)=O(Kt^2)
\tag{7.9}
\]

before energy-exchange cancellations are resummed.

This does not by itself prove exponential growth of the exact Galerkin
semigroup; the \(O(K)\) operator may contain skew or normal-form
structure.  It does prove that the first-step bound (2.4) cannot be
iterated and that no tame forward estimate follows solely from the
common-annulus geometry.  Since Euler has no dispersive time phase, any
successful repair must exhibit an explicit algebraic/energy cancellation
for chains such as (7.3)--(7.8).

---

## 8. What would complete or refute the repair

The next analytic gate should not be stated merely as control of
"unmatched products."  The explicit chain above shows what must be
controlled.

A valid forward-chain theorem needs:

1. a phase-resolved endpoint map for the **full** state, including every
   \(q_a-q_b\) low mode and every first high wake;
2. the Jacobian at the nonzero six-sideband pump, not at the
   sideband-free Beltrami bath;
3. a normal-form, symmetrizer, or invariant energy estimate that controls
   the \(O(K)\) second-generation couplings;
4. preservation of the rank-eight point-jet projection after those
   corrections; and
5. a proof that the outgoing low modes define a recurrent localized pump,
   rather than only the correct gradient at one instant and one point.

Until those statements are proved, the star is a worthwhile principal
ansatz but not a transition theorem and not evidence of a Navier--Stokes
breakthrough.

---

## 9. Important escape boundary: the first six modes form a cleaner single-parent strain star

The \(O(K)\) chain in Section 7 uses the second central direction
\(w_2=(1,-1,0)\).  It does **not** apply to the smaller subsystem
consisting only of pairs \(1,\ldots,6\), all of which share

\[
 w=(1,1,0),\qquad w\cdot q_a=-1.
\tag{9.1}
\]

This distinction is material if the return construction needs only a
symmetric trace-free strain, rather than an arbitrary element of
\(\mathfrak{sl}(3)\).

### 9.1 Exact symmetric rank

Let

\[
 S_a=\operatorname{sym}
 \left(
   (\mathbb P_{q_a}w)\otimes q_a
 \right),
 \qquad 1\le a\le6.
\tag{9.2}
\]

In the coordinates

\[
 (M_{11},M_{22},M_{12},M_{13},M_{23})
\tag{9.3}
\]

on \(\operatorname{Sym}_0(3)\), the five columns with zero-based indices

\[
 0,2,3,4,5
\tag{9.4}
\]

have determinant

\[
 \boxed{-\frac18.}
\tag{9.5}
\]

Hence the six single-parent children have exact normalized symmetric rank
five.  They also obey

\[
 \sum_{a=1}^6
 \left(-2,-\frac52,1,1,\frac54,\frac54\right)_a S_a
 =
 \operatorname{diag}\left(-1,-\frac54,\frac94\right).
\tag{9.6}
\]

For all six physical homochiral symbols, (3.7) gives the same scalar

\[
 g_a=\frac12.
\tag{9.7}
\]

Therefore the actual leading symmetric minor is

\[
 \left(\frac12\right)^5\left(-\frac18\right)
 =
 \boxed{-\frac1{256}},
\tag{9.8}
\]

still uniformly nonzero.  With unit common-parent amplitude, the physical
partner coefficients are simply twice those in (9.6).

Thus the second parent direction and the seventh and eighth sidebands are
not needed to synthesize the desired symmetric pump.

### 9.2 All unmatched high--high products are smaller

For the first six partners there are only three leading wave classes:

\[
 Kw,\qquad -Kw+q_a,\qquad Kw-q_a,
\tag{9.9}
\]

and every \(q_a\) has the same \(w\)-projection.

The intended near-opposite products satisfy

\[
 Kw+(-Kw+q_a)=q_a
\tag{9.10}
\]

and have an \(O(1)\) nonzero symbol by (3.5).

Every other initial high--high product is smaller.

* For two modes with the same leading sign, the helical vectors differ by
  \(O(K^{-1})\), so their cross product in (2.2) is \(O(K^{-1})\).
* For two partners with opposite leading signs,
  \[
    w\cdot(q_a-q_b)=0,
  \tag{9.11}
  \]
  so their radial difference is \(O(K^{-1})\).

Consequently,

\[
 \boxed{
 \text{all unmatched initial high--high symbols in the six-star}
 =O(K^{-1}),
 }
\tag{9.12}
\]

with many vanishing exactly because of the two exact partner-radius
groups in (4.3).  In particular, there is no analogue of the
\(k_2+l_1\) cross-angle source used in Section 7.

This is a genuine improvement over the full eight-column
\(\mathfrak{sl}(3)\) star.

### 9.3 What still remains

Once the \(O(1)\) low children have been produced, their interaction with
the high carrier has size \(O(K)\):

\[
 {\cal B}_{Kw+q}
 (u_q,h_s(Kw))
 =
 \mathbb P_{Kw+q}
 \left(
   (u_q\cdot Kw)h_s(Kw)
   +(h_s(Kw)\cdot q)u_q
 \right).
\tag{9.13}
\]

Generically \(u_q\cdot w\ne0\).  This is not the unwanted cross-family
mechanism of Section 7; it is the expected leading transport/eikonal
coupling between the new low pump and its high carrier.  In material
coordinates a divergence-free low flow may make this coupling
energy-neutral and resum it without an \(e^{CK}\) norm.  No contradiction
to such a tame transport reduction was found in this audit.

The remaining exact requirements are:

1. solve the coupled creation-and-transport problem while the children
   grow from zero;
2. use the full low field, not only its gradient at the core, in the
   material phase equation;
3. prove a Gevrey-tame endpoint map for the resulting one-parent charge
   chain; and
4. restore global helicity cancellation by a disjoint
   orientation-reversed copy, without overlapping opposite helicities.

The revised adversarial verdict is therefore:

> The full eight-star still has the explicit second-generation
> cross-family \(O(K)\) chain.  The first-six single-parent restriction
> avoids that chain, retains the exact rank-five strain chart, and is a
> materially more plausible principal transition ansatz.  Its remaining
> \(O(K)\) term is structured low--high transport, not an immediate
> algebraic obstruction.
