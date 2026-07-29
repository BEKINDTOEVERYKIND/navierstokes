# Audit of the claimed reduction of 3D Navier--Stokes to axisymmetric swirl

Date: 2026-07-29

## Verdict

The claimed reduction in R. Shahmurov,
[*Hypothetical Singularity of 3D Navier-Stokes in Clay Institute set up
Reduces to Axisymmetric with Swirl class*](https://arxiv.org/abs/2606.07875v1),
is not established by version 1.

There is an exact first logical failure: the proof declares the descendant
scale order to be well founded even though the dyadic parabolic scales have
the infinite chain

\[
 R,\quad R/2,\quad R/4,\quad R/8,\ldots .
\]

This is not a cosmetic omission.  Excluding such an infinite descent is the
regularity problem itself.

The scalar vorticity-amplitude identity used at the start of the paper is
correct.  The audit therefore does not reject the manuscript by association
with the author's separate axisymmetric-profile preprint; it isolates the
failure inside this manuscript's own terminal-selection argument.

## 1. The circular terminal selection

Definition 9.1 calls a descendant strict if it improves, in lexicographic
order,

\[
 \text{scale}\prec\text{terminal time}\prec
 \text{spatial lattice position}\prec\text{active-hull rank},
\]

and states that the scale coordinate is dyadic and that these coordinates
are well founded.  It then selects a terminal sequence for which no strict
descendant remains.

If improvement in scale means passage to the smaller child packet, as it
must throughout the descendant construction, the assertion is false.
Writing \(R_n=2^{-n}R_0\), every \(R_n\) has the strict child \(R_{n+1}\).
There is no minimal dyadic scale.  Reversing the order would make it well
founded only by declaring the smaller packets not to be improvements, in
which case it cannot eliminate the cascade used later in the proof.

Theorem 9.8 uses precisely the missing conclusion.  Its proof says that if
an active descendant has score at least the singularity-forcing threshold,
the well-founded order selects a strict singularity-forcing descendant
"contrary to terminality."  But a first singular point supplies
singularity-forcing packets at arbitrarily small radii.  A diagonal
subsequence can extract limits at successive scales; it cannot manufacture
a final smallest scale.

Thus

\[
 \text{arbitrarily small active descendants}
 \quad\Longrightarrow\quad
 \text{a terminal packet with no descendant}
\]

is assumed, not proved.

## 2. The later finite rank does not repair the scale descent

Definition 10.9 introduces four selected-output ranks and calls a strict
descendant a rank-zero exit.  Theorem 10.18 then observes that the outer
rank list

\[
 4>3>2>1>0
\]

is finite.

That proves no termination of the inner process.  An iteration may exit
rank four to "strict descendant," rescale, and restart at rank four on the
child packet indefinitely.  Treating the whole infinite dyadic chain as
rank zero changes its label, not its well-foundedness.  Lemma 12.3 again
refers back to the same terminal selection instead of supplying a monotone
quantity whose positive decrement has a finite total budget.

## 3. Classification of failures is not exclusion of failures

The symmetry part has a second independent gap.  Definition 9.9 declares
the only coherent terminal classes to be

1. a fixed-vorticity-direction class, and
2. a one-axis-equivariant class.

Definition 9.10 then calls active mass outside a neighborhood of one of
those two classes "fragmentation."  Consequently, Theorem 9.17's statement
that a zero-fragmentation terminal measure lies in one of the two classes
is essentially built into the definitions.

What is needed is a quantitative theorem excluding persistent positive
fragmentation.  Lemma 9.13 instead argues that any positive
axis-equivariance defect is a named output, fragmentation, or a strict
descendant.  Theorem 10.7 says that repeated outputs either spend a finite
currency, create a descendant, or yield another ancient output profile;
its proof defines positive measures but does not establish a global
summability inequality for all of them.  Theorem 10.8 then routes those
profiles by a prose priority list.  The strict-descendant branch returns to
the non-well-founded scale chain above.

The same issue appears transparently in Lemma L.6.  The "shape output" is
defined as the tail oscillation in a topology containing the desired
strong \(L^3_{\rm loc}\) velocity and \(L^2_{\rm loc}\) vorticity
distances.  Its vanishing is therefore exactly the desired Cauchy
compactness.  Calling failure of Cauchy compactness an output does not prove
that the failure cannot occur.

## 4. What remains valid

For \(\omega=A\xi\), the regularized limiting identity

\[
 (\partial_t+u\cdot\nabla-\nu\Delta)A
 =
 A\bigl(\xi\cdot S\xi-\nu|\nabla\xi|^2\bigr)
\]

away from the zero set, with the usual inequality/regularization at
\(A=0\), follows correctly from the vorticity equation.  It is a useful
starting identity, but it does not impose an axis of symmetry.

The actual missing theorem would have to supply at least one of:

- a genuinely well-founded scale functional;
- a positive descendant cost bounded by a finite critical norm;
- a compactness theorem that excludes recurrent loss at successively
  smaller scales; or
- a rigidity theorem deriving rotational equivariance from the PDE rather
  than defining all non-equivariant data as an output.

None is proved in version 1.  Therefore this preprint does not provide a
shortcut to either global regularity or the constructive singularity route.

