# Reconstructed C114--C117: the terminal \(A_2\) triad, first-sideband cancellation, and the second-Picard leak

Date: 2026-08-05
Status: **reconstructed SELF / exact Fourier algebra; bare-gate closure refuted**
Checker: [`checks/terminal_triad_hexagon_c114_c117.py`](../checks/terminal_triad_hexagon_c114_c117.py)

## 0. Claim boundary

The old registry says only that the missing in-session sequence contained
"terminal triads C114--C117."  The original prose is unavailable.  This
file conservatively reconstructs the exact terminal-triad calculation that
feeds the later \(A_2\) hexagon ladder.  The numbering is reconstructed;
the identities are derived afresh and checked exactly.

The outcome has one positive and one negative half.  The integer \(A_2\)
terminal triad has exact wavevector closure, and two helicity channels can cancel the
reality-difference sideband while retaining the terminal output.  Thus the
desired child is isolated at first Picard order.  Once that child is active,
however, at least one named off-hexagon mode is forced at the next Picard
order.  A bare terminal triad is not a one-cell closure; a ladder, an
envelope cancellation, or spatial export is still required.

## 1. Fourier and helical conventions

For nonzero \(p,q\in\mathbb Z^3\), transverse polarizations
\(a\cdot p=b\cdot q=0\), and \(k=p+q\ne0\), define the symmetrized projected
Euler interaction

\[
 \mathcal N_k(p,a;q,b)
 :=-iP_k\{(a\cdot q)b+(b\cdot p)a\},
 \qquad
 P_kv=v-k{k\cdot v\over |k|^2}. \tag{1.1}
\]

This is the coefficient at \(k\) produced by the two Fourier modes, up to
the global sign convention for the evolution equation.

Take

\[
 \begin{aligned}
 k_1&=(1,-1,0),\\
 k_2&=(0,1,-1),\\
 k_3&=(-1,0,1),\\
 n&=(1,1,1).
 \end{aligned} \tag{1.2}
\]

For \(s\in\{-1,1\}\), use the unnormalised helical polarization

\[
 H_s(k):=n+i\,s\,{k\times n\over\sqrt2}. \tag{1.3}
\]

For the six roots \(\pm k_j\),

\[
 k\cdot H_s(k)=0,
 \qquad
 i k\times H_s(k)=s\sqrt2\,H_s(k),
 \qquad
 |H_s(k)|^2=6. \tag{1.4}
\]

Thus \(h_s(k)=H_s(k)/\sqrt6\) is unit normalized.  The convention also has
\(h_s(-k)=\overline{h_s(k)}\), exactly the Fourier reality relation.

## 2. C114 -- exact integer \(A_2\) hexagon

The vectors in (1.2) satisfy

\[
 k_1+k_2+k_3=0,
 \qquad |k_i|^2=2,
 \qquad k_i\cdot k_j=-1\quad(i\ne j). \tag{2.1}
\]

Therefore

\[
 \mathcal H=\{\pm k_1,\pm k_2,\pm k_3\}
\]

is an exact integer \(A_2\) root hexagon on the shell \(|k|=\sqrt2\).
Every sum of two distinct positive roots lands at the negative third root.
This is an exact terminal wavevector identity; there is no approximate
integer rounding.  No separate time-frequency resonance is asserted.

The three differences are off the terminal shell:

\[
 d_{12}=k_1-k_2=(1,-2,1),
 \quad d_{23}=k_2-k_3=(1,1,-2),
 \quad d_{31}=k_3-k_1=(-2,1,1), \tag{2.2}
\]

and \(|d_{ij}|^2=6\).

## 3. C115 -- exact helicity selection and terminal coefficient

Direct substitution in (1.1) gives, for either the terminal sum
\(k_1+k_2=-k_3\) or its cyclic rotations,

\[
 \mathcal N_{-k_3}
   (k_1,H_s(k_1);k_2,H_s(k_2))=0. \tag{3.1}
\]

This is the pairwise form of homochiral Beltrami darkness.

For opposite helicities the interaction is nonzero and especially simple:

\[
 \begin{aligned}
 \mathcal N_{-k_3}
  (k_1,H_{-}(k_1);k_2,H_{+}(k_2))&=+3\sqrt2\,n,\\
 \mathcal N_{-k_3}
  (k_1,H_{+}(k_1);k_2,H_{-}(k_2))&=-3\sqrt2\,n.
 \end{aligned} \tag{3.2}
\]

After normalizing both inputs, the output vector has squared norm

\[
 {54\over6^2}={3\over2}. \tag{3.3}
\]

Its projection onto each output helicity \(h_\pm(-k_3)\) has squared
magnitude \(3/4\), hence coefficient magnitude \(\sqrt3/2\).  A single
mixed-helicity pair excites both output helicities equally; it is not a
one-helicity terminal selector.

## 4. C116 -- two helicity channels cancel the first reality sideband

A real field containing the \(k_2\) mode also contains its conjugate at
\(-k_2\), with the same helical sign in convention (1.3).  A sign-pure
mixed-helicity pair \((k_1,-k_2)\) produces the named off-hexagon mode

\[
 d_{12}=k_1-k_2=(1,-2,1). \tag{4.1}
\]

For the same mixed-helicity assignments as (3.2), exact calculation gives

\[
 \begin{aligned}
 \mathcal N_{d_{12}}
  (k_1,H_{-}(k_1);-k_2,H_{+}(-k_2))&=-3\sqrt2\,n,\\
 \mathcal N_{d_{12}}
  (k_1,H_{+}(k_1);-k_2,H_{-}(-k_2))&=+3\sqrt2\,n.
 \end{aligned} \tag{4.2}
\]

Thus a **single** mixed-helicity channel has normalized sideband norm squared
\(3/2\), exactly the terminal norm in (3.3).  This is the sign-pure leakage
gate, but it is not yet an absolute no-go because each root may carry both
helicities.

Write the polarizations at the two positive roots as

\[
 a_i=x_iH_+(k_i)+y_iH_-(k_i),\qquad i=1,2. \tag{4.3}
\]

By (3.2)--(4.2), the desired terminal coefficient and the unwanted
difference coefficient are proportional respectively to

\[
 T_{12}=y_1x_2-x_1y_2,
 \qquad
 D_{12}=x_1\overline{y_2}-y_1\overline{x_2}. \tag{4.4}
\]

These two bilinear forms are not the same.  The exact choice

\[
 y_1=y_2=1,\qquad x_1=i,\qquad x_2=-i \tag{4.5}
\]

gives

\[
 D_{12}=0,qquad T_{12}=-2i\ne0. \tag{4.6}
\]

Self-products at \(2k_i\) vanish by transversality, opposite-root products
at zero vanish, and the conjugate difference is also cancelled by reality.
Consequently the real initial support \(\{\pm k_1,\pm k_2\}\) with (4.5)
has quadratic Euler output only at \(\pm(k_1+k_2)=\mp k_3\).  This is a
genuine first-order terminal gate, not an exact invariant subsystem.

There is no second unordered **wavevector** pair in \(\mathcal H\) summing
to \(d_{12}\); the cancellation in (4.6) uses the second helicity channel
on the same pair.  That distinction matters for the later polarization
and depletion ledgers.

## 5. C117 -- terminal isolation lasts only one Picard order

The output in (3.2) is parallel to \(n\).  Since

\[
 H_+(k)+H_-(k)=2n, \tag{5.1}
\]

the generated child at \(k_c=k_1+k_2=-k_3\) has equal helical amplitudes:
\(x_c=y_c\ne0\).  Its interaction with parent \(i\) emits the off-hexagon
mode

\[
 e_i=k_i+k_c=k_i-k_3,
 \qquad |e_i|^2=6, \tag{5.2}
\]

with coefficient proportional to

\[
 y_ix_c-x_iy_c=(y_i-x_i)x_c. \tag{5.3}
\]

Both secondary outputs \(e_1,e_2\) vanish only if
\(x_1=y_1\) and \(x_2=y_2\), but those equalities force
\(T_{12}=0\).  Therefore any nonzero terminal output forces at least one of

\[
 e_1=k_1-k_3=(2,-1,-1),
 \qquad e_2=k_2-k_3=(1,1,-2) \tag{5.4}
\]

when the child is fed back into the quadratic equation.  These are the
first unavoidable off-hexagon products relevant to the one-cell leakage
gate.  Within the six roots, \(e_1\) is produced only by
\(\{k_1,-k_3\}\) and \(e_2\) only by \(\{k_2,-k_3\}\), so there is no
second wavevector pair available for cancellation at this order.

In a short-time Taylor/Picard ledger, the desired child is \(O(t)\), while
the parent--child products in (5.3) integrate to \(O(t^2)\).  Viscosity is
diagonal in Fourier space and does not alter this support ordering.  Thus
the paired terminal gate delays mode leakage by one Picard order; it does
not cancel it to all orders.

This does **not** refute the later \(A_2\) ladder.  It explains why the
ladder is necessary: it must absorb (5.4), control their descendants, and
arrange terminal extraction only after the reality-symmetric evolution.
Nor does it rule out cancellation by additional modes, a spatial envelope,
or separated/unfolded gates.  Each is a larger PDE construction and must
carry its own energy, pressure, viscosity, and wake ledger.

For the one-cell theorem, the exact remaining obstruction is now named:

> **second-Picard terminal leakage:** the paired-helicity gate cancels the
> first difference sideband, but its equal-amplitude two-helicity child necessarily emits
> at least one of \(e_1,e_2\) on the next parent--child interaction.

## 6. Claim ledger

| ID | Reconstructed claim | Status |
|---|---|---|
| C114 | (1.2) is an exact integer equal-shell \(A_2\) terminal hexagon; its differences have squared length 6. | EXACT |
| C115 | Same-helicity terminal pairs are dark; opposite-helicity pairs have normalized output norm squared \(3/2\) and equal \(\sqrt3/2\) projections onto both output helicities. | EXACT |
| C116 | Two helicity channels can cancel the first reality difference mode while retaining a nonzero terminal output; (4.5) is explicit. | EXACT |
| C117 | The generated child has equal helicities and necessarily forces at least one named length-\(\sqrt6\) mode at the next Picard order; terminal leakage is delayed from \(O(t)\) to \(O(t^2)\), not removed. | EXACT formal support/Taylor ordering; ladder/envelope/export PDE closure OPEN |
