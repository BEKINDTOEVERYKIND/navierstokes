# C173: a paired multipole Piola collar cancels the leading active moments

**Date:** 2026-08-06
**Status:** exact paired-core and moment identities; exact finite-shell
Fourier kernel for the affine coefficient and its formal first-curvature
variation; conditional low-frequency active-fiber estimate; the broad-band
MCKC/BAFL theorem remains open
**Checker:**
[checks/paired_multipole_piola_collar_c173.py](../checks/paired_multipole_piola_collar_c173.py)

## 0. Claim boundary

This note stays on the C107/C143/C171 one-cell material-collar closure.
It tests whether two co-moving compact affine curl collars can cancel the
leading parent-cross wake without canceling the affine strain seen by the
intended child.

There is a genuine finite-dimensional kernel mismatch.

* Two concentric copies at radii \(r\) and \(2r\), with coefficients
  \(128/127\) and \(-1/127\), still equal the prescribed affine field in
  the inner core. Their local strains therefore add to the original
  strain. On the other hand, their complete third velocity-moment tensor
  cancels. For a symmetric affine core the paired field vanishes to at
  least Fourier order five at the origin, rather than at least order three.
* This gives two additional powers of \(r|k|\). On the parent-scale
  \(A_2\) fibers, where \(r|k|\asymp q^{-1}\), a fixed Fourier-coordinate
  parent-cross or wake-to-active block gains \(q^{-2}\) at fixed local-core
  normalization. Consequently
  \(b^2q^{-2}\log h=o(b^3)\) on the existing schedule.
* A fixed scalar differential notch can additionally kill any fixed
  finite collection of child-scale shells in material-label Fourier
  coordinates. A triple shell zero leaves
  the affine core unchanged and kills the constant-affine coefficient and
  its complete first curvature variation. For the named \(A_2\) child
  and wakes the multiplier has spatial order only twelve and a fixed
  coefficient sum \(343/64\).

These statements do not close MCKC. The paired construction does not make
the full \(L^2\) residual smaller; it only improves specified
low-frequency moments/fibers. The C144 endpoint needs a three-dimensional
broad band. An underlying-jet-independent scalar symbol that robustly
notches even the \(q\) axial members of that packet through their first two
Fourier jets requires differential order at least \(3q\), far beyond the
existing growing-order Gevrey budget. Thus a
fixed finite collar multipole can pay the named finite \(A_2\) return, but
it cannot be promoted to BAFL by declaring the entire broad active packet
dark. A non-polynomial dynamic slaving or export estimate remains open.

No unforced stage map, BAFL estimate, cascade, singularity, or Millennium
conclusion is claimed.

## 1. The compact affine curl family

Let \(\chi\) be a fixed smooth radial function, equal to one on \(B_1\)
and supported in \(B_2\). Let \(S=S^T\) be trace free. For \(R>0\), set

\[
 {\cal A}_R(x)=\chi(x/R){(Sx)\times x\over3},
 \qquad w_R=\nabla\times{\cal A}_R.                 \tag{1.1}
\]

As in C143/C171,

\[
 \nabla\cdot w_R=0,\qquad
 \operatorname {supp}w_R\subset B_{2R},\qquad
 w_R(x)=Sx\quad (|x|<R).                           \tag{1.2}
\]

If \(\|S\|\asymp\lambda\), scaling gives

\[
 \|\nabla^m w_R\|_2\le C_{\chi,m}\lambda
 R^{5/2-m}.                                        \tag{1.3}
\]

The symmetry assumption is the one used by the C142 affine selector. It
is also exactly what removes the first velocity moment below.

## 2. An exact two-radius moment cancellation

Define

\[
 \boxed{
 w_{\rm pair}={128\over127}w_r-{1\over127}w_{2r}.}
                                                            \tag{2.1}
\]

Both summands equal \(Sx\) on \(B_r\), and

\[
 {128\over127}-{1\over127}=1.
\]

Therefore

\[
 \boxed{w_{\rm pair}=Sx\quad\hbox{on }B_r.}        \tag{2.2}
\]

The negative outer coefficient cancels a collar moment, not the local
strain. The field is still real, compactly supported in \(B_{4r}\), and
divergence free.

Here is the exact moment calculation. Since \({\cal A}_R\) is even,
\(w_R\) is odd. Thus its velocity moments of even total degree vanish.
For the degree-one moment, integration by parts gives

\[
 \int x_\ell (w_R)_j\,dx
 =-\epsilon_{j\ell n}\int({\cal A}_R)_n\,dx.       \tag{2.3}
\]

Radiality gives \(\int\chi(x/R)x_a x_b\,dx=c_R\delta_{ab}\), so

\[
 \int({\cal A}_R)_n\,dx
 ={c_R\over3}\epsilon_{npq}S_{pq}=0               \tag{2.4}
\]

because \(S\) is symmetric. Hence every moment through degree two
vanishes.

Writing \(w_R(x)=R w_1(x/R)\), with \(S\) held fixed, shows that every
degree-\(d\) moment scales as \(R^{d+4}\). In particular, the full
degree-three tensor of (2.1) has coefficient

\[
 {128\over127}r^7-{1\over127}(2r)^7=0.             \tag{2.5}
\]

Parity also kills degree four. Consequently

\[
 \boxed{\widehat w_{\rm pair}(\xi)=O(\lambda r^4
                    (r|\xi|)^5)\quad(r|\xi|\le1).} \tag{2.6}
\]

The corresponding single-collar bound starts no earlier than
\(O(\lambda r^4(r|\xi|)^3)\). More generally, for
\(|\alpha|\le5\),

\[
 |\partial_\xi^\alpha\widehat w_{\rm pair}(\xi)|
 \le C_{\chi,\alpha}\lambda r^{4+|\alpha|}
 (r|\xi|)^{\max(5-|\alpha|,0)}.                   \tag{2.7}
\]

The norm condition of the pair is harmless. Before changing the fixed
cutoff constant, (1.3) gives

\[
 \|\nabla^m w_{\rm pair}\|_2
 \le {128+2^{5/2-m}\over127}
 C_{\chi,m}\lambda r^{5/2-m}.                     \tag{2.8}
\]

## 3. Affine Piola transport, curvature, and pressure moments

For an affine volume-preserving map \(x=c+Fa\), Piola transport gives

\[
 U(x)=F w_{\rm pair}(a),\qquad
 \widehat U(k)=e^{-ik\cdot c}F\widehat w_{\rm pair}(F^Tk). \tag{3.1}
\]

Thus (2.2), all moment cancellations, and their reality partners survive
an arbitrary common affine chart. If \(A=\nabla V\) is constant, the
leading C171 residual is

\[
 G_{pc}=2AU,\qquad \widehat{G}_{pc}(k)=2AF
 e^{-ik\cdot c}\widehat w_{\rm pair}(F^Tk).       \tag{3.2}
\]

Leray projection cannot undo a zero or enlarge a fixed Fourier
coefficient. Conditional on a uniformly nondegenerate affine chart, the
parent-scale \(A_2\) fibers obey \( |F^Tk|\asymp K\) and
\(rK=q^{-1}\). Equations (2.6)--(3.2) therefore gain \(q^{-2}\) relative
to the first nonzero moment of one symmetric collar. This statement is
uniform in polarization; in particular it applies simultaneously to the
real partners of

\[
 k_c=(1,0,-1),\qquad e_1=(2,-1,-1),\qquad
 e_2=(1,1,-2),                                    \tag{3.3}
\]

whose squared lengths are \(2,6,6\).

The same two powers survive the first spatial curvature jet. To see the
mechanism without hiding a cancellation, use \(a=ry\), put
\(\varepsilon=r/\ell\), and write, with dimensionless tensors bounded
independently of the stage,

\[
 \begin{aligned}
 X_\varepsilon(ry)&=c+r\left(Fy+
              {\varepsilon\over2}H[y,y]\right)+O(r\varepsilon^2),\\
 D_aX_\varepsilon(ry)&=F+\varepsilon H[y,\cdot]+O(\varepsilon^2),\\
 A_\varepsilon(X_\varepsilon(ry))
   &=A_0+\varepsilon A_1[y]+O(\varepsilon^2).
 \end{aligned}                                    \tag{3.4}
\]

Put \(\omega(y)=r^{-1}w_{\rm pair}(ry)\). For a fixed affine Kelvin
covector \(k_0\), put \(\zeta=rF^Tk_0\). Suppressing the common powers of
\(r\), the coefficient is

\[
 I_\varepsilon=\int e^{-ik_0\cdot X_\varepsilon(ry)}
 A_\varepsilon(X_\varepsilon(ry))D_aX_\varepsilon(ry)
 \omega(y)\,dy.
\]

Differentiation at zero gives exactly

\[
 \begin{aligned}
 I'_0=\int e^{-i\zeta\cdot y}\big\{&
 A_1[y]F\omega
 +A_0H[y,\cdot]\omega\\
 &-{i\over2}(rk_0\cdot H[y,y])A_0F\omega\big\}\,dy,
 \end{aligned}                                    \tag{3.5}
\]

after absorbing the normalization of \(A_1,H\) into the tensors. A first
variation of the selected covector adds only another linear-label term.
Thus (3.5) is a sum of terms containing

\[
 \widehat\omega,\qquad \partial_\zeta\widehat\omega,\qquad
 \zeta\,\partial_\zeta^2\widehat\omega.            \tag{3.6}
\]

At parent scale, \(|\zeta|\asymp q^{-1}\) and the explicit curvature
parameter also satisfies \(\varepsilon=q^{-1}\). Equation (2.7) then
makes \(\varepsilon\partial_\zeta\widehat\omega\) and
\(\varepsilon\zeta\partial_\zeta^2\widehat\omega\) order
\(|\zeta|^5\), while \(\varepsilon\widehat\omega\) is smaller. Thus the
\(q^{-2}\) improvement is not destroyed by the first curvature jet under
the bounded tensors and affine-chart comparability assumed in (3.4).
This is only the derivative at \(\varepsilon=0\); a uniform Taylor
remainder through the full gain window still belongs to MCKC(i).

The moment gain also improves, but does not remove, the global pressure
tail in fixed Euclidean coordinates. In the constant-affine case (3.2) is
\(O(|k|^5)\) at \(k=0\), so
\(\widehat{\nabla\cdot G}_{pc}=O(|k|^6)\). The resulting exterior
pressure gradient is \(O(d^{-8})\), versus the general \(O(d^{-6})\)
bound for one symmetric radial collar. Terms in the first
curvature variation contain at most one label factor, or a quadratic
label factor accompanied by one Fourier \(k\); they are \(O(|k|^4)\).
Their pressure gradient is therefore \(O(d^{-7})\), versus the general
\(O(d^{-5})\) single-collar bound. These exponents require the displayed compact source,
fixed-coordinate multipole expansion, and bounded affine/curvature
coefficients. They are exterior orders on \(\mathbb R^3\), not local
\(L^2\) gains, backward non-normal estimates, or a periodic-image
summation theorem.

## 4. A fixed finite-shell kernel through the formal first curvature variation

The paired moments handle parent-scale modes. A separate fixed-order
notch can make selected child-scale Fourier fibers exactly dark. For a
finite set of positive squared normalized radii \({\cal R}\), put

\[
 Q_{\cal R}(z)=\prod_{\rho^2\in{\cal R}}
       \left(1-{z\over\rho^2}\right)^3,
 \qquad
 w^\sharp=Q_{\cal R}(-r^2\Delta)w_{\rm pair}.      \tag{4.1}
\]

This is still a compact divergence-free curl field. Since
\(\Delta(Sx)=0\), \(Q_{\cal R}(0)=1\), and (2.2) holds on an open core,

\[
                 w^\sharp=Sx\quad\hbox{on }B_r.   \tag{4.2}
\]

At every \(\zeta\) with \(|\zeta|^2\in{\cal R}\), the Fourier multiplier
\(Q_{\cal R}(r^2|\xi|^2)\) and all of its derivatives through order two
vanish at \(\xi=\zeta/r\). It follows that

\[
 \partial_\xi^\alpha\widehat w^\sharp(\zeta/r)=0,
 \qquad |\alpha|\le2.                             \tag{4.3}
\]

Equation (4.3) kills (3.2) exactly under an affine Piola map. It also
kills every term in the first variation (3.5): the amplitude corrections
use at most one label moment and the quadratic phase correction uses at
most two. Hence the selected active Fourier coefficient, including its
Leray pressure part, is zero through first curvature order. Real data are
automatic because (4.1) has real even symbol, so the same assertion holds
at \(-\zeta/r\).

The shell is a **material-label** shell. Under the affine map, a physical
Fourier coordinate \(k\) is killed only when

\[
                  rF^Tk=\zeta,\qquad |\zeta|^2\in{\cal R}. \tag{4.3a}
\]

Thus the notch follows a prescribed Kelvin covector through the common
chart. It does not simultaneously kill an arbitrary fixed Eulerian shell
when \(F(t)\) changes, and tracking a moving broad projection remains part
of MCKC/BAFL.

For the finite \(A_2\) set (3.3), take

\[
 Q_{A_2}(z)=\left(1-{z\over2}\right)^3
             \left(1-{z\over6}\right)^3
 =1-2z+{19\over12}z^2-{17\over27}z^3
   +{19\over144}z^4-{1\over72}z^5+{1\over1728}z^6.
                                                            \tag{4.4}
\]

Its spatial order is twelve and

\[
                 \sum_{m=0}^6|[z^m]Q_{A_2}|={343\over64}.  \tag{4.5}
\]

Consequently, after enlarging the fixed cutoff seminorm,

\[
 \|\nabla^m w^\sharp\|_2
 \le C_{\chi,m+12}\lambda r^{5/2-m},              \tag{4.6}
\]

with a constant independent of \(q,b,K,\lambda\). The notch therefore
does not alter C171's raw \(L^2\) parent-cross ledger.

## 5. The \(L^2\), backward, and \(b\)-power ledgers

Neither pairing nor the fixed notch reduces the full residual norm.
Using (4.6) in C171 gives the same time-integrated \(L^2\) ratio

\[
 C M_Fq^{-5/2}\log h,                              \tag{5.1}
\]

and, after the full backward focus \(h=q^{3/2}\),

\[
 C M_Fq^{-1}\log h=12CM_Fn^{-8}\log n.            \tag{5.2}
\]

Thus the construction is not an \(L^2\)-operator proof of BAFL.

If the retained \(b^2\) wake is realized in this paired collar profile and
the active projection is one of the fixed parent-scale \(A_2\) Fourier
coordinates, the paired moment gain changes its direct coefficient scale
from \(b^2\log h\) to

\[
 b^2q^{-2}\log h.                                  \tag{5.3}
\]

On

\[
 q=n^8,\qquad b=n^{-2},\qquad h=n^{12},            \tag{5.4}
\]

its ratio to the required active allowance \(b^3\) is

\[
 {b^2q^{-2}\log h\over b^3}
 =12n^{-14}\log n\longrightarrow0.                \tag{5.5}
\]

The normalization in (5.3) is load-bearing: \(b^2\) multiplies a profile
whose local core/collar scale is fixed as in Sections 1--2. If instead one
insists that a particular suppressed low-frequency Fourier coefficient of
the paired field itself remain exactly \(b^2\), rescaling the profile by
\(q^2\) erases the claimed gain. Connecting the actual C140 Fourier wake to
the fixed-core paired normalization is therefore an additional realization
obligation.

This is not a proof that the actual C140 retained wake has the paired
profile. The finite-shell notch makes the affine and first-curvature coefficients
exactly zero rather than merely small. With a uniformly bounded
second-curvature Taylor remainder, its first possible targeted term has
the same favorable schematic scale as (5.3), or better. Neither (5.3)
nor that conditional remainder estimate controls a moving broad active
projection followed by the non-normal BAFL propagator.

## 6. Why the broad endpoint defeats a fixed multipole notch

The C144 packet contains the \(q\) axial normalized frequencies

\[
 \zeta_a=\left(4+{a\over q}\right)N,
 \qquad 0\le a<q.                                  \tag{6.1}
\]

Let \(P_q(D)\) be a scalar constant-coefficient differential notch designed
to work independently of the underlying collar Fourier jets. Preserving the
affine slope requires at least \(P_q(0)=1\) (and may impose further
conditions on the linear part). If its symbol and all first and second
Fourier derivatives vanish at every point in (6.1), as in the robust
product-rule kernel (4.3), the univariate restriction

\[
 p_q(t)=P_q(tN)                                    \tag{6.2}
\]

has \(q\) distinct roots, each of multiplicity at least three, while
\(p_q(0)=1\). Therefore

\[
                   \deg P_q\ge3q.                 \tag{6.3}
\]

A radial polynomial in \(|D|^2\) with the same two-jet shell zeros pays
spatial order at least \(6q\).
Since \(q=n^8\), this is far above the existing growing-order Gevrey
window \(m_j\asymp j^2/\log n\), and its cutoff seminorm constant is not
stage uniform. The obstruction already appears on the axial subset; the
other \(q^3-q\) modes can only add conditions. This degree lower bound does
not apply to a notch tailored to accidental zeros of one particular collar
transform, nor to a non-polynomial or time-dependent cancellation.

There is also an analytic boundary. A compactly supported field has an
entire Fourier transform, so forcing the full vector transform to vanish
on an open three-dimensional band forces the field to vanish identically,
contradicting (4.2) when \(S\ne0\). For an invertible constant parent
gradient \(A\), the same conclusion follows from

\[
 k\cdot\widehat w=0,\qquad P_kA\widehat w=0
\]

on an open band: away from the quadratic characteristic set,
\(\widehat w=0\), and analyticity finishes the argument. For a fixed
rank-two matrix, the analogous argument applies only on an open set where
the two-dimensional map
\(P_kA|_{k^\perp}:k^\perp\to k^\perp\) is injective; no uniform assertion
about the full moving \(A_2\) gradient family is made here. This does not
rule out an integrated time-dependent slaving graph; it rules out replacing
that theorem by an exact static open-band collar kernel under the stated
injectivity hypothesis.

## 7. Exact surviving obstruction

C173 gives a bounded, explicit repair for the finite-dimensional part of
MCKC(ii): two material affine curl collars cancel the first nonzero
pressure/active moment while retaining the full local strain, and a
fixed-order triple notch kills the named finite \(A_2\) fibers through the
first curvature jet. The scheduled \(b\)-power then closes on those fixed
material-label coordinates, or on physical coordinates tracked by (4.3a).

What survives is precise:

> **Broad-band material-collar response.** Control the paired collar's
> full Duhamel response on the actual \(q^2\)-to-\(q^3\) moving active
> bundle, with a stage-uniform chart and pressure/heat tails, without a
> differential order growing like \(q\). Equivalently, prove a dynamic
> slaving/export cancellation rather than an exact static Fourier notch.

This is the same MCKC/BAFL obstruction, now with its finite-moment part
removed. The full \(L^2\) residual, self/viscous terms inherited from
C143/C171, material-chart bound, periodic pressure images, and unforced
stage remain open.
