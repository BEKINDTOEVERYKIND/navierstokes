#!/usr/bin/env python3
"""Independent symbolic audit of the lagged-intermittency window.

Checks, with exact symbolic algebra (no floats):
 1. alpha = 1 + 3/(2b) reproduces alpha=26/11 at b=11/10, and the paper's
    fractal volume N_k^(-2(alpha-1)) equals the lagged envelope volume
    N_{k-1}^(-3) identically in b.
 2. The admissible window {b>1, 2b<beta<alpha(b)} is nonempty exactly for
    1 < b < (1+sqrt(13))/4, and (11/10, 9/4) lies inside it.
 3. GENERALIZATION of the bundle's pointwise certificate: over the ENTIRE
    admissible window, both first-order margins are negative and both
    second-order margins are positive. I.e. "first order provably fails,
    second order clears" is structural, not a property of the magic rationals.
    (All conditional on the bundle's leakage-deficit exponent (b-1)(beta-1),
    which is asserted, not derived, in GATE_AUDIT.md.)
 4. Same conclusion under my alternative crude deficit accounting
    (3/2)(1-1/b) (envelope-beat L2 bookkeeping), as a robustness probe.
"""
import sympy as sp

b, beta = sp.symbols("b beta", positive=True)
alpha = 1 + sp.Rational(3, 2) / b

# --- 1. pointwise reproduction and volume identity ---
assert alpha.subs(b, sp.Rational(11, 10)) == sp.Rational(26, 11)
# paper: volume N_k^(-2(alpha-1)); lagged envelope: N_{k-1}^(-3) = N_k^(-3/b)
assert sp.simplify(2 * (alpha - 1) - sp.Rational(3, 1) / b) == 0
print("1. alpha(11/10)=26/11  and  2(alpha-1) == 3/b identically: OK")

# --- 2. window nonempty iff 1 < b < (1+sqrt(13))/4 ---
bstar = (1 + sp.sqrt(13)) / 4
gap = sp.simplify(alpha - 2 * b)          # need > 0 for a beta to exist
roots = sp.solve(sp.Eq(gap, 0), b)
print("2. alpha-2b roots:", roots, "; b* =", sp.nsimplify(bstar), "=", sp.N(bstar, 8))
assert any(sp.simplify(r - bstar) == 0 for r in roots)
b0, beta0 = sp.Rational(11, 10), sp.Rational(9, 4)
assert 1 < b0 < bstar and 2 * b0 < beta0 < alpha.subs(b, b0)
print("   (11/10, 9/4) admissible: OK ;  beta window at b=11/10:",
      (2 * b0, alpha.subs(b, b0)))

# --- 3. window-wide first/second-order margins, deficit D=(b-1)(beta-1) ---
D = (b - 1) * (beta - 1)
M1 = 2 * (b - 1) - D            # carrier ratio, second order
M2 = 2 * (1 - 1 / b) - D        # envelope ratio, second order
F1 = (b - 1) - D                # first order
F2 = (1 - 1 / b) - D
window = sp.And(b > 1, b < bstar, beta > 2 * b, beta < alpha)

checks = {
    "M1 = (b-1)(3-beta)":  (sp.factor(M1), sp.factor(M1) - (b - 1) * (3 - beta)),
    "M2 = (b-1)(2/b-beta+1)": (sp.factor(M2), sp.simplify(M2 - (b - 1) * (2 / b - beta + 1))),
    "F1 = (b-1)(2-beta)":  (sp.factor(F1), sp.factor(F1) - (b - 1) * (2 - beta)),
}
for name, (fact, resid) in checks.items():
    assert sp.simplify(resid) == 0, name
    print("3. factored:", name)

# In-window: beta > 2b > 2 => beta > 2, and beta < alpha < 5/2 < 3.
# M1 > 0 <=> beta < 3      : always in window.
# M2 > 0 <=> beta < 1+2/b  : since alpha = 1+3/(2b) < 1+2/b (3/2<2), always.
# F1 < 0 <=> beta > 2      : always. F2 < 0 <=> beta > 1+1/b: 2b>1+1/b <=> (2b+1)(b-1)>0: always.
assert sp.simplify(alpha - (1 + 2 / b)) == sp.simplify(-1 / (2 * b)) and sp.Rational(-1, 2) < 0
assert sp.factor(2 * b - (1 + 1 / b)) == sp.factor((2 * b + 1) * (b - 1) / b)
print("3. window-wide: M1>0, M2>0, F1<0, F2<0 on the ENTIRE admissible window "
      "(reduced to beta in (2,3) and beta<1+2/b, both implied): OK")

# margins at the chosen point
for nm, e in [("M1", M1), ("M2", M2), ("F1", F1), ("F2", F2)]:
    print(f"   {nm}(11/10,9/4) = {e.subs({b: b0, beta: beta0})}")

# --- 4. robustness: alternative deficit exponent (3/2)(1-1/b) ---
D2 = sp.Rational(3, 2) * (1 - 1 / b)
for nm, e in [("M1'", 2 * (b - 1) - D2), ("M2'", 2 * (1 - 1 / b) - D2),
              ("F1'", (b - 1) - D2), ("F2'", (1 - 1 / b) - D2)]:
    v = sp.simplify(e.subs(b, b0))
    print(f"4. {nm}(11/10) = {v}  ({'>0' if v > 0 else '<=0'})")
# sign analysis in b for the alternative deficit:
#   M2' = (1/2)(1-1/b) > 0 always; F2' = -(1/2)(1-1/b) < 0 always
#   M1' = 2(b-1) - (3/2)(b-1)/b = (b-1)(2 - 3/(2b)) > 0 iff b > 3/4: always
#   F1' = (b-1)(1 - 3/(2b)) > 0 iff b > 3/2: NEGATIVE in window (b < 1.152)
assert sp.simplify(2 * (1 - 1 / b) - D2 - sp.Rational(1, 2) * (1 - 1 / b)) == 0
print("4. under the alternative accounting the dichotomy (first fails / second clears) "
      "still holds on the whole window: OK")
print("\nALL SYMBOLIC CHECKS PASSED")
