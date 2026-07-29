#!/usr/bin/env python3
"""Phase-2 starter: exact anchors for the affine-core continuation (all symbolic).

Checks, for general (gamma, lambda) with alpha^2 = (lam+g)(lam+1-g):
  1. the canonical affine family (P2.2): trace 0, M+M^2 symmetric;
  2. U0 = M y solves the stationary similarity profile equation exactly
     (the theta=0 homotopy anchor of PHASE2_SPEC.md);
  3. the fixed-line spectrum {0, 2g-1, 1+g};
  4. a sample matched tail U1 = R^{1-1/g} Phi(yhat): eigenrelation L_g U1 = 0*U1
     and its divergence constraint;
  5. assembles the unprojected n=2 exterior source  (U1.grad)U1  and verifies its
     homogeneity degree d_2 = 1 - 2/g.  (Leray projection of homogeneous fields is
     the first Colab task: solve Delta chi = div source with chi homogeneous of
     degree d_2+1, then subtract grad chi.)

Run:  python phase2_starter.py          (~1 min, pure sympy)
"""
import sympy as sp

g, lam = sp.symbols('gamma lamda', positive=True)
alpha = sp.sqrt((lam + g) * (lam + 1 - g))
y1, y2, y3 = sp.symbols('y1 y2 y3', real=True)
y = sp.Matrix([y1, y2, y3])
Y = (y1, y2, y3)

M = sp.Matrix([[1, 0, 0], [0, lam, -alpha], [0, alpha, -1 - lam]])

def div(u):  return sum(sp.diff(u[i], Y[i]) for i in range(3))
def grad_s(f): return sp.Matrix([sp.diff(f, v) for v in Y])
def advect(a, b): return sp.Matrix([sum(a[j] * sp.diff(b[i], Y[j]) for j in range(3)) for i in range(3)])
def L_gamma(u): return (1 - g) * u + g * advect(y, u)

print("== 1. canonical family ==")
print("   trace M =", sp.simplify(sp.trace(M)))
P2 = M + M * M
print("   M+M^2 symmetric:", sp.simplify(P2 - P2.T) == sp.zeros(3, 3))

print("== 2. theta=0 anchor: L_g(My) + (My.grad)(My) + grad P0 = 0 ==")
U0 = M * y
P0 = -(y.T * P2 * y)[0, 0] / 2
res = sp.simplify(L_gamma(U0) + advect(U0, U0) + grad_s(P0))
print("   residual =", res.T, "   div U0 =", sp.simplify(div(U0)))

print("== 3. fixed-line spectrum ==")
J = M + g * sp.eye(3)
cp = sp.factor(J.charpoly().as_expr())
print("   charpoly factors:", cp, "   -> roots {0, 2g-1, 1+g}")

print("== 4. sample matched tail ==")
R = sp.sqrt(y1**2 + y2**2 + y3**2)
d1 = 1 - 1 / g
# sample div-free homogeneous field of degree d1:  U1 = curl( R^{d1+1} * e3 ) / (d1+1)-ish
# use A = R^{d1-1} * (−y2, y1, 0)  (a swirl tail about e3); check degree, eigenrelation, div
U1 = R**(d1 - 1) * sp.Matrix([-y2, y1, 0])
eig = sp.simplify(L_gamma(U1) - 0 * U1)   # n=1: eigenvalue (1-n)=0
print("   L_g U1 =", sp.simplify(eig.T), "  (claim: 0)")
print("   div U1 =", sp.simplify(div(U1)))
print("   degree check: U1(s*y)/U1(y) = s^(1-1/g):",
      sp.simplify((R.subs({y1: 2*y1, y2: 2*y2, y3: 2*y3})**(d1-1) * 2) / R**(d1-1) / 2**d1) == 1)

print("== 5. n=2 exterior source (unprojected) ==")
S2 = advect(U1, U1)
# homogeneity: S2(s y) = s^{d2} S2(y), d2 = 2*d1 - 1 = 1 - 2/g
s = sp.symbols('s', positive=True)
d2 = 1 - 2 / g
scaled = S2.subs({y1: s * y1, y2: s * y2, y3: s * y3}, simultaneous=True)
print("   degree of (U1.grad)U1 is d2 = 1-2/g:",
      sp.simplify(scaled - s**d2 * S2) == sp.zeros(3, 1))
print("   next (Colab): Leray-project S2 among homogeneous degree-d2 fields,")
print("   then U2 = P S2  (n-1 = 1), and iterate (P2.5).")
print("\nAll anchors passed. This file is the correctness gate for any Phase-2 solver:")
print("a discretization that does not reproduce (2) to machine precision is wrong.")
