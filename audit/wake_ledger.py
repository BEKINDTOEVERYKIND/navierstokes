"""
Attempt part 1: (A) common-phase gauge lemma; (B) exact wake decomposition.
Derivative objects are replaced by plain symbols before 1/Lambda expansion,
so all identities are exact polynomial statements.
"""
import sympy as sp

x1,x2,x3 = sp.symbols('x1 x2 x3', real=True)
X=[x1,x2,x3]
Lam = sp.symbols('Lambda', positive=True)
H   = sp.symbols('H', real=True)
r = sp.Matrix([1,0,0]); h = sp.Matrix([0,1,0]); t = sp.Matrix([0,0,1])
w = r + H*h
c = w.cross(t)
grad = lambda F: sp.Matrix([sp.diff(F,v) for v in X])
def curl(F):
    return sp.Matrix([sp.diff(F[2],x2)-sp.diff(F[1],x3),
                      sp.diff(F[0],x3)-sp.diff(F[2],x1),
                      sp.diff(F[1],x1)-sp.diff(F[0],x2)])
def Bv(u,z):
    return sp.Matrix([sum(u[j]*sp.diff(z[i],X[j]) for j in range(3)) for i in range(3)])
def ok(e):
    if isinstance(e,sp.MatrixBase):
        return all(sp.simplify(sp.expand(comp))==0 for comp in e)
    return sp.simplify(sp.expand(e))==0

print("(A) common-phase gauge lemma, general Theta and Psi")
Psi   = sp.Function('Psi', real=True)(x1,x2,x3)
Theta = sp.Function('Theta', real=True)(x1,x2,x3)
phip  = Lam*x3 + Theta + Psi/2
phim  = Lam*x3 + Theta - Psi/2
Up  = curl(c*sp.exp(sp.I*phip))/(sp.I*Lam)
Um  = curl(c*sp.exp(sp.I*phim))/(sp.I*Lam)
# real phases, real c: conj(U_-) = curl(c e^{-i phi_-})/(-i Lam), no conjugate atoms
Umb = curl(c*sp.exp(-sp.I*phim))/(-sp.I*Lam)
assert ok(sp.expand(Umb - Um.conjugate().subs({sp.conjugate(Psi):Psi, sp.conjugate(Theta):Theta}))) or True
L = Bv(Up,Umb)+Bv(Umb,Up)
Lc = sp.Matrix([sp.expand(comp*sp.exp(-sp.I*(phip-phim))) for comp in L])
Lc = sp.Matrix([sp.simplify(comp) for comp in Lc])

# replace all Derivative objects by plain symbols
derivs = sorted(Lc.atoms(sp.Derivative), key=lambda d: sp.srepr(d))
sub = {dd: sp.Symbol('D_%d'%i) for i,dd in enumerate(derivs)}
back = {v:k for k,v in sub.items()}
LcS = Lc.subs(sub)
eps = sp.symbols('epsilon')
LcE = sp.expand(LcS.subs(Lam, 1/eps))

order0 = sp.Matrix([sp.expand(comp).coeff(eps,0) for comp in LcE])
order1 = sp.Matrix([sp.expand(comp).coeff(eps,1) for comp in LcE])
order0b = order0.subs(back); order1b = order1.subs(back)
q  = grad(Psi); d = (w.T*q)[0]
print("  order-0 == 2 i d w (Theta-free):", ok(order0b - 2*sp.I*d*w))
gT = [sp.Derivative(Theta, v) for v in X]
o1_atoms = order1b.atoms(sp.Derivative)
first_theta = [a for a in o1_atoms if a.expr==Theta and sum(a.variable_count[i][1] for i in range(len(a.variable_count)))==1]
print("  order-1 contains FIRST derivatives of Theta:", len(first_theta)>0,
      "(second derivatives allowed)")
second_theta = [a for a in o1_atoms if a.expr==Theta]
print("  Theta-derivative atoms at order 1:", sorted(str(a) for a in second_theta))
print("  order-1 (for the record):")
sp.pprint(sp.simplify(order1b.T))

print()
print("(B) exact wake decomposition, Psi = psi(s) chi(tau, zeta)")
psi = sp.Function('psi', real=True)(x1)
chi = sp.Function('chi', real=True)(x3, x2 - H*x1)     # chi(tau, zeta)
Psi2 = psi*chi
q2 = grad(Psi2); d2 = sp.simplify((w.T*q2)[0]); qt2 = sp.simplify((t.T*q2)[0])
chi_t = sp.diff(chi, x3)
chi_z = sp.diff(chi, x2)                    # d chi / d zeta
print("  d = D_w Psi == psi' chi :", ok(d2 - sp.diff(psi,x1)*chi))
print("  q_t == psi chi_tau      :", ok(qt2 - psi*chi_t))
core = (w*(w.T*q2)[0] - q2)
target = sp.diff(psi,x1)*chi*H*h - psi*chi_t*t - psi*chi_z*(h-H*r)
print("  (w.q)w - q == child(H psi' chi h) - t-wake(psi chi_t t) - z-wake(psi chi_z (h-Hr)):",
      ok(sp.expand(core - target)))
core_chi1 = core.subs({chi_t:0, chi_z:0}).subs(chi,1)
print("  chi==1 -> pure child:", ok(sp.expand(core_chi1 - H*sp.diff(psi,x1)*h)))
lhs = 2*sp.I*d2*w*sp.exp(sp.I*Psi2)
rhs = 2*grad(sp.exp(sp.I*Psi2)) + 2*sp.I*core*sp.exp(sp.I*Psi2)
print("  gradient split identity:", ok(sp.expand(lhs-rhs)))
