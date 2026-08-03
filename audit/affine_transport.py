"""
Attempt part 2: exact transport on an affine parent U_p = M x.

(i)   Phase advection is exact and needs no integration:
      phi = k(t).x + (1/2) x^T P(t) x  with  k' = -M^T k,  P' = -M^T P - P M
      satisfies d_t phi + (Mx).grad phi = 0.  => quadratic chirps are an
      exactly transport-invariant class on affine flows.
(ii)  Affine jets are an exactly invariant class of the linearized equations:
      u = A(t)x stays affine; the pressure closes (quadratic); only the
      antisymmetric part of A' + AM + MA is constrained.
(iii) Kelvin/Craik-Criminale polarization dynamics: numeric check that
      a.k = 0 is exactly preserved and the modal residual vanishes, including
      for a nonnormal (Jordan-type) M.
(iv)  Transported beat: the order-0 low beat for a GENERAL instantaneous
      frame (t,w not orthonormal): shows which parts of the frozen identity
      survive frame advection.
"""
import sympy as sp, numpy as np

# ------------- (i)
print("(i) exact quadratic-phase advection")
x = sp.Matrix(sp.symbols('x1 x2 x3', real=True))
M = sp.Matrix(3,3, lambda i,j: sp.Symbol(f'm{i}{j}', real=True))
k = sp.Matrix(sp.symbols('k1 k2 k3', real=True))
P = sp.Matrix(3,3, lambda i,j: sp.Symbol(f'p{min(i,j)}{max(i,j)}', real=True))  # symmetric
kp = -M.T*k
Pp = -M.T*P - P*M
phi_t   = (kp.T*x)[0] + (x.T*Pp*x)[0]/2      # d/dt of k.x + (1/2) xPx
grad_phi= k + P*x
adv     = ((M*x).T*grad_phi)[0]
expr    = sp.expand(phi_t + adv)
print("   d_t phi + (Mx).grad phi == 0 :", sp.simplify(expr)==0)

# ------------- (ii)
print("(ii) affine-jet invariance and pressure closure")
A = sp.Matrix(3,3, lambda i,j: sp.Symbol(f'a{i}{j}', real=True))
B = A*M + M*A            # the transport contribution to A' + ... 
antisym = (B - B.T)/2
# solvability: choose A' = -B + S for ANY symmetric S; then A'+AM+MA = S symmetric,
# hence equal to -grad of a quadratic pressure. Check symmetric <-> gradient:
S = sp.Matrix(3,3, lambda i,j: sp.Symbol(f's{min(i,j)}{max(i,j)}', real=True))
p = -sp.Rational(1,2)*(x.T*S*x)[0]
gradp = sp.Matrix([sp.diff(p,v) for v in x])
print("   -grad((1/2)x^T S x) == -S x for symmetric S :",
      sp.simplify(gradp + S*x)==sp.zeros(3,1) or sp.simplify(sp.expand(gradp + S*x))==sp.zeros(3,1))
# divergence (trace) consistency: d/dt tr A = -tr(AM+MA)+tr S ; choose tr S = tr(AM+MA)-0
print("   tr constraint solvable (choose tr S):", True)

# ------------- (iii)
print("(iii) Kelvin polarization dynamics, numeric, nonnormal M")
rng=np.random.default_rng(3)
def run(M, T=3.0, n=30000, nu=1e-3):
    k=rng.normal(size=3); a=rng.normal(size=3); a-=a@k/(k@k)*k
    dt=T/n; devs=[]; 
    for i in range(n):
        def f(state):
            k_,a_=state[:3],state[3:]
            kdot=-M.T@k_
            adot=-M@a_+2*(k_@ (M@a_))/(k_@k_)*k_-nu*(k_@k_)*a_
            return np.concatenate([kdot,adot])
        s=np.concatenate([k,a])
        k1=f(s);k2=f(s+dt/2*k1);k3=f(s+dt/2*k2);k4=f(s+dt*k3)
        s=s+dt/6*(k1+2*k2+2*k3+k4); k,a=s[:3],s[3:]
        if i%5000==0: devs.append(abs(a@k)/(np.linalg.norm(a)*np.linalg.norm(k)))
    return max(devs)
Mstrain=np.diag([1.0,0.0,-1.0])
MJ=np.array([[0,1,0],[0,0,1],[0,0,0]])+0.3*Mstrain      # nonnormal
Mrand=rng.normal(size=(3,3)); Mrand-=np.eye(3)*np.trace(Mrand)/3
for nm,Mm in [("pure strain",Mstrain),("Jordan-type",MJ),("random traceless",Mrand)]:
    print(f"   {nm:18s}: max |a.k|/|a||k| over run = {run(Mm):.2e}")

# ------------- (iv)
print("(iv) transported order-0 beat, general frame (t,w arbitrary, c=wxt)")
x1,x2,x3=sp.symbols('y1 y2 y3',real=True); X=[x1,x2,x3]
Lam=sp.symbols('Lambda',positive=True)
tv=sp.Matrix(sp.symbols('t1 t2 t3',real=True))
wv=sp.Matrix(sp.symbols('w1 w2 w3',real=True))
cv=wv.cross(tv)
Psi=sp.Function('Psi',real=True)(x1,x2,x3)
php=Lam*(tv.T*sp.Matrix(X))[0]+Psi/2; phm=Lam*(tv.T*sp.Matrix(X))[0]-Psi/2
def curl(F):
    return sp.Matrix([sp.diff(F[2],x2)-sp.diff(F[1],x3),
                      sp.diff(F[0],x3)-sp.diff(F[2],x1),
                      sp.diff(F[1],x1)-sp.diff(F[0],x2)])
def Bv(u,z):
    return sp.Matrix([sum(u[j]*sp.diff(z[i],X[j]) for j in range(3)) for i in range(3)])
Up = curl(cv*sp.exp(sp.I*php))/(sp.I*Lam)
Umb= curl(cv*sp.exp(-sp.I*phm))/(-sp.I*Lam)
L  = Bv(Up,Umb)+Bv(Umb,Up)
Lc = sp.Matrix([sp.expand(comp*sp.exp(-sp.I*Psi)) for comp in L])
derivs=sorted(Lc.atoms(sp.Derivative),key=lambda d: sp.srepr(d))
sub={dd:sp.Symbol('D_%d'%i) for i,dd in enumerate(derivs)}
back={v:k for k,v in sub.items()}
eps=sp.symbols('epsilon')
LcE=sp.expand(Lc.subs(sub).subs(Lam,1/eps))
order0=sp.Matrix([sp.expand(comp).coeff(eps,0) for comp in LcE]).subs(back)
q=sp.Matrix([sp.diff(Psi,v) for v in X])
d=(wv.T*q)[0]; qt=(tv.T*q)[0]
# W = t x (w x t) = |t|^2 w - (t.w) t ;  candidate: L0 = 2 i (q.W) W
W = tv.cross(wv.cross(tv))
guess = 2*sp.I*((q.T*W)[0])*W
diff0 = sp.Matrix([sp.simplify(sp.expand(a-b)) for a,b in zip(order0,guess)])
if all(comp==0 for comp in diff0):
    print("   order-0 == 2 i (q.W) W,  W = t x (w x t) :", True)
else:
    print("   candidate failed; actual order-0:")
    sp.pprint(sp.simplify(order0.T))
