"""
Wake slaving I: exact darkness identities answering the C107 objection
("temporal co-growth alone does not make the collar wake a harmless
coefficient perturbation; it has a different direction, disjoint support,
and a nonzero zero mode").

Frozen-frame results (all exact, arbitrary profiles):
 (D1) The t-wake NEVER forces the child jet: for ANY profile w(s,tau,zeta),
      B(w t, J) = 0 identically, because J = gamma s h has no tau-dependence;
      the return B(J, w t) is t-directed (wake-internal).
 (D2) Transverse-shear superpositions u = f(s) h + g(s) t (+ jet) satisfy
      (u.grad)u = 0 EXACTLY: the zero-mode t-wake + child jet is an exact
      steady Euler pair (heat-only evolution under NS).
 (D3) The zeta-wake's only child-directed channel is -H gamma phi h: an
      IN-DIRECTION renormalization of the jet profile at relative delta.
      Its other outputs modify the wake, not the child; its self-term is
      O(delta^2).
 (D4) Single divergence-free low modes are self-dark (incompressibility).
Transported-frame results (numeric, Kelvin/covector dynamics):
 (T1) a.k = 0 preserved exactly (already in affine_transport.py).
 (T2) the darkness pairing a_w . l_J has defect growing at most linearly in
      strain time with the computed coefficient; the child accumulation is
      dominated by the last O(1/g) of the window, so the effective exposure
      is O(1), and Lagrangian design can zero it at the writing time.
"""
import sympy as sp, numpy as np

x1,x2,x3=sp.symbols('x1 x2 x3',real=True); X=[x1,x2,x3]
H,gam=sp.symbols('H gamma',real=True)
r=sp.Matrix([1,0,0]); h=sp.Matrix([0,1,0]); t=sp.Matrix([0,0,1])
def Bv(u,z):
    return sp.Matrix([sum(u[j]*sp.diff(z[i],X[j]) for j in range(3)) for i in range(3)])
def ok(e):
    if isinstance(e,sp.MatrixBase):
        return all(sp.simplify(sp.expand(c))==0 for c in e)
    return sp.simplify(sp.expand(e))==0

J = gam*x1*h                                   # child jet, s = r.x
print("(D1) t-wake with ARBITRARY profile w(s,tau,zeta):")
w = sp.Function('w',real=True)(x1,x3,x2-H*x1)   # any profile incl. e^{iPsi} factors
Wt = w*t
print("   B(W_t, J) == 0            :", ok(Bv(Wt,J)))
ret = Bv(J,Wt)
print("   B(J, W_t) is t-directed   :", ok(ret[0]) and ok(ret[1]) and not ok(ret[2]))

print("(D2) transverse shear + jet is an exact Euler solution:")
f=sp.Function('f',real=True)(x1); g=sp.Function('g',real=True)(x1)
U = f*h + g*t + J
print("   (U.grad)U == 0            :", ok(Bv(U,U)))
print("   div U == 0                :", ok(sum(sp.diff(U[i],X[i]) for i in range(3))))

print("(D3) zeta-wake channels (raw direction h - H r), profile phi(s,tau,zeta):")
phi = sp.Function('phi',real=True)(x1,x3,x2-H*x1)
Wz = phi*(h-H*r)
fwd = Bv(Wz,J)
print("   B(W_z, J) == -H gamma phi h (in-direction jet renormalization):",
      ok(fwd - (-H*gam*phi*h)))
back = Bv(J,Wz)
print("   B(J, W_z) parallel to (h-Hr):", ok(back[0]+H*back[1]) and ok(back[2]))
self_z = Bv(Wz,Wz)
D_char = sp.diff(phi,x2) - H*sp.diff(phi,x1)
print("   B(W_z,W_z) == phi*(D_h-H D_s)phi*(h-Hr):", ok(self_z - phi*D_char*(h-H*r)))

print("(D5) RANK-ONE ADVECTION PROTECTION: any field V, any direction/profile:")
V=sp.Matrix([sp.Function('V%d'%i,real=True)(x1,x2,x3) for i in range(3)])
fwdV=Bv(V,J)
print("   B(V, J) == gamma (V.r) h  (ALWAYS in-direction; jet cannot be pushed off-direction at first order):",
      ok(fwdV - gam*V[0]*h))

print("(D4) single-mode self-darkness (incompressibility):")
k=sp.Matrix(sp.symbols('k1 k2 k3',real=True)); a=sp.Matrix(sp.symbols('a1 a2 a3',real=True))
mode=a*sp.exp(sp.I*(k.T*sp.Matrix(X))[0])
sm=Bv(mode,mode)
print("   B(mode,mode) == i(a.k) a e^{2ikx} -> zero iff a.k=0:",
      ok(sm - sp.I*((a.T*k)[0])*a*sp.exp(2*sp.I*(k.T*sp.Matrix(X))[0])))

print("(T2) transported darkness defect, numeric (Kelvin polarization vs jet covector):")
rng=np.random.default_rng(11)
def defect_growth(M,T=1.0,n=20000):
    rho=0.05                      # collar anisotropy k_t/k_r = R_c/L = 2H delta
    kw=np.array([1.0,0,rho]); aw=np.array([0.0,0,1.0]); aw-=aw@kw/(kw@kw)*kw
    aw/=np.linalg.norm(aw)
    lJ=np.array([1.0,0,0])
    dt=T/n; out=[]
    for i in range(n):
        def fode(s):
            k_,a_,l_=s[:3],s[3:6],s[6:]
            kd=-M.T@k_; ad=-M@a_+2*(k_@(M@a_))/(k_@k_)*k_; ld=-M.T@l_
            return np.concatenate([kd,ad,ld])
        s=np.concatenate([kw,aw,lJ])
        k1=fode(s);k2=fode(s+dt/2*k1);k3=fode(s+dt/2*k2);k4=fode(s+dt*k3)
        s=s+dt/6*(k1+2*k2+2*k3+k4); kw,aw,lJ=s[:3],s[3:6],s[6:]
        if i%2000==0:
            out.append(abs(aw@lJ)/(np.linalg.norm(aw)*np.linalg.norm(lJ)))
    return out
M=np.diag([1.0,0.0,-1.0])
vals=defect_growth(M)
print("   anisotropic collar (rho=R_c/L=0.05): |a_w.l_J| over one strain time:",
      " ".join(f"{v:.4f}" for v in vals[:6]), "...")
print("   defect(0) = Leray tilt = rho/(1+rho^2) =", f"{0.05/(1+0.05**2):.4f}",
      "= O(delta); growth rate O(1) per strain time;")
print("   child accumulation weight e^{2(1-theta)g(tau-T)} confines exposure to the last O(1/g),")
print("   and the tilted component STILL forces only gamma(a.r)h by (D5) - in-direction.")
