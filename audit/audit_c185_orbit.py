"""Independent audit of C185 (+C159 premise spot-check) at 0c5bf3b.
Own implementation from the printed equations only: phase orbit, gamma
quadrature/shooting for beta, covector reconstruction, cooperative frame,
B(s) entries, the four certificate inequalities, and a fully independent
3D Kelvin monodromy check of M w >= e^{1/5} w, det M = 1, tr M > 2."""
import numpy as np

d = 4/5
N = np.array([1.,1.,1.]); n = N/np.sqrt(3)
r1 = np.array([1.,-1.,0.]); r2 = np.array([0.,1.,-1.])

def fab(a,b): return np.cos(a)+np.cos(b)+d*np.cos(a+b)
def phase_rhs(a,b):
    s = np.sin(a+b)
    return 3*(np.sin(b)+d*s), -3*(np.sin(a)+d*s)
def g_vec(a,b):
    s = np.sin(a+b)
    fa = -(np.sin(a)+d*s); fb = -(np.sin(b)+d*s)
    return fa*r1 + fb*r2
def H_mat(a,b):
    ca, cb, cab = np.cos(a), np.cos(b), np.cos(a+b)
    faa = -(ca + d*cab); fbb = -(cb + d*cab); fab_ = -d*cab
    return faa*np.outer(r1,r1) + fab_*(np.outer(r1,r2)+np.outer(r2,r1)) + fbb*np.outer(r2,r2)
Ncross = np.array([[0,-1,1],[1,0,-1],[-1,1,0]], float)   # [N]_x
def U_vec(a,b): return Ncross@g_vec(a,b) - np.sqrt(2)*fab(a,b)*N
def A_mat(a,b): return Ncross@H_mat(a,b) - np.sqrt(2)*np.outer(N, g_vec(a,b))

# basepoint
a0 = np.arctan2(np.sqrt(21)/5, -2/5); b0 = -a0
c_inv = 3*(2*( ( (np.sin(a0)+d*np.sin(a0+b0))**2 + (np.sin(b0)+d*np.sin(a0+b0))**2
                + (np.sin(a0)+d*np.sin(a0+b0))*(np.sin(b0)+d*np.sin(a0+b0))*(-1)*(-1) ) ))
h0 = np.dot(g_vec(a0,b0), g_vec(a0,b0))
c = 3*h0
print(f"h(X0) = {h0:.10f} (claim 126/25 = {126/25}), c = {c:.10f} (claim 378/25 = {378/25})")
print(f"f(X0) = {fab(a0,b0):.2e};  |U(X0)|^2 = {np.dot(U_vec(a0,b0),U_vec(a0,b0)):.10f}")

# integrate orbit, find period (return to (a0,b0))
def rk4(y, dt, rhs):
    k1 = rhs(y); k2 = rhs(y+dt/2*k1); k3 = rhs(y+dt/2*k2); k4 = rhs(y+dt*k3)
    return y + dt/6*(k1+2*k2+2*k3+k4)
def orbit_rhs(y):
    da, db = phase_rhs(y[0], y[1]); return np.array([da, db])
dt = 2e-6
y = np.array([a0,b0]); t = 0.0; traj=[(0.0,a0,b0)]
prev = None
T = None
for i in range(6_000_000):
    y = rk4(y, dt, orbit_rhs); t += dt
    if i % 25 == 0: traj.append((t, y[0], y[1]))
    if t > 0.05:
        dd = (y[0]-a0)**2 + (y[1]-b0)**2
        if prev is not None and dd > prev and prev < 1e-8:
            T = t - dt; break
        prev = dd
print(f"period T = {T:.8f}   (f-drift {abs(fab(y[0],y[1])):.2e})")

# beta via exact linear shooting on gamma' = 3 sqrt2 beta + (2c/3) U^T H U / h^2
def UHU_h2(a,b):
    U = U_vec(a,b); h = np.dot(g_vec(a,b), g_vec(a,b))
    return (U@H_mat(a,b)@U)/h**2
# integral of UHU/h^2 over period via the stored fine trajectory (re-integrate finely)
y = np.array([a0,b0]); I1 = 0.0; nstep = int(T/dt)
for i in range(nstep):
    v1 = UHU_h2(y[0],y[1]); ymid = rk4(y, dt/2, orbit_rhs); v2 = UHU_h2(ymid[0],ymid[1])
    y2 = rk4(y, dt, orbit_rhs); v3 = UHU_h2(y2[0],y2[1])
    I1 += dt*(v1+4*v2+v3)/6; y = y2
beta = -(2*c/3)*I1/(3*np.sqrt(2)*T)
m = np.sqrt(3)*beta
print(f"beta = {beta:.8f},  m = k.n = {m:.8f}")

# covector k(t) and frame; verify k' = -A^T k and invariants; build B and Kelvin monodromy
def k_of(a,b,gam):
    h = np.dot(g_vec(a,b), g_vec(a,b))
    return (c/(3*h))*U_vec(a,b) + beta*N + gam*g_vec(a,b)
def full_rhs(z):
    a,b,gam = z[0], z[1], z[2]
    da, db = phase_rhs(a,b)
    dgam = 3*np.sqrt(2)*beta + (2*c/3)*UHU_h2(a,b)
    return np.array([da, db, dgam])
z = np.array([a0,b0,0.0])
kv0 = k_of(*z)
# combined 12-dim system: (a,b,gamma, V flattened) with one RK4, dt=2e-5
def big_rhs(Y):
    a,b,gam = Y[0],Y[1],Y[2]; Vm = Y[3:].reshape(3,3)
    da, db = phase_rhs(a,b)
    dgam = 3*np.sqrt(2)*beta + (2*c/3)*UHU_h2(a,b)
    kv = k_of(a,b,gam); q = kv@kv; A = A_mat(a,b)
    AV = A@Vm
    dV = -AV + 2*np.outer(kv, (kv@AV))/q
    return np.concatenate([[da,db,dgam], dV.ravel()])
DT = 2e-5; nstep = int(T/DT)
Y = np.concatenate([z,[0]*0, np.eye(3).ravel()])
Bviol = {"B12":1e9,"B21":1e9,"row1":1e9,"row2":1e9}; maxinv=0.0; samp=0
Pn = np.eye(3) - np.outer(n,n)
for i in range(nstep):
    a,b,gam = Y[0],Y[1],Y[2]
    if i % 200 == 0:
        kv = k_of(a,b,gam); q = kv@kv; A = A_mat(a,b)
        p = kv - m*n; D = p@p
        lv = A.T@n
        Sw = lambda w: Pn@(A@w)
        E2v = np.cross(kv, n)
        B11 = T*(m*(lv@p))/D
        B22 = T*(2*(p@Sw(p)) + m*(lv@p))/D
        B21 = T*(m*m*np.sqrt(2)*c)/(q*D)
        B12 = T*(2*m*(p@Sw(E2v)) + np.sqrt(2)*c*(D-m*m))/D
        Bviol["B12"]=min(Bviol["B12"],B12); Bviol["B21"]=min(Bviol["B21"],B21)
        Bviol["row1"]=min(Bviol["row1"],B11+0.15*B12)
        Bviol["row2"]=min(Bviol["row2"],(20/3)*B21+B22)
        samp+=1
        # covector ODE residual (finite difference)
        eps=1e-7
        z3=np.array([a,b,gam]); z3e=rk4(z3,eps,full_rhs)
        resid=(k_of(*z3e)-kv)/eps + A.T@kv
        maxinv=max(maxinv,np.linalg.norm(resid)/np.linalg.norm(A.T@kv))
    k1=big_rhs(Y); k2=big_rhs(Y+DT/2*k1); k3=big_rhs(Y+DT/2*k2); k4=big_rhs(Y+DT*k3)
    Y=Y+DT/6*(k1+2*k2+2*k3+k4)
z=Y[:3]; V=Y[3:].reshape(3,3)
print(f"gamma(T) = {z[2]:.2e};  covector-ODE max rel residual = {maxinv:.2e}")
print(f"k(T)-k(0) rel = {np.linalg.norm(k_of(*z)-kv0)/np.linalg.norm(kv0):.2e}")
print(f"min over orbit ({samp} samples): B12={Bviol['B12']:.3f} (>32?)  B21={Bviol['B21']:.4f} (>9/10?)")
print(f"  B11+(3/20)B12={Bviol['row1']:.4f} (>7/10?)  (20/3)B21+B22={Bviol['row2']:.4f} (>1/5?)")
kv=kv0; q=kv@kv
E1 = n-(m/q)*kv; E2=np.cross(kv,n)
E1n,E2n = E1/np.linalg.norm(E1), E2/np.linalg.norm(E2)
M = np.array([[E1n@(V@E1n), E1n@(V@E2n)],[E2n@(V@E1n), E2n@(V@E2n)]])
w=np.array([1.,0.15]); Mw=M@w
print(f"M = {M.round(6).tolist()}  det={np.linalg.det(M):.6f} tr={np.trace(M):.6f}")
print(f"Mw = {Mw.round(5).tolist()} vs e^0.2 w = {(np.e**0.2*w).round(5).tolist()};  Mw>=e^0.2 w:",bool(np.all(Mw>=np.e**0.2*w)))
print(f"rho(M) = {max(abs(np.linalg.eigvals(M))):.6f} > e^0.2 = {np.e**0.2:.6f}")
