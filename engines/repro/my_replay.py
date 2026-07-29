#!/usr/bin/env python3
"""Independent double-precision replay of the shipped seeds on MY solver.

Conventions per their README: amplitude set so u_rms = nu*N0*r0 with
N0 = C_S0/E_S0 from the seed; nu=0.05; tau = t * (nu*N0^2*r0); shells
S0=[4,8), S1=[8,16), S2=[16,kcut]; C_S = sum_S |k||u_hat|^2 (coefficient units).
Solver: my rotational-form RK4 2/3-dealiased pseudospectral code (numpy, complex128)
-- an implementation independent of both of their validators.
"""
import numpy as np, sys, json, time
from scipy import fft as sfft

seedfile = sys.argv[1]
M = int(sys.argv[2]) if len(sys.argv) > 2 else 64
r0 = float(sys.argv[3]) if len(sys.argv) > 3 else 30.0
tau_max = float(sys.argv[4]) if len(sys.argv) > 4 else 1.0
nu, dtau = 0.05, 0.0025

d = np.load(seedfile)
u32 = d['u0_hat'].astype(np.complex128)
Ms = u32.shape[1]

k1 = np.fft.fftfreq(M, d=1.0/M)
K = np.array(np.meshgrid(k1, k1, k1, indexing='ij'))
K2 = (K**2).sum(0); K2z = K2.copy(); K2z[0,0,0] = 1.0
Kmag = np.sqrt(K2)
kcut = M//3
mask = (np.abs(K) <= kcut).all(0)

# ---- embed 32^3 coefficient array into M^3 numpy-fft spectral array ----
uh = np.zeros((3, M, M, M), complex)
ks = np.fft.fftfreq(Ms, d=1.0/Ms).astype(int)
for a in range(Ms):
    for b in range(Ms):
        for c in range(Ms):
            ka, kb, kc = ks[a], ks[b], ks[c]
            uh[:, ka % M, kb % M, kc % M] = u32[:, a, b, c] * M**3

# Hermitian + Leray projection (report sizes)
idx = np.ix_(range(3), (-np.arange(M)) % M, (-np.arange(M)) % M, (-np.arange(M)) % M)
uh_h = 0.5*(uh + np.conj(uh[idx]))
herm = np.abs(uh_h - uh).max()/np.abs(uh).max()
uh = uh_h
dv = (K*uh).sum(0)/K2z
uh = uh - K*dv
print(f"projections: hermitian {herm:.2e}, leray {np.abs(K*dv).max()/np.abs(uh).max():.2e}")

def E_of(uh): return float((np.abs(uh)**2).sum())/M**6
def C_shell(uh, lo, hi):
    sh = (Kmag >= lo) & (Kmag < hi)
    return float((Kmag*(np.abs(uh)**2).sum(0))[sh].sum())/M**6

E0 = E_of(uh); C0 = C_shell(uh, 4, 8)
N0 = C0/E0
target_urms = nu*N0*r0
uh *= target_urms/np.sqrt(E0)
E0 = E_of(uh); C0 = C_shell(uh, 4, 8); N0 = C0/E0
t_nl = 1.0/(nu*N0**2*r0)
dt = dtau*t_nl
print(f"seed N0 = {N0:.6f}, u_rms = {np.sqrt(E0):.6f} (target {target_urms:.6f}), t_nl = {t_nl:.6e}")

def fftv(u):  return np.array([sfft.fftn(c, workers=-1) for c in u])
def ifftv(uh):return np.array([np.real(sfft.ifftn(c, workers=-1)) for c in uh])
def rhs(uh):
    uhm = uh*mask
    u = ifftv(uhm)
    oh = 1j*np.cross(K, uhm, axisa=0, axisb=0, axisc=0)
    w = ifftv(oh)
    ch = fftv(np.cross(u, w, axisa=0, axisb=0, axisc=0))*mask
    dvv = (K*ch).sum(0)/K2z
    return ch - K*dvv - nu*K2*uh

steps = int(round(tau_max/dtau))
hist = []
umax0 = np.abs(ifftv(uh*mask)).max()
t0 = time.time()
for n in range(steps+1):
    C1 = C_shell(uh, 8, 16); C2 = C_shell(uh, 16, kcut + 1e-9)
    hist.append((n*dtau, C_shell(uh, 4, 8), C1, C2))
    if n == steps: break
    a1 = rhs(uh); a2 = rhs(uh+0.5*dt*a1); a3 = rhs(uh+0.5*dt*a2); a4 = rhs(uh+dt*a3)
    uh = uh + (dt/6)*(a1+2*a2+2*a3+a4)

hist = np.array(hist)
i1 = int(np.argmax(hist[:,2]))
print(json.dumps({
  "solver": "independent numpy RK4 rotational (mine)", "M": M, "r0": r0,
  "precision": "complex128", "dtau": dtau, "steps": steps,
  "max_C1_over_C0": hist[i1,2]/C0, "tau_at_peak": hist[i1,0],
  "q1_over_q0": float(np.sqrt(hist[i1,2]/C0)),
  "C2_max_over_C0": float(hist[:,3].max()/C0),
  "walltime_s": round(time.time()-t0,1)}, indent=1))
np.save(f"replay_hist_M{M}_r{r0}.npy", hist)
