import numpy as np, sys, json, time
from scipy import fft as sfft

M   = int(sys.argv[1]) if len(sys.argv)>1 else 64
r0  = float(sys.argv[2]) if len(sys.argv)>2 else 10.0
nu  = float(sys.argv[3]) if len(sys.argv)>3 else 0.05
seed= int(sys.argv[4]) if len(sys.argv)>4 else 1
mode = sys.argv[5] if len(sys.argv)>5 else "pure"
inviscid_test = (mode=="inviscid")

k1 = np.fft.fftfreq(M, d=1.0/M)
K  = np.array(np.meshgrid(k1,k1,k1,indexing='ij'), dtype=np.float64)
K2 = (K**2).sum(0); K2z = K2.copy(); K2z[0,0,0]=1.0
Kmag = np.sqrt(K2)
kmax23 = M//3
mask = (np.abs(K)<=kmax23).all(0)

def fftv(u):  return np.array([sfft.fftn(c, workers=-1) for c in u])
def ifftv(uh):return np.array([np.real(sfft.ifftn(c, workers=-1)) for c in uh])
def leray(uh):
    d=(K*uh).sum(0)/K2z; return uh-K*d

def rhs(uh):
    uhm = uh*mask
    u = ifftv(uhm)
    oh = 1j*np.cross(K, uhm, axisa=0, axisb=0, axisc=0)
    w = ifftv(oh)
    ch = fftv(np.cross(u, w, axisa=0, axisb=0, axisc=0))*mask
    return leray(ch) - nu*K2*uh

# helical basis with h_s(-k) = conj(h_s(k))
ref = np.array([0.618,0.785,1.0])[:,None,None,None]
e1 = np.cross(K, np.broadcast_to(ref,K.shape), axisa=0,axisb=0,axisc=0)
sgn = np.sign(K[0] + 1e-9*np.sign(K[1] + 1e-9*np.sign(K[2]))); sgn[K2==0]=1.0
ref2 = np.array([1.0,0.0,0.0])[:,None,None,None]
e1b = np.cross(K, np.broadcast_to(ref2,K.shape), axisa=0,axisb=0,axisc=0)
bad = np.linalg.norm(e1,axis=0)<1e-12
for c in range(3): e1[c][bad]=e1b[c][bad]
n = np.linalg.norm(e1,axis=0); n[n<1e-12]=1.0
e1 = e1/n*sgn
khat = K/np.sqrt(K2z)
e2 = np.cross(khat, e1, axisa=0,axisb=0,axisc=0)
hp = (e1+1j*e2)/np.sqrt(2); hm = (e1-1j*e2)/np.sqrt(2)

def hel_split(uh):
    ap = (np.conj(hp)*uh).sum(0); am = (np.conj(hm)*uh).sum(0)
    return ap, am

# ---- initial data: pure + chirality, shell S0: 4<=|k|<8 ----
rng = np.random.default_rng(seed)
w = rng.normal(size=(3,M,M,M))
uh = leray(fftv(w))
ap, am = hel_split(uh)
shell0 = (Kmag>=4)&(Kmag<8)
ap = np.where(shell0, ap, 0)
if mode == "mixed":
    am = np.where(shell0, am, 0)
else:
    am = np.zeros_like(am)
uh = ap*hp + am*hm
# reality check
uu = np.array([sfft.ifftn(c, workers=-1) for c in uh])
assert np.abs(uu.imag).max() < 1e-10*max(1e-30,np.abs(uu.real).max()), "field not real"
N0 = 5.5
E0_target = (nu*N0*r0)**2      # space-filling shell: critical strength u_rms = nu*N0*r0
def Enorm(uh):
    ap, am = hel_split(uh)
    return ((np.abs(ap)**2+np.abs(am)**2)).sum()/M**6
uh *= np.sqrt(E0_target/Enorm(uh))

t_nl = 1.0/(nu*N0**2*r0)
umax = np.abs(ifftv(uh*mask)).max()
dt = min(0.25*(2*np.pi/M)/max(umax,1e-9), 2.0/(nu*(kmax23**2)+1e-30), t_nl/60)
T_end = 10*t_nl
if inviscid_test: nu=0.0; T_end = 20*dt

shells = [(4,8),(8,16),(16,kmax23+0.001)]
def diags(uh):
    ap, am = hel_split(uh)
    E  = ((np.abs(ap)**2+np.abs(am)**2)).sum()/M**6
    H  = ((Kmag)*(np.abs(ap)**2-np.abs(am)**2)).sum()/M**6
    out = {"E":E, "H":H}
    for j,(a,b) in enumerate(shells):
        sh = (Kmag>=a)&(Kmag<b)
        Cp = (Kmag*np.abs(ap)**2)[sh].sum()/M**6
        Cm = (Kmag*np.abs(am)**2)[sh].sum()/M**6
        out[f"C{j}"]=Cp+Cm; out[f"Cm{j}"]=Cm
    return out

t=0.0; hist=[]
d0=diags(uh); hist.append((t,d0))
steps=0; t0=time.time()
while t < T_end:
    k1_=rhs(uh); k2_=rhs(uh+0.5*dt*k1_); k3_=rhs(uh+0.5*dt*k2_); k4_=rhs(uh+dt*k3_)
    uh = uh + (dt/6)*(k1_+2*k2_+2*k3_+k4_)
    t += dt; steps+=1
    if steps%4==0: hist.append((t,diags(uh)))
d0 = hist[0][1]
res = {"M":M,"r0":r0,"nu":nu,"seed":seed,"dt":dt,"t_nl":t_nl,"steps":steps,
       "walltime":time.time()-t0,
       "E0":d0["E"],"H0":d0["H"],"C0_0":d0["C0"]}
# extract module metrics
ts   = np.array([h[0] for h in hist])
C1   = np.array([h[1]["C1"] for h in hist]); Cm1 = np.array([h[1]["Cm1"] for h in hist])
C2   = np.array([h[1]["C2"] for h in hist])
Ctot = np.array([h[1]["C0"]+h[1]["C1"]+h[1]["C2"] for h in hist])
Es   = np.array([h[1]["E"] for h in hist]); Hs = np.array([h[1]["H"] for h in hist])
i1 = int(np.argmax(C1))
res.update({
  "f_transfer_S1": float(C1[i1]/d0["C0"]),
  "t_peak_over_tnl": float(ts[i1]/t_nl),
  "r1_over_r0": float((np.sqrt(C1[i1]/11.0)/(nu*11.0)) / (np.sqrt(d0["C0"]/N0)/(nu*N0))),
  "minus_share_S1_at_peak": float(Cm1[i1]/max(C1[i1],1e-300)),
  "C1_over_Ctot_at_peak": float(C1[i1]/max(Ctot[i1],1e-300)),
  "f_transfer_S2": float(C2.max()/d0["C0"]),
  "Ctot_max_over_C0": float(Ctot.max()/Ctot[0]),
  "E_drift": float(abs(Es[-1]-Es[0])/Es[0]),
  "H_drift_rel_C": float(abs(Hs[-1]-Hs[0])/Ctot[0]),
})
print(json.dumps(res, indent=1))
np.save(f"hist_M{M}_r{r0}_s{seed}_{mode}.npy", np.array([(h[0],h[1]["C0"],h[1]["C1"],h[1]["C2"],h[1]["Cm1"],h[1]["E"],h[1]["H"]) for h in hist]))
