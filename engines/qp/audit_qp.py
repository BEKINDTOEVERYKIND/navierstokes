#!/usr/bin/env python3
"""Independent audit of Q_p-optimized seeds: full ledger + gaming-vector checks."""
import sys, json, math
import numpy as np, torch

path, M = sys.argv[1], int(sys.argv[2])
r0, nu, dtau, tau_max = 30.0, 0.05, 0.005, 0.85
DT, CDT = torch.float64, torch.complex128
k1 = torch.fft.fftfreq(M, d=1.0/M).to(DT)
KX,KY,KZ = torch.meshgrid(k1,k1,k1,indexing='ij')
K = torch.stack([KX,KY,KZ]); K2=(K**2).sum(0); K2z=K2.clone(); K2z[0,0,0]=1.0
Kmag = torch.sqrt(K2); kcut = M//3
mask = ((K.abs()<=kcut).all(0)).to(DT)
S0 = (Kmag>=4)&(Kmag<8); S1 = (Kmag>=8)&(Kmag<16); S2 = (Kmag>=16)&(Kmag<=kcut+1e-9)

def fftv(u): return torch.fft.fftn(u,dim=(1,2,3))
def ifftv(uh): return torch.fft.ifftn(uh,dim=(1,2,3)).real
def leray(uh):
    d=(K*uh).sum(0)/K2z; return uh-K.to(CDT)*d
def cross(a,b):
    return torch.stack([a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]])
def rhs(uh):
    um=uh*mask; u=ifftv(um); w=ifftv(1j*cross(K.to(CDT),um))
    return leray(fftv(cross(u,w)+0j)*mask)-nu*K2*uh
def rk4(x,dt):
    a1=rhs(x);a2=rhs(x+0.5*dt*a1);a3=rhs(x+0.5*dt*a2);a4=rhs(x+dt*a3)
    return x+(dt/6)*(a1+2*a2+2*a3+a4)
def E_of(uh): return (uh.abs()**2).sum()/M**6
def C_sh(uh,sh): return (Kmag*(uh.abs()**2).sum(0))[sh].sum()/M**6

def Qp(uh, sh, pw):
    proj = uh*sh.to(DT)
    u = ifftv(proj)
    E = E_of(proj); C = C_sh(proj,sh); N = (C/(E+1e-300))
    lp = (u.abs()**pw).mean()**(1.0/pw)
    return (lp/N**(1.0-3.0/pw)).item(), N.item(), lp.item()

d = np.load(path); u32 = d['u0_hat']; Ms = 32
uh = torch.zeros((3,M,M,M),dtype=CDT)
ks = np.fft.fftfreq(Ms,d=1.0/Ms).astype(int)
for a in range(Ms):
    for b in range(Ms):
        for c in range(Ms):
            uh[:,ks[a]%M,ks[b]%M,ks[c]%M]=torch.tensor(u32[:,a,b,c],dtype=CDT)*M**3
idxn=((-torch.arange(M))%M)
uh=0.5*(uh+torch.conj(uh[:,idxn][:,:,idxn][:,:,:,idxn]))
uh=leray(uh*S0.to(DT))
E=E_of(uh); C=C_sh(uh,S0); N0=C/E
uh=uh*(nu*N0*r0)/torch.sqrt(E)
C0=C_sh(uh,S0).item()
u0=ifftv(uh); urms0=torch.sqrt((u0**2).mean())
seed_stats={"conc":(u0.abs().max()/urms0).item(),"i4":((u0**4).mean()/urms0**4).item(),"N0":N0.item()}
Q30,_,_ = Qp(uh,S0,3.0); Q80,_,_ = Qp(uh,S0,8.0)

t_nl=1.0/(nu*N0.item()**2*r0); dt=dtau*t_nl
rows=[]
x=uh.clone()
with torch.no_grad():
    for n in range(int(round(tau_max/dtau))):
        x=rk4(x,dt); tau=(n+1)*dtau
        if abs(tau*20-round(tau*20))<1e-9:
            Q31,N1,_=Qp(x,S1,3.0); Q81,N1b,_=Qp(x,S1,8.0)
            u1=ifftv(x*S1.to(DT)); u1rms=torch.sqrt((u1**2).mean())
            rows.append({"tau":round(tau,3),"Q3":Q31/Q30,"Q8":Q81/Q80,
                "C1_over_C0":C_sh(x,S1).item()/C0,"N1":N1,
                "child_conc":(u1.abs().max()/u1rms).item(),
                "rho_inf":((u1.abs().max()/(nu*N1))/(u0.abs().max()/(nu*N0))).item()})
pk3=max(rows,key=lambda r:r["Q3"]); pk8=max(rows,key=lambda r:r["Q8"])
print(json.dumps({"file":path.split('/')[-1],"M":M,"seed":seed_stats,
  "peak_Q3":pk3,"peak_Q8":pk8,
  "C1_peak":max(r["C1_over_C0"] for r in rows)},indent=1))
