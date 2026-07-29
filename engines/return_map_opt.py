#!/usr/bin/env python3
"""Differentiable renormalization-return optimization for 3D Navier-Stokes (novel experiment N1/N2).

Searches over smooth, band-limited, divergence-free seeds for approximate one-step
fixed points of the true NS dyadic renormalization map: evolve one turnover, rescale
the first-octave child back onto the parent lattice, and demand BOTH critical-content
recurrence (C1/C0 >= 1, already achieved by prior seeds) AND shape fidelity
(the current wall: best known dyadic return correlation ~0.20).

Objectives (--obj):
  c1      max C1(tau*)/C0(0)                        [baseline / sanity]
  return  max [C1/C0] * F^eta                       [FLAGSHIP: F = dyadic return fidelity in
                                                     the critical metric, with learnable shift]
  qp      max scale-normalized ||P_S1 u||_p ratio   [--p 3 or 8; anti-cheat penalties]

Modes:
  --smoke               tiny CPU self-test (gradients flow, physics sane)
  --crosscheck SEED     roll a stored seed forward, print C-ledger (validates engine)
  --replay SEED         full diagnostic ledger for any saved seed (also for >=160^3 replays)
  (default)             optimize; saves best seed as npz {u0_hat (3,32,32,32), meta}

Engine: rotational-form RK4 2/3-dealiased pseudospectral, gradient-checkpointed,
float64 on CPU / complex64 or complex128 on GPU (--dtype). Conventions match the
audited reproduction package: u_rms = nu*N0*r0, N0 = C_S0/E_S0, tau = t*nu*N0^2*r0,
shells S0=[4,8), S1=[8,16), S2=[16,kcut], C_S = sum_S |k||u_hat|^2.
"""
import argparse, json, math, time
import numpy as np
import torch
import torch.utils.checkpoint as cp

p = argparse.ArgumentParser()
p.add_argument('--M', type=int, default=48)
p.add_argument('--r0', type=float, default=30.0)
p.add_argument('--nu', type=float, default=0.05)
p.add_argument('--tau', type=float, default=0.65)
p.add_argument('--dtau', type=float, default=0.005)
p.add_argument('--obj', default='return', choices=['c1', 'return', 'qp'])
p.add_argument('--p', type=float, default=3.0)
p.add_argument('--eta', type=float, default=1.0)
p.add_argument('--match_conc', type=float, default=0.0, help='penalty weight for child/parent concentration mismatch (qp objective)')
p.add_argument('--iters', type=int, default=200)
p.add_argument('--lr', type=float, default=0.03)
p.add_argument('--seed', type=int, default=0)
p.add_argument('--restarts', type=int, default=1)
p.add_argument('--init', default=None, help='npz seed to warm-start from')
p.add_argument('--dtype', default='auto', choices=['auto', 'c64', 'c128'])
p.add_argument('--out', default='best_seed.npz')
p.add_argument('--smoke', action='store_true')
p.add_argument('--crosscheck', default=None)
p.add_argument('--replay', default=None)
args = p.parse_args()

dev = 'cuda' if torch.cuda.is_available() else 'cpu'
if args.dtype == 'auto':
    cdt = torch.complex64 if dev == 'cuda' else torch.complex128
else:
    cdt = torch.complex64 if args.dtype == 'c64' else torch.complex128
rdt = torch.float32 if cdt == torch.complex64 else torch.float64
if args.smoke:
    args.M, args.iters, args.tau, args.dtau = 24, 8, 0.05, 0.01

M, nu = args.M, args.nu
kcut = M // 3
k1 = torch.fft.fftfreq(M, d=1.0 / M).to(rdt).to(dev)   # integer frequencies
KX, KY, KZ = torch.meshgrid(k1, k1, k1, indexing='ij')
K = torch.stack([KX, KY, KZ])
K2 = (K ** 2).sum(0); K2z = K2.clone(); K2z[0, 0, 0] = 1.0
Kmag = torch.sqrt(K2)
mask = ((K.abs() <= kcut).all(0)).to(rdt)
S0 = (Kmag >= 4) & (Kmag < 8)
S1 = (Kmag >= 8) & (Kmag < 16)
S2 = (Kmag >= 16) & (Kmag <= kcut + 1e-9)

def fftv(u):  return torch.fft.fftn(u, dim=(1, 2, 3))
def ifftv(uh): return torch.fft.ifftn(uh, dim=(1, 2, 3)).real

def leray(uh):
    d = (K * uh).sum(0) / K2z
    return uh - K.to(cdt) * d

def cross(a, b):
    return torch.stack([a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]])

def rhs(uh):
    um = uh * mask
    u = ifftv(um)
    w = ifftv(1j * cross(K.to(cdt), um))
    ch = fftv(cross(u, w) + 0j) * mask
    return leray(ch) - nu * K2 * uh

def rk4(uh, dt):
    a1 = rhs(uh); a2 = rhs(uh + 0.5*dt*a1); a3 = rhs(uh + 0.5*dt*a2); a4 = rhs(uh + dt*a3)
    return uh + (dt/6)*(a1 + 2*a2 + 2*a3 + a4)

def E_of(uh): return (uh.abs()**2).sum() / M**6
def C_sh(uh, sh): return (Kmag*(uh.abs()**2).sum(0))[sh].sum() / M**6

def make_seed(wraw):
    """real tensor -> band-limited, div-free, amplitude-normalized spectral field."""
    uh = fftv(wraw + 0j)
    uh = uh * S0.to(rdt)
    uh = leray(uh)
    E = E_of(uh); C = C_sh(uh, S0)
    N0 = C / E
    uh = uh * (nu * N0 * args.r0) / torch.sqrt(E)
    return uh, N0

def embed_npz(path):
    d = np.load(path); u32 = d['u0_hat']; Ms = u32.shape[1]
    uh = torch.zeros((3, M, M, M), dtype=torch.complex128)
    ks = np.fft.fftfreq(Ms, d=1.0/Ms).astype(int)
    for a in range(Ms):
        for b in range(Ms):
            for c in range(Ms):
                uh[:, ks[a] % M, ks[b] % M, ks[c] % M] = torch.tensor(u32[:, a, b, c]) * M**3
    idx = [torch.arange(3)[:,None,None,None],
           ((-torch.arange(M)) % M)[None,:,None,None].expand(3,M,M,M),
           ((-torch.arange(M)) % M)[None,None,:,None].expand(3,M,M,M),
           ((-torch.arange(M)) % M)[None,None,None,:].expand(3,M,M,M)]
    uh = 0.5*(uh + torch.conj(uh[idx[0], idx[1], idx[2], idx[3]]))
    uh = uh.to(cdt).to(dev)
    uh = leray(uh * S0.to(rdt))          # enforce band+div exactly
    E = E_of(uh); C = C_sh(uh, S0); N0 = C/E
    return uh * (nu*N0*args.r0)/torch.sqrt(E), N0

# dyadic index map: for k in S0, child bin at 2k
idxs = torch.nonzero(S0, as_tuple=False)
kvec = torch.stack([k1[idxs[:,0]], k1[idxs[:,1]], k1[idxs[:,2]]], 1)
child = ((2*kvec).round().long() % M)
ci = (child[:,0], child[:,1], child[:,2])
pi = (idxs[:,0], idxs[:,1], idxs[:,2])
kmag_par = Kmag[pi]; kmag_ch = Kmag[ci]

# neighborhood offsets for the coarse-grained dyadic pullback (sees ALL parities)
_deltas = torch.tensor([[dx,dy,dz] for dx in (0,1) for dy in (0,1) for dz in (0,1)], device=dev)
_cn = []
for _d in _deltas:
    _cb = ((2*kvec + _d[None,:].to(rdt)).round().long() % M)
    _cn.append((_cb[:,0], _cb[:,1], _cb[:,2]))

def fidelity(u0h, uth, shift):
    """Corrected return fidelity: parent vs neighborhood-averaged child spectrum around 2k,
    normalized by the FULL child-band critical content -> sliver concentration cannot game it."""
    A = (kmag_par.sqrt()[None,:] + 0j) * u0h[:, pi[0], pi[1], pi[2]]
    B = 0
    for _d, cix in zip(_deltas, _cn):
        kk = 2*kvec.to(rdt) + _d[None,:].to(rdt)
        phase = torch.exp(1j * (kk @ shift.to(rdt)))
        B = B + uth[:, cix[0], cix[1], cix[2]] * phase[None,:]
    B = (kmag_ch.sqrt()[None,:] + 0j) * B / 8.0
    acc2 = 0
    for cix in _cn:
        acc2 = acc2 + (uth[:, cix[0], cix[1], cix[2]].abs()**2)
    C1raw = (Kmag*(uth.abs()**2).sum(0))[S1].sum()
    cov = (kmag_ch[None,:]*acc2/8.0).sum() / (C1raw + 1e-300)      # child content seen by the map
    cs = ((A.conj()*B).sum().abs()**2) / ((A.abs()**2).sum()*(B.abs()**2).sum() + 1e-300)
    return cs * cov                                                 # in [0,1]; sliver-proof (M>=48)

def qp_ratio(u0h, uth, pw):
    # v2.3: STANDARD vector-magnitude L^p (pointwise |u|), matching the other model's
    # validator -- our previous componentwise pooling overstated Q3 by ~7% on designed seeds
    def Q(uh, sh):
        proj = uh * sh.to(rdt)
        u = ifftv(proj)
        speed = torch.sqrt((u**2).sum(0) + 1e-300)
        E = E_of(proj); C = C_sh(proj, sh); N = C/(E + 1e-300)
        lp = (speed**pw).mean()**(1.0/pw)
        return lp / N**(1.0 - 3.0/pw)
    return Q(uth, S1) / (Q(u0h, S0) + 1e-300)

def anticheat(u0h):
    u = ifftv(u0h)
    urms = torch.sqrt((u**2).mean())
    conc = u.abs().max()/urms
    i4 = (u**4).mean()/urms**4
    pen = torch.relu(conc - 22.5)**2 + torch.relu(i4 - 137.0)**2
    subs = []
    for lo in (4,5,6,7):
        sb = (Kmag >= lo) & (Kmag < lo+1)
        subs.append((u0h.abs()**2).sum(0)[sb].sum())
    fr = torch.stack(subs); fr = fr/fr.sum()
    pen = pen + torch.relu(fr - 0.55).sum()**2 * 100
    return pen, conc, i4

def rollout(uh, steps, dt, want_hist=False):
    hist = []
    for n in range(steps):
        uh = cp.checkpoint(rk4, uh, torch.tensor(dt), use_reentrant=False) if uh.requires_grad else rk4(uh, dt)
        if want_hist and (n % 4 == 0): hist.append((n+1, C_sh(uh, S0).item(), C_sh(uh, S1).item(), C_sh(uh, S2).item()))
    return uh, hist

def ledger(u0h, N0):
    t_nl = 1.0/(nu*N0.item()*N0.item()*args.r0)
    dt = args.dtau * t_nl
    steps = int(round(args.tau/args.dtau))
    C0 = C_sh(u0h, S0).item()
    uh = u0h.clone(); best = (0.0, -1.0, u0h.clone())
    hist=[]
    with torch.no_grad():
        for n in range(steps):
            uh = rk4(uh, dt)
            c1 = C_sh(uh, S1).item()
            hist.append((n+1)*args.dtau)
            if c1 > best[1]: best = ((n+1)*args.dtau, c1, uh.clone())
    tau_pk, c1pk, upk = best
    u1 = ifftv(upk * S1.to(rdt)); u0 = ifftv(u0h)
    E1 = E_of(upk*S1.to(rdt)).item(); N1 = C_sh(upk,S1).item()/max(E1,1e-300)
    out = {
      "M": M, "r0": args.r0, "C1_over_C0": c1pk/C0, "tau_peak": tau_pk,
      "q1_over_q0": math.sqrt(c1pk/C0), "C2_over_C0": C_sh(upk, S2).item()/C0,
      "L3_ratio": (((u1.abs()**3).mean()**(1/3)) / ((u0.abs()**3).mean()**(1/3))).item(),
      "rho_inf": ((u1.abs().max()/(nu*N1)) / (u0.abs().max()/(nu*N0))).item(),
      "fidelity_id_shift": fidelity(u0h, upk, torch.zeros(3, device=dev)).item(),
      "fidelity_shiftopt": fidelity_best(u0h, upk),
      "even_sublattice_share_C1": (( Kmag*(upk.abs()**2).sum(0))[S1 & EVEN].sum() /
                                   ((Kmag*(upk.abs()**2).sum(0))[S1].sum() + 1e-300)).item(),
    }
    return out

EVEN = ((KX.round().long()%2==0) & (KY.round().long()%2==0) & (KZ.round().long()%2==0))

def fidelity_best(u0h, uth):
    best, barg = -1.0, None
    gr = torch.linspace(0, math.pi, 8, device=dev)[:-1]
    with torch.no_grad():
        for ax in gr:
            for ay in gr:
                for az in gr:
                    f = fidelity(u0h, uth, torch.stack([ax,ay,az])).item()
                    if f > best: best, barg = f, torch.stack([ax,ay,az])
    sh = barg.clone().detach().requires_grad_(True)
    o = torch.optim.Adam([sh], lr=0.03)
    for _ in range(50):
        o.zero_grad(); fv = fidelity(u0h, uth, sh); (-fv).backward(); o.step()
    return fv.item()

# ---------------- modes ----------------
if args.crosscheck or args.replay:
    path = args.crosscheck or args.replay
    uh, N0 = embed_npz(path)
    print(json.dumps({"mode": "replay", "device": dev, "dtype": str(cdt)} | ledger(uh, N0), indent=1))
    raise SystemExit

torch.manual_seed(args.seed)
best_global = None
for rs in range(args.restarts):
    if args.init:
        u_init, _ = embed_npz(args.init)
        wraw = torch.nn.Parameter(ifftv(u_init).to(rdt))
    else:
        wraw = torch.nn.Parameter(torch.randn(3, M, M, M, dtype=rdt, device=dev))
    shift = torch.nn.Parameter(torch.zeros(3, dtype=rdt, device=dev))
    opt = torch.optim.Adam([wraw, shift], lr=args.lr)
    for it in range(args.iters):
        opt.zero_grad()
        u0h, N0 = make_seed(wraw)
        t_nl = 1.0/(nu*N0*N0*args.r0)
        steps = int(round(args.tau/args.dtau))
        uth, _ = rollout(u0h, steps, (args.dtau*t_nl).item() if torch.is_tensor(t_nl) else args.dtau*t_nl)
        C0 = C_sh(u0h, S0); C1 = C_sh(uth, S1)
        pen, conc, i4 = anticheat(u0h)
        if args.obj == 'c1':
            objv = C1/C0
        elif args.obj == 'return':
            F = fidelity(u0h, uth, shift)
            objv = (C1/C0) * F**args.eta
        else:
            objv = qp_ratio(u0h, uth, args.p)
            if args.match_conc > 0:
                u0r = ifftv(u0h); u1r = ifftv(uth * S1.to(rdt))
                c0 = u0r.abs().max()/torch.sqrt((u0r**2).mean())
                c1x = u1r.abs().max()/torch.sqrt((u1r**2).mean())
                objv = objv - args.match_conc * (torch.log(c1x/c0))**2
        loss = -objv + pen
        loss.backward()
        opt.step()
        if it % 10 == 0 or it == args.iters-1:
            F_ = fidelity(u0h, uth, shift).item()
            print(f"[rs{rs} it{it:4d}] obj={objv.item():.6f}  C1/C0={(C1/C0).item():.6f}  F={F_:.4f}  conc={conc.item():.1f}  pen={pen.item():.2e}")
    with torch.no_grad():
        u0h, N0 = make_seed(wraw)
        score = objv.item()
        if best_global is None or score > best_global[0]:
            best_global = (score, u0h.detach().clone(), N0.item())

score, u0h, N0 = best_global
# save at native 32^3 band (all energy in |k|<8): unit-energy coefficient array
u32 = np.zeros((3, 32, 32, 32), complex)
ks32 = np.fft.fftfreq(32, d=1/32).astype(int)
uh_np = u0h.cpu().numpy() / M**3                     # coefficient units
uh_np = uh_np / np.sqrt((np.abs(uh_np)**2).sum())    # unit coefficient energy
for a in range(32):
    for b in range(32):
        for c in range(32):
            ka, kb, kc = ks32[a], ks32[b], ks32[c]
            if abs(ka) < 16 and abs(kb) < 16 and abs(kc) < 16:
                u32[:, a, b, c] = uh_np[:, ka % M, kb % M, kc % M]
np.savez(args.out, u0_hat=u32.astype(np.complex64), r0=args.r0, obj=args.obj, score=score, shift=shift.detach().cpu().numpy())
print(json.dumps({"saved": args.out, "objective": args.obj, "score": score,
                  "final_ledger": ledger(u0h, torch.tensor(N0))}, indent=1))
