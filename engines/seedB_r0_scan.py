#!/usr/bin/env python
"""Reverse r0-scan of the Stage-5 crossing seed (seed B) on the independent engine.

Purpose (pre-registered, Stage 6 item 3): locate the critical strength R1* below
which the one-octave scale-fair endpoint crossing (Q3_oct >= 1) dies.

Engine: the independent complex128 pseudospectral solver used for all prior
cross-audits (rotational form P[u x omega], classical RK4, componentwise 2/3
dealiasing |k_i| <= M/3, Leray projection). Shares no code with the other
model's optimizer or validator. Conventions identical to the Stage-5 replay
that confirmed the crossing at 96^3:
  shells S0=[4,8), S1=[8,16); nu=0.05; u_rms = nu*N0*r0 with N0=C_S0/E_S0;
  tau = t * nu * N0^2 * r0; Q3 = ||u_S1(t)||_L3 / ||u_S0(0)||_L3 (pointwise
  vector magnitude |u|; at p=3 the scale prefactor (N1/N0)^(1-3/p) is 1);
  lambda = (C1/E1)/N0 at the peak; Q3_oct = Q3^(ln2/ln lambda).
Both bookkeepings are reported: time-max Q3 and Q3 at the C1-peak time.

Status: finite-resolution numerical diagnostic; not a proof and not evidence
of a Navier-Stokes singularity.

Usage (Colab A100, after uploading this file and the seed npz):
  python seedB_r0_scan.py --M 96 --dtau 0.0025 \
      --r0list 30,100,320,1000,3200,10000,100000,1000000
  python seedB_r0_scan.py --M 128 --dtau 0.0025 --r0list 1000000   # N-trend point
"""
import argparse, glob, json, math, os, sys
import numpy as np
import torch

ap = argparse.ArgumentParser()
ap.add_argument("--npz", default=None, help="seed npz (key u0_hat); auto-globbed if omitted")
ap.add_argument("--M", type=int, default=96)
ap.add_argument("--dtau", type=float, default=0.0025)
ap.add_argument("--taumax", type=float, default=1.5)
ap.add_argument("--r0list", default="30,100,320,1000,3200,10000,100000,1000000")
ap.add_argument("--nu", type=float, default=0.05)
ap.add_argument("--out", default="seedB_r0_scan.json")
args = ap.parse_args()

path = args.npz
if path is None:
    for pat in ("q3_rate_near_euler_B/optimized_q3_peak_u0_hat.npz",
                "optimized_q3_peak_u0_hat.npz", "*optimized_q3_peak*u0_hat*.npz",
                "*seed*B*.npz"):
        hits = sorted(glob.glob(pat))
        if hits:
            path = hits[0]; break
if path is None or not os.path.exists(path):
    sys.exit("seed npz not found — upload optimized_q3_peak_u0_hat.npz (seed B) "
             "next to this script, or pass --npz PATH")

dev = "cuda" if torch.cuda.is_available() else "cpu"
M, nu = args.M, args.nu
DT, CDT = torch.float64, torch.complex128
print(json.dumps({"engine": "independent-c128", "device": dev, "M": M,
                  "dtau": args.dtau, "taumax": args.taumax, "seed": path}))

k1 = torch.fft.fftfreq(M, d=1.0 / M).to(DT).to(dev)
KX, KY, KZ = torch.meshgrid(k1, k1, k1, indexing="ij")
K = torch.stack([KX, KY, KZ]); K2 = (K ** 2).sum(0)
K2z = K2.clone(); K2z[0, 0, 0] = 1.0
Kmag = torch.sqrt(K2); kcut = M // 3
mask = ((K.abs() <= kcut).all(0)).to(DT)
S0 = (Kmag >= 4) & (Kmag < 8); S1 = (Kmag >= 8) & (Kmag < 16)

def fftv(u):  return torch.fft.fftn(u, dim=(1, 2, 3))
def ifftv(uh): return torch.fft.ifftn(uh, dim=(1, 2, 3)).real
def leray(uh):
    d = (K * uh).sum(0) / K2z
    return uh - K.to(CDT) * d
def cross(a, b):
    return torch.stack([a[1] * b[2] - a[2] * b[1],
                        a[2] * b[0] - a[0] * b[2],
                        a[0] * b[1] - a[1] * b[0]])
def rhs(uh):
    um = uh * mask
    u = ifftv(um); w = ifftv(1j * cross(K.to(CDT), um))
    return leray(fftv(cross(u, w) + 0j) * mask) - nu * K2 * uh
def rk4(x, dt):
    a1 = rhs(x); a2 = rhs(x + 0.5 * dt * a1)
    a3 = rhs(x + 0.5 * dt * a2); a4 = rhs(x + dt * a3)
    return x + (dt / 6) * (a1 + 2 * a2 + 2 * a3 + a4)
def L3mag(p):
    u = ifftv(p)
    return ((torch.sqrt((u ** 2).sum(0)) ** 3).mean()) ** (1.0 / 3.0)

# --- embed seed once (identical to the confirming Stage-5 replay) ---
dd = np.load(path)
key = "u0_hat" if "u0_hat" in dd else list(dd.keys())[0]
u32 = dd[key]; Ms = u32.shape[1]
uh = torch.zeros((3, M, M, M), dtype=CDT, device=dev)
ks = np.fft.fftfreq(Ms, d=1.0 / Ms).astype(int)
src = torch.tensor(u32, dtype=CDT, device=dev) * M ** 3
for a in range(Ms):
    for b in range(Ms):
        for c in range(Ms):
            if abs(ks[a]) < M // 2 and abs(ks[b]) < M // 2 and abs(ks[c]) < M // 2:
                uh[:, ks[a] % M, ks[b] % M, ks[c] % M] = src[:, a, b, c]
idxn = ((-torch.arange(M, device=dev)) % M)
uh = 0.5 * (uh + torch.conj(uh[:, idxn][:, :, idxn][:, :, :, idxn]))
uh = leray(uh * S0.to(DT))
E = (uh.abs() ** 2).sum() / M ** 6
C = (Kmag * (uh.abs() ** 2).sum(0))[S0].sum() / M ** 6
N0 = C / E
uh_unit = uh / torch.sqrt(E)          # unit-energy template; amplitude set per r0

results = []
for r0 in [float(s) for s in args.r0list.split(",")]:
    x = uh_unit * (nu * N0 * r0)
    base = L3mag(x)
    t_nl = 1.0 / (nu * N0.item() ** 2 * r0)
    dt = args.dtau * t_nl
    nsteps = int(round(args.taumax / args.dtau))
    best = (0.0, 0.0, 0.0)             # (Q3, tau, lambda) at time-max
    cbest = (0.0, 0.0)                 # (C1, Q3 at that time)
    with torch.no_grad():
        for n in range(nsteps):
            x = rk4(x, dt); tau = (n + 1) * args.dtau
            p = x * S1.to(DT)
            q3 = (L3mag(p) / base).item()
            E1 = (p.abs() ** 2).sum() / M ** 6
            C1 = ((Kmag * (p.abs() ** 2).sum(0))[S1].sum() / M ** 6).item()
            if q3 > best[0]:
                best = (q3, tau, (C1 / E1 / N0).item())
            if C1 > cbest[0]:
                cbest = (C1, q3)
    q3, tau, lam = best
    q3o = math.exp(math.log(2) * math.log(q3) / math.log(lam)) if q3 > 0 and lam > 1 else float("nan")
    rec = {"r0": r0, "M": M, "Q3_timemax": round(q3, 6), "tau_peak": round(tau, 4),
           "lambda": round(lam, 5), "Q3_per_octave": round(q3o, 6),
           "Q3_at_C1peak": round(cbest[1], 6)}
    results.append(rec)
    print(json.dumps(rec), flush=True)

json.dump({"engine": "independent-c128", "device": dev, "M": M, "dtau": args.dtau,
           "taumax": args.taumax, "seed": path, "results": results},
          open(args.out, "w"), indent=1)
print("wrote", args.out)
