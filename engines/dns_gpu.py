#!/usr/bin/env python3
"""GPU (CuPy) port of the forward-cascade module-stage measurement.

Measures, on the true incompressible Navier-Stokes equations on T^3, how much
critical (H^{1/2}) content one octave hands to the next, at what arriving
strength, and how the answer depends on packet strength r0, chirality
composition, and *designed* (aligned-seeded) versus generic data.

Usage (Colab, after `pip install cupy-cuda12x`):
    python dns_gpu.py M R0 NU SEED MODE [ETA] [OUT.json]
      M     grid (128 / 192 / 256; 256^3 complex64 fits a 16 GB T4, use A100 for 384+)
      R0    packet strength u_rms/(nu*N0)   (try 10, 30, 100)
      NU    viscosity (0.05 is fine; r0 is what matters)
      SEED  rng seed
      MODE  pure | mixed | inviscid | aligned
      ETA   (aligned mode) seed amplitude in units of critical, e.g. 0.1; negative flips sign
Examples:
    python dns_gpu.py 256 10 0.05 1 mixed
    python dns_gpu.py 256 30 0.05 1 pure
    python dns_gpu.py 256 10 0.05 1 aligned 0.1
Validation first:  python dns_gpu.py 128 3 0.05 1 inviscid   (expect E,H drift ~ 1e-6 single / 1e-13 double)
"""
import sys, json, time

try:
    import cupy as xp
    GPU = True
except Exception:
    import numpy as xp
    GPU = False

DTYPE = xp.complex64 if GPU else xp.complex128   # switch to complex128 on A100 if desired

M    = int(sys.argv[1]) if len(sys.argv) > 1 else 128
r0   = float(sys.argv[2]) if len(sys.argv) > 2 else 10.0
nu   = float(sys.argv[3]) if len(sys.argv) > 3 else 0.05
seed = int(sys.argv[4]) if len(sys.argv) > 4 else 1
mode = sys.argv[5] if len(sys.argv) > 5 else "mixed"
eta  = float(sys.argv[6]) if len(sys.argv) > 6 else 0.1
out  = sys.argv[7] if len(sys.argv) > 7 else f"gpu_M{M}_r{r0}_{mode}_s{seed}.json"

k1 = xp.fft.fftfreq(M, d=1.0 / M).astype(xp.float32 if GPU else xp.float64)
K = xp.stack(xp.meshgrid(k1, k1, k1, indexing="ij"))
K2 = (K ** 2).sum(0); K2z = K2.copy(); K2z[0, 0, 0] = 1.0
Kmag = xp.sqrt(K2)
kcut = M // 3
mask = (xp.abs(K) <= kcut).all(0)

def fftv(u):   return xp.stack([xp.fft.fftn(c) for c in u]).astype(DTYPE)
def ifftv(uh): return xp.stack([xp.real(xp.fft.ifftn(c)) for c in uh])
def leray(uh):
    d = (K * uh).sum(0) / K2z
    return uh - K * d

def rhs(uh):
    uhm = uh * mask
    u = ifftv(uhm)
    oh = 1j * xp.cross(K, uhm, axisa=0, axisb=0, axisc=0)
    w = ifftv(oh)
    ch = fftv(xp.cross(u, w, axisa=0, axisb=0, axisc=0)) * mask
    return leray(ch) - nu * K2 * uh

# helical basis with h_s(-k) = conj(h_s(k))
ref = xp.asarray([0.618, 0.785, 1.0])[:, None, None, None]
e1 = xp.cross(K, xp.broadcast_to(ref, K.shape), axisa=0, axisb=0, axisc=0)
sgn = xp.sign(K[0] + 1e-9 * xp.sign(K[1] + 1e-9 * xp.sign(K[2]))); sgn[K2 == 0] = 1.0
ref2 = xp.asarray([1.0, 0.0, 0.0])[:, None, None, None]
e1b = xp.cross(K, xp.broadcast_to(ref2, K.shape), axisa=0, axisb=0, axisc=0)
bad = xp.linalg.norm(e1, axis=0) < 1e-12
for c in range(3): e1[c][bad] = e1b[c][bad]
n = xp.linalg.norm(e1, axis=0); n[n < 1e-12] = 1.0
e1 = e1 / n * sgn
khat = K / xp.sqrt(K2z)
e2 = xp.cross(khat, e1, axisa=0, axisb=0, axisc=0)
hp = ((e1 + 1j * e2) / xp.sqrt(2.0)).astype(DTYPE)
hm = ((e1 - 1j * e2) / xp.sqrt(2.0)).astype(DTYPE)

def hel_split(uh):
    return (xp.conj(hp) * uh).sum(0), (xp.conj(hm) * uh).sum(0)

# ------------------ initial data ------------------
N0 = 5.5
rng = xp.random.default_rng(seed) if GPU else __import__("numpy").random.default_rng(seed)
w = xp.asarray(rng.standard_normal(size=(3, M, M, M), dtype=None) if GPU
               else rng.normal(size=(3, M, M, M)))
uh = leray(fftv(w))
ap, am = hel_split(uh)
shell0 = (Kmag >= 4) & (Kmag < 8)
ap = xp.where(shell0, ap, 0)
am = xp.where(shell0, am, 0) if mode in ("mixed", "inviscid", "aligned") else xp.zeros_like(am)
uh = ap * hp + am * hm

E0_target = (nu * N0 * r0) ** 2
def Enorm(uh):
    a, b = hel_split(uh)
    return float(((xp.abs(a) ** 2 + xp.abs(b) ** 2)).sum() / M ** 6)
uh *= (E0_target / Enorm(uh)) ** 0.5

if mode == "aligned":
    # designed data: seed the receiving octave with the transfer direction itself
    nh = rhs(uh) + nu * K2 * uh          # nonlinear part only
    shell1 = (Kmag >= 8) & (Kmag < 16)
    sh = xp.where(shell1, 1.0, 0.0)
    seed_h = leray(nh * sh)
    a, b = hel_split(seed_h)
    Cs = float((Kmag * (xp.abs(a) ** 2 + xp.abs(b) ** 2)).sum() / M ** 6)
    if Cs > 0:
        seed_h *= eta * (nu * 11.0 * 11.0 ** 0.5) / Cs ** 0.5   # eta x critical C at shell 1
        uh = uh + seed_h

t_nl = 1.0 / (nu * N0 ** 2 * r0)
umax = float(xp.abs(ifftv(uh * mask)).max())
dt = min(0.25 * (6.2832 / M) / max(umax, 1e-9), 2.0 / (nu * kcut ** 2 + 1e-30), t_nl / 60)
T_end = 10 * t_nl
if mode == "inviscid":
    nu = 0.0; T_end = 20 * dt

shells = [(4, 8), (8, 16), (16, 32), (32, kcut + 0.001)]
def diags(uh):
    a, b = hel_split(uh)
    E = float((xp.abs(a) ** 2 + xp.abs(b) ** 2).sum() / M ** 6)
    H = float((Kmag * (xp.abs(a) ** 2 - xp.abs(b) ** 2)).sum() / M ** 6)
    o = {"E": E, "H": H}
    for j, (lo, hi) in enumerate(shells):
        sh = (Kmag >= lo) & (Kmag < hi)
        Cp = float((Kmag * xp.abs(a) ** 2)[sh].sum() / M ** 6)
        Cm = float((Kmag * xp.abs(b) ** 2)[sh].sum() / M ** 6)
        o[f"C{j}"] = Cp + Cm; o[f"Cm{j}"] = Cm
    return o

t = 0.0; hist = [(0.0, diags(uh))]; steps = 0; t0 = time.time()
while t < T_end:
    a1 = rhs(uh); a2 = rhs(uh + 0.5 * dt * a1); a3 = rhs(uh + 0.5 * dt * a2); a4 = rhs(uh + dt * a3)
    uh = uh + (dt / 6) * (a1 + 2 * a2 + 2 * a3 + a4)
    t += dt; steps += 1
    if steps % 4 == 0: hist.append((t, diags(uh)))

import numpy as np
d0 = hist[0][1]
ts = np.array([h[0] for h in hist])
C1 = np.array([h[1]["C1"] for h in hist]); Cm1 = np.array([h[1]["Cm1"] for h in hist])
C2 = np.array([h[1]["C2"] for h in hist])
Es = np.array([h[1]["E"] for h in hist]); Hs = np.array([h[1]["H"] for h in hist])
i1 = int(np.argmax(C1))
res = {
    "M": M, "r0": r0, "nu": nu, "seed": seed, "mode": mode, "eta": eta if mode == "aligned" else None,
    "gpu": GPU, "steps": steps, "walltime_s": round(time.time() - t0, 1), "dt": dt, "t_nl": t_nl,
    "f_transfer_S1": float(C1[i1] / d0["C0"]),
    "t_peak_over_tnl": float(ts[i1] / t_nl),
    "r1_over_r0": float((np.sqrt(C1[i1] / 11.0) / (nu * 11.0)) / (np.sqrt(d0["C0"] / N0) / (nu * N0))),
    "minus_share_S1_at_peak": float(Cm1[i1] / max(C1[i1], 1e-300)),
    "f_transfer_S2": float(C2.max() / d0["C0"]),
    "E_drift": float(abs(Es[-1] - Es[0]) / Es[0]),
    "H_drift_rel": float(abs(Hs[-1] - Hs[0]) / max(abs(Hs[0]), d0["C0"])),
}
print(json.dumps(res, indent=1))
open(out, "w").write(json.dumps({"res": res,
    "hist": [[h[0]] + [h[1][k] for k in ("C0", "C1", "C2", "Cm1", "E", "H")] for h in hist]}))
print("saved", out)
