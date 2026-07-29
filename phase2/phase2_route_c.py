#!/usr/bin/env python3
"""Phase-2 Route C: PINN least-squares search for the affine-core continuation (novel experiment N3).

Solves the stationary Type-II Euler similarity profile equation
    (1-gamma) U + gamma (y.grad)U + (U.grad)U + grad P = 0,   div U = 0
with the structured ansatz of PHASE2_SPEC.md:
    U = chi(R) M y + (1-chi) R^{d1} N_Phi(yhat) + g(R) N_c(y),   d1 = 1-1/gamma,
    P = chi P_aff + (1-chi) R^{2 d1} N_Psi(yhat) + h(R) N_p(y),
inner anchor M(lambda, alpha) with alpha^2=(lam+g)(lam+1-g)  [fixed-line spectrum
{0, 2g-1, 1+g} automatic], envelopes g = R/(1+R)^{1/gamma}, h = R^2/(1+R)^{2/gamma}.

Modes:
  --anchor      correctness gate: networks off, chi==1 -> residual must be ~1e-12
  (default)     optimize; prints region-wise residual ledger + JSON summary
Usage (A100): python phase2_route_c.py --gamma 0.45 --lam 0.0 --iters 20000 --pts 8192
"""
import argparse, json, math, time
import torch
from torch.func import jacrev, vmap

ap = argparse.ArgumentParser()
ap.add_argument('--gamma', type=float, default=0.45)
ap.add_argument('--lam', type=float, default=0.0)
ap.add_argument('--iters', type=int, default=3000)
ap.add_argument('--pts', type=int, default=4096)
ap.add_argument('--width', type=int, default=64)
ap.add_argument('--lr', type=float, default=2e-3)
ap.add_argument('--Rb', type=float, default=3.0)
ap.add_argument('--theta', type=float, default=1.0, help='homotopy: 0=pure affine BC, 1=full tail')
ap.add_argument('--seed', type=int, default=0)
ap.add_argument('--anchor', action='store_true')
ap.add_argument('--out', default='phase2_summary.json')
args = ap.parse_args()

torch.manual_seed(args.seed)
dev = 'cuda' if torch.cuda.is_available() else 'cpu'
DT = torch.float64
g, lam = args.gamma, args.lam
alpha = math.sqrt((lam + g) * (lam + 1 - g))
d1 = 1.0 - 1.0 / g
M = torch.tensor([[1., 0., 0.], [0., lam, -alpha], [0., alpha, -1. - lam]], dtype=DT, device=dev)
P2m = M + M @ M
ell = torch.tensor([0., alpha, lam + g], dtype=DT, device=dev)
ell = ell / ell.norm()

class Net(torch.nn.Module):
    def __init__(self, nin, nout, w):
        super().__init__()
        self.f = torch.nn.Sequential(
            torch.nn.Linear(nin, w), torch.nn.Tanh(),
            torch.nn.Linear(w, w), torch.nn.Tanh(),
            torch.nn.Linear(w, nout))
        for m in self.f:
            if isinstance(m, torch.nn.Linear):
                torch.nn.init.xavier_uniform_(m.weight, gain=0.3); torch.nn.init.zeros_(m.bias)
    def forward(self, x): return self.f(x)

W = args.width
NPhi = Net(3, 3, W).to(dev).to(DT)   # tail angular profile
NPsi = Net(3, 1, W).to(dev).to(DT)   # tail pressure profile
Nc   = Net(4, 3, W).to(dev).to(DT)   # correction (rho, yhat)
Np   = Net(4, 1, W).to(dev).to(DT)   # pressure correction
ANCHOR = args.anchor

def UP(y):
    """y: (3,) -> (U (3,), P (1,)) with structured ansatz."""
    R = torch.sqrt((y * y).sum() + 1e-30)
    yh = y / R
    chi = torch.tensor(1.0, dtype=DT, device=dev) if ANCHOR else 0.5 * (1 - torch.tanh((R - args.Rb) / 1.0))
    Uaff = M @ y
    Paff = -0.5 * (y @ (P2m @ y))
    if ANCHOR:
        return Uaff, Paff.reshape(1)
    tailU = args.theta * R**d1 * NPhi(yh)
    tailP = args.theta * R**(2 * d1) * NPsi(yh).squeeze(-1)
    ge = R / (1 + R)**(1.0 / g)
    he = R * R / (1 + R)**(2.0 / g)
    feat = torch.cat([ (R / (1 + R)).reshape(1), yh ])
    U = chi * Uaff + (1 - chi) * tailU + ge * Nc(feat)
    P = chi * Paff + (1 - chi) * tailP + he * Np(feat).squeeze(-1)
    return U, P.reshape(1)

def UP_flat(y):
    U, P = UP(y)
    return torch.cat([U, P])

Jf = jacrev(UP_flat)

def residuals(ys):
    out = vmap(UP_flat)(ys)                     # (B,4)
    J = vmap(Jf)(ys)                            # (B,4,3)
    U = out[:, :3]; dU = J[:, :3, :]; dP = J[:, 3, :]
    adv = torch.einsum('bj,bij->bi', U, dU)
    ydotU = torch.einsum('bj,bij->bi', ys, dU)
    res = (1 - g) * U + g * ydotU + adv + dP
    div = dU[:, 0, 0] + dU[:, 1, 1] + dU[:, 2, 2]
    return res, div

def sample(n):
    rho = torch.rand(n, dtype=DT, device=dev) * 0.97 + 0.005
    R = rho / (1 - rho)
    v = torch.randn(n, 3, dtype=DT, device=dev); v = v / v.norm(dim=1, keepdim=True)
    ys = R[:, None] * v
    nline = n // 8
    s = (torch.rand(nline, dtype=DT, device=dev) * 2 - 1) * 1.5
    off = 0.05 * torch.randn(nline, 3, dtype=DT, device=dev)
    ys2 = s[:, None] * ell[None, :] + off
    return torch.cat([ys, ys2]), R

def wgt(ys):
    R = ys.norm(dim=1)
    scale = R + R.clamp(min=1e-6)**(2 * d1 - 1) + 0.1
    return 1.0 / scale**2, R

if ANCHOR:
    ys, _ = sample(2048)
    res, div = residuals(ys)
    print(json.dumps({"mode": "anchor", "gamma": g, "lam": lam,
                      "max_residual": res.abs().max().item(),
                      "max_div": div.abs().max().item()}, indent=1))
    raise SystemExit

params = list(NPhi.parameters()) + list(NPsi.parameters()) + list(Nc.parameters()) + list(Np.parameters())
opt = torch.optim.Adam(params, lr=args.lr)
t0 = time.time()
for it in range(args.iters):
    opt.zero_grad()
    ys, _ = sample(args.pts)
    res, div = residuals(ys)
    w, R = wgt(ys)
    jetmask = (R < 1.0).to(DT)
    Uaffs = (M[None, :, :] @ ys[:, :, None]).squeeze(-1)
    Uv = vmap(UP_flat)(ys)[:, :3]
    L = (w * (res**2).sum(1)).mean() + (w * div**2).mean() \
        + 0.3 * (jetmask * ((Uv - Uaffs)**2).sum(1)).mean()
    L.backward()
    opt.step()
    if it % max(1, args.iters // 20) == 0:
        print(f"[{it:6d}] loss={L.item():.6e}  t={time.time()-t0:.0f}s")

# region-wise ledger
with torch.no_grad():
    ys, _ = sample(16384)
    res, div = residuals(ys)
    R = ys.norm(dim=1)
    def reg(lo, hi):
        m = (R >= lo) & (R < hi)
        if m.sum() == 0: return None
        return {"rms_res": (res[m]**2).sum(1).mean().sqrt().item(),
                "rms_div": (div[m]**2).mean().sqrt().item(), "n": int(m.sum())}
    summary = {"gamma": g, "lam": lam, "alpha": alpha, "theta": args.theta,
               "iters": args.iters, "final_loss": L.item(),
               "regions": {"inner R<1": reg(0, 1), "blend 1-6": reg(1, 6),
                           "outer 6-30": reg(6, 30), "far >30": reg(30, 1e9)},
               "walltime_s": round(time.time() - t0, 1)}
    print(json.dumps(summary, indent=1))
    open(args.out, 'w').write(json.dumps(summary, indent=1))
