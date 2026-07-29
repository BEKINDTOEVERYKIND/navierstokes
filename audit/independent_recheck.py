#!/usr/bin/env python3
"""Independent (numpy-only) recomputation of the gate's decision metrics.

1. Loads best_downshifter.npz, recomputes every decision metric from the raw
   child_hat with an independent implementation; compares to their JSON.
2. Decomposes the winning child by helicity sign and radius (detects the
   quasi-Beltrami same-helicity mechanism).
3. NULL BASELINE: evaluates the same metrics for (a) their helical seed and
   (b) nrand random divergence-free band fields, at a given geometry, with NO
   optimization. This is the control the go/no-go needs: rho2_scaled must be
   compared against its unoptimized value, and the null's rho-scaling tells us
   whether "flat in rho" is already true before optimization.
"""
import argparse, json, math
import numpy as np

def make_grid(n):
    k1 = np.fft.fftfreq(n, d=1.0 / n)
    KX, KY, KZ = np.meshgrid(k1, k1, k1, indexing="ij")
    K = np.stack([KX, KY, KZ]); K2 = (K ** 2).sum(0)
    K2z = np.where(K2 > 0, K2, 1.0)
    return K, K2, K2z, np.sqrt(K2)

def leray(fh, K, K2z):
    d = (K * fh).sum(0) / K2z
    out = fh - K * d
    out[:, 0, 0, 0] = 0.0
    return out

def l2_hat(fh, n):
    return math.sqrt(float((np.abs(fh) ** 2).sum()) / n ** 3)

def metrics(child_hat, K, K2z, Kmag, n, carrier, parent_scale,
            plow=0.5, phigh=2.0):
    child = np.fft.ifftn(child_hat, axes=(1, 2, 3), norm="ortho").real
    adv = np.zeros_like(child)
    for ax in range(3):
        d = np.fft.ifftn(1j * K[ax] * child_hat, axes=(1, 2, 3), norm="ortho").real
        adv += child[ax][None] * d
    fh = leray(np.fft.fftn(adv, axes=(1, 2, 3), norm="ortho"), K, K2z)
    pmask = ((Kmag >= plow * parent_scale) & (Kmag <= phigh * parent_scale))
    pn = l2_hat(fh * pmask, n); ln = l2_hat(fh * (~pmask), n)
    ratio = ln / pn if pn > 0 else float("inf")
    return {"parent_norm": pn, "leakage_norm": ln, "leakage_ratio": ratio,
            "parent_strength": pn / parent_scale ** 2.5,
            "order2_scaled_ratio": ratio * (carrier / parent_scale) ** 2,
            "parent_fraction": pn ** 2 / (pn ** 2 + ln ** 2)}

def helicity_split(child_hat, K, Kmag, n):
    """Energy share by helicity sign: project each k onto h_pm(k)."""
    e_plus = e_minus = 0.0
    # build h_pm per k via a reference vector; vectorized
    kx, ky, kz = K
    ref = np.zeros_like(K); ref[2] = 1.0
    para = np.abs(kx * ref[0] + ky * ref[1] + kz * ref[2]) / np.maximum(Kmag, 1e-12)
    ref2 = np.zeros_like(K); ref2[0] = 1.0
    use2 = para > 0.9
    for i in range(3):
        ref[i] = np.where(use2, ref2[i], ref[i])
    e1 = np.cross(ref, K, axis=0)
    e1n = np.sqrt((e1 ** 2).sum(0)); e1 = e1 / np.maximum(e1n, 1e-300)
    e2 = np.cross(K / np.maximum(Kmag, 1e-300), e1, axis=0)
    hp = (e1 + 1j * e2) / math.sqrt(2.0)
    hm = (e1 - 1j * e2) / math.sqrt(2.0)
    cp = (child_hat * np.conj(hp)).sum(0)
    cm = (child_hat * np.conj(hm)).sum(0)
    e_plus = float((np.abs(cp) ** 2).sum()); e_minus = float((np.abs(cm) ** 2).sum())
    tot = e_plus + e_minus
    return e_plus / tot, e_minus / tot

def radius_profile(child_hat, Kmag, lo, hi):
    E = (np.abs(child_hat) ** 2).sum(0)
    out = {}
    for r in range(int(lo), int(hi) + 1):
        m = (Kmag >= r - 0.5) & (Kmag < r + 0.5)
        out[r] = float(E[m].sum())
    tot = sum(out.values())
    return {r: v / tot for r, v in out.items() if v / tot > 1e-4}

def rand_child(n, K, K2z, Kmag, carrier, hw, rng):
    raw = rng.standard_normal((3, n, n, n))
    fh = np.fft.fftn(raw, axes=(1, 2, 3), norm="ortho")
    cmask = (Kmag >= carrier - hw) & (Kmag <= carrier + hw)
    fh = leray(fh * cmask, K, K2z)
    fh /= l2_hat(fh, n)
    return fh

def helical_seed(n, K, K2z, Kmag, carrier, hw, order=1):
    x = 2 * np.pi * np.arange(n) / n
    XX, YY, ZZ = np.meshgrid(x, x, x, indexing="ij")
    def diri(c):
        r = np.ones_like(c)
        for m in range(1, order + 1):
            r += 2 * np.cos(m * c)
        return r / math.sqrt(2 * order + 1)
    env = diri(XX) * diri(YY) * diri(ZZ)
    f = np.zeros((3, n, n, n))
    f[1] = env * np.cos(carrier * XX); f[2] = -env * np.sin(carrier * XX)
    fh = np.fft.fftn(f, axes=(1, 2, 3), norm="ortho")
    cmask = (Kmag >= carrier - hw) & (Kmag <= carrier + hw)
    fh = leray(fh * cmask, K, K2z)
    nn = l2_hat(fh, n)
    return fh / nn if nn > 0 else fh

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--npz", default=None)
    ap.add_argument("--null-geometry", default=None,
                    help="n,carrier,parent,halfwidth e.g. 128,20,5,4")
    ap.add_argument("--nrand", type=int, default=6)
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()

    if a.npz:
        d = np.load(a.npz)
        theirs = json.loads(str(d["metrics_json"]))
        cfg = json.loads(str(d["args_json"]))
        ch = d["child_hat"].astype(np.complex128)
        n = ch.shape[1]; carrier = float(cfg["carrier"])
        P = float(cfg["parent_scale"]); hw = float(cfg["child_halfwidth"])
        K, K2, K2z, Kmag = make_grid(n)
        div = np.abs((K * ch).sum(0)).max()
        mine = metrics(ch, K, K2z, Kmag, n, carrier, P)
        print(json.dumps({"file": a.npz, "n": n, "max_div_residual": float(div),
                          "child_l2": l2_hat(ch, n)}))
        for key in ("parent_norm", "leakage_norm", "leakage_ratio",
                    "parent_strength", "order2_scaled_ratio"):
            print(f"  {key:>22}: theirs={theirs[key]:.6g}  mine={mine[key]:.6g}  "
                  f"rel={abs(mine[key]-theirs[key])/max(abs(theirs[key]),1e-300):.2e}")
        hp, hm = helicity_split(ch, K, Kmag, n)
        print(json.dumps({"helicity_share_plus": round(hp, 4),
                          "helicity_share_minus": round(hm, 4),
                          "radius_profile": radius_profile(ch, Kmag, carrier - hw, carrier + hw)}))

    if a.null_geometry:
        n, carrier, P, hw = [float(v) for v in a.null_geometry.split(",")]
        n = int(n)
        K, K2, K2z, Kmag = make_grid(n)
        rng = np.random.default_rng(a.seed)
        rows = []
        h = helical_seed(n, K, K2z, Kmag, carrier, hw)
        m = metrics(h, K, K2z, Kmag, n, carrier, P)
        rows.append({"init": "helical", **{k: round(v, 6) for k, v in m.items()}})
        for i in range(a.nrand):
            ch = rand_child(n, K, K2z, Kmag, carrier, hw, rng)
            m = metrics(ch, K, K2z, Kmag, n, carrier, P)
            rows.append({"init": f"rand{i}", **{k: round(v, 6) for k, v in m.items()}})
        lr = [r["leakage_ratio"] for r in rows[1:]]
        o2 = [r["order2_scaled_ratio"] for r in rows[1:]]
        print(json.dumps({"geometry": {"n": n, "carrier": carrier, "parent": P,
                                       "halfwidth": hw, "rho": carrier / P},
                          "rows": rows,
                          "rand_median_leakage_ratio": sorted(lr)[len(lr)//2],
                          "rand_median_order2_scaled": sorted(o2)[len(o2)//2]},
                         indent=1))
