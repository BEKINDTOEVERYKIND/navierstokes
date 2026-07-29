# Round 2b — corrected cells (upload `warmstart_bundle.zip` to the Colab first)

Audit verdict on your four corrected-metric seeds: the honest shape-return frontier now has
its first real number. Best seed (s2): **F = 0.109 at 48³ AND 0.109 at 64³** (resolution-stable,
unlike s1 whose F collapsed cross-grid), at C₁/C₀ = 0.80–0.89, with the gaming detector at
0.126 ≈ the random parity baseline (clean). Consistent scores across all four restarts.
The warm-start failures were my path assumption — the bundle fixes that.

**Cell 0 — setup (after uploading warmstart_bundle.zip):**
```python
import subprocess, glob, json, os, numpy as np
subprocess.run(["unzip","-o","warmstart_bundle.zip"], check=True)
def run(cmd): print(" ".join(cmd)); subprocess.run(cmd, check=True)
def maybe(path): return path if os.path.exists(path) else None
```

**Cell A — the decisive fidelity run: refine s2 directly at 64³ (~2 h):**
```python
run(["python","return_map_opt.py","--M","64","--r0","30","--tau","0.65","--dtau","0.0025",
     "--obj","return","--iters","600","--lr","0.008","--init","ret2_s2.npz","--out","ret_64ref.npz"])
run(["python","return_map_opt.py","--replay","ret_64ref.npz","--M","96","--r0","30",
     "--tau","0.85","--dtau","0.0025","--dtype","c128"])
```
Question: does training *at* 64³ push F past ~0.11, and does the result hold at 96³?

**Cell B — content-vs-shape Pareto front (~2 h):**
```python
for eta in ["0.3","1.0","3.0"]:
    init = maybe("ret_64ref.npz") or "ret2_s2.npz"
    run(["python","return_map_opt.py","--M","48","--r0","30","--obj","return","--eta",eta,
         "--iters","400","--lr","0.02","--init",init,"--out",f"pareto_eta{eta}.npz"])
```

**Cell C — the L³ push (highest single-result value; ~2 h):**
```python
for mc in ["0.0","0.05","0.2"]:
    run(["python","return_map_opt.py","--M","64","--r0","30","--obj","qp","--p","3",
         "--iters","500","--lr","0.01","--match_conc",mc,"--init","qp3.npz",
         "--out",f"qp3_push_mc{mc}.npz"])
```
Pre-registered: Q₃ ≥ 1 with child/parent concentration ratio ≲ 2, stable at 96³ = the clean
novel crossing ("one turnover hands the next octave undiminished endpoint-critical mass").
A parity-penalized plateau at ~0.99 = the endpoint depletion constant. Either is a result.

**Cell D — audit ledgers for everything (send these back + the npz files):**
```python
for f in sorted(set(glob.glob("ret_64ref.npz")+glob.glob("pareto_*.npz")+glob.glob("qp3_push_*.npz"))):
    for Mg in ["64","96"]:
        run(["python","return_map_opt.py","--replay",f,"--M",Mg,"--r0","30",
             "--tau","0.85","--dtau","0.0025","--dtype","c128"])
```
