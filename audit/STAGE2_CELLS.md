# Round-2 Colab cells (corrected metric v2.1 — use the NEW return_map_opt.py)

Stage-1 verdict from your four seeds, in one line: the A100 pipeline works end-to-end
(s0's objective replayed on an independent CPU solver to 3 digits), and the optimizer was
strong enough to find a hole in my fidelity metric — F=0.92 was a sliver artifact (the child
put 0.03% of its content on the exact dyadic bins and aligned only that). The metric is now
sliver-proof: F = (full-band coarse-grained correlation) x (content coverage), in [0,1],
and the ledger prints the parity share so any future gaming is visible instantly.
Old-metric seeds score ~0.0005 under the corrected metric — round 2 asks the real question.

**Cell 1 — setup** (adjust ZIPDIR to wherever you unzipped the reproduction bundle):
```python
import subprocess, glob, json, numpy as np
ZIPDIR = "ns_one_stage_crossing_reproduction"   # contains seeds/seed_r30.npz
def run(cmd): print(" ".join(cmd)); subprocess.run(cmd, check=True)
```

**Cell 2 — stage 1 v2, corrected objective (~1.5–2 h):**
```python
for s in range(4):
    run(["python","return_map_opt.py","--M","48","--r0","30","--tau","0.65","--dtau","0.005",
         "--obj","return","--iters","400","--lr","0.03","--seed",str(s),"--out",f"ret2_s{s}.npz"])
run(["python","return_map_opt.py","--M","48","--r0","30","--obj","return","--iters","400",
     "--lr","0.02","--init","ad574bbc-ret48_s0.npz" if glob.glob("ad574bbc*") else "ret48_s0.npz",
     "--out","ret2_warm_s0.npz"])
run(["python","return_map_opt.py","--M","48","--r0","30","--obj","return","--iters","400",
     "--lr","0.02","--init",f"{ZIPDIR}/seeds/seed_r30.npz","--out","ret2_warm_r30.npz"])
```

**Cell 3 — auto-pick best, refine at 64³ (~1–2 h):**
```python
cands = {f: float(np.load(f)["score"]) for f in glob.glob("ret2_*.npz")}
best = max(cands, key=cands.get); print("scores:", cands, "\nbest:", best)
run(["python","return_map_opt.py","--M","64","--r0","30","--tau","0.65","--dtau","0.0025",
     "--obj","return","--iters","600","--lr","0.01","--init",best,"--out","ret2_64.npz"])
```

**Cell 4 — the content-vs-shape Pareto sweep (novel object; ~2 h):**
```python
for eta in ["0.3","1.0","3.0"]:
    run(["python","return_map_opt.py","--M","48","--r0","30","--obj","return","--eta",eta,
         "--iters","400","--init","ret2_64.npz" if glob.glob("ret2_64.npz") else best,
         "--out",f"pareto_eta{eta}.npz"])
```

**Cell 5 — audit ledgers for everything you'll send back:**
```python
for f in sorted(glob.glob("ret2_*.npz")+glob.glob("pareto_*.npz")):
    run(["python","return_map_opt.py","--replay",f,"--M","64","--r0","30",
         "--tau","0.85","--dtau","0.0025","--dtype","c128"])
```
Send back: all printed JSON + the npz files. I re-audit on the independent solver as before.

Interpretation, pre-registered: corrected-F plateau after honest optimization is the novel
depletion constant for dyadic shape return; ledger lines to watch are `fidelity_shiftopt`
(the real number), `even_sublattice_share_C1` (gaming detector — should sit near ~0.1,
not ~0 or ~1), and the C1/C0-vs-F trade in Cell 4.

## Cell 6 — the L³ push (added after the Q_p audit; highest-value run in the batch)

Your Q_p results audited clean and resolution-robust on the independent engine
(64³/96³, c128): **Q₃ = 0.991–0.992** (prior best 0.76 — one percent from unity at the
endpoint norm), **Q₈ = 1.62–1.63** (a crossing, but via child spikes: concentration 70 vs
parent 12.5 — intermittency lever, not shape return), and **ρ∞ = 1.17–2.14** (pointwise
local-Reynolds ratio above one for the first time; generic value 0.50). Next: push Q₃
over the line *honestly* — warm-start from your qp3 seed with the new concentration-parity
penalty, so a crossing can't come from spike-morphing:

```python
for mc in ["0.0","0.05","0.2"]:
    run(["python","return_map_opt.py","--M","64","--r0","30","--obj","qp","--p","3",
         "--iters","500","--lr","0.01","--match_conc",mc,
         "--init","0c8d306b-qp3.npz" if glob.glob("0c8d306b*") else "qp3.npz",
         "--out",f"qp3_push_mc{mc}.npz"])
for f in sorted(glob.glob("qp3_push_*.npz")):
    run(["python","return_map_opt.py","--replay",f,"--M","96","--r0","30",
         "--tau","0.85","--dtau","0.0025","--dtype","c128"])
```

Pre-registered: Q₃ ≥ 1 with child/parent concentration ratio ≤ ~2 and stability at 96³
would be the cleanest novel crossing of this whole program — "one turnover of true NS
hands the next octave undiminished scale-normalized L³ mass." A plateau at ~0.99 under
the parity penalty is equally reportable: the endpoint-critical depletion constant.
