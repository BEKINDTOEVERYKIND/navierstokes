# Colab cells — quadratic downshifter gate

This run is deliberately different from the earlier direct-DNS optimizations.
It tests one make-or-break local lemma for a scale-separated shell embedding:
can a high-frequency divergence-free cell send its quadratic Euler interaction
mostly into one much lower parent band, with leakage decaying quadratically in
the scale ratio?

Every output is written directly to Google Drive. A recycled runtime therefore
does not erase completed checkpoints.

## Cell 1 — mount Drive and upload the bundle

```python
from google.colab import drive, files
drive.mount("/content/drive")

uploaded = files.upload()  # choose ns_breakthrough_gate.zip
!rm -rf /content/ns_breakthrough_gate
!unzip -oq /content/ns_breakthrough_gate.zip -d /content
%cd /content/ns_breakthrough_gate

import json, pathlib, subprocess, sys, torch
print(torch.cuda.get_device_name(0))
!python -m py_compile ns_downshifter_gate.py parameter_certificate.py
!python parameter_certificate.py

OUT = pathlib.Path("/content/drive/MyDrive/ns_breakthrough_gate_results")
OUT.mkdir(parents=True, exist_ok=True)
print("persistent output:", OUT)
```

## Cell 2 — smoke test

This should finish quickly and verifies the code and persistent write path.

```python
subprocess.run([
    sys.executable, "ns_downshifter_gate.py",
    "--resolution", "48",
    "--carrier", "9",
    "--child-halfwidth", "2",
    "--parent-scale", "2",
    "--iterations", "20",
    "--restarts", "1",
    "--report-every", "5",
    "--checkpoint-every", "5",
    "--output-dir", str(OUT / "smoke_rho4p5"),
], check=True)
```

## Cell 3 — decisive ratio sweep

Run all three. The grids are chosen so the high-frequency quadratic leakage is
represented without aliasing. On an A100, start in single precision.

```python
runs = [
    # resolution, carrier K, parent P, child half-width, iterations
    (128, 20, 5, 4, 350),   # rho = 4
    (160, 25, 5, 4, 350),   # rho = 5
    (192, 32, 4, 3, 450),   # rho = 8
]

for n, carrier, parent, halfwidth, iterations in runs:
    name = f"gate_N{n}_K{carrier}_P{parent}"
    subprocess.run([
        sys.executable, "ns_downshifter_gate.py",
        "--resolution", str(n),
        "--carrier", str(carrier),
        "--child-halfwidth", str(halfwidth),
        "--parent-scale", str(parent),
        "--iterations", str(iterations),
        "--restarts", "3",
        "--learning-rate", "0.015",
        "--min-parent-strength", "0.01",
        "--strength-penalty", "8",
        "--report-every", "10",
        "--checkpoint-every", "10",
        "--output-dir", str(OUT / name),
    ], check=True)
```

## Cell 4 — summarize the scaling decision

The key column is `rho2_scaled`. If leakage is genuinely second order,
`leakage_ratio * rho**2` should remain roughly bounded as `rho` increases,
while `parent_strength` remains comparable and nontrivial.

```python
rows = []
for path in sorted(OUT.glob("gate_*/summary.json")):
    data = json.loads(path.read_text())
    d = data["decision_metrics"]
    rows.append({
        "run": path.parent.name,
        "rho": d["carrier_parent_ratio"],
        "parent_strength": d["parent_strength"],
        "leakage_ratio": d["leakage_ratio"],
        "rho2_scaled": d["order2_scaled_ratio"],
    })

for row in rows:
    print(row)

(OUT / "ratio_sweep_summary.json").write_text(
    json.dumps(rows, indent=2, sort_keys=True) + "\n"
)
```

## Cell 5 — zip results from Drive

```python
archive = "/content/ns_breakthrough_gate_results.zip"
subprocess.run([
    "zip", "-qr", archive, OUT.name
], cwd=str(OUT.parent), check=True)
files.download(archive)
```

## Interpretation

- Strong positive signal: parent strength stays comparable and
  `rho2_scaled` is flat or decreasing from ratio 4 to 8.
- Negative signal: the optimizer can reduce leakage only by collapsing the
  parent strength, or `rho2_scaled` grows clearly with ratio.
- Ambiguous signal: a single restart wins, the ratio-8 run is resolution
  sensitive, or float32 and float64 disagree. In that case validate only the
  best saved candidate in float64; do not launch another broad DNS search.

This experiment does not prove blow-up even if positive. It tests whether the
specific second-order cell lemma left by the analytic audit is numerically
plausible.
