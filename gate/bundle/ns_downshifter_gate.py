#!/usr/bin/env python3
"""Search for a scale-separated Navier--Stokes quadratic downshifter.

This is a deliberately narrow go/no-go experiment for the remaining
Palasek-embedding route.  It searches over real, divergence-free vector fields
W whose Fourier support lies in a band around a high carrier K.  For

    F(W) = P div(W tensor W),

the desired part is the projection of F to a prescribed low parent band near
P, while every other Fourier component is counted as leakage.

The optimization target is therefore

    ||F_leak||_2 / ||F_parent||_2,

subject to a lower bound on the dimensionless parent interaction strength
||F_parent||_2 / P^(5/2).  A successful second-order WKB/normal-form cell
should exhibit leakage proportional to (P/K)^2 (or better) at fixed
dimensionless parent strength.

This numerical search is not a proof.  Its purpose is to decide whether it is
worth attempting the much harder analytic cell lemma.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np
import torch


@dataclass
class Metrics:
    iteration: int
    restart: int
    loss: float
    parent_norm: float
    leakage_norm: float
    total_norm: float
    leakage_ratio: float
    parent_fraction: float
    parent_strength: float
    order2_scaled_ratio: float
    child_l4: float
    child_linf: float
    parent_linf: float
    parent_mean_wavenumber: float
    wall_seconds: float


class SpectralGrid:
    def __init__(self, n: int, device: torch.device, dtype: torch.dtype):
        self.n = n
        self.device = device
        self.dtype = dtype
        self.cdtype = torch.complex128 if dtype == torch.float64 else torch.complex64
        modes = torch.fft.fftfreq(n, d=1.0 / n, device=device, dtype=dtype)
        self.kx, self.ky, self.kz = torch.meshgrid(
            modes, modes, modes, indexing="ij"
        )
        self.k = torch.stack((self.kx, self.ky, self.kz), dim=0)
        self.k2 = self.kx.square() + self.ky.square() + self.kz.square()
        self.kabs = torch.sqrt(self.k2)
        self.inv_k2 = torch.where(
            self.k2 > 0, self.k2.reciprocal(), torch.zeros_like(self.k2)
        )
        self.zero_mask = self.k2 == 0
        self.points = float(n**3)

    def fft(self, field: torch.Tensor) -> torch.Tensor:
        return torch.fft.fftn(field, dim=(-3, -2, -1), norm="ortho")

    def ifft(self, field_hat: torch.Tensor) -> torch.Tensor:
        return torch.fft.ifftn(field_hat, dim=(-3, -2, -1), norm="ortho").real

    def leray(self, field_hat: torch.Tensor) -> torch.Tensor:
        dot = (field_hat * self.k).sum(dim=0)
        projected = field_hat - self.k * (dot * self.inv_k2)
        return torch.where(self.zero_mask.unsqueeze(0), 0.0, projected)

    def l2_hat(self, field_hat: torch.Tensor) -> torch.Tensor:
        return torch.sqrt(
            field_hat.abs().square().sum() / self.points + torch.finfo(self.dtype).tiny
        )

    def normalize_field(
        self, raw: torch.Tensor, mask: torch.Tensor
    ) -> tuple[torch.Tensor, torch.Tensor]:
        field_hat = self.leray(self.fft(raw) * mask)
        norm = self.l2_hat(field_hat)
        field_hat = field_hat / norm
        return self.ifft(field_hat), field_hat

    def quadratic_force(
        self, field: torch.Tensor, field_hat: torch.Tensor
    ) -> torch.Tensor:
        adv = torch.zeros_like(field)
        for axis, kval in enumerate((self.kx, self.ky, self.kz)):
            derivative = self.ifft((1j * kval).unsqueeze(0) * field_hat)
            adv = adv + field[axis].unsqueeze(0) * derivative
        return self.leray(self.fft(adv))


def radial_mask(grid: SpectralGrid, lower: float, upper: float) -> torch.Tensor:
    return ((grid.kabs >= lower) & (grid.kabs <= upper)).to(grid.dtype)


def helical_initialization(
    grid: SpectralGrid, carrier: int, envelope_order: int
) -> torch.Tensor:
    """A real intermittent helical wave, subsequently projected to the search band."""
    x = (
        2.0
        * math.pi
        * torch.arange(grid.n, device=grid.device, dtype=grid.dtype)
        / grid.n
    )
    xx, yy, zz = torch.meshgrid(x, x, x, indexing="ij")

    def dirichlet(coord: torch.Tensor) -> torch.Tensor:
        result = torch.ones_like(coord)
        for mode in range(1, envelope_order + 1):
            result = result + 2.0 * torch.cos(mode * coord)
        return result / math.sqrt(2 * envelope_order + 1)

    envelope = dirichlet(xx) * dirichlet(yy) * dirichlet(zz)
    phase = carrier * xx
    field = torch.zeros(
        (3, grid.n, grid.n, grid.n), device=grid.device, dtype=grid.dtype
    )
    field[1] = envelope * torch.cos(phase)
    field[2] = -envelope * torch.sin(phase)
    return field


def random_initialization(
    grid: SpectralGrid, generator: torch.Generator
) -> torch.Tensor:
    return torch.randn(
        (3, grid.n, grid.n, grid.n),
        generator=generator,
        device=grid.device,
        dtype=grid.dtype,
    )


def scalar(value: torch.Tensor) -> float:
    return float(value.detach().cpu())


def jsonable_args(args: argparse.Namespace) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in vars(args).items():
        result[key] = str(value) if isinstance(value, Path) else value
    return result


def evaluate(
    grid: SpectralGrid,
    raw: torch.Tensor,
    child_mask: torch.Tensor,
    parent_mask: torch.Tensor,
    parent_scale: float,
    carrier: float,
    min_parent_strength: float,
    strength_penalty: float,
) -> tuple[torch.Tensor, dict[str, torch.Tensor], torch.Tensor, torch.Tensor]:
    child, child_hat = grid.normalize_field(raw, child_mask)
    force_hat = grid.quadratic_force(child, child_hat)
    parent_hat = force_hat * parent_mask
    leakage_hat = force_hat * (1.0 - parent_mask)

    parent_norm = grid.l2_hat(parent_hat)
    leakage_norm = grid.l2_hat(leakage_hat)
    total_norm = grid.l2_hat(force_hat)
    leakage_ratio = leakage_norm / parent_norm
    parent_strength = parent_norm / (parent_scale**2.5)

    # The logarithmic ratio is scale-free.  The one-sided penalty prevents the
    # optimizer from approaching a Beltrami field for which both numerator and
    # denominator vanish.
    shortfall = torch.relu(
        torch.log(
            torch.as_tensor(
                min_parent_strength, device=grid.device, dtype=grid.dtype
            )
            / parent_strength
        )
    )
    loss = torch.log(leakage_ratio) + strength_penalty * shortfall.square()

    parent_unit_hat = -parent_hat / parent_norm
    parent = grid.ifft(parent_unit_hat)
    child_speed = torch.sqrt(child.square().sum(dim=0))
    parent_speed = torch.sqrt(parent.square().sum(dim=0))
    l4 = torch.mean(child_speed**4).pow(0.25)
    mean_k = torch.sqrt(
        (
            parent_unit_hat.abs().square()
            * grid.kabs.square().unsqueeze(0)
        ).sum()
        / (parent_unit_hat.abs().square().sum() + torch.finfo(grid.dtype).tiny)
    )

    diagnostics = {
        "parent_norm": parent_norm,
        "leakage_norm": leakage_norm,
        "total_norm": total_norm,
        "leakage_ratio": leakage_ratio,
        "parent_fraction": parent_norm.square() / total_norm.square(),
        "parent_strength": parent_strength,
        "order2_scaled_ratio": leakage_ratio * (carrier / parent_scale) ** 2,
        "child_l4": l4,
        "child_linf": child_speed.max(),
        "parent_linf": parent_speed.max(),
        "parent_mean_wavenumber": mean_k,
    }
    return loss, diagnostics, child_hat, parent_unit_hat


def write_history(path: Path, rows: list[Metrics]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(rows[0])))
        writer.writeheader()
        writer.writerows(asdict(row) for row in rows)


def save_candidate(
    output_dir: Path,
    child_hat: torch.Tensor,
    parent_hat: torch.Tensor,
    metrics: Metrics,
    args: argparse.Namespace,
) -> None:
    np.savez_compressed(
        output_dir / "best_downshifter.npz",
        child_hat=child_hat.detach().cpu().numpy(),
        parent_hat=parent_hat.detach().cpu().numpy(),
        metrics_json=np.asarray(json.dumps(asdict(metrics), sort_keys=True)),
        args_json=np.asarray(json.dumps(jsonable_args(args), sort_keys=True)),
    )
    checkpoint = {
        "best": asdict(metrics),
        "configuration": jsonable_args(args),
        "interpretation": {
            "candidate_second_order_scaling": "leakage_ratio ~ (parent/carrier)^2",
            "success_indicator": (
                "order2_scaled_ratio stays bounded or falls as carrier/parent grows, "
                "at comparable parent_strength"
            ),
        },
    }
    (output_dir / "checkpoint_summary.json").write_text(
        json.dumps(checkpoint, indent=2, sort_keys=True) + "\n"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resolution", type=int, default=64)
    parser.add_argument("--carrier", type=int, default=12)
    parser.add_argument("--child-halfwidth", type=float, default=3.0)
    parser.add_argument("--parent-scale", type=float, default=3.0)
    parser.add_argument("--parent-low-factor", type=float, default=0.5)
    parser.add_argument("--parent-high-factor", type=float, default=2.0)
    parser.add_argument("--iterations", type=int, default=300)
    parser.add_argument("--restarts", type=int, default=3)
    parser.add_argument("--learning-rate", type=float, default=0.02)
    parser.add_argument("--min-parent-strength", type=float, default=0.01)
    parser.add_argument("--strength-penalty", type=float, default=8.0)
    parser.add_argument("--envelope-order", type=int, default=1)
    parser.add_argument("--seed", type=int, default=20260724)
    parser.add_argument("--report-every", type=int, default=10)
    parser.add_argument("--checkpoint-every", type=int, default=25)
    parser.add_argument("--precision", choices=("single", "double"), default="single")
    parser.add_argument("--device", choices=("auto", "cpu", "cuda"), default="auto")
    parser.add_argument("--output-dir", type=Path, default=Path("ns_downshifter_gate"))
    args = parser.parse_args()

    if args.resolution % 2:
        raise ValueError("resolution must be even")
    if args.carrier + args.child_halfwidth >= args.resolution / 2:
        raise ValueError("child band reaches the Nyquist frequency")
    if 2 * (args.carrier + args.child_halfwidth) >= args.resolution / 2:
        raise ValueError(
            "resolution is too small to represent the quadratic high-frequency leakage "
            "without aliasing; require 2*(carrier+child-halfwidth) < resolution/2"
        )
    if args.parent_scale <= 0 or args.carrier <= args.parent_scale:
        raise ValueError("require 0 < parent-scale < carrier")

    if args.device == "auto":
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    else:
        device = torch.device(args.device)
    if device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA requested but unavailable")
    dtype = torch.float64 if args.precision == "double" else torch.float32
    torch.manual_seed(args.seed)
    if device.type == "cuda":
        torch.cuda.manual_seed_all(args.seed)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    grid = SpectralGrid(args.resolution, device, dtype)
    child_mask = radial_mask(
        grid,
        args.carrier - args.child_halfwidth,
        args.carrier + args.child_halfwidth,
    )
    parent_mask = radial_mask(
        grid,
        args.parent_low_factor * args.parent_scale,
        args.parent_high_factor * args.parent_scale,
    )
    if torch.count_nonzero(child_mask) == 0 or torch.count_nonzero(parent_mask) == 0:
        raise ValueError("one of the spectral masks is empty")

    config = {
        **jsonable_args(args),
        "device_resolved": str(device),
        "dtype_resolved": str(dtype),
        "child_mode_count": int(torch.count_nonzero(child_mask).cpu()),
        "parent_mode_count": int(torch.count_nonzero(parent_mask).cpu()),
        "carrier_parent_ratio": args.carrier / args.parent_scale,
    }
    (args.output_dir / "configuration.json").write_text(
        json.dumps(config, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(config, indent=2, sort_keys=True), flush=True)

    history: list[Metrics] = []
    best: Metrics | None = None
    best_child: torch.Tensor | None = None
    best_parent: torch.Tensor | None = None
    started = time.time()

    for restart in range(args.restarts):
        generator = torch.Generator(device=device)
        generator.manual_seed(args.seed + restart)
        if restart == 0:
            initial = helical_initialization(
                grid, args.carrier, max(1, args.envelope_order)
            )
            initial = initial + 0.02 * random_initialization(grid, generator)
        else:
            initial = random_initialization(grid, generator)
        raw = torch.nn.Parameter(initial)
        optimizer = torch.optim.Adam([raw], lr=args.learning_rate)

        for iteration in range(args.iterations + 1):
            optimizer.zero_grad(set_to_none=True)
            loss, diagnostics, child_hat, parent_hat = evaluate(
                grid,
                raw,
                child_mask,
                parent_mask,
                args.parent_scale,
                args.carrier,
                args.min_parent_strength,
                args.strength_penalty,
            )
            record = Metrics(
                iteration=iteration,
                restart=restart,
                loss=scalar(loss),
                parent_norm=scalar(diagnostics["parent_norm"]),
                leakage_norm=scalar(diagnostics["leakage_norm"]),
                total_norm=scalar(diagnostics["total_norm"]),
                leakage_ratio=scalar(diagnostics["leakage_ratio"]),
                parent_fraction=scalar(diagnostics["parent_fraction"]),
                parent_strength=scalar(diagnostics["parent_strength"]),
                order2_scaled_ratio=scalar(diagnostics["order2_scaled_ratio"]),
                child_l4=scalar(diagnostics["child_l4"]),
                child_linf=scalar(diagnostics["child_linf"]),
                parent_linf=scalar(diagnostics["parent_linf"]),
                parent_mean_wavenumber=scalar(
                    diagnostics["parent_mean_wavenumber"]
                ),
                wall_seconds=time.time() - started,
            )
            history.append(record)

            feasible = record.parent_strength >= 0.8 * args.min_parent_strength
            if feasible and (
                best is None or record.leakage_ratio < best.leakage_ratio
            ):
                best = record
                best_child = child_hat.detach().clone()
                best_parent = parent_hat.detach().clone()

            if (
                iteration % args.report_every == 0
                or iteration == args.iterations
            ):
                print(
                    f"restart={restart} iter={iteration:04d} "
                    f"parent={record.parent_norm:.6g} "
                    f"strength={record.parent_strength:.6g} "
                    f"leak/parent={record.leakage_ratio:.6g} "
                    f"rho2_scaled={record.order2_scaled_ratio:.6g} "
                    f"parent_frac={record.parent_fraction:.6g}",
                    flush=True,
                )

            if (
                best is not None
                and best_child is not None
                and best_parent is not None
                and (
                    iteration % args.checkpoint_every == 0
                    or iteration == args.iterations
                )
            ):
                write_history(args.output_dir / "optimization_history.csv", history)
                save_candidate(
                    args.output_dir, best_child, best_parent, best, args
                )

            if iteration == args.iterations:
                break
            loss.backward()
            torch.nn.utils.clip_grad_norm_([raw], max_norm=10.0)
            optimizer.step()
            # Remove directions that never influence the projected field.  This
            # also prevents unconstrained low modes in the raw parameter from
            # accumulating optimizer state.
            with torch.no_grad():
                projected, _ = grid.normalize_field(raw, child_mask)
                raw.copy_(projected)

    write_history(args.output_dir / "optimization_history.csv", history)
    if best is None or best_child is None or best_parent is None:
        raise RuntimeError(
            "No candidate met 80% of min-parent-strength. Lower the threshold "
            "only after inspecting the reported parent strengths."
        )
    save_candidate(args.output_dir, best_child, best_parent, best, args)
    summary: dict[str, Any] = {
        "best": asdict(best),
        "configuration": config,
        "decision_metrics": {
            "carrier_parent_ratio": args.carrier / args.parent_scale,
            "leakage_ratio": best.leakage_ratio,
            "order2_scaled_ratio": best.order2_scaled_ratio,
            "parent_strength": best.parent_strength,
        },
        "warning": (
            "This is a finite-resolution numerical gate, not a Navier--Stokes "
            "regularity or blow-up proof."
        ),
    }
    (args.output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    print("BEST", json.dumps(asdict(best), indent=2, sort_keys=True), flush=True)
    print(f"wrote {args.output_dir.resolve()}", flush=True)


if __name__ == "__main__":
    main()
