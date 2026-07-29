#!/usr/bin/env python3
"""Exact finite-Fourier counterexample to a K^{-1} energy-flux bound."""

from __future__ import annotations

from collections import defaultdict


Vector = tuple[complex, complex, complex]
Mode = tuple[int, int, int]


def dot(a: Vector, b: tuple[int, int, int] | Vector) -> complex:
    return sum(x * y for x, y in zip(a, b))


def add(a: Vector, b: Vector) -> Vector:
    return tuple(x + y for x, y in zip(a, b))  # type: ignore[return-value]


def scale(c: complex, a: Vector) -> Vector:
    return tuple(c * x for x in a)  # type: ignore[return-value]


def leray(n: Mode, a: Vector) -> Vector:
    nn = sum(x * x for x in n)
    coefficient = dot(a, n) / nn
    return tuple(a[i] - coefficient * n[i] for i in range(3))  # type: ignore[return-value]


def neg(n: Mode) -> Mode:
    return (-n[0], -n[1], -n[2])


def nonlinear_at(n: Mode, modes: dict[Mode, Vector]) -> Vector:
    raw: Vector = (0j, 0j, 0j)
    for r, ur in modes.items():
        s = (n[0] - r[0], n[1] - r[1], n[2] - r[2])
        us = modes.get(s)
        if us is not None:
            raw = add(raw, scale(dot(ur, s), us))
    return scale(1j, leray(n, raw))


def norm2(a: Vector) -> float:
    return float(sum(abs(x) ** 2 for x in a))


def main() -> None:
    for nscale in (1, 2, 4, 8, 16):
        amplitude = 3.0
        p = (nscale, 0, 0)
        q = (0, nscale, 0)
        k = (nscale, nscale, 0)
        ep: Vector = (0j, amplitude + 0j, 0j)
        eq: Vector = (0j, 0j, amplitude + 0j)
        ek: Vector = (0j, 0j, -1j * amplitude)
        modes = {
            p: ep,
            neg(p): tuple(x.conjugate() for x in ep),
            q: eq,
            neg(q): tuple(x.conjugate() for x in eq),
            k: ek,
            neg(k): tuple(x.conjugate() for x in ek),
        }

        for mode, value in modes.items():
            assert abs(dot(value, mode)) < 1e-12
            conjugate = tuple(x.conjugate() for x in value)
            assert modes[neg(mode)] == conjugate

        flux = 0.0
        for mode, value in modes.items():
            if sum(x * x for x in mode) > nscale * nscale:
                rhs = scale(-1, nonlinear_at(mode, modes))
                flux += 2 * dot(
                    tuple(x.conjugate() for x in value), rhs
                ).real

        energy = sum(norm2(value) for value in modes.values())
        high_energy = sum(
            norm2(value)
            for mode, value in modes.items()
            if sum(x * x for x in mode) > nscale * nscale
        )
        expected = 4 * nscale * amplitude**3
        assert abs(flux - expected) < 1e-9
        ratio_without_c = (
            abs(flux) * nscale / (high_energy * energy) ** 0.5
        )
        print(
            f"N={nscale:2d} flux={flux:8.1f} "
            f"E={energy:.1f} E_hi={high_energy:.1f} "
            f"required_C={ratio_without_c:.6g}"
        )

    print("exact triad flux scales as N A^3; advertised RHS scales as A^2/N")


if __name__ == "__main__":
    main()
