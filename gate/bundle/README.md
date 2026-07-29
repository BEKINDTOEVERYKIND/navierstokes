# Navier–Stokes quadratic downshifter gate

This bundle contains:

- `GATE_AUDIT.md`: the mathematical reduction and exact parameter window;
- `parameter_certificate.py`: exact rational checks of every exponent margin;
- `ns_downshifter_gate.py`: the GPU search;
- `COLAB_CELLS.md`: copy/paste Colab cells with direct-to-Drive checkpoints.

Run the exact certificate locally with:

```bash
python parameter_certificate.py
```

For Colab, follow `COLAB_CELLS.md`. The output is written directly under
`MyDrive/ns_breakthrough_gate_results`, so a recycled runtime does not destroy
the checkpoints.

This is a research gate, not a claimed proof of Navier–Stokes blow-up or
regularity.
