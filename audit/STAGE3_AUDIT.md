# Stage 3 numerical audit

Status: finite-resolution numerical evidence only. This is not a proof, a
computer-assisted proof, or evidence of a Navier--Stokes singularity.

## Recovered archive

`stage3_results.zip` was recovered intact:

- byte size: `3,097,033`
- SHA-256: `2fa51b2566016b90a3c4a6a6c4f76bfc8825e42d87f8aaface09ac8fe435e5ac`
- all 59 ZIP members passed the integrity test

## Inviscid-limit scan

At `N=192`, float64, the first-octave Q3 values were:

| R0 | baseline | direct-S2 | robust-p16 |
|---:|---:|---:|---:|
| 30 | 0.767305 | 0.758792 | 0.769280 |
| 80 | 0.787340 | 0.778037 | 0.789455 |
| 160 | 0.793365 | 0.783850 | 0.795508 |
| 320 | 0.796371 | 0.786755 | 0.798523 |
| 640 | 0.797870 | 0.788206 | 0.800026 |
| 1280 | 0.798619 | 0.788931 | 0.800776 |

A linear fit in `1/R0` over the final four points gives extrapolated limits
0.79937, 0.78966, and 0.80153 respectively. Thus viscosity explains only
about 0.03 of the Q3 deficit at R0=30 for these seeds. The dominant loss is
geometric/dephasing rather than viscous.

Critical Fourier content itself is not enough: at R0=1280 the baseline has
`max C1/C0 = 1.09991` while `max Q3 = 0.79862`.

## Robust-p16 through octave 2

At N=256:

| octave | max Cj/C0 | max Q3 | max Q8 | max Q16 | max Rinf |
|---:|---:|---:|---:|---:|---:|
| 1 | 0.963604 | 0.769421 | 0.598392 | 0.606511 | 0.649038 |
| 2 | 0.475738 | 0.388609 | 0.193198 | 0.179382 | 0.182493 |

The first-octave pointwise advantage does not persist through a second
transfer.

## First-octave adversarial Q3 seed

The global N=64 search found Q3=0.88224. N=96 polishing raised its internal
score to 0.89484. Independent float64 validation gives:

| resolution | max Q3 | peak tau |
|---:|---:|---:|
| 128 | 0.896704691 | 0.7850 |
| 192 | 0.897466292 | 0.7875 |
| 256 | 0.897641323 | 0.7875 |

The N=192 to N=256 change is 0.0195%. At N=256, halving `dtau` from 0.0025
to 0.00125 changes max Q3 by only `2.7e-10`. This is a well-resolved
finite-dimensional effect.

At the N=256 Q3 peak:

- critical fraction: 0.789957
- scale ratio: 1.507461
- Q8: 0.806187
- Q16: 0.795724
- Rinf: 0.814739

Relative to robust-p16, max Q3 improves by 16.67%.

## Scale-step loophole

The optimizer moved 70.89% of the source's critical weight into `[7,8)`,
versus 37.52% for robust-p16, and approached the allowed minimum scale ratio.
Raw Q3 therefore overstates its progress toward a full octave.

For fair comparison define

`Q3_oct = exp(log(2) * log(Q3) / log(lambda))`,

the multiplicative retention corresponding to one full frequency doubling.
At N=256:

| seed/target | Q3 | lambda | Q3 per octave |
|---|---:|---:|---:|
| robust-p16, shell 1 | 0.769421 | 1.676654 | 0.703590 |
| Q3-polished, shell 1 | 0.897641 | 1.507461 | 0.833293 |

The improvement remains substantial after correction, but it is not a
near-recurrence.

## Second-octave adversary

The N=96 optimizer predicted Q3=0.53380 in `[16,32)`. Independent validation
gives 0.43996, 0.43816, and 0.44558 at N=128, 192, and 256. The N=256 value
corresponds to Q3 per octave 0.61376.

The second shell is not fully converged: from N=192 to N=256, Q3 changes
1.69%, Q16 changes 13.36%, and Rinf changes 19.00%. A 384-grid validation is
appropriate, although the result is already far below recurrence.

## Next tests

1. Scan the Q3-polished seed over R0 up to at least 1280.
2. Optimize Q3 retention per octave instead of raw Q3.
3. Optimize directly at high R0 and validate any crossing at 128/192/256 with
   a half time step.
4. Run N=384 on the second-octave candidate.
5. Independently validate the other search's reported Q3 approximately 0.99
   seed as soon as its NPZ is available.

