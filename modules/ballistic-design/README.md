# Ballistic Design and Performance

Solid rocket motor pre-design for the Stark-100k: Vieille-law fit from
strand-burner data, BATES grain and nozzle sizing, transient diagnosis driven
by NASA CEA in the loop, Monte Carlo uncertainty, and erosive-burning
assessment. Full derivation and results are in
[`stark-100k-report.pdf`](../../report/stark-100k-report.pdf).

## Key results

### Propellant and ballistic fit

| Quantity | Value |
|---|---|
| Vieille law `r = a·P^n` | a = 1.617743 mm/s/bar^n (95% CI ± 0.064570), n = 0.381551 (± 0.010621) |
| R² | 0.990027 |
| Burn rate at design Pc | 7.954788 mm/s |
| Formulation | AP/HTPB, O/F = 6.0, 85.71% AP / 14.29% HTPB, ρ = 1685.22 kg/m³ |

### Grain and nozzle (best BATES configuration)

| Quantity | Value |
|---|---|
| Segments | N = 1 |
| Outer / port diameter | 0.8714 m / 0.4738 m |
| Segment length | 1.5441 m |
| Kn range | 304.551 → 322.625 |
| Propellant mass | 1093.22 kg (0.03% from target) |
| Pressure flatness | 9.18% (peak-to-peak 5.96 bar) |
| Nozzle | conical, d_t = 0.1145 m, d_e = 0.3408 m, L = 0.8006 m |

### Transient diagnosis (real motor performance)

| Quantity | Value |
|---|---|
| Total impulse | 2544634.90 N·s (+1.79% vs target) |
| Burn time | 25.54 s |
| MEOP | 66.06 bar |
| Mean thrust | 99596.86 N (−0.40%) |
| Isp (theoretical / real) | 245.47 s / 228.53 s |

Assumed non-ideal corrections: η_c* = 0.950, η_Cf = 0.980.

### Monte Carlo uncertainty (2000 samples)

| Quantity | Mean ± std |
|---|---|
| MEOP | 61.80 ± 3.40 bar |
| Burn time | 25.83 ± 1.31 s |
| Total impulse | 2380535 ± 33658 N·s |

The Monte Carlo draws are reproducible: the script seeds the RNG
(`rng(0)`) before sampling. Values differ slightly from the report
(62.06 ± 6.29 bar, 25.91 ± 2.55 s, 2.380564e6 ± 3.362e4 N·s) because the
report used the same uncertainty model without a fixed seed.

### Erosive burning

The design operates deeply in the non-erosive regime: core Mach 0.0348 at the
aft end and mass flux 174.7 kg/(m²·s) at ignition. The Rogers thresholds at the
design chamber pressure (65 bar), as reported in the report, are
**G ≤ 703 kg/(m²·s)** (non-erosive) and **G ≈ 1757 kg/(m²·s)** (max recommended
erosivity).

## Run order

MATLAB R2026a. Run from this folder:

```matlab
run('design.m');                % fit, BATES grain, nozzle; writes results/*.png
run('diagnosis.m');             % transient + Monte Carlo (long: CEA in the loop)
run('erosive_burning_rodgers.m');  % erosive-burning regime and parametric study
```

`diagnosis.m` expects the workspace variables produced by `design.m` (run the
two in the same session, in this order). `erosive_burning_rodgers.m` is
standalone.

## Outputs

`design.m` writes `results/vieille_law_fit.png`, `results/2d_motor_sketch.png`
and `results/quasi_steady_evolution.png`. `diagnosis.m` writes
`results/diagnosis_results.mat` and five figures (`temporal_histories`,
`chamber_nozzle_geometry`, `cea_nozzle_snapshots`, `real_motor_performance`,
`monte_carlo_uncertainty`). `erosive_burning_rodgers.m` writes
`results/erosive_burning_parametric.png`.

Key expected console outputs: `MEOP = 66.06 bar`,
`Total impulse = 2544634.90 N*s`, `Avg Isp real = 228.53 s`, and the Monte
Carlo summary above.

## Notes

- The Rogers thresholds are hard-coded to the report values
  (703 / 1757 kg/(m²·s)). The delivered code interpolated the original Rogers
  table at 943 psia, giving 1272.19 / 2109.21 kg/(m²·s); the report values are
  authoritative.
- `diagnosis.m` runs the NASA CEA executable in the loop. The `cea/` folder
  ships the FCEA2 binaries for Windows, Linux and macOS together with
  `thermo.lib`, `trans.lib` and the reference `.inp`/`.out` runs. NASA CEA is
  open software distributed by NASA; on macOS the Gatekeeper quarantine must be
  removed before the first run (`xattr -d com.apple.quarantine cea/FCEA2_mac`).

## Folder contents

```
ballistic-design/
├── design.m                    # Vieille fit, BATES grain, nozzle geometry
├── diagnosis.m                 # transient diagnosis + Monte Carlo
├── erosive_burning_rodgers.m   # erosive-burning regime study
├── engine/
│   ├── run_cea_transient.m     # CEA-in-the-loop nozzle cases
│   └── Uncertainty.m           # Vieille-law regression (a, n, R², CI)
├── cea/                        # FCEA2 binaries + thermo/trans lib + runs
└── results/                    # generated figures and data
```