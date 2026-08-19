# Stark-100k SRM

Pre-design of the **Stark-100k**, a solid-propellant motor delivering 100 kN
of thrust at sea level over a total impulse of 2.5 MN·s, with a near-constant
thrust profile from BATES-style radial and axial burning. The primary
engineering objective is the **quantification of the heat load at the nozzle
throat**: the throat carries no insert and is instrumented with a straight,
co-flowing water cooling jacket spanning area ratio 2 (subsonic) to 2
(supersonic), whose water temperature rise yields the transferred heat. This
repository documents the ballistic design and performance of the motor, its
thermal protections, and the cooling jacket together with the supporting tank
sizing.

The single authoritative source is
[`stark-100k-report.pdf`](report/stark-100k-report.pdf): every module folder
mirrors one section of the report, and every number in these READMEs traces
back to it.

## Contributors

- [Matteo Mafrici](https://github.com/matteomafrici)
- [Mario Guida](https://github.com/marioguida27)
- Juan Álvaro Cobos Franco
- Martina Lucia Magarelli
- Camilla Martino
- Margherita Palitta
- Jules Jean Laurence Simon

## Objective

The Stark-100k is a solid-propellant unit providing 100 kN of thrust at sea
level over a total impulse of 2.5 MN·s. The motor burns a standard
non-aluminized AP/HTPB propellant at about 70 bar with a quasi-constant
thrust profile: the BATES grain burns radially and axially so that the
initial and final chamber pressures coincide. The nozzle is optimal at sea
level.

The engineering focus is the throat heat load. To monitor the steady-state
heat transfer, the throat does not use a conventional insert: a straight
co-flow water jacket (pressurized, tap-water inlet at 18 °C) covers the
region between area ratio 2 (subsonic) and 2 (supersonic), and the
transferred heat is derived from the water temperature rise. A thermal
barrier coating is considered, but its use must be justified with numbers —
which is exactly what the thermal-protection module provides. The jacket is
reusable; the rest of the motor is expendable, so the engine is not
performance-optimized beyond its primary instrumentation role.

### Design targets

| Quantity | Target | Achieved (design point) |
|---|---|---|
| Mean thrust | 100 kN | 99596.86 N (−0.40%) |
| Total impulse | 2.5 MN·s | 2544634.90 N·s (+1.79%) |
| Burn time | ~25 s | 25.54 s |
| Chamber pressure | ≤ 70 bar | MEOP 66.06 bar |
| Propellant mass | 1.09 t (AP/HTPB, O/F = 6) | 1093.22 kg (0.03% error) |

### Motor overview

The motor is 0.8714 m in outer diameter with a port diameter of 0.4738 m,
a single segment 1.5441 m long and a conical nozzle of 0.1145 m throat
diameter. The propellant burning rate follows a Vieille law
`r = 1.617743 · P^0.381551 mm/s` (R² = 0.990) fitted from strand-burner data.
Ignition is deep in the non-erosive regime: core Mach number 0.0348 and mass
flux 174.7 kg/(m²·s) against a non-erosive threshold of 703 kg/(m²·s) at
65 bar.

The chamber wall carries a RSZ thermal barrier coating (0.7242 mm) to keep the
metal/ceramic interface at 1300 K; the nozzle throat is cooled by a water
jacket sized on the plane, cylindrical and Bartz heat-transfer models (TBC
39 / 37 / 20 µm, coolant velocity 14.94–15.05 m/s against a 15 m/s limit), and
the blowdown pressurant tank (400 L, 16.45 bar) feeds the jacket over the
burn.

## Module overview

Each module folder is a self-contained design package mirroring the
corresponding section of the report:

| Module | Folder | Scope |
|---|---|---|
| Ballistic design & performance | [`modules/ballistic-design`](modules/ballistic-design) | Vieille-law fit, BATES grain and nozzle sizing, CEA-in-the-loop transient diagnosis, Monte Carlo uncertainty, erosive-burning assessment |
| Thermal protections | [`modules/thermal-protection`](modules/thermal-protection) | RSZ/YSZ TBC sizing (transient FEM + resistance models), C/C configurations, passive vs active cooling at the convergent |
| Cooling jacket & tank | [`modules/cooling-jacket-tank`](modules/cooling-jacket-tank) | Throat water jacket (thermal + hydraulic), coolant velocity loop, preliminary blowdown tank sizing |

## Key results

| Quantity | Value |
|---|---|
| MEOP | 66.06 bar |
| Total impulse | 2544634.90 N·s |
| Burn time | 25.54 s |
| Isp (theoretical / real) | 245.47 s / 228.53 s |
| Chamber-wall TBC (RSZ) | 0.7242 mm |
| Convergent RSZ (throat h_g) | 1.1828 mm |
| Convergent active cooling | 0.293 mm (CuCrZr) / 0.197 mm (Inconel 718) |
| Jacket TBC (YSZ, plane) | 39 µm |
| Jacket water mass flow | 2.30 kg/s |
| Tank (blowdown) | 400 L, P_in = 16.45 bar |

Monte Carlo uncertainty over the ballistic parameters (2000 samples, fixed
seed) gives MEOP 61.80 ± 3.40 bar and total impulse 2.38 ± 0.03 MN·s; the
report values (62.06 ± 6.29 bar) come from the same model without a fixed
seed.

## What you can clone and use

The repository is designed to be cloned and explored without external
dependencies beyond the free tools listed below:

- **Read-only usage** — every `modules/*/README.md` is a complete summary; no
  tool is needed to read the design.
- **MATLAB** — runnable pipelines live in each `modules/*/` folder
  (MATLAB R2026a). Output data and figures land in the corresponding
  `results/` folders.
- **NASA CEA** — the CEA executable ships inside the ballistic module
  (`modules/ballistic-design/cea/`) so the transient diagnosis runs
  out-of-the-box on Windows, Linux and macOS.
- **Report** — [`stark-100k-report.pdf`](report/stark-100k-report.pdf) is the
  single source for every file in this repo.

## Reproducibility

Everything in this repository is reproducible from the sources below. Each
module README documents the exact run order and expected outputs for its own
pipelines.

- **Simulations** — the MATLAB pipelines in `modules/*/` are self-contained;
  their inputs are tracked in the repo, so re-running them regenerates the
  same figures and numbers. The Monte Carlo script uses a fixed RNG seed
  (`rng(0)`) for exact reproducibility.
- **Report** — every number in these READMEs traces back to
  [`stark-100k-report.pdf`](report/stark-100k-report.pdf), the single
  authoritative source for the design.

## License

The design work in this repository (code, data, text, figures) is released
under the [MIT License](LICENSE). NASA CEA, used by the ballistic module, is
open software distributed by NASA and subject to its own terms.

## Acronyms

| Acronym | Meaning |
|---|---|
| AP / HTPB | Ammonium perchlorate / hydroxyl-terminated polybutadiene |
| BATES | Ballistic Analysis and Technical Evaluation System (grain shape) |
| CEA | Chemical Equilibrium with Applications (NASA) |
| Kn | Ratio of burning area to throat area |
| MEOP | Maximum expected operating pressure |
| RSZ / YSZ | Rare-earth-stabilized / yttria-stabilized zirconia |
| SRM | Solid rocket motor |
| TBC | Thermal barrier coating |