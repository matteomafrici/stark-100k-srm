# Cooling Jacket and Tank Design

Thermal and hydraulic design of the nozzle-throat water cooling jacket for the
Stark-100k: gas-side heat transfer, YSZ TBC sizing on the coolant side,
required vs actual coolant velocity and tube count, and the preliminary
blowdown tank sizing. Full derivation and results are in
[`stark-100k-report.pdf`](../../report/stark-100k-report.pdf).

## Key results

### Gas properties at the throat

| Quantity | Value |
|---|---|
| T_throat | 2589.2 K |
| p_throat | 37.1521 bar |
| ρ_throat | 4.3146 kg/m³ |
| v_throat | 1001.14 m/s |
| Pr_gas | 0.6088 |

### Nozzle geometry and jacket

| Quantity | Value |
|---|---|
| Throat diameter D_t | 0.1146 m |
| Radius at ε = 2 | 0.0810 m |
| L_convergent / L_divergent | 0.0237 m / 0.0885 m |
| Jacket length | 0.1123 m |
| Wetted area A_jacket | 0.05439 m² |

### Coolant-side design (three heat-transfer models)

Optimal YSZ thickness on the coolant side:

| Model | TBC thickness | Required v | Real v | Water m_dot |
|---|---|---|---|---|
| Plane | 39 µm | 14.8033 m/s | 14.9361 m/s | 2.3041 kg/s |
| Cylindrical | 37 µm | 14.9020 m/s | 15.0460 m/s | 2.3164 kg/s |
| Bartz | 20 µm | 14.8416 m/s | 14.9670 m/s | 2.3089 kg/s |

The real velocity always exceeds the required one (v_real ≥ v_req, the
wall-temperature constraint), so the cooling requirement is met. The
TBC/gas interface temperature is 2500.7 K (plane), 2499.2 K (cylindrical) and
2038.3 K (Bartz) — the Bartz model gives the most conservative interface
temperature; the report adopts the plane model as reference (39 µm).

Water side: p_coolant = 9.75 bar, T_sat = 455.1 K, T_out design = 323.1 K
(50 °C, 132 K margin), Pr_water = 3.5707.

### Tank sizing (blowdown, from the report)

The 400 L pressurant tank with initial pressure P_in = 16.45 bar provides the
required outflow over the burn. This result comes from the report (Section 5.4)
and is not recomputed by this module.

## Run order

MATLAB R2026a. Single script:

```matlab
run('throat_cooling_jacket.m');
```

It prints the gas properties, geometry, adiabatic wall temperature, coolant
properties, the three-model TBC and velocity results, and produces eleven
figures (h_g and T_aw profiles, wall temperatures vs TBC thickness, nozzle
profile, Mach distribution, heat flux, convergence of v_req vs thickness and
of v_real vs tube count), exported to `results/`.

## Notes

- The wall-velocity constraint (v_req < 15 m/s) forces the small TBC thickness
  (tens of µm); the standalone wall-temperature constraint alone would allow
  ~350 µm. Both constraints are evaluated in the report.
- Propulsion inputs used here come from `diagnosis.m` of the ballistic module
  (m_dot = 43.377 kg/s, c* = 1544.3 m/s).

## Folder contents

```
cooling-jacket-tank/
├── throat_cooling_jacket.m   # throat jacket thermal/hydraulic design + figures
└── results/                  # generated figures
```