# Thermal Protections

Thermal protection system for the Stark-100k motor: RSZ/YSZ thermal barrier
coating (TBC) sizing on the Inconel 718 chamber wall and on the nozzle
convergent, with transient FEM-like analysis, series-resistance models and
alternative C/C configurations. Full derivation and results are in
[`stark-100k-report.pdf`](../../report/stark-100k-report.pdf).

## Key results

### Combustion-chamber wall (RSZ on Inconel 718)

| Quantity | Value |
|---|---|
| Required RSZ thickness | 0.7242 mm |
| RSZ surface temperature | 2144.02 K |
| RSZ/metal interface (target) | 1300.04 K |
| Alloy external face | 1210.33 K |
| Gas-side coefficient h_g | 964.74 W/(m²·K) |

The transient analysis (`convergent_tbc_sizing.m`) sizes the convergent-section
protection with the throat heat-transfer coefficient: RSZ thickness
**1.1828 mm**, h_g = 7312.99 W/(m²·K), Bi = 14.42. At 29 s the RSZ hot wall
reaches 2734.42 K.

### Nozzle convergent: passive vs active cooling (`convergent_passive_vs_active_cooling.m`)

Air-only (passive) and water-jacket (active) cases at A/A_t = 2:

| Case | Material | Required TBC thickness |
|---|---|---|
| Air only (passive) | CuCrZr | 8.668 mm |
| Air only (passive) | Inconel 718 | 1.972 mm |
| Water jacket (active) | CuCrZr | 0.293 mm |
| Water jacket (active) | Inconel 718 | 0.197 mm |

The air-only thicknesses are outside the normally used TBC range, which is why
the active cooling case is required (see the cooling jacket module). The
verification blocks report the TBC gas-side surface temperature: it exceeds
the 1600 K design limit in all cases (e.g. 2276.9 K for the active CuCrZr
case), consistent with the report, which notes that the gas–TBC interface
temperature is well above recommended operating conditions but below the YSZ
melting point (~2700 °C).

### C/C and bilayer alternatives

- `carbon_carbon_single_layer_tbc.m` — parametric sweep of C/C thermal
  conductivity (2–233 W/(m·K)) with a single TBC layer: required thickness
  drops from 3.206 mm (k = 2) to 0.328 mm (k = 233); interface held at 2200 K
  (C/C limit).
- `bilayer_tbc_inconel.m` — RSZ top coat (1.905 mm) over a fixed 0.100 mm YSZ
  bond coat on Inconel 718: total 2.005 mm, heat flux 0.449 MW/m².
- `bilayer_tbc_carbon_carbon.m` — same bilayer on C/C, resistance model: total
  thickness 3.40 mm (k = 2) down to 0.40 mm (k = 233).

## Run order

MATLAB R2026a. All scripts are standalone; run from this folder in any order:

```matlab
run('chamber_wall_tbc_sizing.m');            % chamber-wall sizing (FEM + bisection)
run('convergent_tbc_sizing.m');              % convergent sizing, transient (FEM + bisection)
run('convergent_passive_vs_active_cooling.m');  % convergent: passive vs active cooling
run('carbon_carbon_single_layer_tbc.m');     % C/C configurations with single TBC layer
run('bilayer_tbc_inconel.m');                % bilayer TBC (RSZ + YSZ) on Inconel 718
run('bilayer_tbc_carbon_carbon.m');          % bilayer TBC (RSZ + YSZ) on C/C
```

## Notes

- `convergent_passive_vs_active_cooling.m` originally computed the
  water-jacket verification with the air-only network (copy-paste error in the
  `fsolve` call). It is fixed in this repository so the active-case thickness
  (0.293 / 0.197 mm) is used consistently in the verification block.
- The scripts print the required thickness and temperature tables shown above;
  they do not export figures.

## Folder contents

```
thermal-protection/
├── chamber_wall_tbc_sizing.m            # chamber wall, RSZ sizing, transient
├── convergent_tbc_sizing.m              # convergent, RSZ sizing, transient tables
├── convergent_passive_vs_active_cooling.m  # passive vs active (water jacket) cooling
├── carbon_carbon_single_layer_tbc.m     # C/C conductivity sweep, single TBC layer
├── bilayer_tbc_inconel.m                # bilayer TBC (RSZ + YSZ) on Inconel 718
└── bilayer_tbc_carbon_carbon.m          # bilayer TBC (RSZ + YSZ) on C/C, resistance model
```