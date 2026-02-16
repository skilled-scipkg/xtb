# xtb documentation map: Inputs and Modeling

Practical references for setting up and validating xtb simulations from CLI inputs.

## Primary manuals
- `man/xtb.1.adoc` | CLI flags, run types, methods, solvent options, warning/strict behavior, output files.

## Runnable input artifacts
- `assets/inputs/coord/caffeine.coord` | quick single-point and optimization starter geometry.
- `assets/inputs/coord/quartz.3d.coord` | periodic-style coordinate example for robustness checks.
- `assets/inputs/vasp/ammonia.vasp` | VASP-format parser sanity check.
- `assets/inputs/xyz/taxol.xyz` | larger XYZ input for parser/performance checks.

## Behavior references from tests
- `test/unit/CMakeLists.txt` | canonical CLI smoke commands (`--strict`, `--gbsa`, method variants) and `XTBPATH` usage.
