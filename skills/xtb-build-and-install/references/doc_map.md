# xtb documentation map: Build and Install

Practical references for configure/build/test/install and post-install runtime smoke checks.

## Build system docs
- `meson/README.adoc` | Meson prerequisites, compiler selection, LAPACK backend flags, test/install flow.
- `cmake/README.adoc` | CMake prerequisites, generator usage, compiler selection, test/install flow.
- `cmake/CMakeLists.txt` | install target wiring for exported headers/modules.
- `test/CMakeLists.txt` | test tree registration (`api`, `unit`).

## Runtime smoke and validation docs
- `test/unit/CMakeLists.txt` | CLI smoke tests (`--version`, `--help`, `--strict`, GBSA) and `XTBPATH` test environment.
- `test/api/CMakeLists.txt` | C API smoke executable and test registration.
- `test/api/c_api_example.c` | End-to-end API sequence for singlepoint, solvent, and Hessian workflows.

## Sample inputs for build verification
- `assets/inputs/coord/caffeine.coord`
- `assets/inputs/coord/quartz.3d.coord`
- `assets/inputs/vasp/ammonia.vasp`
- `assets/inputs/xyz/taxol.xyz`
