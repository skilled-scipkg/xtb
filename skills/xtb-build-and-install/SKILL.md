---
name: xtb-build-and-install
description: Use this skill for xtb build, compile, test, and install workflows (Meson/CMake), including toolchain selection, LAPACK backend choices, install prefixes, and post-install simulation smoke checks.
---

# xtb: Build and Install

## High-Signal Playbook
### Route conditions
- Use this skill for compiler/toolchain selection, configure/build/test/install flow, install prefix setup, and build flags (`meson/README.adoc`, `cmake/README.adoc`).
- Route runtime input/modeling questions to `xtb-inputs-and-modeling`.
- Route `xcontrol` instruction syntax and parser behavior to `xtb-xcontrol-7`.

### Triage questions
- Which build system is required: Meson+Ninja or CMake (+Ninja/Make)?
- Which compiler family is available: Intel (`ifort`/`icc`) or GCC (`gfortran`/`gcc`)?
- Which LAPACK backend is expected (MKL/OpenBLAS/netlib/custom)?
- Is this a local developer build, CI build, or production install with a custom prefix?
- Are tests mandatory before install?
- Are there platform constraints (Darwin requires real GCC, not clang wrappers)?

### Canonical workflow
1. Pick one build system and one compiler family (`meson/README.adoc`, `cmake/README.adoc`).
2. Export Fortran/C compilers before configure.
3. Configure release build (`--buildtype release --optimization 2` for Meson, `-DCMAKE_BUILD_TYPE=Release` for CMake).
4. Build using the same generator used during configure.
5. Run tests before installation.
6. Set install prefix (`$HOME/.local` for user install) before install.
7. Run a post-install xtb smoke command with a known input.

### Minimal working examples
```bash
# Meson + GCC + OpenBLAS
export FC=gfortran CC=gcc
meson setup build --buildtype release --optimization 2 -Dlapack=openblas
ninja -C build
ninja -C build test
meson configure build --prefix "$HOME/.local"
ninja -C build install
```

```bash
# CMake + GCC + Ninja
export FC=gfortran CC=gcc
cmake -B _build -S . -GNinja -DCMAKE_BUILD_TYPE=Release
ninja -C _build
ctest --test-dir _build --output-on-failure
cmake -B _build -S . -DCMAKE_INSTALL_PREFIX="$HOME/.local"
ninja -C _build install
```

```bash
# Post-install simulation smoke check
export XTBPATH="$PWD"
xtb assets/inputs/coord/caffeine.coord --gfn 2 --sp --strict --namespace smoke1
```

### Pitfalls and fixes
- GCC 7 build/runtime issues: use GCC 8+ or apply feedstock patches (`meson/README.adoc`, `cmake/README.adoc`).
- Darwin compiler confusion: ensure `gcc` is GNU GCC, not clang alias (`meson/README.adoc`, `cmake/README.adoc`).
- Numerical instability with aggressive optimization: keep `-O2` unless validated (`meson/README.adoc`).
- Mixed build directory names (`build` vs `_build`): use one consistently.
- Install permission failures under `/usr/local`: use user prefix or privileged install flow.
- Conda MKL edge case: when using conda-forge MKL, use netlib backend (`meson/README.adoc`).

### Validation checkpoints
- Tool minimums pass before configure (`meson >= 0.51`, `ninja >= 1.7`; `cmake >= 3.9`, `ninja >= 1.10`).
- Build completes without linker errors for the chosen backend.
- Tests pass (`ninja -C build test` or `ctest --test-dir _build --output-on-failure`).
- Runtime smoke command exits cleanly and produces an energy summary.
- `XTBPATH` resolves parameter files during smoke tests (`test/unit/CMakeLists.txt`, `test/api/CMakeLists.txt`).

### Source drill-down (when docs stop)
- Build graph wiring: `src/CMakeLists.txt`, `src/meson.build`
- CLI startup and option handling: `src/prog/main.F90` (`xtbMain`, `parseArguments`)
- Environment and strict/error behavior: `src/type/environment.f90` (`initEnvironment`, `warning`, `error`)
- Parameter/input file lookup: `src/readin.f90` (`xfind`, `get*Value`)
- Namespace-compatible output naming: `src/setparam.f90` (`get_namespace`)
- C API simulation path: `include/xtb.h`, `src/api/interface.f90`, `src/api/calculator.f90`, `src/api/results.f90`
- Executable behavior checks: `test/unit/CMakeLists.txt`, `test/api/CMakeLists.txt`, `test/api/c_api_example.c`

## Scope
- Handle questions about build, installation, compilation, and environment setup.
- Keep responses architectural first; only go function-level when needed to resolve ambiguity.

## Primary documentation references
- `meson/README.adoc`
- `cmake/README.adoc`
- `test/CMakeLists.txt`
- `cmake/CMakeLists.txt`
- `test/unit/CMakeLists.txt`
- `test/api/CMakeLists.txt`

## Workflow
- Start with the primary references above.
- If details are missing, inspect `references/doc_map.md` for more runnable docs/examples.
- Use tests as behavior/regression references for smoke and install checks.
- If ambiguity remains after docs, inspect `references/source_map.md` and start with the ranked source entry points.
- Cite exact file paths when answering.

## Tutorials and examples
- `assets/inputs`

## Test references
- `test`

## Optional deeper inspection
- `include`
- `src`
- `test`

## Source entry points for unresolved issues
- `src/prog/main.F90` (`xtbMain`, `parseArguments`)
- `src/type/environment.f90` (`initEnvironment`, `warning`, `error`)
- `src/readin.f90` (`xfind`, `getIntValue`, `getRealValue`, `getBoolValue`)
- `src/setparam.f90` (`get_namespace`)
- `src/api/interface.f90` (`singlepoint_api`, `hessian_api`, `cpcmx_calc_api`)
- `src/api/calculator.f90` (`loadGFN2xTB_api`, `setSolvent_api`, `setAccuracy_api`, `setMaxIter_api`)
- `src/api/results.f90` (`getEnergy_api`, `getGradient_api`, `getCharges_api`)
- `include/xtb.h`
- `src/CMakeLists.txt`
- `src/meson.build`
- `test/unit/CMakeLists.txt`
- `test/api/CMakeLists.txt`
- `test/api/c_api_example.c`
- Prefer targeted source search: `rg -n "<symbol_or_keyword>" include src test`.
