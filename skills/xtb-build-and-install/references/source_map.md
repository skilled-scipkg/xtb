# xtb source map: Build and Install

Generated from source roots:
- `include`
- `src`
- `test`

Use this map only after exhausting topic docs in `references/doc_map.md`.

## Topic query tokens
- `build`
- `cmake`
- `meson`
- `install`
- `test`
- `xtbpath`
- `strict`
- `api`

## Fast source navigation
- `rg -n "parseArguments|xtbMain|--strict|--define|--input" src/prog/main.F90`
- `rg -n "initEnvironment|XTBPATH|XTBHOME|warning|error" src/type/environment.f90`
- `rg -n "xfind|get(Int|Real|Bool|List)Value" src/readin.f90`
- `rg -n "singlepoint_api|hessian_api|setSolvent_api|getEnergy_api" src/api/*.f90 include/xtb.h`
- `rg -n "add_test|XTBPATH|tests_capi" test/unit/CMakeLists.txt test/api/CMakeLists.txt`

## Suggested source entry points
- `src/prog/main.F90` | `xtbMain`, `parseArguments`; runtime startup and CLI behavior after install.
- `src/type/environment.f90` | `initEnvironment`, `warning`, `error`; strict-mode and environment variable handling.
- `src/readin.f90` | `xfind`, `getIntValue`, `getRealValue`, `getBoolValue`; parameter lookup/parsing path.
- `src/setparam.f90` | `get_namespace`; namespaced output behavior used in smoke and regression runs.
- `src/api/interface.f90` | `singlepoint_api`, `cpcmx_calc_api`, `hessian_api`; high-level C API execution bridge.
- `src/api/calculator.f90` | `loadGFN2xTB_api`, `setSolvent_api`, `setAccuracy_api`, `setMaxIter_api`; calculator configuration path.
- `src/api/results.f90` | `getEnergy_api`, `getGradient_api`, `getCharges_api`; API result extraction behavior.
- `include/xtb.h` | exported C API entry points used by downstream simulation code.
- `test/api/c_api_example.c` | complete behavior check for API lifecycle from init to cleanup.
- `test/unit/CMakeLists.txt` | CLI smoke invocations and `XTBPATH` test environment contract.
- `test/api/CMakeLists.txt` | C API smoke target registration and test execution wiring.
- `src/CMakeLists.txt` | CMake subdirectory graph for compile target composition.
- `src/meson.build` | Meson subdirectory graph for compile target composition.
