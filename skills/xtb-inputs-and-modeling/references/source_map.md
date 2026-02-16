# xtb source map: Inputs and Modeling

Generated from source roots:
- `include`
- `src`
- `test`

Use this map only after exhausting topic docs in `references/doc_map.md`.

## Topic query tokens
- `input`
- `geometry`
- `run type`
- `gfn`
- `scc`
- `charge`
- `spin`
- `solvent`
- `hessian`
- `xtbpath`

## Fast source navigation
- `rg -n "parseArguments|--sp|--opt|--hess|--md|--gfn|--gbsa|--strict|--define" src/prog/main.F90`
- `rg -n "set_(chrg|spin|gfn|scc|opt|md|hess|define)|rdcontrol" src/set_module.f90`
- `rg -n "initEnvironment|XTBPATH|XTBHOME" src/type/environment.f90`
- `rg -n "xfind|get(Int|Real|Bool|List)Value" src/readin.f90`
- `rg -n "initSolvModel|normalizeSolventName|getParamFile|pgiWrapper" src/solv/model.f90 src/solv/input.F90`
- `rg -n "scc\(|eeq_chrgeq_|new_charge_model_|mh_lindh|mh_eeq" src/scc_core.f90 src/eeq_model.f90 src/charge_model.f90 src/model_hessian.f90`

## Suggested source entry points
- `src/prog/main.F90` | `xtbMain`, `parseArguments`; run-type and flag precedence.
- `src/prog/argparser.f90` | argument iteration and extraction (`nextFlag`, `nextArg`, `findArg`).
- `src/set_module.f90` | setters for CLI/xcontrol-facing controls (`set_chrg`, `set_spin`, `set_gfn`, `set_scc`, `set_opt`, `set_md`, `set_hess`).
- `src/setparam.f90` | global defaults and namespace handling via `get_namespace`.
- `src/type/environment.f90` | `initEnvironment`; `HOME`, `XTBHOME`, `XTBPATH` resolution and strict behavior.
- `src/readin.f90` | parameter and list parsing (`xfind`, `get*Value`).
- `src/main/geometry.f90` | geometry summaries and derived metrics (`calc_distances`, `calc_angles`, `check_cold_fusion`).
- `src/scc_core.f90` | SCC solver core (`scc`, `fermismear`, `occ`).
- `src/charge_model.f90` | charge-model initialization routines.
- `src/eeq_model.f90` | EEQ charge-equilibration core and GBSA variants.
- `src/solv/input.F90` | solvent input object construction (`pgiWrapper`).
- `src/solv/model.f90` | solvation model construction and solvent-name normalization.
- `src/model_hessian.f90` | model Hessian generation paths used by `--hess` workflows.
- `test/unit/test_gfn2.f90` | function-level checks around GFN2 and solvent-related behavior.
- `test/unit/test_gfnff.f90` | force-field model behavior checks.
- `test/unit/test_eeq.f90` | EEQ and charge/solvation-related regression checks.
- `test/unit/CMakeLists.txt` | executable smoke commands used as baseline behavior checks.
