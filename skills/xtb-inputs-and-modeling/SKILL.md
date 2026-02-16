---
name: xtb-inputs-and-modeling
description: "Use this skill for xtb simulation setup from input files and CLI flags: geometry formats, run-type/method selection, SCC controls, charge/spin handling, solvation models, and modeling behavior validation."
---

# xtb: Inputs and Modeling

## High-Signal Playbook
### Route conditions
- Use this skill for `xtb(1)` geometry input formats, run-type selection, SCC controls, charge/spin handling, solvent flags, and output artifact expectations (`man/xtb.1.adoc`).
- Route build/toolchain/install questions to `xtb-build-and-install`.
- Route detailed `xcontrol(7)` block syntax (`$constrain`, `$scan`, `$wall`, etc.) to `xtb-xcontrol-7`.

### Triage questions
- Which geometry format is used (`coord`, `xyz`, `sdf/mol`, `pdb`, `vasp`, `json`)?
- Which run type is intended (`--sp`, `--opt`, `--hess`, `--md`)?
- Which method is required (`--gfn INT`, `--gfnff`, `--tblite`, `--ptb`)?
- What are charge/spin, and are `.CHRG`/`.UHF` files present?
- Is implicit solvation needed (`--alpb`, `--gbsa`, `--cosmo`, `--cpcmx`)?
- Is preflight checking needed (`--define`) before production (`--strict`)?

### Canonical workflow
1. Validate geometry readability with `xtb info` or a quick single-point.
2. Set `XTBPATH` explicitly to the parameter-file root.
3. Pick one run type (only the first run type is honored).
4. Set method (`--gfn`/`--gfnff`) and SCC controls (`--acc`, `--iterations`, `--etemp`).
5. Set charge/spin explicitly (`--chrg`, `--uhf`) when reproducibility matters.
6. Add solvent flags if needed.
7. Run `--define` first when checking complex setups.
8. Run production with `--strict` and inspect warning summary/exit code.

### Minimal working examples
```bash
# Input/readability preflight
export XTBPATH="$PWD"
xtb info assets/inputs/coord/caffeine.coord
```

```bash
# Single-point in implicit solvent
xtb assets/inputs/coord/caffeine.coord \
  --gfn 2 --sp --gbsa h2o \
  --acc 1.0 --iterations 250 --etemp 300 \
  --namespace sp1 > sp1.out
```

```bash
# Geometry optimization with explicit preflight and strict production
xtb assets/inputs/coord/caffeine.coord --opt tight --cycles 200 --define
xtb assets/inputs/coord/caffeine.coord --opt tight --cycles 200 --strict --namespace opt1 > opt1.out
```

### Pitfalls and fixes
- Multiple run types in one command: only the first is used; split into separate runs.
- `XTBPATH` unset/misconfigured: parameter lookup falls back through deprecated paths; set `XTBPATH` explicitly.
- `.CHRG`/`.UHF` surprises: CLI values override file/xcontrol values.
- Solvent naming mistakes: use documented solvent names/states from `man/xtb.1.adoc`.
- Point-charge embedding + parallel Hessian can fail due to I/O; avoid that combination.
- Hessian on poorly optimized structures can trigger reliability warnings.
- Non-zero exit (`128`) means failure and requires triage before trusting outputs.

### Validation checkpoints
- Preflight: `xtb info <geometry>` succeeds for selected format.
- Production runs return zero and include an energy summary in stdout.
- Optimization runs produce expected artifacts (`xtbopt.xyz`/`xtbopt.log`, namespaced when `--namespace` is used).
- Warning summary is clean under `--strict` for production jobs.
- Solvation sanity: compare gas-phase vs solvent energies before scaling to long campaigns.

### Source drill-down (when docs stop)
- CLI parsing and option precedence: `src/prog/main.F90`, `src/prog/argparser.f90`
- Runtime defaults and option setters: `src/set_module.f90`, `src/setparam.f90`
- Environment/parameter resolution: `src/type/environment.f90`, `src/readin.f90`
- Geometry analysis and checks: `src/main/geometry.f90`
- SCC/charge internals: `src/scc_core.f90`, `src/charge_model.f90`, `src/eeq_model.f90`
- Solvation parsing/modeling: `src/solv/input.F90`, `src/solv/model.f90`
- Hessian/model behavior: `src/model_hessian.f90`

## Scope
- Handle questions about inputs, system setup, models, and physical parameterization.
- Keep responses architectural first; use function-level detail only for unresolved behavior.

## Primary documentation references
- `man/xtb.1.adoc`

## Workflow
- Start with the primary references above.
- If details are missing, inspect `references/doc_map.md` for additional docs and runnable inputs.
- Use tests as behavior references when modeling behavior is ambiguous.
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
- `src/prog/argparser.f90` (`nextFlag`, `nextArg`, `findArg`)
- `src/set_module.f90` (`set_chrg`, `set_spin`, `set_gfn`, `set_scc`, `set_opt`, `set_md`, `set_hess`, `set_define`)
- `src/setparam.f90` (`get_namespace`)
- `src/type/environment.f90` (`initEnvironment`)
- `src/readin.f90` (`xfind`, `getIntValue`, `getRealValue`, `getBoolValue`)
- `src/main/geometry.f90` (`main_geometry`, `calc_distances`, `check_cold_fusion`)
- `src/scc_core.f90` (`scc`, `fermismear`, `occ`)
- `src/charge_model.f90` (`new_charge_model_2019`, `gfn0_charge_model`)
- `src/eeq_model.f90` (`eeq_chrgeq_core`, `eeq_chrgeq_gbsa`)
- `src/solv/input.F90` (`pgiWrapper`)
- `src/solv/model.f90` (`initSolvModel`, `normalizeSolventName`, `getParamFile`)
- `src/model_hessian.f90` (`mh_lindh`, `mh_eeq`)
- Prefer targeted source search: `rg -n "<symbol_or_keyword>" include src test`.
