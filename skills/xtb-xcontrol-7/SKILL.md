---
name: xtb-xcontrol-7
description: Use this skill for xtb xcontrol instruction files and parser behavior, including $constrain/$fix/$scan/$wall/$scc/$opt/$md/$write groups, CLI integration with -I/--copy/--define/--strict, and source-level debugging of xcontrol behavior.
---

# xtb: Xcontrol 7

## High-Signal Playbook
### Route conditions
- Use this skill for `xcontrol(7)` instruction files and parser behavior (`man/xcontrol.7.adoc`).
- Route CLI-only method/run-type/solvent activation questions to `xtb-inputs-and-modeling`.
- Route toolchain/build/install questions to `xtb-build-and-install`.

### Triage questions
- Is the user changing constraints, scans, walls, MD/opt parameters, or output control groups?
- Should atoms be softly constrained (`$constrain`) or exactly fixed (`$fix`)?
- For scans, should mode be `sequential` or `concerted`?
- Is GBSA only set in `$gbsa` metadata, or also activated on CLI (`--gbsa`)?
- Are failures parser errors (syntax/indexing) or runtime convergence issues?
- Is this preflight (`--define`) or production (`--strict`)?

### Canonical workflow
1. Start from a minimal `xcontrol` file with only required groups.
2. Keep run type/method activation on CLI; use `-I FILE` only for grouped instructions.
3. Use `$constrain` for potential-based bias; use `$fix` for hard freezing.
4. Add `$scan` only after constraints parse cleanly.
5. Run `--define` first to expose warnings and parser issues.
6. Run with `--strict` for production so warnings become hard failures.
7. If behavior remains unclear, inspect parser and potential routines in source map.

### Minimal working example
```bash
cat > xtb.inp <<'XTBINPUT'
$chrg 0
$spin 0
$constrain
  force constant=0.5
  distance: 1,2,1.40
$opt
  optlevel=tight
$end
XTBINPUT

xtb -I xtb.inp assets/inputs/coord/caffeine.coord --opt tight --define
xtb -I xtb.inp assets/inputs/coord/caffeine.coord --opt tight --strict --namespace xc1
```

### Pitfalls and fixes
- Group flags must begin in column 1 with `$`; otherwise parser routing fails.
- `$scan` is strict and aborts on malformed instructions; fix syntax before tuning parameters.
- `$constrain` may warn/skip entries; warnings do not imply constraints were applied.
- `$fix` can be deactivated in unstable MD paths; use MD-appropriate constraints.
- `$gbsa` group sets metadata only; activate solvent on CLI as needed.
- Legacy set-blocks remain for compatibility but are slated for deprecation.
- `--copy` rewrites settings and strips comments; keep authoritative files in version control.

### Validation checkpoints
- `--define` run prints no parser errors for the selected groups.
- `--strict` production run exits cleanly for validated inputs.
- For scan workflows, `xtbscan.log` (or namespaced equivalent) is generated.
- Constraint summary in output reflects requested atom indices/targets.
- Warning block is empty for production settings.

### Source drill-down (when docs stop)
- Main xcontrol parser and dispatch: `src/set_module.f90`
- User-data parser for constraints/scans/walls: `src/constrain_param.f90`
- Constraint potential implementation: `src/constrain_pot.f90`
- Scan setup/runtime: `src/scanparam.f90`, `src/scan_driver.f90`
- Namespace naming for scan logs: `src/setparam.f90`
- CLI integration for `-I`, `--copy`, `--define`, `--strict`: `src/prog/main.F90`

## Scope
- Handle questions about `xcontrol(7)` instruction syntax, parser behavior, and xcontrol-driven runtime controls.
- Keep responses architectural first; go function-level only when required to resolve ambiguity.

## Primary documentation references
- `man/xcontrol.7.adoc`
- `man/xtb.1.adoc`

## Workflow
- Start with the primary references above.
- If details are missing, inspect `references/doc_map.md` for additional runnable anchors.
- If ambiguity remains after docs, inspect `references/source_map.md` and start with parser-oriented entry points.
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
- `src/set_module.f90` (`rdcontrol`, `rdblock`, `set_define`, `set_fix`, `set_constr`, `set_scan`, `set_wall`, `set_gbsa`)
- `src/constrain_param.f90` (`read_userdata`, `rdblock`, `set_fix`, `set_constr`, `set_scan`, `set_wall`)
- `src/constrain_pot.f90` (`constrain_dist`, `constrain_angle`, `constrain_dihedral`, `constrain_pot`, `constrain_hess`)
- `src/scanparam.f90` (`init_scan`, `setup_constrain_pot`, `pot_info`)
- `src/scan_driver.f90` (`relaxed_scan`)
- `src/setparam.f90` (`get_namespace`)
- `src/prog/main.F90` (`parseArguments`)
- `src/readin.f90` (`get*Value` helpers used by parser paths)
- `src/type/environment.f90` (`initEnvironment`, strict warning handling)
- Prefer targeted source search: `rg -n "<symbol_or_keyword>" include src test`.
