# xtb source map: Xcontrol 7

Generated from source roots:
- `include`
- `src`
- `test`

Use this map only after exhausting topic docs in `references/doc_map.md`.

## Topic query tokens
- `xcontrol`
- `constrain`
- `fix`
- `scan`
- `wall`
- `define`
- `strict`
- `copy`

## Fast source navigation
- `rg -n "rdcontrol|rdsetbl|rdblock|set_(define|fix|constr|scan|wall|gbsa)" src/set_module.f90`
- `rg -n "read_userdata|set_(fix|constr|scan|wall)|warning\(|error\(" src/constrain_param.f90`
- `rg -n "constrain_(dist|angle|dihedral|pot|hess)" src/constrain_pot.f90`
- `rg -n "init_scan|setup_constrain_pot|pot_info|relaxed_scan|xtbscan\.log" src/scanparam.f90 src/scan_driver.f90`
- `rg -n -e "--input|--copy|--define|--strict|parseArguments" src/prog/main.F90`

## Suggested source entry points
- `src/set_module.f90` | top-level xcontrol parser and group dispatch (`rdcontrol`, `rdblock`, `set_*`).
- `src/constrain_param.f90` | constraint/scan/wall parser semantics with warning/error behavior.
- `src/constrain_pot.f90` | force/energy/hessian implementation for constrained coordinates.
- `src/scanparam.f90` | scan state setup and constraint-potential initialization.
- `src/scan_driver.f90` | scan execution loop and `xtbscan.log` write path.
- `src/setparam.f90` | `get_namespace`; namespaced scan/log artifact naming.
- `src/prog/main.F90` | CLI options controlling xcontrol ingestion (`-I`, `--copy`, `--define`, `--strict`).
- `src/readin.f90` | typed value readers (`getIntValue`, `getRealValue`, `getBoolValue`, `getListValue`) used by parser routines.
- `src/type/environment.f90` | strict-mode warning-to-error behavior used during xcontrol validation.
- `test/unit/CMakeLists.txt` | strict CLI behavior baseline and environment setup pattern for checks.
