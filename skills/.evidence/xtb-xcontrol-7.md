# Evidence: xtb-xcontrol-7

## Primary docs
- `man/xcontrol.7.adoc`

## Primary source entry points
- `skills/xtb-xcontrol-7/references/doc_map.md`
- `src/scc_core.f90`
- `src/model_hessian.f90`
- `src/grad_core.f90`
- `src/eeq_model.f90`
- `src/charge_model.f90`
- `include/xtb.h`
- `src/zmatpr.f90`
- `src/xhelp.f90`
- `src/wrmodef.f90`
- `src/wrmo.f90`
- `src/wrgbw.f90`
- `src/wrbas.f90`
- `src/vertical.f90`
- `src/topology.f90`
- `src/timing.f90`
- `src/thermo.f90`
- `src/surfac.f`
- `src/stm.f`
- `src/splitparam.f90`

## Extracted headings
- (none extracted)

## Executable command hints
- $fit
- $samerand
- $chrg 'int'
- $spin 'int'
- $cma
- $constrain
- $cube
- $embedding
- $external
- $fix
- $gbsa
- $gfn

## Warnings and pitfalls
- convergence thresholds for the `ancopt(3)`:
- average the energy and gradient before checking for convergence to accelerate
- NOTE: the scan parser will always terminate in error if the instruction could
- instructions with wrong input by raising a warning.
- damping for the Broyden convergence accelerator
