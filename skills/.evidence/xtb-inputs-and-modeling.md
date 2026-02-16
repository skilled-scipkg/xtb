# Evidence: xtb-inputs-and-modeling

## Primary docs
- `man/xtb.1.adoc`

## Primary source entry points
- `skills/xtb-inputs-and-modeling/references/doc_map.md`
- `src/main/geometry.f90`
- `src/model_hessian.f90`
- `src/eeq_model.f90`
- `src/charge_model.f90`
- `src/type/dispersion_model.f90`
- `src/solv/model.f90`
- `src/solv/input.F90`
- `src/mctc/boundaryconditions.f90`
- `src/ptb/corebasis.F90`
- `src/xtb/basisset.f90`
- `src/type/basisset.f90`

## Extracted headings
- (none extracted)

## Executable command hints
- (none extracted)

## Warnings and pitfalls
- `xtb(1)` will issue a warning if `XTBHOME` is not part of the `XTBPATH`
- `xtb(1)` can generate the two types of warnings, the first warning section
- summary of each warning with its respective string representation in the
- This warning will be issued twice, once before the Hessian,
- this this warning could be detected) and in the warning block
- in the end. The warning will be generated if the gradient norm
- Failure (termination via error stop generates 128 as return value)
