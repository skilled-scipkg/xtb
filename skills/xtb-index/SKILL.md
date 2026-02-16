---
name: xtb-index
description: Use this skill to route xtb requests to the correct topic skill before deep inspection. It is the start point for build/install, simulation input/modeling, and xcontrol instruction workflows.
---

# xtb Skills Index

## Start here
- Work from repository root so all relative paths resolve.
- Run a quick preflight before routing:

```bash
test -f assets/inputs/coord/caffeine.coord
test -f man/xtb.1.adoc
(xtb --version || ./_build/xtb --version || ./build/xtb --version)
```

- If no runnable binary is available, start with `xtb-build-and-install`.

## Route the request
- Build/install/toolchain/test/install-prefix questions: `xtb-build-and-install`
- CLI input formats, run types, method flags, charge/spin, solvent, and model behavior: `xtb-inputs-and-modeling`
- `xcontrol(7)` blocks and parser behavior (`$constrain`, `$fix`, `$scan`, `$wall`, `$scc`, `$opt`, `$md`, `$write`): `xtb-xcontrol-7`

## Escalate in this order
1. Start from the target skill `SKILL.md` playbook and primary docs.
2. Use the target skill `references/doc_map.md` for additional docs and runnable examples.
3. Use the target skill `references/source_map.md` for function-level source checks.
4. Run targeted symbol search only after docs and maps: `rg -n "<symbol_or_keyword>" include src test`.

## Source directories for deep inspection
- `include`
- `src`
- `test`
