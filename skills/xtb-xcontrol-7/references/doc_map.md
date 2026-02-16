# xtb documentation map: Xcontrol 7

Practical references for xcontrol instruction syntax, parser interaction, and simulation preflight.

## Primary manuals
- `man/xcontrol.7.adoc` | xcontrol group syntax and parser behavior (`$constrain`, `$fix`, `$scan`, `$wall`, `$scc`, `$opt`, `$md`, `$write`).
- `man/xtb.1.adoc` | CLI integration for `-I/--input`, `--copy`, `--define`, `--strict`, and warning semantics.

## Runnable artifacts
- `assets/inputs/coord/caffeine.coord` | baseline geometry for xcontrol preflight and strict runs.

## Validation anchors
- `test/unit/CMakeLists.txt` | strict CLI test style and `XTBPATH` test environment conventions.
