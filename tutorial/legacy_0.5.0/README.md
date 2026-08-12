# Legacy material — OtsuFire 0.5.x

**These scripts do not run against OtsuFire 2.0.0.** They are kept as a record
of how the 0.5.x experiments were produced, not as usable documentation.

For the current workflow see [`../OtsuFire_2.0.0/`](../OtsuFire_2.0.0/).

## Why they no longer run

The 2.0.0 rebuild removed 17 of the 19 functions exported by 0.1.x and
reorganised the rest behind config objects. Two changes in particular break
every script in this folder:

* The supervised stage is driven by `build_supervised_burned_config()` +
  `run_oneyear_supervised_pipeline()`. The older entry points these scripts
  call are gone or are now internal.
* The `scenario` argument of the supervised config was renamed to `run_label`
  in 0.11.0. Anything here that passes `scenario =` fails immediately.

Several of them also describe package behaviour that has since changed —
negative-pool policies that no longer exist, a registry that was removed on
2026-06-05, and a feature recipe that was unified across the OOF, final and
scoring stages.

## What is in here

| File | What it was |
|---|---|
| `SUPERVISED_ONE_YEAR_Otsu_0.5.0.R` | One-year supervised baseline run under 0.5.0. |
| `TUTORIAL_SUPERVISED_2005_balanced_VERSION_B_pragmatic.R` | Step-by-step walkthrough of the 2005 balanced combo, calling internals directly. |
| `TUTORIAL_PASO_A_PASO/` | The same walkthrough split into blocks, starting from a full setup block. |
| `Section_1.*.R`, `Section_2.*.R` | Analyses written for specific paper sections. |
| `experimento_1.R` | Experiment E1: the model without the 13 hotspot features. |
| `copia_seguridad_*.R` | Snapshot scripts that froze the E1/E2/E4 and baseline outputs. |

The comments were translated from Spanish to English in 2026-08 along with the
rest of the repository; the code itself was not modified, so these files still
reflect exactly what was run at the time.
