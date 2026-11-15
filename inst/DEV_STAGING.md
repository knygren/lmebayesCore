# Staged features (lmebayesCore -> glmbayesCore)

lmebayesCore holds features under development until they are ready for
glmbayesCore (CRAN). Staged code is written in glmbayesCore's final shape:
same file names (with an `_indexed` suffix where glmbayesCore already has a
file), same C++ namespaces, same `.two_block_*_cpp` wrapper names. Migrating a
feature is a copy plus the mechanical changes listed below.

Every staged file is either **migrating** or **dev only**. Dev-only symbols
carry a `dev_` / `.dev_` prefix (C++ namespace `lmebayes::dev`), so

    rg -l "dev_|lmebayes::dev" R src

lists exactly what must not be copied.

Status levels: `experimental` -> `parity-verified` -> `distribution-verified`
-> `ready`.

## Indexed two-block sweep (groups / population)

### Checks

1. **Parity (bitwise)**: `lmebayes/data-raw/check_port_steps.R`. Every R step
   and its C++ twin give identical output on the same inputs; the R and C++
   loops give identical state under `use_mode = TRUE`; a recorded R run
   replayed through both loops reproduces every logged input, the final
   state, and the prior sums exactly; name-based references confirm chain and
   level alignment. Ends with "All steps match."
2. **Distribution**: `lmebayes/data-raw/compare_sweep_indexed.R` (n = 1000).
   Welch comparison of `.rGLMM_reg_run_pipeline_cpp` and `_indexed` against
   the current `rGLMM_sweep` pipeline; no systematic `|z| > 2` pattern.

### In use in lmebayesCore

`rGLMM_reg_known_vcov` (non-Gaussian, every population component
`dNormal`; the `glmm_known` route under `rglmerb` / `glmerb`) runs its pilot
and main stages on `.rGLMM_sweep_cpp` and reports
`draw_engine = "rGLMM_sweep_cpp"`. ING (`rGLMM_reg_estimated_vcov`) and
Gaussian models stay on `rGLMM_sweep`. `.rGLMM_reg_run_pipeline` keeps
`rGLMM_sweep` as its default, so `compare_sweep_indexed.R` still compares
against the R sweep.

Both indexed loops return `sweep_mean` / `sweep_cov` (cross-chain summary of
the stacked population coefficients after each sweep); `.two_block_finish`
turns them into the same `sweep_history` object `rGLMM_sweep` stores.

### Migrating files

| lmebayesCore file | glmbayesCore destination | Status |
|---|---|---|
| `src/simfuncs_twoblock.h` | merge declarations into `src/simfuncs_twoblock.h` (`glmbayes::sim::twoblock`) | parity-verified |
| `src/two_block_index_spec.cpp` | `src/two_block_index_spec.cpp` | parity-verified |
| `src/two_block_groups_indexed.cpp` | `src/two_block_groups_indexed.cpp`, or merge into `two_block_groups.cpp` | parity-verified |
| `src/two_block_population_indexed.cpp` | `src/two_block_population_indexed.cpp`, or merge into `two_block_population.cpp` | parity-verified |
| `src/two_block_sweep.cpp` | `src/two_block_sweep.cpp` | parity-verified |
| `src/export_wrappers_twoblock.cpp` | two-block tier of `src/export_wrappers.cpp` | parity-verified |
| `src/rng_utils.h` (`glmbayes::rng::RSyncScope`) | `src/rng_utils.h` | parity-verified |
| `R/two_block_index_map.R` | `R/two_block_index_map.R` | parity-verified |
| `R/two_block_groups_indexed.R` | `R/two_block_groups_indexed.R` | parity-verified |
| `R/two_block_population_indexed.R` | `R/two_block_population_indexed.R` | parity-verified |
| `R/two_block_sweep_indexed.R` | `R/two_block_sweep_indexed.R` | parity-verified |
| two-block indexed block of `R/rcpp_wrappers_two_block.R` | `R/rcpp_wrappers_two_block.R` (`.rNormalGLMGroups_cpp` already exists in `R/rcpp_wrappers.R`; do not copy it) | parity-verified |

The C++ loop passed the distribution check before the staging refactor
(C++ vs current: 5/94 `|z| > 2`, max population `|z|` 0.94). Re-run
`compare_sweep_indexed.R` before promoting to `distribution-verified`.

### Staging-only shims (drop on migration)

| File | Why |
|---|---|
| `src/simfuncs_groups.h` | forwards `glmbayes::sim::group::rNormalGLMGroups` to lmebayesCore's `glmbayes::sim::rNormalGLMBlocks`; glmbayesCore has the real function |
| `.rNormalGLMGroups_cpp` in `R/rcpp_wrappers_two_block.R` | calls `_lmebayesCore_rNormalGLMBlocks_cpp_export`; glmbayesCore's wrapper calls `_glmbayesCore_rNormalGLMGroups_cpp_export` |

### Dev-only files (never migrate)

| File | Contents |
|---|---|
| `src/dev_twoblock.h`, `src/dev_two_block_replay.cpp` | `lmebayes::dev::twoblock`: `ReplayObserver` (tape replay, record log, prior sums), `log_new` |
| `src/dev_exports.cpp` | `dev_two_block_sweep_cpp_export`, `dev_two_block_log_new_cpp_export` |
| `R/dev_two_block_replay.R` | `.dev_two_block_log_new`, `.dev_two_block_observer`, `.dev_two_block_sweep_replay`, `.dev_two_block_prior_track`, `.dev_rGLMM_sweep_indexed_tracked`, `.dev_rGLMM_sweep_cpp_tracked` |
| `R/dev_rcpp_wrappers.R` | `.dev_two_block_sweep_cpp`, `.dev_two_block_log_new_cpp` |
| `.rGLMM_reg_run_pipeline_indexed` / `_cpp` in `R/rGLMM_reg.R` | development entry points; they use the tracked sweeps so results carry `prior_track` |

Replay, record, and tracking attach to the production loops through their
observer hooks (`SweepObserver` in C++, `observer =` in `.two_block_sweep_r`),
so the parity checks exercise the loops that ship. With no observer the
loops do no extra work.

### What migration changes

1. Copy the migrating `src/` files; merge `simfuncs_twoblock.h` declarations
   into glmbayesCore's header; add `RSyncScope` to its `rng_utils.h`.
2. Drop `src/simfuncs_groups.h`; include glmbayesCore's `simfuncs_groups.h`
   in `two_block_groups_indexed.cpp`.
3. Move the exports into glmbayesCore's `export_wrappers.cpp` two-block tier.
4. Copy the migrating R files and the two-block indexed wrappers; change
   `_lmebayesCore_` to `_glmbayesCore_` in the `.Call` symbols; skip
   `.rNormalGLMGroups_cpp`.
5. Do not copy any `dev_` file or symbol.
6. Route glmbayesCore's `rGLMM_sweep` through `.rGLMM_sweep_cpp` only as a
   separate reviewed change, keeping the exported signature. Gaussian
   families and non-`dNormal` population components stay on the current
   sweep.
7. Separately review glmbayesCore's `twoBlockGibbs.cpp` loops, which call R-RNG
   kernels with R callbacks inside one `.Call` and may need `RSyncScope`.
