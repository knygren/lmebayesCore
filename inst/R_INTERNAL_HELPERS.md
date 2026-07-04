# `R/` — internal helpers (undocumented in `man/`)

Functions and symbols in **`R/`** with **`@noRd`**, **`@keywords internal`**, or
no roxygen, without a dedicated help page. Intended for **`glmbayesCore:::`**
or in-package use — not part of the exported API unless promoted.

**Columns:** *File* is the defining source; *Called from* lists direct callers
in `R/` (comma-separated). Helpers with no callers are marked *(unused)*.

Companion: [R_EXPORTED_AND_DOCUMENTED.md](R_EXPORTED_AND_DOCUMENTED.md).

**lmebayes** resolves these via `getFromNamespace` / `:::`: `.two_block_as_staged_names`, `.two_block_tau2_ref_from_pfamily`, `.mrglmb_normalize_pfamily_lists`, `.validate_pfamily_for_rlmb`, `.lmebayes_priors_from_pfamily_list`, `.lmebayes_block2_icm_labels`, `.lmebayes_mer_optional_args`, `extract_mer_variance_components`.

---

## Mixed-model matrix samplers (`mixed_rmerb_helpers.R`, `rlmerb.R`, `rglmerb.R`)

| Function | File | Role | Called from |
|----------|------|------|-------------|
| `.lmebayes_resolve_dispersion_ranef()` | `mixed_rmerb_helpers.R` | Map `dispersion_ranef` / legacy `dispersion` to per-RE weights | `rlmerb()`, `rglmerb()`, `.lmebayes_priors_from_pfamily_list()` |
| `.lmebayes_validate_dispersion_ranef()` | `mixed_rmerb_helpers.R` | Validate dispersion-ranef vector | *(unused)* |
| `.lmebayes_priors_from_pfamily_list()` | `mixed_rmerb_helpers.R` | Normalize `pfamily_list` → sampler `prior` list | `rlmerb()`, `rglmerb()`; **lmebayes** `lmerb()`, `glmerb()` |
| `.lmebayes_run_lmm_engine()` | `mixed_rmerb_helpers.R` | Route Gaussian GLMM to `rLMMNormal_reg` / ING path | `rglmerb()` (Gaussian), `rlmerb()` |
| `.lmebayes_block1_prior_list()` | `mixed_rmerb_helpers.R` | Assemble Block~1 `prior_list` from design + pfamily | `rlmerb()`, `rglmerb()` |
| `.lmebayes_add_fixef_summaries()` | `mixed_rmerb_helpers.R` | Attach fixef summary slots to sampler output | `rlmerb()`, `rglmerb()` |
| `.lmebayes_block2_icm_labels()` | `mixed_rmerb_helpers.R` | ICM column labels for Block~2 fixef table | `rlmerb()`, `rglmerb()`; **lmebayes** `lmerb()`, `glmerb()` |
| `.lmebayes_print_icm_fixef_table()` | `mixed_rmerb_helpers.R` | Console ICM fixef table | `rlmerb()`, `rglmerb()` |
| `.lmebayes_print_ranef_mode_reference()` | `mixed_rmerb_helpers.R` | Ranef mode reference line (GLMM) | `rglmerb()` |
| `.lmebayes_print_fixef_init()` | `mixed_rmerb_helpers.R` | Fixef init banner | `rlmerb()`, `rglmerb()` |

---

## Model setup (`model_setup.R`)

| Function | File | Role | Called from |
|----------|------|------|-------------|
| `.lmebayes_normalize_family()` | `model_setup.R` | Coerce `family` argument | `model_setup()` |
| `.lmebayes_mer_convergence_issues()` | `model_setup.R` | Collect lme4 convergence warnings | `Prior_Setup_lmebayes()` |
| `.lmebayes_mer_optional_args()` | `model_setup.R` | Optional `lmer` / `glmer` call args | `model_setup()`, **lmebayes** `glmerb()` |

---

## Mixed-model prior setup (`Prior_Setup_lmebayes.R`)

| Function | File | Role | Called from |
|----------|------|------|-------------|
| `.lmebayes_resolve_pwt()` | `Prior_Setup_lmebayes.R` | Resolve `pwt` / per-RE weights | `Prior_Setup_lmebayes()` |
| `.lmebayes_resolve_disp_prior()` | `Prior_Setup_lmebayes.R` | Resolve dispersion hyperprior fields | `Prior_Setup_lmebayes()` |
| `.lmebayes_block_glm_estimable()` | `Prior_Setup_lmebayes.R` | Block-GLM estimability check for binomial | `Prior_Setup_lmebayes()` |

---

## lme4 design chain (`lme4_design_utilities.R`)

| Function | File | Role | Called from |
|----------|------|------|-------------|
| `.lmebayes_validate_uncorrelated_re_formula()` | `lme4_design_utilities.R` | Require uncorrelated RE formula | `extract_re_hyper_matrices()` |
| `.lme4_Z_random_column_map()` | `lme4_design_utilities.R` | Map Z columns to RE terms | `.lme4_Z_random_colnames()`, `get_lme4_components()` |
| `.lme4_Z_random_colnames()` | `lme4_design_utilities.R` | Column names for sparse Z | `.lme4_label_Z_random_sparse()` |
| `.lme4_Z_random_rownames()` | `lme4_design_utilities.R` | Row names for sparse Z | `.lme4_Z_random_row_map()`, `.lme4_label_Z_random_sparse()` |
| `.lme4_Z_random_row_map()` | `lme4_design_utilities.R` | Map Z rows to obs × group | `get_lme4_components()` |
| `.lme4_label_Z_random_sparse()` | `lme4_design_utilities.R` | Dimnames on sparse Z | `get_lme4_components()` |

Exported entry points that reach the Z-label chain: `model_setup()` →
`extract_re_hyper_matrices()` → `get_lme4_components()`.

---

## Two-block rate / TV (`two_block_rate.R`, `two_block_tv_bound.R`)

| Function | File | Role | Called from |
|----------|------|------|-------------|
| `.two_block_rate_inputs` | `two_block_rate.R` | — | *(unused)* |
| `.two_block_S_P11` | `two_block_rate.R` | — | *(unused)* |
| `.two_block_gen_eigen` | `two_block_rate.R` | — | *(unused)* |
| `.two_block_erfn` | `two_block_tv_bound.R` | — | *(unused)* |
| `.two_block_tv_bound_one` | `two_block_tv_bound.R` | — | *(unused)* |

---

## Two-block pilot / GLMM (`two_block_glmm_pilot_helpers.R`, `two_block_pilot_cost.R`)

| Function | File | Role | Called from |
|----------|------|------|-------------|
| `.two_block_print_pilot_stage_diagnostics` | `two_block_glmm_pilot_helpers.R` | Print pilot-stage diagnostics between pilot and main sampling (UB path) | rGLMM.R, two_block_lmm_staged_sweep_outer.R |
| `.two_block_as_staged_names` | `two_block_glmm_pilot_helpers.R` | — | rGLMM.R, rLMMNormal_reg.R, two_block_lmm_staged_sweep_outer.R |
| `.two_block_pilot_chisq_test` | `two_block_glmm_pilot_helpers.R` | Hotelling chi-squared test: pilot fixef mean vs start | rGLMM.R, two_block_lmm_staged_sweep_outer.R |
| `.two_block_fixef_colmeans` | `two_block_glmm_pilot_helpers.R` | Column means of fixef draws with names copied from fixef mode | rGLMM.R, two_block_lmm_staged_sweep_outer.R |
| `.two_block_pilot_eigenvalue_ub` | `two_block_glmm_pilot_helpers.R` | Post-pilot eigenvalue upper bounds (per-draw rate maxima) | *(unused)* |
| `.two_block_pilot_ub_from_coefficients` | `two_block_glmm_pilot_helpers.R` | — | rGLMM.R, two_block_lmm_staged_sweep_outer.R |
| `.two_block_resolve_n_pilot` | `two_block_measurement_prior.R` | — | *(unused)* |
| `.two_block_pilot_will_run` | `two_block_pilot_cost.R` | still run (Gaussian LMM with ING Block~2 components). | rGLMM.R, two_block_lmm_staged_sweep_outer.R |
| `.two_block_resolve_pilot_plan` | `two_block_pilot_cost.R` | — | rGLMM.R, two_block_lmm_staged_sweep_outer.R |
| `.two_block_print_pilot_plan` | `two_block_pilot_cost.R` | Print resolved pilot / main sampling plan (before pilot stage) | rGLMM.R, two_block_lmm_staged_sweep_outer.R |
| `.two_block_print_pilot_cost_opt` | `two_block_pilot_cost.R` | Print pilot cost optimization advisory (before pilot stage) | *(unused)* |

---

## Two-block sweep history (`two_block_sweep_history.R`)

| Function | File | Role | Called from |
|----------|------|------|-------------|
| `.two_block_print_sweep_boundary` | `two_block_batch_gibbs.R` | Print sweep/block enter or exit line with wall-clock timestamp | rGLMM_sweep.R |
| `.two_block_print_sweep_early_diagnostics` | `two_block_batch_gibbs.R` | — | rGLMM_sweep.R |
| `.two_block_build_sweep_history` | `two_block_sweep_history.R` | Build structured two-block sweep history from per-sweep fixef snapshots | rGLMM_sweep.R |
| `.two_block_filter_sweep_history_table` | `two_block_sweep_history.R` | Filter sweep-history table rows for printing | *(unused)* |
| `.two_block_sweep_history_header_n` | `two_block_sweep_history.R` | Sweep count shown in the print header | *(unused)* |
| `.two_block_print_sweep_history_body` | `two_block_sweep_history.R` | Print one Block~2 sweep-history table (mode + optional sweep rows) | *(unused)* |
| `.two_block_print_sweep_history_tables` | `two_block_sweep_history.R` | Print stage-end table via structured sweep history (legacy helper) | *(unused)* |

---

## Two-block measurement / tau2 (`two_block_measurement_prior.R`, `two_block_tau2_ref.R`)

| Function | File | Role | Called from |
|----------|------|------|-------------|
| `.two_block_tau2_start_from_pfamily` | `two_block_batch_gibbs.R` | Starting tau2 vector from pfamily prior fields (plug-in dispersions) | rGLMM.R, rGLMM_sweep.R, two_block_lmm_staged_sweep_outer.R |
| `.two_block_measurement_prior_list` | `two_block_measurement_prior.R` | — | *(unused)* |
| `.two_block_icm_at_start` | `two_block_measurement_prior.R` | — | rGLMM.R, rLMMNormal_reg.R |
| `.two_block_validate_gap_tol` | `two_block_measurement_prior.R` | Validate gap tolerance for pilot chain count derivation | rGLMM.R, two_block_lmm_staged_sweep_outer.R, two_block_pilot_cost.R |
| `.two_block_tau2_ref_from_pfamily` | `two_block_tau2_ref.R` | — | two_block_batch_gibbs.R, two_block_measurement_prior.R |
| `.two_block_tau2_ref_vector` | `two_block_tau2_ref.R` | Named plug-in tau^2 vector from a validated pfamily_list | *(unused)* |
| `.two_block_tau2_start_from_dispersion_draws` | `two_block_tau2_ref.R` | — | rGLMM.R, two_block_lmm_staged_sweep_outer.R |

---

## Two-block drivers (`two_block_rNormal_reg*.R`)

| Function | File | Role | Called from |
|----------|------|------|-------------|
| `.two_block_block1_prior_with_tau2_cpp` | `rcpp_wrappers.R` | — | two_block_batch_gibbs.R |
| `.two_block_block1_iters_mean_cpp` | `rcpp_wrappers.R` | — | two_block_batch_gibbs.R |
| `.two_block_block1_one_chain_draw_cpp` | `rcpp_wrappers.R` | — | *(unused)* |
| `.two_block_block1_one_chain_cpp` | `rcpp_wrappers.R` | — | *(unused)* |
| `.two_block_block1_one_chain_from_mu_P_cpp` | `rcpp_wrappers.R` | — | two_block_batch_gibbs.R |
| `.two_block_block1_one_chain_v2_cpp` | `rcpp_wrappers.R` | — | *(unused)* |
| `.two_block_block1_all_chains_v2_internal_cpp` | `rcpp_wrappers.R` | — | *(unused)* |
| `.two_block_block1_all_chains_v2_internal_loop_cpp` | `rcpp_wrappers.R` | — | *(unused)* |
| `.two_block_block1_all_chains_cpp` | `rcpp_wrappers.R` | — | two_block_batch_gibbs.R |
| `.two_block_block1_prior_with_tau2_r` | `two_block_batch_gibbs.R` | Refresh Block 1 prior precision for ING components (R reference; v5 any_ing) | *(unused)* |
| `.two_block_block1_prior_with_tau2` | `two_block_batch_gibbs.R` | Refresh Block 1 prior precision for ING components (mirrors C++ twoBlockGibbs) | *(unused)* |
| `.two_block_block1_prep_one_chain` | `two_block_batch_gibbs.R` | One-chain Block 1 prep: fixef -> mu_all -> prior_list (no sampling) | *(unused)* |
| `.two_block_block1_prep_all_chains` | `two_block_batch_gibbs.R` | All-chain Block 1 prep: mu_all and prior_list for every chain | *(unused)* |
| `.two_block_block1_iters_mean_r` | `two_block_batch_gibbs.R` | Mean envelope candidates per group from a Block~1 draw (R reference) | *(unused)* |
| `.two_block_block1_iters_mean` | `two_block_batch_gibbs.R` | Mean envelope candidates per group from a Block~1 draw | *(unused)* |
| `.two_block_block1_reorder_b_r` | `two_block_batch_gibbs.R` | — | *(unused)* |
| `.two_block_block1_reorder_b` | `two_block_batch_gibbs.R` | — | *(unused)* |
| `.two_block_block1_draw_block` | `two_block_batch_gibbs.R` | — | *(unused)* |
| `.two_block_block1_draw_reorder` | `two_block_batch_gibbs.R` | Reorder Block~1 draw coefficients and summarize envelope iterations | *(unused)* |
| `.two_block_block1_draw_one_chain` | `two_block_batch_gibbs.R` | One-chain Block 1 draw given a prepared prior_list | *(unused)* |
| `.two_block_block1_one_chain` | `two_block_batch_gibbs.R` | — | *(unused)* |
| `.two_block_block1_one_chain_from_batch` | `two_block_batch_gibbs.R` | — | *(unused)* |
| `.two_block_block1_draw_all_chains` | `two_block_batch_gibbs.R` | All-chain Block 1 draw from prepared prior_lists | *(unused)* |
| `.two_block_block1_one_chain_r` | `two_block_batch_gibbs.R` | Block 1 one-chain prep + draw (R reference; all piece flags FALSE) | *(unused)* |
| `.two_block_block1_glmbfamfunc` | `two_block_batch_gibbs.R` | — | *(unused)* |
| `.two_block_block1_all_chains` | `two_block_batch_gibbs.R` | — | rGLMM_sweep.R |
| `.two_block_block1_all_chains_v2` | `two_block_batch_gibbs.R` | — | *(unused)* |
| `.two_block_block1_all_chains_via_cpp` | `two_block_batch_gibbs.R` | — | *(unused)* |
| `.two_block_normalize_family` | `two_block_rNormal_reg.R` | — | rGLMM.R, two_block_measurement_prior.R, two_block_mode_weights.R, two_block_rate.R, two_block_rNormal_reg.R |
| `.two_block_validate_block1_prior` | `two_block_rNormal_reg.R` | — | rGLMM.R, rLMMNormal_reg.R, two_block_lmm_staged_sweep_outer.R, two_block_rNormal_reg.R |
| `.two_block_block1_prior_list` | `two_block_rNormal_reg.R` | — | *(unused)* |
| `.two_block_mu_all` | `two_block_rNormal_reg.R` | — | *(unused)* |
| `.two_block_format_cpp_out` | `two_block_rNormal_reg.R` | Format raw C++ output into a `two_block_rNormal_reg` object | *(internal)* |
| `.two_block_validate_pfamily_list` | `two_block_rNormal_reg.R` | — | rGLMM.R, rLMMNormal_reg.R |

---

## Two-block batch Gibbs (`two_block_batch_gibbs.R`)

| Function | File | Role | Called from |
|----------|------|------|-------------|
| `.two_block_progbar_prefix` | `two_block_batch_gibbs.R` | Build a progress-bar prefix for sweep-outer Gibbs stages | rGLMM_sweep.R |
| `.two_block_progress_bar` | `two_block_batch_gibbs.R` | Text progress bar matching glmbayesCore C++ style | rGLMM_sweep.R |
| `.two_block_progress_bar_finish` | `two_block_batch_gibbs.R` | — | rGLMM_sweep.R |
| `.two_block_print_block1_phase` | `two_block_batch_gibbs.R` | Print Block 1 prep/draw sub-phase boundary with wall-clock timestamp | *(unused)* |
| `.two_block_fixef_mode_at` | `two_block_batch_gibbs.R` | Lookup ICM mode for one fixef (re_name, covariate) pair | two_block_sweep_history.R |
| `.two_block_rglmb_iter_count` | `two_block_batch_gibbs.R` | — | *(unused)* |
| `.two_block_snapshot_fixef_stats` | `two_block_batch_gibbs.R` | Snapshot chain colMeans and SDs of Block 2 fixef after one sweep | rGLMM_sweep.R |
| `.two_block_sweep_ranef_chain_mean` | `two_block_batch_gibbs.R` | — | *(unused)* |
| `.two_block_sweep_mu_all_chain_mean` | `two_block_batch_gibbs.R` | — | *(unused)* |
| `.two_block_print_block_diag` | `two_block_batch_gibbs.R` | Print per-sweep block diagnostics (fixef table across chains; b vs mode optional) | *(unused)* |
| `.rGLMM_sweep_initialize` | `two_block_batch_gibbs.R` | — | rGLMM_sweep.R |
| `.two_block_ensure_batch_b_dimnames` | `two_block_batch_gibbs.R` | — | *(unused)* |
| `.two_block_batch_b_array_to_master_matrix` | `two_block_batch_gibbs.R` | — | *(unused)* |
| `.two_block_batch_b_master_matrix_to_array` | `two_block_batch_gibbs.R` | — | *(unused)* |
| `.two_block_batch_fixef_chain` | `two_block_batch_gibbs.R` | Extract chain-i fixef list from batch state | *(unused)* |
| `.two_block_lapply_chains` | `two_block_batch_gibbs.R` | Apply FUN to each chain index, optionally in parallel (Unix/macOS only) | *(unused)* |
| `.two_block_batch_tau2_chain_row_r` | `two_block_batch_gibbs.R` | Extract chain-i tau2 row from batch state (R reference) | *(unused)* |
| `.two_block_batch_tau2_chain_row` | `two_block_batch_gibbs.R` | — | *(unused)* |
| `.two_block_batch_b_assign_slice_r` | `two_block_batch_gibbs.R` | — | *(unused)* |
| `.two_block_batch_b_assign_slice` | `two_block_batch_gibbs.R` | uses the R reference subassignment. | *(unused)* |
| `.two_block_batch_iters_ranef_add_r` | `two_block_batch_gibbs.R` | — | *(unused)* |
| `.two_block_batch_iters_ranef_add` | `two_block_batch_gibbs.R` | uses the R reference. | *(unused)* |
| `.two_block_block2_all_chains` | `two_block_batch_gibbs.R` | — | rGLMM_sweep.R |
| `.rGLMM_sweep_save` | `two_block_batch_gibbs.R` | — | rGLMM_sweep.R |

---

## LMM engines (`rLMMNormal_reg.R`, `two_block_lmm_staged_sweep_outer.R`)

| Function | File | Role | Called from |
|----------|------|------|-------------|
| `.rLMM_validate_matrix_inputs` | `rLMMNormal_reg.R` | Shared matrix-level validation for LMM replicate-chain engines | *(unused)* |
| `.rLMM_validate_P` | `rLMMNormal_reg.R` | — | *(unused)* |
| `.rLMM_validate_fixed_dispersion_prior_list` | `rLMMNormal_reg.R` | — | *(unused)* |
| `.rLMM_validate_dGamma_dispersion_prior_list` | `rLMMNormal_reg.R` | — | *(unused)* |
| `.rLMM_observation_mu` | `rLMMNormal_reg.R` | Observation-level linear predictor from group random effects | *(unused)* |
| `.rLMM_b_matrix_from_coefficients` | `rLMMNormal_reg.R` | — | *(unused)* |
| `.rLMM_icm_at_start` | `rLMMNormal_reg.R` | — | two_block_lmm_staged_sweep_outer.R |
| `.rLMM_calibrate_m_convergence` | `rLMMNormal_reg.R` | — | *(unused)* |
| `.rLMMNormal_reg_run` | `rLMMNormal_reg.R` | — | *(unused)* |
| `.rLMM_format_v2_out` | `rLMMNormal_reg.R` | — | *(unused)* |
| `.rLMM_format_sweep_out` | `two_block_lmm_staged_sweep_outer.R` | — | *(unused)* |
| `.rLMMNormal_reg_run_with_pilot` | `two_block_lmm_staged_sweep_outer.R` | ING LMM replicate chains: pilot then main via sweep-outer v6 (mirrors rGLMM) | rLMMNormal_reg.R |

---

## GLMM sweep (`rGLMM.R`, `two_block_batch_gibbs.R`)

| Function | File | Role | Called from |
|----------|------|------|-------------|
| `.rGLMM_rate_at_mode` | `rGLMM.R` | Local-Gaussian rate at the ICM mode | *(unused)* |
| `.rGLMM_format_v6_out` | `rGLMM.R` | — | *(unused)* |

---

## Multi-response / pfamily validation (`multi_rlmb.R`, `multi_rNormal_reg.R`)

| Function | File | Role | Called from |
|----------|------|------|-------------|
| `.mrglmb_check_inputs` | `multi_rlmb.R` | — | multi_rNormal_reg.R |
| `.mrglmb_n_draw` | `multi_rlmb.R` | — | multi_rNormal_reg.R |
| `.mrglmb_normalize_prior_lists` | `multi_rlmb.R` | — | multi_rNormal_reg.R |
| `.mrglmb_assemble` | `multi_rlmb.R` | — | multi_rNormal_reg.R |
| `.mrglmb_normalize_pfamily_lists` | `multi_rlmb.R` | — | *(unused)* |
| `.validate_pfamily_for_rlmb` | `multi_rlmb.R` | — | *(unused)* |
| `.validate_rindep_prior_list` | `multi_rlmb.R` | — | *(unused)* |
| `.validate_normal_gamma_prior_list` | `multi_rlmb.R` | — | *(unused)* |
| `.check_symmetric_pd` | `multi_rlmb.R` | — | multi_rNormal_reg.R |
| `.validate_normal_prior_list` | `multi_rNormal_reg.R` | — | *(unused)* |

---

## Block / simfunction utils (`simfunction_block_utils.R`)

| Function | File | Role | Called from |
|----------|------|------|-------------|
| `.prior_payload_for_rNormalGLMBlocks_cpp` | `simfunction_block_utils.R` | — | *(unused)* |
| `.prior_list_to_P_Sigma` | `simfunction_block_utils.R` | — | *(unused)* |
| `.check_P_pd` | `simfunction_block_utils.R` | — | *(unused)* |
| `normalize_prior_for_blocks` | `simfunction_block_utils.R` | — | *(unused)* |

---

## lmerb / build_mu (`build_mu_all.R`, `lmerb_posterior_mean.R`)

| Function | File | Role | Called from |
|----------|------|------|-------------|
| `build_mu_all_r` | `build_mu_all.R` | Build per-group random-effect prior means (R reference implementation) | *(unused)* |
| `.lmerb_validate_design` | `build_mu_all.R` | — | glmerb_posterior_mode.R, lmerb_posterior_mean.R |
| `.lmerb_validate_measurement_prior_list` | `lmerb_posterior_mean.R` | — | glmerb_posterior_mode.R |

---

## ING guard (`ing_prior_guard.R`)

| Function | File | Role | Called from |
|----------|------|------|-------------|
| `.ing_n_prior_from_shape` | `ing_prior_guard.R` | — | *(unused)* |
| `.ing_stop_if_prior_exceeds_data` | `ing_prior_guard.R` | — | simfunction.R, two_block_rNormal_reg.R |

---

## C++ R wrappers (`rcpp_wrappers.R`)

| Function | File | Role | Called from |
|----------|------|------|-------------|
| `.rNormalGLM_cpp` | `rcpp_wrappers.R` | — | simfunction.R, simulationpipeline.R |
| `.rNormalGLMBlocks_cpp` | `rcpp_wrappers.R` | — | *(unused)* |
| `.rNormalRegBlocks_cpp` | `rcpp_wrappers.R` | — | *(unused)* |
| `.block_rNormalReg_cpp` | `rcpp_wrappers.R` | — | simfunction_block.R |
| `.block_rNormalGLM_cpp` | `rcpp_wrappers.R` | — | simfunction_block.R |
| `.two_block_rNormal_reg_cpp` | `rcpp_wrappers.R` | — | two_block_rNormal_reg.R |
| `.two_block_rNormal_reg_staged_cpp` | `rcpp_wrappers.R` | — | *(unused)* |
| `.two_block_build_mu_all_cpp` | `rcpp_wrappers.R` | — | build_mu_all.R |
| `.two_block_batch_fixef_chain_cpp` | `rcpp_wrappers.R` | — | *(unused)* |
| `.two_block_batch_tau2_chain_row_cpp` | `rcpp_wrappers.R` | — | two_block_batch_gibbs.R |
| `.two_block_batch_b_assign_slice_cpp` | `rcpp_wrappers.R` | — | two_block_batch_gibbs.R |
| `.two_block_batch_iters_ranef_add_cpp` | `rcpp_wrappers.R` | — | two_block_batch_gibbs.R |
| `.two_block_reorder_b_to_group_levels_cpp` | `rcpp_wrappers.R` | — | two_block_batch_gibbs.R |
| `.two_block_align_b_to_xhyper_cpp` | `rcpp_wrappers.R` | — | two_block_batch_gibbs.R |
| `.two_block_block2_one_chain_cpp` | `rcpp_wrappers.R` | — | two_block_batch_gibbs.R |
| `.rNormalReg_cpp` | `rcpp_wrappers.R` | — | simfunction.R |
| `.rIndepNormalGammaReg_cpp` | `rcpp_wrappers.R` | — | simfunction.R |
| `.rNormalGammaReg_cpp` | `rcpp_wrappers.R` | — | simfunction.R |
| `.rGammaGaussian_cpp` | `rcpp_wrappers.R` | — | simfunction.R |
| `.rGammaGamma_cpp` | `rcpp_wrappers.R` | — | simfunction.R |
| `.rNormalGLM_std_cpp` | `rcpp_wrappers.R` | — | envelopeorchestrator.R, simulationpipeline.R |
| `.rIndepNormalGammaReg_std_cpp` | `rcpp_wrappers.R` | — | envelopeorchestrator.R, simulationpipeline.R |
| `.rIndepNormalGammaReg_std_parallel_cpp` | `rcpp_wrappers.R` | — | envelopeorchestrator.R |
| `.EnvelopeCentering_cpp` | `rcpp_wrappers.R` | — | envelopeorchestrator.R |
| `.EnvelopeSize_cpp` | `rcpp_wrappers.R` | — | simulationpipeline.R |
| `.EnvelopeBuild_cpp` | `rcpp_wrappers.R` | — | simulationpipeline.R |
| `.EnvelopeBuild_Ind_Normal_Gamma_cpp` | `rcpp_wrappers.R` | — | simulationpipeline.R |
| `.EnvelopeEval_cpp` | `rcpp_wrappers.R` | — | simulationpipeline.R |
| `.EnvelopeDispersionBuild_cpp` | `rcpp_wrappers.R` | — | simulationpipeline.R |
| `.EnvelopeOrchestrator_cpp` | `rcpp_wrappers.R` | — | envelopeorchestrator.R |
| `.EnvelopeSet_Grid_cpp` | `rcpp_wrappers.R` | — | simulationpipeline.R |
| `.EnvelopeSet_LogP_cpp` | `rcpp_wrappers.R` | — | simulationpipeline.R |
| `.glmb_Standardize_Model_cpp` | `rcpp_wrappers.R` | — | simulationpipeline.R |
| `.glmbayesCore_has_opencl_cpp` | `rcpp_wrappers.R` | — | gpu_diagnostics.R |
| `.gpu_names_cpp` | `rcpp_wrappers.R` | — | *(unused)* |

---

## Build, attach, misc.

| Function | File | Role | Called from |
|----------|------|------|-------------|
| `use_RcppParallel` | `internal_rcppparallel.R` | — | *(unused)* |
| `DIC_Info` | `dic_info.R` | — | summary.rglmb.R |
| `dpois2` | `simulationpipeline.R` | — | *(unused)* |
| `simfunction.default` | `simfunction.R` | — | *(unused)* |

---

## Other internals

| Function | File | Role | Called from |
|----------|------|------|-------------|
| `.two_block_summarize_pfamily_list` | `two_block_rNormal_reg.R` | — | rGLMM.R, rLMMNormal_reg.R |

---

## Review checklist (internals)

| Priority | Item |
|----------|------|
| 1 | Avoid new `@noRd` helpers unless tied to exported behavior. |
| 2 | Keep `.mrglmb_normalize_pfamily_lists` / `.validate_pfamily_for_rlmb` stable for **lmebayes** `block_core_pfamily.R`. |
| 3 | Remove or wire up *(unused)* helpers when touching related code. |

