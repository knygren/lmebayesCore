## Scratch: exercise every model class planned for the Chapter-C series using
## the ordinary front door (Prior_Setup_GLMM -> pfamily_list / dGamma_list ->
## rlmerb / rglmerb).  Modest n; the shipped artifacts use a large n.
## Run: Rscript data-raw/_scratch_probe_routes.R

suppressMessages(library(lmebayesCore))

con <- file("c:/Rpackages/_routes.txt", open = "wt")
on.exit(close(con), add = TRUE)
say <- function(...) writeLines(paste0(...), con)

N <- 50L

report <- function(label, expr) {
  t0 <- proc.time()[["elapsed"]]
  res <- tryCatch(force(expr), error = function(e) e)
  el <- round(proc.time()[["elapsed"]] - t0, 1)
  say("")
  say("=== ", label, "  [", el, "s]")
  if (inherits(res, "error")) {
    say("  FAILED: ", conditionMessage(res))
    return(invisible(NULL))
  }
  ci <- res$convergence_info
  say("  m_convergence   : ", paste(res$m_convergence, collapse = ","))
  say("  lambda_star     : ", if (is.null(ci$lambda_star)) "NULL" else
    round(ci$lambda_star, 6))
  for (nm in c("lambda_star_upper", "lambda_star_marginal",
               "lambda_star_combined", "sim_method_used", "draw_engine")) {
    if (!is.null(ci[[nm]])) say("  ", nm, " : ",
                                paste(round_if_num(ci[[nm]]), collapse = ", "))
  }
  sh <- res$sweep_history
  say("  sweep_history   : ", if (is.null(sh)) "NULL" else
    paste0("stage=", sh$stage, " n_sweeps=", sh$n_sweeps,
           " table_rows=", nrow(sh$table)))
  say("  pilot           : ", if (is.null(res$pilot)) "NULL" else
    paste0("n=", res$pilot$n))
  invisible(res)
}

round_if_num <- function(x) if (is.numeric(x)) round(x, 6) else x

## ------------------------------------------------------------- Gaussian
data(sleepstudy, package = "lme4", envir = environment())
dat <- sleepstudy
dat$Days_c <- dat$Days - mean(dat$Days)
form <- Reaction ~ Days_c + (1 + Days_c || Subject)

ps <- Prior_Setup_GLMM(form, data = dat, pop.pwt = 0.01)

say("ps$design present : ", !is.null(ps$design))
say("ps$design class   : ", paste(class(ps$design), collapse = "/"))

design <- ps$design
say("J = ", nlevels(design$group),
    ", groupef.names = ", paste(design$groupef.names, collapse = ", "))

set.seed(1)

## C01 - exact Gaussian iid sampler
report("C01 exact iid (sim_method = DEFAULT)", rlmerb(
  n = N, design = design,
  pfamily_list   = pfamily_list(ps),
  dispprior_list = list(dispersion = ps$group.dispersion),
  sim_method = "DEFAULT",
  progbar = FALSE, verbose = FALSE, print_icm_table = FALSE
))

## C03 - Gaussian two-block, known vcov
report("C03 two-block Gaussian, known vcov", rlmerb(
  n = N, design = design,
  pfamily_list   = pfamily_list(ps),
  dispprior_list = list(dispersion = ps$group.dispersion),
  sim_method = "TWO_BLOCK_GIBBS",
  progbar = FALSE, verbose = FALSE, print_icm_table = FALSE
))

## C07 - Gaussian, fixed dispersion, estimated vcov (ING Block 2)
report("C07 Gaussian, estimated vcov (ING)", rlmerb(
  n = N, design = design,
  pfamily_list   = pfamily_list(ps, ptypes = "dIndependent_Normal_Gamma"),
  dispprior_list = list(dispersion = ps$group.dispersion),
  progbar = FALSE, verbose = FALSE, print_icm_table = FALSE
))

## --------------------------------- per-group measurement dispersion routes
ps_g <- tryCatch(
  Prior_Setup_GLMM(form, data = dat, pop.pwt = 0.01, dispformula = ~Subject),
  error = function(e) e
)
if (inherits(ps_g, "error")) {
  say("")
  say("=== Prior_Setup_GLMM(dispformula = ~Subject) FAILED: ",
      conditionMessage(ps_g))
} else {
  dl <- dGamma_list(ps_g, warn_asymmetric = FALSE)
  say("")
  say("=== dGamma_list length ", length(dl))

  ## C06 - Gaussian, estimated group dispersion, known vcov
  report("C06 Gaussian, estimated group dispersion, known vcov", rlmerb(
    n = N, design = ps_g$design,
    pfamily_list   = pfamily_list(ps_g),
    dispprior_list = dl,
    progbar = FALSE, verbose = FALSE, print_icm_table = FALSE
  ))

  ## C09 - Gaussian, estimated dispersion AND estimated vcov
  report("C09 Gaussian, estimated dispersion AND estimated vcov", rlmerb(
    n = N, design = ps_g$design,
    pfamily_list   = pfamily_list(ps_g, ptypes = "dIndependent_Normal_Gamma"),
    dispprior_list = dl,
    progbar = FALSE, verbose = FALSE, print_icm_table = FALSE
  ))
}

## ------------------------------------------------------------- Poisson
if (requireNamespace("bayesrules", quietly = TRUE)) {
  data(airbnb, package = "bayesrules", envir = environment())
  ad <- airbnb
  ad$rating_c <- ad$rating - mean(ad$rating, na.rm = TRUE)
  ad <- ad[complete.cases(ad[, c("reviews", "rating_c", "neighborhood")]), ]
  gform <- reviews ~ rating_c + (rating_c || neighborhood)

  ps_p <- Prior_Setup_GLMM(gform, data = ad, family = poisson(),
                           pop.pwt = 0.01)
  say("")
  say("Poisson: J = ", nlevels(ps_p$design$group), ", n = ", nrow(ad))

  ## C05 - log-concave likelihood, known vcov
  report("C05 Poisson GLMM, known vcov", rglmerb(
    n = N, design = ps_p$design,
    pfamily_list = pfamily_list(ps_p),
    family = poisson(),
    progbar = FALSE, verbose = FALSE
  ))

  ## C08 - log-concave likelihood, estimated vcov
  report("C08 Poisson GLMM, estimated vcov (ING)", rglmerb(
    n = N, design = ps_p$design,
    pfamily_list = pfamily_list(ps_p, ptypes = "dIndependent_Normal_Gamma"),
    family = poisson(),
    progbar = FALSE, verbose = FALSE
  ))
} else {
  say("")
  say("bayesrules not installed - Poisson routes skipped")
}

say("")
say("DONE")
