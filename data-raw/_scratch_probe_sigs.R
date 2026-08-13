## Scratch: dump exported signatures needed for the Chapter-C vignette build.
con <- file("c:/Rpackages/_sig.txt", open = "wt")
on.exit(close(con), add = TRUE)

say <- function(...) writeLines(paste0(...), con)

say("R: ", R.version.string)
suppressMessages(library(lmebayesCore))
say("lmebayesCore loaded OK, version ",
    as.character(utils::packageVersion("lmebayesCore")))

fns <- c(
  "rLMMNormal_reg_known_vcov_iid", "rLMMNormal_reg_known_vcov_two_bg",
  "rlmerb", "rglmerb", "two_block_rate", "two_block_l_for_tv",
  "two_block_tv_bound", "two_block_rate_from_pfamily_list",
  "dGamma_list", "pfamily_list", "model_setup", "Prior_Setup_GLMM",
  "plot_var_convergence", "plot_mean_convergence"
)

for (f in fns) {
  if (!exists(f)) {
    say(f, ": NOT EXPORTED")
    next
  }
  fo <- formals(get(f))
  parts <- vapply(names(fo), function(nm) {
    d <- paste(deparse(fo[[nm]]), collapse = " ")
    if (identical(d, "")) nm else paste0(nm, " = ", d)
  }, character(1))
  say(f, "(", paste(parts, collapse = ", "), ")")
  say("")
}

say("DONE")
