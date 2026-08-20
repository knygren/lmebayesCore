devtools::load_all(".", quiet = TRUE)
source("data-raw/_dbg_gamma_eps_compare.R")

cat("\n\n=== epsilon* exact vs MC simulation (same data) ===\n")

mode_exact <- mode_c05
eps_exact_cl <- epsilon_star(mode_exact, method = "closure")$eps_star
eps_exact_opt <- epsilon_optimize(mode_exact)$eps_star

run_mc <- function(n) {
  mode_mc <- population_mode(
    design, pf, gaussian(), dispprior_list,
    estep = "mc", n = n, mc_seed = 42L,
    icm_init = TRUE, maxit = 200L
  )
  eps_cl <- epsilon_star(mode_mc, method = "closure")$eps_star
  eps_opt <- epsilon_optimize(mode_mc, n = n, mc_seed = 42L)$eps_star
  c(
    max_gamma_gap = max(abs(unlist(mode_mc$fixef) - unlist(mode_exact$fixef))),
    eps_closure = eps_cl,
    eps_optimize = eps_opt
  )
}

r200 <- run_mc(200L)
r10k <- run_mc(10000L)

tab <- data.frame(
  route = c("exact (closure)", "exact (optimize)", "MC n=200 closure", "MC n=200 optimize",
            "MC n=10k closure", "MC n=10k optimize"),
  eps_star = c(eps_exact_cl, eps_exact_opt, r200["eps_closure"], r200["eps_optimize"],
               r10k["eps_closure"], r10k["eps_optimize"])
)
print(tab, digits = 6)

cat("\nmax |gamma_MC200 - gamma_exact|:", signif(r200["max_gamma_gap"], 4), "\n")
cat("max |gamma_MC10k - gamma_exact|:", signif(r10k["max_gamma_gap"], 4), "\n")
