## Convergence plot methods for rGLMM_reg / rglmerb (glmerb only).
## Must collate after plot_convergence_methods.R (aliases the rLMMNormal_reg methods).

#' @rdname plot_var_convergence
#' @method plot_var_convergence rGLMM_reg
#' @export
plot_var_convergence.rGLMM_reg <- plot_var_convergence.rLMMNormal_reg

#' @rdname plot_var_convergence
#' @method plot_var_convergence rglmerb
#' @export
plot_var_convergence.rglmerb <- plot_var_convergence.rLMMNormal_reg

#' @rdname plot_mean_convergence
#' @method plot_mean_convergence rGLMM_reg
#' @export
plot_mean_convergence.rGLMM_reg <- plot_mean_convergence.rLMMNormal_reg

#' @rdname plot_mean_convergence
#' @method plot_mean_convergence rglmerb
#' @export
plot_mean_convergence.rglmerb <- plot_mean_convergence.rLMMNormal_reg
