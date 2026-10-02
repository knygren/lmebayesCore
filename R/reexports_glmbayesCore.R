## MER setup and prior symbols (Tier 1) — implemented in glmbayesCore.

#' @inherit glmbayesCore::model_setup title description details params return seealso
#' @export
model_setup <- glmbayesCore::model_setup

#' @inherit glmbayesCore::Prior_Setup_GLMM title description details params return seealso
#' @export
Prior_Setup_GLMM <- glmbayesCore::Prior_Setup_GLMM

#' @inherit glmbayesCore::check_identifiability title description details params return seealso
#' @export
check_identifiability <- glmbayesCore::check_identifiability

#' @inherit glmbayesCore::Prior_SetupGroup title description details params return seealso
#' @export
Prior_SetupGroup <- glmbayesCore::Prior_SetupGroup

#' @inherit glmbayesCore::dGamma_list title description details params return seealso
#' @export
dGamma_list <- glmbayesCore::dGamma_list

#' @inherit glmbayesCore::pfamily_list title description details params return seealso examples
#' @importFrom glmbayesCore pfamily_list
#' @export
pfamily_list <- glmbayesCore::pfamily_list

## Row-block group samplers (Tier 2) — implemented in glmbayesCore.

#' @inherit glmbayesCore::simfuncs_group title description details params return seealso examples
#' @seealso \code{\link{rNormal_reg_group_safe}}, \code{\link{rNormalGLM_reg_group_safe}}
#' @name simfuncs_group
#' @aliases rNormalGLM_reg_group rNormal_reg_group
NULL

#' @rdname simfuncs_group
#' @export
rNormalGLM_reg_group <- glmbayesCore::rNormalGLM_reg_group

#' @rdname simfuncs_group
#' @export
rNormal_reg_group <- glmbayesCore::rNormal_reg_group
