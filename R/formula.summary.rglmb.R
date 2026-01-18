#' Model Formulae
#'
#' This function is a method function for the \code{"summary.rglmb"} class used to 
#' Extract a formulae for the objective and the family
#' @param x an object of class \code{summary.rglmb}, typically the result of a call to \link{summary.glmb}
#' @param ... further arguments to or from other methods
#' @return The function returns model formulae
#' @export
#' @method formula summary.rglmb



formula.summary.rglmb<-function(x,...){
  
  z=x
  y=z$y
  x=z$x
  return(formula(glm(y~x-1,family=family(z))))
  
}