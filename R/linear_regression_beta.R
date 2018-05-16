#' A functions to extract beta coefficient from linear regressions
#'
#' Extract beta coefficient from linear regression
#'
#' @param x numerical vector, explanatory variable
#' @param y numerical vector, depentent variable
#' @return beta from linear regressions
#' @export
#'
#' @examples
#' #linear regression beta coefficient
#' beta <- linearRegcov_effect(y = mtcars$mpg,x=mtcars$cyl)
#' #linear regression beta from data frame
#' effect_size <- sapply(x=df, linearRegcov_effect, y=df$y, data=covariable)
linearRegcov_effect <- function(y,...) stats::coef(base::summary(stats::glm(y ~ . , data.frame(y,...) ,family="gaussian")))[2,1]
