#' A function to extract standard error from linear regressions
#'
#' Extract standard error from linear regression
#' #'
#' @param x numerical vector, explanatory variable
#' @param y numerical vector, depentent variable
#' @return standard error from regressions
#' @export
#'
#' @examples
#' #linear regression standard error
#' standard_error <- linearRegcov_stderr(y = mtcars$mpg,x=mtcars$cyl)
#' #linear regression pvalue to a data frame
#' standard_error <- sapply(x=df, linearRegcov_stderr, y=traitdf$y, data=covariable)
linearRegcov_stderr <- function(y,...) stats::coef(base::summary(stats::glm(y ~ . , data.frame(y,...) ,family="gaussian")))[2,2]
