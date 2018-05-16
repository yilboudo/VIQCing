#' A functions to extract p-value from linear regressions
#'
#' Extract p-value from linear regression
#'
#' @param x numerical vector, explanatory variable
#' @param y numerical vector, depentent variable
#' @return p-value from linear regressions
#' @export
#'
#' @examples
#' #linear regression pvalue
#' pval <- linearRegcov_p(y = mtcars$mpg,x=mtcars$cyl)
#' #linear regression pvalue from data frame
#' observed_pval <- sapply(x=df, linearRegcov_p, y=trait, data=covariable)
linearRegcov_p <- function(y,...) stats::coef(base::summary(stats::glm(y ~ . , data.frame(y,...) ,family="gaussian")))[2,4]
