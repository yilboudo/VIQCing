#' A function to extract p-value from logistic regressions
#'
#' Extract p-value from logistic regression
#'
#' @param x numerical vector, explanatory variable
#' @param y numerical vector, depentent variable
#' @return p-value from regressions
#' @export
#'
#' @examples
#' #logistics regression p-value from logistic regression
#' pval <- logisticRegcov_p(y = mtcars$am ,x=mtcars$mpg)
#' #linear regression p-value
#' pval <- sapply(x=df, logisticRegcov_p, y=df$y, data=covariable)
logisticRegcov_p <- function(y,...) stats::coef(base::summary(stats::glm(y ~ . , data.frame(y,...) ,family="binomial")))[2,4]
