#' A function to extract odds ratio from logistic regressions
#'
#' Extract odds ratio from logistic regression
#'
#' @param x numerical vector, explanatory variable
#' @param y numerical vector, depentent variable
#' @return odds ratio from regressions
#' @export
#'
#' @examples
#' #logistics regression odds ratio from logistic regression
#' confidence_interval_upper <- logisticRegcov_effect(y = mtcars$am ,x=mtcars$mpg)
#' #linear regression odds ratio
#' odds_ratio <- sapply(x=df, logisticRegcov_effect, y=df$y, data=covariable)
logisticRegcov_effect <- function(y,...) base::exp(stats::coef(base::summary(stats::glm(y ~ . , data.frame(y,...),family="binomial")))[2,1])
