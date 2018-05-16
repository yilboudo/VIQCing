#' A functions to extract confidence interval upper bound from logistic regressions
#'
#' Extract CI-U from logistic regression
#'
#' @param x numerical vector, explanatory variable
#' @param y numerical vector, depentent variable
#' @return confidence interval upper regressions
#' @export
#'
#' @examples
#' #logistics regression confidence interval upper from logistic regression
#' confidence_interval_upper <- logisticRegcov_confint_upper(y = mtcars$am ,x=mtcars$mpg)
#' #linear regression confidence interval upper bounds
#' confidence_interval_upper <- sapply(x=df, logisticRegcov_confint_upper, y=df$y, data=covariable)
logisticRegcov_confint_upper   <- function(y,...) base::exp(stats::coef(base::summary(stats::glm(y ~ . , data.frame(y,...),family="binomial")))[2,1] + (1.96*stats::coef(base::summary(stats::glm(y ~ . , data.frame(y,...),family="binomial")))[2,2]))
