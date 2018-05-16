#' A functions to extract confidence interval lower bound from logistic regressions
#'
#' Extract CI-L from logistic regression
#'
#' @param x numerical vector, explanatory variable
#' @param y numerical vector, depentent variable
#' @return confidence interval lower regressions
#' @export
#'
#' @examples
#' #logistics regression confidence interval lower from logistic regression
#' confidence_interval_lower <- logisticRegcov_confint_lower(y = mtcars$am ,x=mtcars$mpg)
#' #linear regression confidence interval lower bounds
#' confidence_interval_lower <- sapply(x=df, logisticRegcov_confint_lower, y=df$y, data=covariable)
logisticRegcov_confint_lower <- function(y,...) base::exp(stats::coef(base::summary(stats::glm(y ~ . , data.frame(y,...),family="binomial")))[2,1] - (1.96*coef(base::summary(stats::glm(y ~ . , data.frame(y,...),family="binomial")))[2,2]))
