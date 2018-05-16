#' A function to extract sample size from logistic regression
#'
#' Extract sample size logistic regression
#'
#' @param x numerical vector, explanatory variable
#' @param y numerical vector, depentent variable
#' @return sample size from logistic regressions
#' @export
#'
#' @examples
#' #logistic regression pvalue
#' pval <- linearRegcov_p(y = mtcars$mpg,x=mtcars$cyl)
#' #logistic regression sample size to a data frame
#' sample_size <- sapply(x=df, logisticRegcov_samplesize, y=df$y, data=covariable)
logisticRegcov_samplesize <- function(y,...) base::nrow(stats::model.frame(stats::glm(y ~ . , data.frame(y,...),family="binomial")))
