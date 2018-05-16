#' A function to extract sample size from linear regressions
#'
#' Extract sample size from linear regression
#'
#' @param x numerical vector, explanatory variable
#' @param y numerical vector, depentent variable
#' @return sample size from regressions
#' @export
#'
#' @examples
#' #linear regression sample size
#' sample_size <- linearRegcov_samplesize(y = mtcars$mpg,x=mtcars$cyl)
#' #linear regression pvalue to a data frame
#' sample_size <- sapply(x=df, linearRegcov_samplesize, y=df$y, data=covariable)
linearRegcov_samplesize <- function(y,...) base::nrow(stats::model.frame(stats::glm(y ~ . , data.frame(y,...) ,family="gaussian")))

