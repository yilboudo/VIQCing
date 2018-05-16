#' A function to inverse normal transform traits
#'
#' Apply inverse normal transformation
#'
#' @param x numerical vector
#'
#' @return object znorm that has been inversed normal transformed
#' @export
#'
#' @examples data= mtcars
#' apply inverse normal transformation to mile per gallon vector
#' data$inv_mpg <- invnorm(mtcars$mpg)

invnorm <- function(x)
{norm <- (x - base::mean(x ))/ stats::sd(x)
znorm<- stats::qnorm((base::rank(norm,na.last="keep")-0.5)/base::sum(!is.na(norm)))
return(znorm)
}

