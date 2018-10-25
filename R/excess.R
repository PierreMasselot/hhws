#' Excess from a baseline
#'
#' \code{excess} computes the relative excesses of \code{x} relatively to a 
#'   baseline. \code{excess_inverse} reverts the process.
#'
#' @param x Numeric vector. Contains original data.
#' @param expected Numeric vector the same length as \code{x}. The baseline 
#'   used in the computation of excesses. If missing in \code{excess}, the  
#'   function \code{\link{baseline}} is called with the additional arguments 
#'   passed in \code{...}. If missing in \code{excess_inverse}, the attribute 
#'   \code{baseline} is used if present.
#' @param fact Numeric. A factor to apply to the relative excess. Default to 100
#'   to express excesses in percentage.
#' @param ... Arguments to be passed to \code{\link{baseline}} if 
#'   \code{expected} is missing. See \code{?baseline}.
#'
#' @details Compute relative excesses from a baseline. Excesses are computed as:
#'   \deqn{excess = (x - expected) / expected * fact}
#'
#' @return \code{excess}: A numeric vector of class \code{excess} the same 
#'   length as \code{x}.In addition, contains the attributes \code{baseline} and
#'   \code{fact} to keep track of the computation.
#'
#' @seealso [baseline()] for computing expected values. 
#'
#' @references
#'   Chebana F., Martel B., Gosselin P., Giroux J.X., Ouarda T.B.M.J., 2013. 
#'     A general and flexible methodology to define thresholds for heat health 
#'     watch and warning systems, applied to the province of Quebec (Canada). 
#'     International journal of biometeorology 57, 631-644. 
#'
#' @examples
#'   data(chicagoNMMAPS)
#'   x <- chicagoNMMAPS$death
#'   dates <- as.POSIXlt(chicagoNMMAPS$date)
#'
#'   om <- excess(x, dates = dates, order = 15)
#'   xrec <- excess_inverse(om)
#'
#'   plot(dates, x)
#'   lines(dates, xrec, col = "red")
excess <- function(x, expected, fact = 100, ...)
{
  if (missing(expected)) expected <- baseline(x, ...)
  om <- fact * (x - expected) / expected
  attributes(om) <- list(baseline = expected, fact = fact)
  return(om)
}

#' @rdname excess
#'
#' @param excess Numeric vector. A series of excesses to be reverted to the 
#'   original scale of \code{x}.
#' @return \code{excess_inverse}: A numeric vector the same length as 
#'   \code{excess}.
excess_inverse <- function(excess, expected, fact)
{  
  if (missing(expected)){
    if (!is.null(attr(excess, "baseline"))){
       expected <- attr(excess, "baseline")
    } else {
       stop("expected must be provided")
    }
  }
  if (missing(fact)){
    if (!is.null(attr(excess, "baseline"))){
       fact <- attr(excess, "fact")
    } else {
       fact <- 100
    }
  }  
  inv <- (excess * expected / fact) + expected
  return(inv)
}