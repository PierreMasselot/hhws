#' Baseline 
#'
#' Compute a baseline from a series through two steps: 1) smoothing and 
#'   2) calendar mean.
#' 
#' @param x Numeric vector. Contains data from which the baseline is computed.
#' @param dates POSIXt vector the same length as \code{x}. Dates corresponding
#'   to each element of \code{x}. Alternatively, a character vector can be
#'   provided as long as it can be coerced to a POSIXt object.
#' @param nyear Integer. Number of years around the current one to consider 
#'   in the calendar mean step.
#' @param center Logical. Is the calendar mean of step 2 centered on current 
#'   year? If FALSE (the default), the past \code{nyear} are used.
#' @param smoothing.fun One of \code{"ma"}, \code{"spline"} or \code{"none"}.
#'   Type of smoothing to perform on step 1. If \code{"ma"} (the default), 
#'   a moving average is performed using the function \code{\link[forecast]{ma}}
#'   from the package \code{forecast}. If \code{"spline"}, 
#'   \code{\link[stats]{smooth.spline}} is used. If \code{"none"}, 
#'   no smoothing is performed.
#' @param ... Additional arguments to be passed to function chosen in 
#'   \code{smoothing.fun}.
#'
#' @details Computes a baseline of expected mortality (or other health issue) 
#'   to use as a reference to later compute a series of excess mortality. 
#'   The baseline is computed through two steps:
#'   \enumerate{
#'     \item Smoothing of the series \code{x} ;
#'     \item Computing of a calendar mean, \emph{i.e.} each day of year of 
#'        the baseline is the mean 
#'        of the smoothed series at the same day of year for \code{nyear}.
#'   }
#'
#' @return A numeric vector the same length as \code{x}.
#'
#' @seealso \code{\link{excess}} to compute excesses from the baseline.
#'
#' @references
#'   Chebana F., Martel B., Gosselin P., Giroux J.X., Ouarda T.B.M.J., 2013. 
#'     A general and flexible methodology to define thresholds for heat health 
#'     watch and warning systems, applied to the province of Quebec (Canada). 
#'     International journal of biometeorology 57, 631-644. 
#'
#' @examples
#'   library(dlnm)
#'   data(chicagoNMMAPS)
#'   x <- chicagoNMMAPS$death
#'   dates <- as.POSIXlt(chicagoNMMAPS$date)
#'
#'   em <- baseline(x, dates, order = 15)
#'   plot(dates, x)
#'   lines(dates, em, col = "red")
#'
#'   em2 <- baseline(x, dates, smoothing.fun = "spline")
#'   plot(dates, x)
#'   lines(dates, em2, col = "red")
#'
#'   em3 <- baseline(x, dates, nyear = 2, order = 15)
#'   plot(dates, x)
#'   lines(dates, em3, col = "red")
#'
#' @export
baseline <- function(x, dates, nyear = Inf, center = FALSE, 
  smoothing.fun = c("ma", "spline", "none"), ...)
{
  n <- length(x)
  if (missing(dates)) dates <- seq.POSIXt(as.POSIXlt("1970-01-01"), by = "day", 
    length.out = n)
  if (is.character(dates)) dates <- as.POSIXlt(dates)
  stopifnot(length(dates) == n)
  x <- x[order(dates)]
  dates <- dates[order(dates)]
  if ((n-1) != diff(range(dates))) {
    complete.seq <- seq(dates[1],dates[n],"day")
    orig.dates <- as.Date(complete.seq) %in% as.Date(dates)
    newy <- rep(NA,length(complete.seq))
    newy[orig.dates] <- x
    x <- newy
    dates <- as.POSIXlt(complete.seq)
    n <- length(x)
  } else {
    orig.dates <- rep(TRUE,n)
  }
  smoothing.fun <- match.arg(smoothing.fun)
  ysm <- switch(smoothing.fun,
    none = x,
    ma = forecast::ma(x, ...),
    spline = stats::smooth.spline(julian(dates, origin = dates[1]),x, ...)$y
  )
  if (nyear > diff(range(dates$year))){
    means <- aggregate(as.vector(ysm), by = list(doy = dates$yday), mean, 
      na.rm=T)
    ms <- means$x
    names(ms) <- means$doy
    ey <- ms[as.character(dates$yday)]
  } else {
    bydoy <- split(ysm, dates$yday)
    means <- lapply(bydoy, function(doy){
      wei <- matrix(0, sum(!is.na(doy)), sum(!is.na(doy)))
      if (center){
        wei[abs(col(wei) - row(wei)) < (nyear / 2)] <- 1
      } else {
        wei[(row(wei) - col(wei)) %in% (1:nyear - 1)] <- 1
      }
      wei <- wei / rowSums(wei)
      res <- rep(NA,length(doy))
      res[!is.na(doy)] <- wei %*% na.omit(doy)
      return(res)
    })
    ey <- vector("numeric",n)
    for (i in unique(dates$yday)) ey[dates$yday == i] <- means[[i+1]]
  }
  return(ey[orig.dates])
}