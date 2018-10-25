#' Extreme episodes
#'
#' Extract episodes of extreme above a given threshold, possibly including 
#'   several days before and after the excess.
#'
#' @param x Numeric vector containing the data from which the extremes are to be
#'   extracted.
#' @param u Numeric value giving the threshold above which values of \code{x}  
#'   are extremes.
#' @param type One of \code{c("absolute","quantile")}. If \code{"absolute"} 
#'   (the default), the threshold is the value of \code{u}. If \code{"quantile"}
#'   the threshold is taken as the quantile of order \code{u}. In the latter
#'   case, \code{u} must be between 0 and 1.
#' @param trend Numeric value or vector indicating an optional trend for the 
#'   threshold. Can be a scalar value for linear trends, or a vector the same 
#'   length as \code{x} for more general trends.
#' @param l Integer. An episode is taken as the day(s) exceeding \code{u} as 
#'   well as the \code{l} days before and after the excesses.
#' @param r Positive integer. Number of consecutive values below threshold  
#'   following an excess to end the episode.
#' @param covariates Numeric matrix. Covariates to for additional contraints on
#'   the extracted episodes.
#' @param uc Numeric vector. Thresholds for covariates. Extremes of \code{x} are
#'   return only if covariates are also above their given thresholds.
#'
#' @return A data.frame object containing the indices, values and 
#'     episode number of all episodes found.
#' 
#'     In addition, contains the attribute
#'     \code{threshold} which gives he threshold for each value of \code{x}.  
#'     Useful for plotting.
#'
#' @references
#'    Chebana F., Martel B., Gosselin P., Giroux J.X., Ouarda T.B.M.J., 2013. 
#'      A general and flexible methodology to define thresholds for heat health 
#'      watch and warning systems, applied to the province of Quebec (Canada). 
#'      International journal of biometeorology 57, 631-644.
#'
#' @examples
#'   data(chicagoNMMAPS)
#'   x <- chicagoNMMAPS$death
#'   dates <- as.POSIXlt(chicagoNMMAPS$date)
#'   n <- nrow(chicagoNMMAPS)
#'
#'   # Compute over-mortality
#'   om <- excess(x, dates = dates, order = 15)
#'   
#'   # Extract all days for which om is above 40%
#'   epis <- episodes(om, u = 40)
#'   # Show the corresponding dates
#'   dates[epis$t]
episodes <- function(x, u, type = c("absolute","quantile"), trend = NULL, l = 0,
  r = 1, covariates = NULL, uc = NULL)
{
  n <- length(x)
  if (!is.null(trend)){
    if (length(trend) == 1) {
      beta0 <- mean(x, na.rm=T) - mean(1:n, na.rm=T) * trend
      trend <- trend * 1:n + beta0
    } else {
      if (length(trend) != n) stop(sprintf("trend must be either of length 1 or 
        of length %i", n))
    }
    x <- x - trend + mean(x, na.rm=T)
  }    
  type <- match.arg(type)
  if (type == "quantile") u <- quantile(x,u)
  inds <- which(x > u)
  if (!is.null(covariates)){
    covariates <- as.matrix(covariates)
    if (nrow(covariates) != n) 
      stop("'covariates' must have the same length as 'y'")
    if (is.null(uc)) 
      stop("'uc' must be provided when 'covariates' is not null.")
    uc <- rep_len(uc, ncol(covariates))
    cond <- mapply(">", as.data.frame(covariates[inds,, drop=F]),uc, 
      SIMPLIFY = "matrix")
    if (is.null(dim(cond))) cond <- matrix(cond, ncol=2)
    totcond <- apply(cond, 1, all)
    inds <- inds[totcond]
  }
  dist_exc <- abs(outer(1:n, inds, "-"))
  inds <- which(apply(dist_exc, 1, function(y) any(y <= l)))
  epis <- c(1, cumsum(diff(inds) > r) + 1)
  exts <- x[inds]
  if (!is.null(trend)){ 
    exts <- exts - mean(x, na.rm=T) + trend[inds]
    u <- u - mean(x, na.rm=T) + trend
  }
  result <- data.frame(t = inds, episode = epis, value = exts)
  attr(result, "threshold") <- rep_len(u, n)
  return(result) 
}