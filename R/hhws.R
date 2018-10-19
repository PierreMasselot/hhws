#!#############################################################################
#!
#!                Heat-Health warning system R functions
#!                     Based on Chebana et al. (2012)
#!
#!#############################################################################

#! Overmortality from mortality baseline

baseline <- function(y, dates, nyear = Inf, side = 1, smoothing.fun = c("ma", "none", "spline"), ...)
# y: vector of health data. 
# dates: a date object
# nyear: number of years before and after for aggregation
# side: side of the aggregation
# smoothing.fun: function to presmooth y before agrgeation for baseline
# ...: arguments for presmoothing fucntion
{ 
  n <- length(y)
baseline <- function(x, dates, nyear = Inf, center = FALSE, smoothing.fun = c("ma", "spline", "none"), ...){
  n <- length(x)
  if (missing(dates)) dates <- as.POSIXlt(1:n)
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
    spline = smooth.spline(julian(dates, origin = dates[1]),x, ...)$y
  )
  if (nyear > diff(range(dates$year))){
    means <- aggregate(as.vector(ysm), by = list(doy = dates$yday), mean, na.rm=T)
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


over_mortality <- function(y, expected, fact = 100)
{
  if (missing(expected)) expected <- baseline(y, smoothing.fun = "ma", order = 15)
  om <- fact * (y - expected) / expected
  return(om)
}

over_mortality_inverse <- function(y, expected, fact = 100)
{
  if (missing(expected)) expected <- baseline(y, smoothing.fun = "ma", order = 15)    
  inv <- (y * expected / fact) + expected
  return(inv)
}

#! Finding overmortality episodes
episodes <- function(y, u, type = c("absolute","quantile"), trend = NULL, l = 0, r = 1, covariates = NULL, uc = NULL)
# x: vector of values from which the excesses are to be extracted
# u : numeric value giving the threshold above which values of x are extremes
# trend: optionnaly, a trend for the threshold. Can be a single value if the trend is linear, or a vector the same length of x for more general trends
# type: type of threshold if it is an absolute value or the quantile
# l: number of days before/after excesses for considering episodes
# r: number of consecutive values below threshold to end the episode ( >= 1)
# covariates: covariates (e.g. tmin, tmax) to filter episodes to keep
# uc: thresholds for covariates
{
  n <- length(y)
  if (!is.null(trend)){
    if (length(trend) == 1) {
      beta0 <- mean(y, na.rm=T) - mean(1:n, na.rm=T)*trend
      trend <- trend * 1:n + beta0
    } else {
      if (length(trend) != n) stop("Inconsistent length for 'trend'")
    }
    y <- y - trend + mean(y, na.rm=T)
  }    
  type <- match.arg(type)
  if (type == "quantile") u <- quantile(y,u)
  inds <- which(y > u)
  if (!is.null(covariates)){
    covariates <- as.matrix(covariates)
    if (nrow(covariates) != n) stop("'covariates' must have the same length as 'y'")
    if (is.null(uc)) stop("'uc' must be provided when 'covariates' is not null.")
    uc <- rep_len(uc, ncol(covariates))
    cond <- mapply(">",as.data.frame(covariates[inds,, drop=F]),uc, SIMPLIFY = "matrix")
    if (is.null(dim(cond))) cond <- matrix(cond,ncol=2)
    totcond <- apply(cond,1,all)
    inds <- inds[totcond]
  }
  dist_exc <- abs(outer(1:n,inds,"-"))
  inds <- which(apply(dist_exc,1,function(x)any(x <= l)))
  epis <- c(1, cumsum(diff(inds) > r) + 1)
  exts <- y[inds]
  if (!is.null(trend)){ 
    exts <- exts - mean(y, na.rm=T) + trend[inds]
    u <- u - mean(y, na.rm=T) + trend
  }
  result <- list(excesses = data.frame(t = inds, episode = epis, value = exts), threshold = rep_len(u, n))
  return(result) 
}

#! Finding threshold
find.threshold <- function(indicators, episodes, u.grid, fixed.alphas = NULL, alpha.step = .1, decreasing.alphas = TRUE, same.alphas = TRUE, order.result = c("None", "Detected", "Missed", "Sensitivity", "False_alarms", "Specificity", "Episodes_found"), order.decreasing = T)
# indicators: an array of indicators
# episodes: episodes to obtain (matrix with ind, episode and value)
# u.grid: list of threshold grids to tests for each indicator in x
# alphas: list of sequence of alphas to test (with the constraint that it sums to one) 
{
  order.result <- match.arg(order.result)
  if (!is.list(indicators)) indicators <- list(indicators)
  indicators <- lapply(indicators, as.matrix)
  p <- length(indicators)
  pvec <- sapply(indicators, ncol)
  inames <- names(indicators)
  nvec <- sapply(indicators, nrow)
  if (length(unique(nvec)) > 1) stop("All indicator matrices must have the same row number")
  n <- unique(nvec)
  if (is.null(inames)) inames <- 1:p
  if (!is.list(u.grid)) u.grid <- list(u.grid)
  u.grid <- rep_len(u.grid, p)
  if (!is.list(fixed.alphas)) fixed.alphas <- list(fixed.alphas)
  fixed.alphas <- rep_len(fixed.alphas, p)
  for (i in 1:p){ # Check that alphas sums are == 1 and that the number of alphas are correct
    al <- fixed.alphas[[i]]
    if (!is.null(al)){
      if ((sum(al) - 1) > .Machine$double.eps) stop("Weigthings in fixed.alphas must sum to one")
      if (length(al) != pvec[i]) stop("Number of weights must be identical to the number of columns of the corresponding indicator")
    }      
  } 
  # Grid of alphas to test if necessary
  alpha.grid <- fixed.alphas
  nulls <- sapply(fixed.alphas,is.null)
  if (any(nulls)){
    for (i in which(nulls)){
      ag <- expand.grid(rep(list(seq(0,1,alpha.step)), pvec[i]))
      ag <- ag[apply(ag, 1, sum) == 1,]
      if (decreasing.alphas) ag <- ag[apply(ag, 1, function(x) all(diff(as.numeric(x))<=0)),]
      alpha.grid[[i]] <- ag
    }
    na <- sapply(alpha.grid[nulls], nrow)
    if (!same.alphas || (length(unique(pvec[nulls])) > 1)) {
      igrid <- expand.grid(lapply(na[nulls], function(x) 1:x))
      for (i in 1:sum(nulls)){
        alpha.grid[[i]] <- alpha.grid[[which(nulls)[i]]][igrid[,i],]
      }
      na <- sapply(alpha.grid[nulls], nrow)
    }
    alpha.grid[which(!nulls)] <- lapply(alpha.grid[which(!nulls)], function(x) matrix(x, max(na), length(x)))
  }
  na <- max(na)
  total.ugrid <- expand.grid(u.grid)
  nu <- nrow(total.ugrid)
  result <- as.data.frame(matrix(NA, nrow = na * nu, ncol = sum(pvec) + p + 6, dimnames = list(NULL, c(unlist(mapply(function(x,y) sprintf("%s_alpha%i", x, 1:y), inames, pvec)), sprintf("s_%s", inames), "Detected", "Missed", "Sensitivity", "False_alarms", "Specificity", "Episodes_found"))))
  pb <- txtProgressBar(min = 0, max = na*nu, style = 3)
  for (i in 1:na){
    indi <- matrix(NA, n, p)
    for (k in 1:p) indi[,k] <- indicators[[k]] %*% unlist(alpha.grid[[k]][i,])
    for (j in 1:nu){
      cond <- paste(sprintf("indi[,%1$i] > total.ugrid[j,%1$i]",1:ncol(total.ugrid)),collapse = " & ")
      found <- eval(parse(text=cond))
      true.positives <- episodes[,1] %in% which(found)
      false.positives <- !which(found) %in% episodes[,1]
      false.negatives <- !episodes[,1] %in% which(found)
      true.negative <- !which(!found) %in% episodes[,1]
      sensitivity <- sum(true.positives) / (sum(true.positives) + sum(false.negatives))
      specificity <- sum(true.negative) / (sum(true.negative) + sum(false.positives))
      episodes.found <- unique(episodes[true.positives,2])
      result[(i-1)*nu + j,] <- c(unlist(sapply(alpha.grid, "[", i, )),total.ugrid[j,], sum(true.positives), sum(false.negatives), sensitivity, sum(false.positives), specificity, length(episodes.found))
      setTxtProgressBar(pb, (i-1)*nu + j)
    }
  }
  close(pb)
  dominated <- sapply(1:(na*nu), function(i) any(((result[-i,"Sensitivity"] >= result[i,"Sensitivity"] & result[-i,"Specificity"] >= result[i,"Specificity"]) & (result[-i,"Sensitivity"] != result[i,"Sensitivity"] | result[-i,"Specificity"] != result[i,"Specificity"])))) #Removing of cases for which both specificity and sensitivity are lower than at least another case
  result <- result[!dominated,]
  if (order.result != "None"){
     ord <- with(result, order(get(order.result), decreasing = order.decreasing))
     result <- result[ord,]
  }
  return(result)
}