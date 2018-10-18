#!#############################################################################
#!
#!                Heat-Health warning system R functions
#!                     Based on Chebana et al. (2012)
#!
#!#############################################################################

#! A FAIRE:
#! find.threshold: problème de format matrice quand il y a une seule variable
#! Ajouter possibilité de lags différents pour les variables
#! Ajouter possibilité de poids alpha différents pour les différentes variables

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
  if (missing(dates)) dates <- as.POSIXlt(1:n)
  if (is.character(dates)) dates <- as.POSIXlt(dates)
  stopifnot(length(dates) == n)
  y <- y[order(dates)]
  dates <- dates[order(dates)]
  if ((n-1) != diff(range(dates))) {
    complete.seq <- seq(dates[1],dates[n],"day")
    orig.dates <- as.Date(complete.seq) %in% as.Date(dates)
    newy <- rep(NA,length(complete.seq))
    newy[orig.dates] <- y
    y <- newy
    dates <- as.POSIXlt(complete.seq)
    n <- length(y)
  } else {
    orig.dates <- rep(TRUE,n)
  }
  smoothing.fun <- match.arg(smoothing.fun)
  ysm <- switch(smoothing.fun,
    none = y,
    ma = forecast::ma(y, ...),
    spline = smooth.spline(julian(dates, origin = dates[1]),y, ...)$y
  )
  if (nyear > diff(range(dates$year))){
    means <- aggregate(as.vector(ysm), by = list(doy = dates$yday), mean, na.rm=T)
    ms <- means$x
    names(ms) <- means$doy
    ey <- ms[as.character(dates$yday)]
  } else {
    bydoy <- split(ysm, dates$yday)
    means <- lapply(bydoy, function(x){
      wei <- matrix(0,sum(!is.na(x)),sum(!is.na(x)))
      wei[abs(col(wei) - row(wei)) < nyear] <- 1
      if (side == 1) wei[(row(wei)-col(wei)) < 0] <- 0
      wei <- wei / rowSums(wei)
      res <- rep(NA,length(x))
      res[!is.na(x)] <- wei %*% na.omit(x)
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
find.threshold <- function(indicators, episodes, u.grid, alphas = seq(0,1,.1), trim = NULL)
# indicators: an array of indicators
# episodes: episodes to obtain (matrix with ind, episode and value)
# u.grid: list of threshold grids to tests for each indicator in x
# alphas: sequence of alphas to test (with the constraint that it sums to one) 
{
  if (length(dim(indicators)) == 2) dim(indicators) <- c(dim(indicators),1) 
  alpha.grid <- expand.grid(rep(list(alphas),dim(indicators)[2]))
  cond1 <- apply(alpha.grid,1,sum) == 1
  cond2 <- apply(alpha.grid,1,function(x) all(diff(as.numeric(x)) <= 0))
  alpha.grid <- alpha.grid[cond1 & cond2,] 
  na <- nrow(alpha.grid)   
  total.ugrid <- expand.grid(u.grid)
  nu <- nrow(total.ugrid)
  result <- as.data.frame(matrix(NA, nrow =na*nu, ncol = sum(dim(indicators)[2:3])+6, dimnames = list(NULL, c(paste("alpha",1:dim(indicators)[2]-1, sep=""), dimnames(indicators)[[3]], "Detected", "Missed", "Sensitivity", "False_alarms", "Specificity", "Episodes_found"))))
  pb <- txtProgressBar(min = 0, max = na*nu, style = 3)
  for (i in 1:na){
    indi <- apply(indicators,3,function(x) x %*% as.numeric(alpha.grid[i,]))
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
      result[(i-1)*nu + j,] <- c(alpha.grid[i,],total.ugrid[j,], sum(true.positives), sum(false.negatives), sensitivity, sum(false.positives), specificity, length(episodes.found))
      setTxtProgressBar(pb, (i-1)*nu + j)
    }
  }
  dominated <- sapply(1:(na*nu), function(i) any((result[-i,"Sensitivity"] > result[i,"Sensitivity"] & result[-i,"Specificity"] > result[i,"Specificity"])))
  result <- result[!dominated,]
  result <- result[order(result[,"Specificity"] + result["Episodes_found"]/max(episodes[,2]), decreasing = T),]
  if (!is.null(trim)) result <- result[1:trim,]
  close(pb)
  return(result)
}