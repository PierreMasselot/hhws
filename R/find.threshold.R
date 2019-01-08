#' Find thresholds
#'
#' Evaluate different combinations of threshold and indicator weightings
#'    according to their sensitivity and specificity.
#'
#' @param indicators List of matrices containing the lagged covariates to 
#'    construct weighted indicators. Can also be a matrix in the case of a
#'    single indicator.
#' @param episodes A matrix or data.frame containing the indices of extreme
#'    episodes. Must contain the time indices in the first column and the 
#'    corresponding episode index in the second column as returned by  
#'    \code{\link{episodes}}.
#' @param u.grid A list of vectors containing the grid of thresholds to be 
#'    considered. The list must have the same number of elements than 
#'    \code{indicators}. If missing, the default is to use all percentiles
#'    from 0.8 to 1.
#' @param fixed.alphas A list of optional prefixed weightings for a subset of
#'    indicators. When provided, must have the same length as \code{indicators}
#'    and contain \code{NULL} for each non-fixed indicator.
#' @param alpha.step Numeric value between 0 and .5. The step of the  
#'    sequence of weightings tested for the indicator. 
#' @param decreasing.alphas Logical. If TRUE (the default), the alpha 
#'    weightings are constrained to decrease with the lag.
#' @param same.alphas Logical. If TRUE, the weightings are constrained to be
#'    identical for each indicator. Note that trying different weightings
#'    (\code{same.alpha = FALSE}) can lead to significantly longer
#'    computation time.
#' @param order.result Character or numeric value indicating a column used
#'    to order the returned table. Can also be a vector.
#' @param order.decreasing Logical. If TRUE (the default), ordering of the
#'    table is made by decreasing order to the column specified by
#'    \code{order.result}.
#' @param keep.all Logical. If FALSE (the default), all combinations having 
#'    both lower specificity and sensitivity than at least another combination
#'    are discarded from the results. Otherwise, all tested combinations
#'    are returned.
#' @param r Positive integer. Number of consecutive values below threshold  
#'   following an excess to end the episode. By default, take the attribute
#'   \code{r} in \code{episodes} or 1 if absent.
#' @param progressBar Logical indicating if a progressBar is displayed during
#'    execution of the function. If TRUE, may greatly increase execution time.
#'
#' @details We consider a warning system as a couple indicator/threshold
#'    used to launch alerts when forecasts of the indicator exceed the
#'    threshold. In the present function, the indicators considered are
#'    linear combinations of all matrix columns in the parameter 
#'    \code{indicators}, with the constraint that, for each indicator, the
#'    weights sum to 1.
#'    
#'    The indicator and threshold are determined by evaluating a large range
#'    of different weightings and threshold (given in \code{u.grid} and 
#'    \code{alpha.step}). For each combination of indicators/thresholds,
#'    the function computes the indices that consist in alerts in the data
#'    and compare them to the actual values given in \code{episodes}. The
#'    function then selects a set of best candidates (unless  
#'    \code{keep.all = TRUE}) based on sensitivity (the proportion of real
#'    episodes detected) and specificity (the proportion of false episode
#'    found). It is left to the user to choose the best combination by a 
#'    trade-off between specificity and sensitivity.   
#'
#' @return A data.frame containing a subset (unless \code{keep.all = TRUE})  
#'    of weightings and thresholds. Weightings correspond to the columns
#'    with a name containing "alpha" and threshold to names beginning with 
#'    "threshold". In addition, several scores are given in each line:
#'    \itemize{
#'      \item{Detected}{The number of indices in \code{episodes} detected by
#'        the combination.}
#'      \item{Missed}{The number of indices in \code{episodes} missed by
#'        the combination.}
#'      \item{Sensitivity}{The proportion of indices in \code{episodes} detected by
#'        the combination, \emph{i.e.} Detected / n.}
#'      \item{False_alarms}{The number of false alarms, \emph{i.e.} of indices
#'        found by the combination which are not in \code{episodes}.}
#'      \item{Specificity}{The proportion of false alarms, \emph{i.e.} 
#'        False_alarms / n.}
#'      \item{Episodes_found}{The number of episodes found, which takes account
#'        of entire episodes instead of only indices.}
#'      \item{False_episodes}{The number of false episodes found, i.e. absent
#'        from the provides \code{episodes}.}
#'    }
#'
#' @seealso [episodes()] for extracting episodes of extreme values and 
#'    [predict_alarms()] for alarms prediction.
#'
#' @references
#'    Chebana F., Martel B., Gosselin P., Giroux J.X., Ouarda T.B.M.J., 2013. 
#'      A general and flexible methodology to define thresholds for heat health 
#'      watch and warning systems, applied to the province of Quebec (Canada). 
#'      International journal of biometeorology 57, 631-644.
#'    Pascal M., Laaidi K., Ledrans M., Baffert E., Caserio-Schonemann C.,
#'      Le Tertre A., Manach J., Medina S., Rudant J., Empereur-Bissonnet P.,
#'      2006. France's heat health watch warning system. International journal
#'      of biometeorology 50, 144-153. 
#'
#' @examples
#'    library(dlnm)
#'    data(chicagoNMMAPS)
#'    x <- chicagoNMMAPS$death
#'    dates <- as.POSIXlt(chicagoNMMAPS$date)
#'    n <- nrow(chicagoNMMAPS)
#'
#'    # Compute over-mortality
#'    om <- excess(x, dates = dates, order = 15)
#'   
#'    # Extract all days for which om is above 40%
#'    epis <- episodes(om, u = 40)
#'
#'    # Prepare indicator based on temperature until lag 2
#'    indic <- matrix(NA, nrow = n, ncol = 3)
#'    indic[,1] <- chicagoNMMAPS$temp  # lag 0
#'    indic[,2] <- c(NA, chicagoNMMAPS$temp[-n]) # Lag 1
#'    indic[,3] <- c(NA, NA, chicagoNMMAPS$temp[1:(n-2)]) # lag 2
#'    # Evaluate different threshold/indicators based on these episodes
#'    find.threshold(indic, epis, u.grid = 20:35)
#'
#' @export
find.threshold <- function(indicators, episodes, u.grid, fixed.alphas = NULL, 
  alpha.step = .1, decreasing.alphas = TRUE, same.alphas = TRUE, 
  order.result = NULL, order.decreasing = TRUE, keep.all = FALSE, r = NULL,
  progressBar = FALSE) 
{
  if (!is.list(indicators)){
    nind <- deparse(substitute(indicators))
    indicators <- list(indicators)
    names(indicators) <- nind
  }
  if (is.null(r)){
    r <- ifelse(is.null(attr(episodes, "r")), 1, attr(episodes, "r"))
  } 
  indicators <- lapply(indicators, as.matrix)
  p <- length(indicators)
  pvec <- sapply(indicators, ncol)
  inames <- names(indicators)
  nvec <- sapply(indicators, nrow)
  if (length(unique(nvec)) > 1) 
    stop("All indicator matrices must have the same row number")
  n <- unique(nvec)
  if (is.null(inames)) inames <- 1:p
  if (missing(u.grid)) lapply(indicators, quantile, probs = seq(.8, 1, by = .01))
  if (!is.list(u.grid)) u.grid <- list(u.grid)
  u.grid <- rep_len(u.grid, p)
  if (!is.list(fixed.alphas)) fixed.alphas <- list(fixed.alphas)
  fixed.alphas <- rep_len(fixed.alphas, p)
  for (i in 1:p){ # Check that alphas sums are == 1 and that the number of alphas are correct
    al <- fixed.alphas[[i]]
    if (!is.null(al)){
      if ((sum(al) - 1) > .Machine$double.eps) 
        stop("Weigthings in fixed.alphas must sum to one")
      if (length(al) != pvec[i]) 
        stop("Number of weights must be identical to the number of columns of the corresponding indicator")
    }      
  } 
  # Grid of alphas to test if necessary
  alpha.grid <- fixed.alphas
  nulls <- sapply(fixed.alphas,is.null)
  if (any(nulls)){
    for (i in which(nulls)){
      ag <- expand.grid(rep(list(seq(0,1,alpha.step)), pvec[i]))
      ag <- ag[apply(ag, 1, sum) == 1,]
      if (decreasing.alphas) ag <- ag[apply(ag, 1, 
        function(x) all(diff(as.numeric(x))<=0)),]
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
    alpha.grid[which(!nulls)] <- lapply(alpha.grid[which(!nulls)], 
      function(x) matrix(x, max(na), length(x)))
  }
  na <- max(na)
  total.ugrid <- expand.grid(u.grid)
  nu <- nrow(total.ugrid)
  result <- as.data.frame(matrix(NA, nrow = na * nu, ncol = sum(pvec) + p + 7, 
    dimnames = list(NULL, c(unlist(mapply(function(x,y) 
      sprintf("%s_alpha%i", x, 1:y - 1), inames, pvec)), 
      sprintf("threshold_%s", inames), 
      "Detected", "Missed", "Sensitivity", "False_alarms", "Specificity",
      "Episodes_found", "False_episodes"))))
  if (progressBar) pb <- txtProgressBar(min = 0, max = na*nu, style = 3)
  count <- 0
  for (i in 1:na){
    indi <- matrix(NA, n, p)
    for (k in 1:p) indi[,k] <- indicators[[k]] %*% unlist(alpha.grid[[k]][i,])
    for (j in 1:nu){
      cond <- paste(sprintf("indi[,%1$i] > total.ugrid[j,%1$i]", 
        1:ncol(total.ugrid)),collapse = " & ")
      found <- eval(parse(text=cond))
      found.ind <- which(found)
#      alarms <- extract_alarms(indicators, alpha = lapply(alpha.grid, "[", i, ), 
#        s = total.ugrid[j,]) # slower
      true.positives <- episodes[,1] %in% found.ind
      false.positives <- !found.ind %in% episodes[,1]
      false.negatives <- !episodes[,1] %in% found.ind
      true.negative <- !which(!found) %in% episodes[,1]
      sensitivity <- sum(true.positives) / (sum(true.positives) + 
        sum(false.negatives))
      specificity <- sum(true.negative) / (sum(true.negative) + 
        sum(false.positives))
      episodes.found <- unique(episodes[true.positives,2])
      distToEpis <- abs(outer(found.ind, episodes[,1], "-"))
      outEpisodes <- apply(distToEpis, 1, function(x) !any(x <= 3))
      false.episodes <- c(1, cumsum(diff(found.ind[outEpisodes]) > r) + 1)
      #Removing of cases for which both specificity and sensitivity are lower than at least another case (if keep.all == FALSE)
      if (count == 0 || keep.all){
        count <- count + 1
        result[count,] <- c(unlist(sapply(alpha.grid, "[", i, )), 
          total.ugrid[j,], sum(true.positives), sum(false.negatives), 
          sensitivity, sum(false.positives), specificity, 
          length(episodes.found), length(unique(false.episodes)))
      } else {
        dominated <- any(((result[1:count, "Sensitivity"] >= sensitivity & 
          result[1:count, "Specificity"] >= specificity) & 
          (result[1:count,"Sensitivity"] != sensitivity | 
          result[1:count,"Specificity"] != specificity)))
        if (!dominated){
          dominator <- ((result[1:count, "Sensitivity"] <= sensitivity & 
            result[1:count, "Specificity"] <= specificity) & 
            (result[1:count,"Sensitivity"] != sensitivity | 
            result[1:count,"Specificity"] != specificity))
          if (any(dominator)){
            result <- result[-which(dominator),]
            count <- count - sum(dominator)
          }
          count <- count + 1
          result[count,] <- c(unlist(sapply(alpha.grid, "[", i, )), 
            total.ugrid[j,], sum(true.positives), sum(false.negatives), 
            sensitivity, sum(false.positives), specificity, 
            length(episodes.found), length(unique(false.episodes)))  
        }        
      }  
      if (progressBar) setTxtProgressBar(pb, (i-1)*nu + j)
    }
  }
  if (progressBar) close(pb)
  result <- result[1:count,]
  if (!is.null(order.result)){
    ord <- do.call(order, c(result[,order.result, drop = F], 
      list(decreasing = order.decreasing)))
    result <- result[ord,]
  }
  return(result)
}

#' Predict alarms
#'
#' Given predefined weights and thresholds, predict which of the indicators
#'    values are alarms.
#'
#' @param indicators List of matrices containing values of the lagged covariates
#'    to test whether they are alarms. Can also be a matrix in the case of a
#'    single indicator. 
#' @param alpha List of vectors giving the weightings used to construct 
#'    indicators. Must have the same length as \code{indicators} and each 
#'    vector must have the same length as \code{ncol} of each matric in 
#'    \code{indicators}.
#' @param s Vector of thresholds. Must have the same length as 
#'    \code{indicators}.
#' @param y Optional vector of corresponding response value. Must have the same
#'    length as the number of rows of each matrix in \code{indicators}.
#' @param r Positive integer. Number of consecutive values below threshold  
#'   following an alarm to end the episode.
#'
#' @details Extracts the indices in \code{indicators} that correspond to an
#'    alarm according to the given weigths and thresholds. If \code{y} is
#'    given, the corresponding values of the response \emph{e.g.}
#'    over-mortality are also returned. In addition, episodes are formed
#'    by consecutive alarms (or by alarms occurring with time differences)
#'    lower or equal than \code{r}.
#'
#' @return A data.frame object containing the indices, episode number, and the 
#'    values of both the response and indicators for all alarms found.
#'
#' @seealso \code{\link{predict_alarms}} to extract the detected alarms given
#'    indicators and thresholds.
#'
#' @references
#'    Chebana F., Martel B., Gosselin P., Giroux J.X., Ouarda T.B.M.J., 2013. 
#'      A general and flexible methodology to define thresholds for heat health 
#'      watch and warning systems, applied to the province of Quebec (Canada). 
#'      International journal of biometeorology 57, 631-644.
#'
#' @examples
#'   library(dlnm)
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
#'
#'   # Prepare indicator based on temperature until lag 2
#'   indic <- matrix(NA, nrow = n, ncol = 3)
#'   indic[,1] <- chicagoNMMAPS$temp  # lag 0
#'   indic[,2] <- c(NA, chicagoNMMAPS$temp[-n]) # Lag 1
#'   indic[,3] <- c(NA, NA, chicagoNMMAPS$temp[1:(n-2)]) # lag 2
#'   # Evaluate different threshold/indicators based on these episodes
#'   tested <- find.threshold(indic, epis, u.grid = 20:35)
#'
#'   # Choose a result and predict
#'   final <- tested[19,]
#'   predict_alarms(indic, final[1:3], s = final[4], y = om)
#'
#' @export
predict_alarms <- function(indicators, alpha, s, y = NA, r = 1)
{
  if (!is.list(indicators)){
    nind <- deparse(substitute(indicators))
    indicators <- list(indicators)
    names(indicators) <- nind
  }
  indicators <- lapply(indicators, as.matrix)
  nvec <- sapply(indicators, nrow)
  n <- unique(nvec)
  stopifnot(length(n) == 1)
  if (!is.list(alpha) || is.data.frame(alpha)) alpha <- list(alpha)
  alpha <- lapply(alpha, unlist)
  stopifnot(length(unique(length(indicators), length(alpha), length(s))) == 1)
  indics <- mapply("%*%", indicators, alpha, SIMPLIFY = FALSE)
  inds.indiv <- mapply(">", indics, s)
  inds <- which(apply(as.matrix(inds.indiv), 1, all))
  epis <- c(1, cumsum(diff(inds) > r) + 1)
  y <- rep_len(y, n)
  result <- data.frame(t = inds, episode = epis, value = y[inds], 
    sapply(indics, "[", inds, ))
  return(result)
}