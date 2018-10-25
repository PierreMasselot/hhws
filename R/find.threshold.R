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
#'    \code{indicators}.
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
#' @param order.result Character of numeric value indicating a column used
#'    to order the returned table.
#' @param order.decreasing Logical. If TRUE (the default), ordering of the
#'    table is made by decreasing order to the column specified by
#'    \code{order.result}.
#' @param keep.all Logical. If FALSE (the default), all combinations having 
#'    both lower specificity and sensitivity than at least another combination
#'    are discarded from the results. Otherwise, all tested combinations
#'    are returned.
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
#'    with a name containing "alpha" and threshold to names beginning with "s_".
#'    In addition, several scores are given in each line:
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
#'    }
#'
#' @seealso [episodes()] for extracting episodes of extreme values. 
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
find.threshold <- function(indicators, episodes, u.grid, fixed.alphas = NULL, 
  alpha.step = .1, decreasing.alphas = TRUE, same.alphas = TRUE, 
  order.result = NULL, order.decreasing = TRUE, keep.all = FALSE) 
{
  if (!is.list(indicators)){
    nind <- deparse(substitute(indicators))
    indicators <- list(indicators)
    names(indicators) <- nind
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
  result <- as.data.frame(matrix(NA, nrow = na * nu, ncol = sum(pvec) + p + 6, 
    dimnames = list(NULL, c(unlist(mapply(function(x,y) 
      sprintf("%s_alpha%i", x, 1:y - 1), inames, pvec)), sprintf("s_%s", inames), 
      "Detected", "Missed", "Sensitivity", "False_alarms", "Specificity",
      "Episodes_found"))))
  pb <- txtProgressBar(min = 0, max = na*nu, style = 3)
  for (i in 1:na){
    indi <- matrix(NA, n, p)
    for (k in 1:p) indi[,k] <- indicators[[k]] %*% unlist(alpha.grid[[k]][i,])
    for (j in 1:nu){
      cond <- paste(sprintf("indi[,%1$i] > total.ugrid[j,%1$i]", 
        1:ncol(total.ugrid)),collapse = " & ")
      found <- eval(parse(text=cond))
      true.positives <- episodes[,1] %in% which(found)
      false.positives <- !which(found) %in% episodes[,1]
      false.negatives <- !episodes[,1] %in% which(found)
      true.negative <- !which(!found) %in% episodes[,1]
      sensitivity <- sum(true.positives) / (sum(true.positives) + 
        sum(false.negatives))
      specificity <- sum(true.negative) / (sum(true.negative) + 
        sum(false.positives))
      episodes.found <- unique(episodes[true.positives,2])
      result[(i-1)*nu + j,] <- c(unlist(sapply(alpha.grid, "[", i, )), 
        total.ugrid[j,], sum(true.positives), sum(false.negatives), sensitivity,
        sum(false.positives), specificity, length(episodes.found))
      setTxtProgressBar(pb, (i-1)*nu + j)
    }
  }
  close(pb)
  if (!keep.all){
    #Removing of cases for which both specificity and sensitivity are lower than at least another case
    dominated <- sapply(1:(na*nu), function(i) any(((result[-i,"Sensitivity"] >= 
      result[i,"Sensitivity"] & result[-i,"Specificity"] >= 
      result[i,"Specificity"]) & (result[-i,"Sensitivity"] != 
      result[i,"Sensitivity"] | result[-i,"Specificity"] != 
      result[i,"Specificity"]))))
    result <- result[!dominated,]
  }
  if (!is.null(order.result)){
     ord <- order(result[,order.result], decreasing = order.decreasing)
     result <- result[ord,]
  }
  return(result)
}