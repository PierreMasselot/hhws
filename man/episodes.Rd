% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/episodes.R
\name{episodes}
\alias{episodes}
\title{Extreme episodes}
\usage{
episodes(x, u, type = c("absolute", "quantile"), trend = NULL, l = 0,
  r = 1, covariates = NULL, uc = NULL)
}
\arguments{
\item{x}{Numeric vector containing the data from which the extremes are to be
extracted.}

\item{u}{Numeric value giving the threshold above which values of \code{x}
are extremes.}

\item{type}{One of \code{c("absolute","quantile")}. If \code{"absolute"} 
(the default), the threshold is the value of \code{u}. If \code{"quantile"}
the threshold is taken as the quantile of order \code{u}. In the latter
case, \code{u} must be between 0 and 1.}

\item{trend}{Numeric value or vector indicating an optional trend for the 
threshold. Can be a scalar value for linear trends, or a vector the same 
length as \code{x} for more general trends.}

\item{l}{Integer. An episode is taken as the day(s) exceeding \code{u} as 
well as the \code{l} days before and after the excesses.}

\item{r}{Positive integer. Number of consecutive values below threshold
following an excess to end the episode.}

\item{covariates}{Numeric matrix. Covariates to for additional contraints on
the extracted episodes.}

\item{uc}{Numeric vector. Thresholds for covariates. Extremes of \code{x} are
return only if covariates are also above their given thresholds.}
}
\value{
A data.frame object containing the indices, values, 
    episode number of all episodes found, and a logical indicating 
    if the day is one of the extremes (useful when \code{l > 0}).

    In addition, contains the attribute
    \code{threshold} which gives he threshold for each value of \code{x}  
    for plotting, as well as \code{r}, \code{l} and \code{uc} that contains 
    values of the corresponding parameters used to extract episodes.
}
\description{
Extract episodes of extreme above a given threshold, possibly including 
  several days before and after the excess.
}
\examples{
  library(dlnm)
  data(chicagoNMMAPS)
  x <- chicagoNMMAPS$death
  dates <- as.POSIXlt(chicagoNMMAPS$date)
  n <- nrow(chicagoNMMAPS)

  # Compute over-mortality
  om <- excess(x, dates = dates, order = 15)
  
  # Extract all days for which om is above 40\%
  epis <- episodes(om, u = 40)
  # Show the corresponding dates
  dates[epis$t]

}
\references{
Chebana F., Martel B., Gosselin P., Giroux J.X., Ouarda T.B.M.J., 2013. 
     A general and flexible methodology to define thresholds for heat health 
     watch and warning systems, applied to the province of Quebec (Canada). 
     International journal of biometeorology 57, 631-644.
}
