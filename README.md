# hhws

R functions for choosing indicator weightings and thresholds for heat-health warning systems.

## Description

Functions for finding indicators and thresholds for an heat-health warning system. These functions can be used for systems similar to those currently implemented in France and in the province of Quebec, Canada. It allows computing excess mortality series, episodes of extreme excesses and find optimal weights and thresholds for the predefined indicators. 

## Content

The package contains several functions, that can be used in the following order: 

* `baseline`: computes the baseline of expected values of a series;
* `excess`, `excess_inverse`: computes the series of excesses relatively to a baseline, or invert the process by computing original values from a series of excesses and a baseline;
* `episodes`: extract episodes of extreme excess values;
* `find.threshold`: evaluates a grid of indicator weights and thresholds to determine the best warning system;
* `predict_alarms`: considering chosen indicators and thresholds, return the corresponding days.
 
See the help to details on how to use each function, *e.g.* `?baseline`. 

## Installation

1. In R, install the package directly from github using the command:
```r
> install_github("PierreMasselot/hhws")
```
2. The package can then be loaded as usual: `library(hhws)`.

## Example

```r

library(dlnm)
data(chicagoNMMAPS)
x <- chicagoNMMAPS$death
dates <- as.POSIXlt(chicagoNMMAPS$date)
n <- nrow(chicagoNMMAPS)

# Compute over-mortality
om <- excess(x, dates = dates, order = 15)
  
# Extract all days for which om is above 40\%
epis <- episodes(om, u = 40)

# Prepare indicator based on temperature until lag 2
indic <- matrix(NA, nrow = n, ncol = 3)
indic[,1] <- chicagoNMMAPS$temp  # lag 0
indic[,2] <- c(NA, chicagoNMMAPS$temp[-n]) # Lag 1
indic[,3] <- c(NA, NA, chicagoNMMAPS$temp[1:(n-2)]) # lag 2

# Evaluate different threshold/indicators based on these episodes
tested <- find.threshold(indic, epis, u.grid = 20:35, thinning = "episodes", order.result = "Episodes_found")

# Choose a result and predict
final <- tested[1,]
predict_alarms(indic, final[1:3], s = final[4], y = om)

```

## References

Chebana, F., Martel, B., Gosselin, P., Giroux, J.-X., Ouarda, T.B., 2013. A general and flexible methodology to define thresholds for heat health watch and warning systems, applied to the province of Québec (Canada). *International journal of biometeorology* **57**, 631–644.

Pascal, M., Laaidi, K., Ledrans, M., Baffert, E., Caserio-Schönemann, C., Le Tertre, A., Manach, J., Medina, S., Rudant, J., Empereur-Bissonnet, P., 2006. France’s heat health watch warning system. *International journal of biometeorology* **50**, 144–153. 
