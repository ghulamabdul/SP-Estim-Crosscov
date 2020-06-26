# Semiparametric Estimation of Cross-covariance Functions for Multivariate Random Fields
This repository provides reproducible code for the manuscript titled *Semiparametric Estimation of Cross-covariance Functions for Multivariate Random Fields* by Ghulam A. Qadir, and Ying Sun. The manuscript describes a new semiparameteric approach to model the multivariate covariance function by flexibly specifying the underlying coherence functions through B-splines. 

## Abstract

The prevalence of spatially referenced multivariate data has impelled researchers to develop procedures for joint modeling of multiple spatial processes.  This ordinarily involves modeling marginal and cross-process dependence for any arbitrary pair of locations using a multivariate spatial covariance function. However, building a flexible multivariate spatial covariance function that is nonnegative definite is challenging. Here, we propose a semiparametric approach for multivariate spatial covariance function estimation with approximate Mat{\'e}rn marginals and highly flexible cross-covariance functions via their spectral representations. The flexibility in our cross-covariance function arises due to B-spline based specification of the underlying coherence functions, which in turn allows us to capture non-trivial cross-spectral features. We then develop a likelihood-based estimation procedure and perform multiple simulation studies to demonstrate the performance of our method, especially on the coherence function estimation. Finally, we analyze particulate matter concentrations (PM_{2.5}) and wind speed data over the West-North-Central climatic region of the United States, where we illustrate that our proposed method outperforms the commonly used full bivariate Mat{\'e}rn model and the linear model of coregionalization for spatial prediction.

## Requirements

The codes are written in R, and reproducing would require installing and loading the following R-packages: `fields`,`sp`,`maps`,`maptools`,`geosphere`,`MASS`,`scoringRules`,`doParallel`,`rgdal`,`ggplot2`,`gridExtra`, `RColorBrewer`, and `viridis`. 

##

