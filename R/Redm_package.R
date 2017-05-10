#' @name rEDM
#' @docType package
#' @title Applications of empirical dynamic modeling from time series.
#' @author Hao Ye
#' @description The rEDM package provides an interface from R to C++ compiled objects 
#' that use time delay embedding to perform state-space reconstruction and nonlinear 
#' forecasting.
#' @keywords package
NULL

#' @name paramecium_didinium
#' @docType data
#' @title Time series for the Paramecium-Didinium laboratory experiment
#' @author Veilleux
#' @description Time series of Paramecium and Didinium abundances (#/mL) from 
#' an experiment by Veilleux (1979)
NULL

#' @name sardine_anchovy_sst
#' @docType data
#' @title Time series for the California Current Anchovy-Sardine-SST system
#' @author ****
#' @description Time series of Pacific sardine landings (CA), Northern anchovy landings (CA), 
#' and sea-surface temperature (3-year average) at the SIO pier and Newport pier
NULL

#' @name tentmap_del
#' @docType data
#' @title Time series for a tent map with mu = 2.
#' @author Hao Ye
#' @description First-differenced time series generated from the tent map recurrence relation with mu = 2.
NULL

#' @name two_species_model
#' @docType data
#' @title Time series for a two-species coupled model.
#' @author Hao Ye
#' @description Time series generated from a discrete-time coupled Lotka-Volterra 
#' model exhibiting chaotic dynamics.
NULL

#' @name block_3sp
#' @docType data
#' @title Time series for a three-species coupled model.
#' @author Hao Ye
#' @description Time series generated from a discrete-time coupled Lotka-Volterra 
#' model exhibiting chaotic dynamics.
NULL

#' @name sockeye_returns
#' @docType data
#' @title Time series for sockeye salmon returns.
#' @author ****
#' @description Time series of sockeye salmon returns from the Fraser River in 
#' British Columbia, Canada.
NULL

#' @name e054_succession
#' @docType data
#' @title Succession data at the Cedar Creek LTER
#' @author ****
#' @description Experiment 054 is a subset of the long-term observational study 
#' of old field succession at the Cedar Creek LTER.
NULL

#' @name e120_biodiversity
#' @docType data
#' @title Biodiversity data at the Cedar Creek LTER
#' @author ****
#' @description Experiment 120, the "Big Biodiversity" experiment at Cedar Creek 
#' LTER. This experiment is the longest running randomized test for the effects 
#' of plant diversity on ecosystem functions.
NULL

#' @name thrips_block
#' @docType data
#' @title Apple-blossom Thrips time series
#' @author ****
#' @description Seasonal outbreaks of Thrips imaginis.
NULL

#' @name compute_stats
#' 
#' @title Compute performance metrics for predictions
#'
#' @description Computes the rho, MAE, RMSE, perc, and p-val performance metrics
#' using the compiled C++ function
#' 
#' @param observed a vector of the observed values
#' @param predicted a vector of the corresponding predicted values
#' 
#' @return A data.frame with components for the various performance metrics:
#' \tabular{ll}{
#'   num_pred \tab number of predictions\cr
#'   rho \tab correlation coefficient between observations and predictions\cr
#'   mae \tab mean absolute error\cr
#'   rmse \tab root mean square error\cr
#'   perc \tab percent correct sign\cr
#'   p_val \tab p-value that rho is significantly greater than 0 using Fisher's 
#' }
#' @examples 
#' compute_stats(rnorm(100), rnorm(100))
#' @export
NULL
