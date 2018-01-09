#' @name rEDM
#' @docType package
#' @title Applications of empirical dynamic modeling from time series.
#' @author \strong{Maintainer}: Hao Ye
#' 
#'   \strong{Authors}: Adam Clark, Ethan Deyle, Steve Munch
#'   
#'   \strong{Contributors}: Oliver Keyes, Jun Cai, Ethan White, Jane Cowles, 
#'     James Stagge, Yair Daon, Andrew Edwards, George Sugihara
#' @description The \pkg{rEDM} package is a new implementation of EDM algorithms 
#'   based on research software previously developed for internal use in the 
#'   Sugihara Lab (UCSD/SIO). Contains C++ compiled objects that use time delay 
#'   embedding to perform state-space reconstruction and nonlinear forecasting 
#'   and an R interface to those objects using \pkg{Rcpp}. It supports both the 
#'   simplex projection method from Sugihara & May (1990) 
#'   <DOI:10.1038/344734a0> and the S-map algorithm in Sugihara (1994) 
#'   <DOI:10.1098/rsta.1994.0106>. In addition, this package implements 
#'   convergent cross mapping as described in Sugihara et al. (2012) 
#'   <DOI:10.1126/science.1227079> and multiview embedding as described in Ye & 
#'   Sugihara (2016) <DOI:10.1126/science.aag0863>.
#' @details This package is divided into a set of main functions to perform 
#'   various analyses, as well as helper functions that perform minor tasks, 
#'   such as generate data, processing output, and wrapper functions.
#' 
#' \strong{Main Functions}: 
#'   \itemize{
#'     \item \code{\link{simplex}} - simplex projection for univariate forecasting
#'     \item \code{\link{s_map}} - S-maps for univariate forecasting
#'     \item \code{\link{block_lnlp}} - simplex or S-map forecasting with a generic reconstructed state-space
#'     \item \code{\link{ccm}} - convergent cross mapping (causal inference)
#'     \item \code{\link{multiview}} - multi-model approach to forecasting
#'     \item \code{\link{tde_gp}} - Gaussian Processes for univariate forecasting
#'     \item \code{\link{block_gp}} - Gaussian Processes with a generic reconstructed state-space
#'   }
#' \strong{Helper Functions}: 
#'   \itemize{
#'     \item \code{\link{compute_stats}} - compute forecast skill metrics
#'     \item \code{\link{ccm_means}} - aggregate output of \code{\link{ccm}} by library size (`lib_size`)
#'     \item \code{\link{make_surrogate_data}} - generate surrogate time series
#'     \item \code{\link{test_nonlinearity}} - test for nonlinearity using surrogate time series
#'   }
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
#' @description Time series of Pacific sardine landings (CA), Northern anchovy 
#'   landings (CA), and sea-surface temperature (3-year average) at the SIO 
#'   pier and Newport pier
NULL

#' @name tentmap_del
#' @docType data
#' @title Time series for a tent map with mu = 2.
#' @author Hao Ye
#' @description First-differenced time series generated from the tent map 
#'   recurrence relation with mu = 2.
NULL

#' @name two_species_model
#' @docType data
#' @title Time series for a two-species coupled model.
#' @author Hao Ye
#' @description Time series generated from a discrete-time coupled 
#'   Lotka-Volterra model exhibiting chaotic dynamics.
NULL

#' @name block_3sp
#' @docType data
#' @title Time series for a three-species coupled model.
#' @author Hao Ye
#' @description Time series generated from a discrete-time coupled 
#'   Lotka-Volterra model exhibiting chaotic dynamics.
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
#' @description Experiment 120, the "Big Biodiversity" experiment at Cedar 
#'   Creek LTER. This experiment is the longest running randomized test for the 
#'   effects of plant diversity on ecosystem functions.
NULL

#' @name thrips_block
#' @docType data
#' @title Apple-blossom Thrips time series
#' @author ****
#' @description Seasonal outbreaks of Thrips imaginis.
NULL

#' @name e120_invnit16
#' @docType data
#' @title Biodiversity data at the Cedar Creek LTER
#' @author ****
#' @description Experiment 120, the "Big Biodiversity" experiment at Cedar 
#'   Creek LTER. This experiment is the longest running randomized test for the 
#'   effects of plant diversity on ecosystem functions.
NULL

#' @name compute_stats
#' 
#' @title Compute performance metrics for predictions
#'
#' @description Computes the rho, MAE, RMSE, perc, and p-val performance 
#'   metrics using the compiled C++ function
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
#' 
NULL
