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

#' @title Time series for the Paramecium-Didinium laboratory experiment
#' @author Veilleux
#' @description Time series of Paramecium and Didinium abundances (#/mL) from 
#'   an experiment by Veilleux (1979)
#' @format
#' \describe{
#'   \item{\code{time}}{time of sampling (in days)}
#'   \item{\code{paramecium}}{paramecium abundance (per mL)}
#'   \item{\code{didinium}}{didinium abundance (per mL)}
#' }
"paramecium_didinium"

#' @title Time series for the California Current Anchovy-Sardine-SST system
#' @author ****
#' @description Time series of Pacific sardine landings (CA), Northern anchovy 
#'   landings (CA), and sea-surface temperature (3-year average) at the SIO 
#'   pier and Newport pier
#' @format
#' \describe{
#'   \item{\code{year}}{year of measurement}
#'   \item{\code{anchovy}}{anchovy landings, scaled to mean = 0, sd = 1}
#'   \item{\code{sardine}}{sardine landings, scaled to mean = 0, sd = 1}
#'   \item{\code{sio_sst}}{3-year running average of sea surface temperature at 
#'     SIO pier, scaled to mean = 0, sd = 1}
#'   \item{\code{np_sst}}{3-year running average of sea surface temperature at 
#'     Newport pier, scaled to mean = 0, sd = 1}
#' }
"sardine_anchovy_sst"

#' @title Time series for a tent map with mu = 2.
#' @author Hao Ye
#' @description First-differenced time series generated from the tent map
#'   recurrence relation with mu = 2.
#' @format A numeric vector with 999 values
"tentmap_del"

#' @title Time series for a two-species coupled model.
#' @author Hao Ye
#' @description Time series generated from a discrete-time coupled 
#'   Lotka-Volterra model exhibiting chaotic dynamics.
#' @format A data frame with 1000 rows and 3 variables: 
#' \describe{
#'   \item{\code{time}}{time index (# of generations)}
#'   \item{\code{x}}{abundance of simulated species x}
#'   \item{\code{y}}{abundance of simulated species y}
#' }
"two_species_model"

#' @title Time series for a three-species coupled model.
#' @author Hao Ye
#' @description Time series generated from a discrete-time coupled 
#'   Lotka-Volterra model exhibiting chaotic dynamics.
#' @format A data frame with 200 rows and 10 variables: 
#' \describe{
#'   \item{\code{time}}{time index (# of generations)}
#'   \item{\code{x_t}}{abundance of simulated species $x$ at time $t$}
#'   \item{\code{x_t-1}}{abundance of simulated species $x$ at time $t-1$}
#'   \item{\code{x_t-2}}{abundance of simulated species $x$ at time $t-2$}
#'   \item{\code{y_t}}{abundance of simulated species $y$ at time $t$}
#'   \item{\code{y_t-1}}{abundance of simulated species $y$ at time $t-1$}
#'   \item{\code{y_t-2}}{abundance of simulated species $y$ at time $t-2$}
#'   \item{\code{z_t}}{abundance of simulated species $z$ at time $t$}
#'   \item{\code{z_t-1}}{abundance of simulated species $z$ at time $t-1$}
#'   \item{\code{z_t-2}}{abundance of simulated species $z$ at time $t-2$}
#' }
"block_3sp"

#' @title Time series for sockeye salmon returns.
#' @author Fisheries and Oceans Canada
#' @description Time series of sockeye salmon returns from the Fraser River in 
#'   British Columbia, Canada.
#' @format A data frame with 55 rows and 10 variables: 
#' \describe{
#'   \item{\code{year}}{year}
#'   \item{\code{Early_Stuart}}{number of sockeye returns for the Early Stuart 
#'     stock (millions of fish)}
#'   \item{\code{Late_Stuart}}{number of sockeye returns for the Late Stuart 
#'     stock (millions of fish)}
#'   \item{\code{Stellako}}{number of sockeye returns for the Stellako 
#'     stock (millions of fish)}
#'   \item{\code{Quesnel}}{number of sockeye returns for the Quesnel 
#'     stock (millions of fish)}
#'   \item{\code{Chilko}}{number of sockeye returns for the Chilko
#'     stock (millions of fish)}
#'   \item{\code{Seymour}}{number of sockeye returns for the Seymour
#'     stock (millions of fish)}
#'   \item{\code{Late_Shuswap}}{number of sockeye returns for the Late Shuswap
#'     stock (millions of fish)}
#'   \item{\code{Birkenhead}}{number of sockeye returns for the Birkenhead
#'     stock (millions of fish)}
#'   \item{\code{Weaver}}{number of sockeye returns for the Weaver
#'     stock (millions of fish)}
#' }
"sockeye_returns"

#' @title Apple-blossom Thrips time series
#' @author ****
#' @description Seasonal outbreaks of Thrips imaginis.
#' @format A data frame with 81 rows and 6 variables: 
#' \describe{
#'   \item{\code{Year}}{year}
#'   \item{\code{Month}}{month}
#'   \item{\code{Thrips_imaginis}}{counts of Thrips imaginis}
#'   \item{\code{maxT_degC}}{mean maximum dail temperature in each month 
#'     (Celsius)}
#'   \item{\code{Rain_mm}}{monthly rainfall (mm)}
#'   \item{\code{Season}}{annual sinusoidal that peaks in December}
#' }
"thrips_block"

#' @title Biodiversity data at the Cedar Creek LTER
#' @author ****
#' @description Experiment 120, the "Big Biodiversity" experiment at Cedar 
#'   Creek LTER. This experiment is the longest running randomized test for the 
#'   effects of plant diversity on ecosystem functions. This is a subset of the 
#'   data from the treatments with 16 planted species.
#' @format A data frame with 238 rows and 7 variables: 
#' \describe{
#'   \item{\code{Exp}}{experiment code (E120)}
#'   \item{\code{Year}}{year}
#'   \item{\code{Plot}}{plot identity}
#'   \item{\code{AbvBioAnnProd}}{annual above-ground productivity of planted 
#'     species (g/m^2)}
#'   \item{\code{noh020tot}}{soil nitrate levels in the top 20cm of soil (ug/kg 
#'     soil)}
#'   \item{\code{invrichness}}{species richness of unplanted species in the 
#'     plot, i.e. weeds}
#'   \item{\code{SummerPrecip.mm}}{annual precipitation from May to August (mm)}
#' }
"e120_invnit16"

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
