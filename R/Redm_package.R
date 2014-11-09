#' @name Redm-package
#' @aliases Redm
#' @docType package
#' @title Applications of empirical dynamic modeling from time series.
#' @author Hao Ye
#' @description The Redm package provides an interface from R to C++ compiled objects 
#' that use time delay embedding to perform state-space reconstruction and nonlinear 
#' forecasting.
#' @keywords package
NULL

#' @name tentmap_1d
#' @docType data
#' @title Time series for a tent map with mu = 2.
#' @author Hao Ye
#' @description Time series generated from the tent map recurrence relation with mu = 2.
NULL

#' @name two_species_model
#' @docType data
#' @title Time series for a two-species coupled model.
#' @author Hao Ye
#' @description Time series generated from a discrete-time coupled Lotka-Volterra 
#' model exhibiting chaotic dynamics.
NULL

#' @name sardine_anchovy_sst
#' @docType data
#' @title Time series for the California Current Anchovy-Sardine-SST system
#' @author Hao Ye
#' @description Time series of Pacific sardine landings (CA), Northern anchovy landings (CA), 
#' and sea-surface temperature (3-year average) at the SIO pier and Newport pier
NULL

#' @name sockeye_data
#' @docType data
#' @title Biological time series for Fraser River sockeye salmon
#' @author Hao Ye
#' @description Time series of recruitment, juveniles, and effective female 
#' spawners for 19 stocks of sockeye salmon in the Fraser River system.
NULL

#' @name environmental_data
#' @docType data
#' @title Physical time series for Fraser River sockeye salmon
#' @author Hao Ye
#' @description Time series of Fraser River discharge, sea-surface temperature at 
#' 2 British Columbia lighthouses, and the Pacific Decadal Oscillation to be used 
#' as environmental covariates of Fraser River sockeye salmon.
NULL

#' @name Rcpp_LNLP-class
#' @docType class
#' @title S4 class for Rcpp compiled object, "LNLP"
NULL

#' @name LNLP
#' @title C++ compiled object for univariate forecasting.
#' @author Hao Ye
#' @seealso \code{\link{simplex}} and \code{\link{s_map}} for R wrappers.
NULL

#' @name Rcpp_BlockLNLP-class
#' @docType class
#' @title S4 class for Rcpp compiled object, "BlockLNLP"
NULL

#' @name BlockLNLP
#' @title C++ compiled object for multivariate forecasting.
#' @author Hao Ye
#' @seealso \code{\link{block_lnlp}} for the R wrapper.
NULL

#' @name Xmap
#' @title C++ compiled object for convergent cross mapping.
#' @author Hao Ye
#' @seealso \code{\link{ccm}} for the R wrapper.
NULL
