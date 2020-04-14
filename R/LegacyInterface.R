
source("R/EDM.R")

#--------------------------------------------------------------------------
# Attempt to mimic the interface of rEDM 0.7
# 
# NOTES
#   make_block() : adds columns parameter     # !!! API Change !!!
#   ccm()        : lib_sizes = c(10, 80, 10), # !!! API Change !!!
#   norm = 2     : others not supported
#   epsilon      : not supported
#   tau          : mathematically correct tau
# 
# Not implemented:
#   block_gp() tde_gp() ccm_means() test_nonlinearity() 
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
# ccm() 
# Wrapper for CCM() : returns data.frame or list if stats_only = FALSE
#--------------------------------------------------------------------------
ccm = function( block,
                lib               = NULL,
                pred              = NULL,
                norm              = 2,
                E                 = 1,
                tau               = -1,
                tp                = 0,
                num_neighbors     = "e+1",
                lib_sizes         = c(10, 75, 5), # !!! API Change !!!
                random_libs       = TRUE,
                num_samples       = 100,
                replace           = FALSE,
                lib_column        = 1,
                target_column     = 2,
                first_column_time = FALSE,
                RNGseed           = NULL,
                exclusion_radius  = NULL,
                epsilon           = NULL,
                stats_only        = TRUE,
                silent            = TRUE )
{
  verbose     = ! silent
  includeData = ! stats_only
  
  if ( norm != 2 ) {
    stop( "ccm(): L2-norm is the only metric currently available." )
  }
  if ( ! is.null( epsilon ) ) {
    stop( "ccm(): epsilon exlcusion not available." )
  }

  #-----------------------------------------------------------------
  # block : either a vector to be used as the time series, or
  #         data.frame or matrix where each column is a time series
  #-----------------------------------------------------------------
  if ( is.null( dim( block ) ) ) {
    # Not a data.frame or matrix, not much sense since E must be 1.
    # Create a time vector and data.frame, play God and make names...
    dataFrame = data.frame( Index = seq( 1:length(block) ),
                            Data  = block )
    columns = "Data"
    target  = "Data"
  }
  else if ( ncol( block ) >= 2 ) {
    if ( first_column_time ) {
      dataFrame = block
    }
    else {
      # Create a time vector and data.frame
      Index     = seq( 1:nrow(block) )
      dataFrame = data.frame( Index = Index, cbind( block ) )
      first_column_time = TRUE # Argh... Bogus
    }
    
    if ( is.numeric( lib_column ) ) {
      lib_column = lib_column + 1
      columns    = names( dataFrame )[ lib_column ]
    }
    else {
      columns = lib_column
    }
    
    if ( is.numeric( target_column ) ) {
      target_column = target_column + 1
      target        = names( dataFrame )[ target_column ]
    }
    else {
      target = target_column
    }
  }
  #-----------------------------------------------------------------

  if ( is.null( lib ) ) {
    # No lib specified.
    lib = c( 1, nrow(dataFrame) )
  }
  if ( is.null( pred ) ) {
    pred = lib
  }
  
  if ( is.null( exclusion_radius ) ) {
    exclusionRadius = 0
  }
  else {
    exclusionRadius = exclusion_radius
  }

  if ( "e+1"   %in% num_neighbors || "E+1"   %in% num_neighbors ||
       "e + 1" %in% num_neighbors || "E + 1" %in% num_neighbors ) {
    knn = 0 # cppEDM Simplex() will set knn = E + 1
  }
  else {
    knn = num_neighbors
  }

  if ( is.null( RNGseed ) ) {
    seed = 0
  }
  else {
    seed = RNGseed
  }

  ccmReturn = CCM( pathIn       = "./",
                   dataFile     = "",
                   dataFrame    = dataFrame,
                   pathOut      = "./",
                   predictFile  = "",
                   E            = E, 
                   Tp           = tp,
                   knn          = knn,
                   tau          = tau,
                   columns      = columns,
                   target       = target,
                   libSizes     = lib_sizes,
                   sample       = num_samples,
                   random       = random_libs,
                   replacement  = replace,
                   seed         = seed,
                   includeData  = includeData,
                   verbose      = verbose,
                   showPlot     = FALSE )
  
  # Add additional fields to the return
  if ( knn == 0 ) { knn = E + 1 }
  if ( includeData ) {
    # ccmReturn is a list
    N = nrow( ccmReturn $ LibMeans )
    ccmReturn $ LibMeans $ E   = rep( E,   N )
    ccmReturn $ LibMeans $ tau = rep( tau, N )
    ccmReturn $ LibMeans $ tp  = rep( tp,  N )
    ccmReturn $ LibMeans $ nn  = rep( knn, N )
    
    # Add lib and target column names to the PredictStats
    ccmReturn $ CCM1_PredictStat $ lib    = rep( columns, N )
    ccmReturn $ CCM1_PredictStat $ target = rep( target,  N )
    ccmReturn $ CCM2_PredictStat $ target = rep( columns, N )
    ccmReturn $ CCM2_PredictStat $ lib    = rep( target,  N )
  }
  else {
    # ccmReturn is data.frame
    N = nrow( ccmReturn )
    ccmReturn $ E   = rep( E,   N )
    ccmReturn $ tau = rep( tau, N )
    ccmReturn $ tp  = rep( tp,  N )
    ccmReturn $ nn  = rep( knn, N )
  }
  
  return( ccmReturn )
}

#--------------------------------------------------------------------------
# block_lnlp()  Presumes the data block is a multivariate state-space
#               No embedding is performed. 
# 
# Note this function has polymorphic return objects
# stats_only TRUE             : data.frame of stats
# stats_only FALSE            : list "stats", "model_output"
# save_smap_coefficients TRUE : list "stats", "model_output", "coef", "cov"
# "model_output" "coef" "cov" are named lists of data.frames for each theta
#--------------------------------------------------------------------------
block_lnlp = function(
  block,
  lib               = NULL,
  pred              = NULL,
  norm              = 2,
  method            = c("simplex", "s-map"),
  tp                = 1,
  num_neighbors     = switch(match.arg(method),
                             simplex = "e+1", `s-map` = 0),
  columns           = NULL,
  target_column     = 1,
  stats_only        = TRUE,
  first_column_time = FALSE,
  exclusion_radius  = NULL,
  epsilon           = NULL,
  theta             = NULL,
  silent            = TRUE,
  save_smap_coefficients = FALSE )
{
  verbose = ! silent
  
  if ( norm != 2 ) {
    stop( "block_lnlp(): L2-norm is the only metric currently available." )
  }
  if ( ! is.null( epsilon ) ) {
    stop( "block_lnlp(): epsilon exlcusion not available." )
  }

  #-----------------------------------------------------------------
  # block : either a vector to be used as the time series, or
  #         data.frame or matrix where each column is a time series
  #
  # Default parameters are set so that passing a vector as the only
  # argument will use that vector to predict itself one time step ahead.
  # If a matrix or data.frame is given as the only argument, the first
  # column will be predicted (one time step ahead), using the remaining
  # columns as the embedding. If the first column is not a time vector,
  # 1:NROW will be used as time values.
  #-----------------------------------------------------------------
  if ( is.null( dim( block ) ) ) {
    # Not a data.frame or matrix, not much sense since E must be 1.
    # Create a time vector and data.frame, play God and make names...
    dataFrame = data.frame( Index = seq( 1:length(block) ),
                            Data  = block )
    columns = "Data"
    target  = "Data"
  }
  else if ( ncol( block ) >= 2 ) {
    if ( first_column_time ) {
      dataFrame = block
    }
    else {
      # Create a time vector and data.frame
      Index     = seq( 1:nrow(block) )
      dataFrame = data.frame( Index = Index, cbind( block ) )
      first_column_time = TRUE # Argh... Bogus
    }
    
    if ( is.numeric( target_column ) ) {
      target_column = target_column + 1 # first col is time/index
      target        = names( dataFrame )[ target_column ]
    }
    else {
      target = target_column # Modern user... provided a column name
    }

    if ( is.null( columns ) ) {
      # presume first data column is the target, all others the columns
      columns = names( dataFrame )[ 3 : ncol(dataFrame) ]
    }
  }
  #-----------------------------------------------------------------

  # Much better to force the user to specify E... So... let's make a guess
  if ( ! is.character( columns ) || length( columns ) > 1 ) {
    # It's not a character string (is it?) length should work
    E = length( columns )
  }
  else {
    # Maybe it's a string with multiple columns... (? who knows...)
    E = length( strsplit( trimws( columns ), "\\s+" )[[1]] )
  }
  
  if ( verbose ) {
    print( paste( 'block_lnlp(): Using target', target,
                  'columns', FlattenToString( columns ), 'E =', E ) )
  }

  if ( is.null( lib ) ) {
    # No lib specified.
    lib = c( 1, nrow(dataFrame) )
  }
  if ( is.null( pred ) ) {
    pred = lib
  }
  
  if ( is.null( exclusion_radius ) ) {
    exclusionRadius = 0
  }
  else {
    exclusionRadius = exclusion_radius
  }

  #--------------------------------------------------------------------
  # Simplex
  #--------------------------------------------------------------------
  if ( 'simplex' %in% method ) {
    if ( "e+1"   %in% num_neighbors || "E+1"   %in% num_neighbors ||
         "e + 1" %in% num_neighbors || "E + 1" %in% num_neighbors ) {
      knn = 0 # cppEDM Simplex() will set knn = E + 1
    }
    else {
      knn = num_neighbors
    }
    
    smplx = Simplex( dataFrame       = dataFrame,
                     pathOut         = "./",
                     predictFile     = "",
                     lib             = lib,
                     pred            = pred,
                     E               = E,
                     Tp              = tp,
                     knn             = knn,
                     tau             = -1,
                     exclusionRadius = exclusionRadius,
                     columns         = columns,
                     target          = target, 
                     embedded        = TRUE,
                     verbose         = verbose,
                     const_pred      = TRUE,
                     showPlot        = FALSE )
    
    if ( knn == 0 ) {
      knn = E + 1 # As set in cppEDM Simplex()
    }

    # First column is cols
    stats = data.frame( cols = FlattenToString( columns ) )
    stats = cbind( stats, ComputeStats( list( smplx ), E, 0, tp, knn, NULL ) )

    if ( stats_only ) {
      return.object = stats
    }
    else {
      # list with "stats_only" and "model_results"
      # model_results is a list of data.frames for each E
      return.object = list( stats = stats, model_output = smplx )
    }
  }
  #--------------------------------------------------------------------
  # SMap
  #--------------------------------------------------------------------
  else if ( 's-map' %in% method ) {

    if ( is.character( num_neighbors ) ) {
      knn = 0
    }
    else {
      knn = num_neighbors
    }

    if ( is.null( theta ) ) {
      theta = c(0, 1e-04, 3e-04, 0.001, 0.003, 0.01, 0.03,
                0.1, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8)
    }
    
    smapList = list()

    for ( i in 1:length(theta) ) {
      theta.i = theta[ i ]
      
      smapList[[ i ]] = SMap( pathIn       = "./",
                              dataFile     = "",
                              dataFrame    = dataFrame,
                              pathOut      = "./",
                              predictFile  = "",
                              lib          = lib,
                              pred         = pred,
                              E            = E, 
                              Tp           = tp,
                              knn          = knn,
                              tau          = -1,
                              theta        = theta.i,
                              exclusionRadius = exclusionRadius,
                              columns      = columns,
                              target       = target,
                              smapFile     = "",
                              jacobians    = "",
                              embedded     = TRUE,
                              const_pred   = TRUE,
                              verbose      = verbose,
                              showPlot     = FALSE )
    }  
  
    names( smapList ) = paste0( "theta", theta )

    smapListPred = lapply( smapList, function(L){ L $ predictions } )
    
    # First column is cols
    stats = data.frame( cols = columns )
    stats = cbind( stats, ComputeStats( smapListPred, E, 0, tp, knn, theta ) )

    if ( stats_only ) {
      return.object = stats
    }
    else {
      # list with "stats_only" and "model_results"
      # model_results is a list of data.frames for each E
      return.object = list( stats = stats, model_output = smapListPred )
    }

    if ( save_smap_coefficients ) {
      smapListCoef = lapply( smapList, function(L){ L $ coefficients } )
      smapListCov  = lapply( smapList,
                             function(L){ cols = ncol( L $ coefficients );
                                          cov( L $ coefficients[,2:cols],
                                               use = 'complete.obs') } )
      return.object[[ "smap_coefficients" ]]             = smapListCoef
      return.object[[ "smap_coefficient_covariances"  ]] = smapListCov
    }
  }
  #---------------------------------------------------------------------------
  else {
    stop( paste( 'block_lnlp(): Invalid method:', method ) )
  }

  return( return.object )
}

#---------------------------------------------------------------------------
# s_map()  Presumes the time_series will be embedded to E.
# 
# Note this function has polymorphic return objects
# stats_only TRUE             : data.frame of stats
# stats_only FALSE            : list "stats", "model_output"
# save_smap_coefficients TRUE : list "stats", "model_output", "coef", "cov"
# "model_output" "coef" "cov" are named lists of data.frames for each theta
#---------------------------------------------------------------------------
s_map = function(
  time_series,
  lib                    = NULL,
  pred                   = NULL,
  norm                   = 2,
  E                      = 1,
  tau                    = -1,
  tp                     = 1,
  num_neighbors          = 0,
  theta                  = NULL,
  stats_only             = TRUE,
  exclusion_radius       = NULL,
  epsilon                = NULL,
  silent                 = TRUE,
  save_smap_coefficients = FALSE )
{
  verbose = ! silent
  
  if ( norm != 2 ) {
    stop( "s_map(): L2-norm is the only metric currently available." )
  }
  if ( ! is.null( epsilon ) ) {
    stop( "s_map(): epsilon exlcusion not available." )
  }

  if ( is.null( theta ) ) {
    theta = c(0, 1e-04, 3e-04, 0.001, 0.003, 0.01, 0.03,
              0.1, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8)
  }

  dataList  = DataFrameFromTimeseries( time_series )
  dataFrame = dataList[[ 'dataFrame' ]]
  columns   = dataList[[ 'columns'   ]]
  target    = dataList[[ 'target'    ]]

  if ( verbose ) {
    print( paste( 's_map(): Using target', target ) )
  }

  if ( is.null( lib ) ) {
    # No lib specified.
    lib = c( 1, nrow(dataFrame) - abs( tau ) * ( E - 1 ) )
  }
  if ( is.null( pred ) ) {
    pred = lib
  }
  
  knn = num_neighbors

  if ( is.null( exclusion_radius ) ) {
    exclusionRadius = 0
  }
  else {
    exclusionRadius = exclusion_radius
  }

  # Would be much more efficient to use PredictNonlinear() to find theta, 
  # then call SMap(), but the legacy returns all (heavily redundant)
  # results. 
  # 
  # Also better to do the embedding once, then call SMap() embedded = TRUE
  # 
  # JP : use mclapply...
  # 
  # Compute SMap for all theta
  smapList = list()

  for ( i in 1:length(theta) ) {
    theta.i = theta[ i ]
    # SMap() : list [[ "predictions",  "coefficients" ]]
    # predictions: "Index" "Observations" "Predictions" "Const_Predictions"
    # coefficients: Index" "C0" "C1"...
    smapList[[ i ]] = SMap( pathIn       = "./",
                            dataFile     = "",
                            dataFrame    = dataFrame,
                            pathOut      = "./",
                            predictFile  = "",
                            lib          = lib,
                            pred         = pred,
                            E            = E, 
                            Tp           = tp,
                            knn          = knn,
                            tau          = tau,
                            theta        = theta.i,
                            exclusionRadius = exclusionRadius,
                            columns      = columns,
                            target       = target,
                            smapFile     = "",
                            jacobians    = "",
                            embedded     = FALSE,
                            const_pred   = TRUE,
                            verbose      = verbose,
                            showPlot     = FALSE )
  }
  
  names( smapList ) = paste0( "theta", theta )

  smapListPred = lapply( smapList, function(L){ L $ predictions } )
  stats        = ComputeStats( smapListPred, E, tau, tp, knn, theta )

  if ( stats_only ) {
    return.object = stats # data.frame
  }
  
  if ( save_smap_coefficients ) {
    # return.object is a list, add coefficients & covariance
    return.object = list( stats = stats )
    
    smapListCoef = lapply( smapList, function(L){ L $ coefficients } )
    smapListCov  = lapply( smapList,
                           function(L){ cols = ncol( L $ coefficients );
                                        cov( L $ coefficients[,2:cols],
                                             use = 'complete.obs') } )
    return.object[[ "smap_coefficients" ]]             = smapListCoef
    return.object[[ "smap_coefficient_covariances"  ]] = smapListCov
  }
  
  if ( ! stats_only ) {
    # ALready a list with "stats"
    # Add model_results as a list of data.frames for each E
    return.object[[ "model_output" ]] = smapListPred
  }

  return( return.object )
}

#------------------------------------------------------------------------
# simplex()  Presumes the time_series will be embedded to E.
# 
# Note this function has polymorphic return objects
# if stats_only TRUE  : data.frame of stats
# if stats_only FALSE : list of "stats", "model_output"
# "model_output" is a named list of data.frames for each E
#------------------------------------------------------------------------
simplex = function(
  time_series,
  lib              = NULL,
  pred             = NULL,
  norm             = 2,
  E                = 1:10,
  tau              = -1,
  tp               = 1,
  num_neighbors    = "e+1",
  stats_only       = TRUE,
  exclusion_radius = NULL,
  epsilon          = NULL,
  silent           = TRUE )
{
  verbose = ! silent
  
  if ( norm != 2 ) {
    stop( "simplex(): L2-norm is the only metric currently available." )
  }
  if ( ! is.null( epsilon ) ) {
    stop( "simplex(): epsilon exlcusion not available." )
  }

  dataList  = DataFrameFromTimeseries( time_series )
  dataFrame = dataList[[ 'dataFrame' ]]
  columns   = dataList[[ 'columns'   ]]
  target    = dataList[[ 'target'    ]]

  if ( verbose ) {
    print( paste( 'simplex(): Using target', target ) )
  }

  if ( is.null( lib ) ) {
    # No lib specified.  Use max(E) for uniform lib size across all E's
    lib = c( 1, nrow(dataFrame) - abs( tau ) * ( max( E ) - 1 ) )
  }
  if ( is.null( pred ) ) {
    pred = lib
  }
  
  if ( "e+1"   %in% num_neighbors || "E+1"   %in% num_neighbors ||
       "e + 1" %in% num_neighbors || "E + 1" %in% num_neighbors ) {
    knn = 0 # cppEDM Simplex() will set knn = E + 1
  }
  else {
    knn = num_neighbors
  }

  if ( is.null( exclusion_radius ) ) {
    exclusionRadius = 0
  }
  else {
    exclusionRadius = exclusion_radius
  }

  # Would be much more efficient to use EmbedDimension() to find E_opt, 
  # then call Simplex(), but the legacy returns all (heavily redundant)
  # simplex results. 
  # 
  # JP : use mclapply...
  # 
  # Compute Simplex for all E
  simplexList = list()
  
  for ( E.i in E ) {
    # data.frame: Time Observations Predictions Const_Predictions
    simplexList[[ E.i ]] = Simplex( dataFrame       = dataFrame,
                                    pathOut         = "./",
                                    predictFile     = "",
                                    lib             = lib,
                                    pred            = pred,
                                    E               = E.i, 
                                    Tp              = tp,
                                    knn             = knn,
                                    tau             = tau,
                                    exclusionRadius = exclusionRadius,
                                    columns         = columns,
                                    target          = target, 
                                    embedded        = FALSE,
                                    verbose         = verbose,
                                    const_pred      = TRUE,
                                    showPlot        = FALSE )
  }
  
  names( simplexList ) = paste0( "E", E )

  if ( knn == 0 ) {
    knn.E = E + 1
  }
  else {
    knn.E = rep( num_neighbors, length( simplexList ) )
  }

  stats = ComputeStats( simplexList, E, tau, tp, knn.E, NULL )

  if ( stats_only ) {
    return.object = stats
  }
  else {
    # list with "stats_only" and "model_results"
    # model_results is a list of data.frames for each E
    return.object = list( stats = stats, model_output = simplexList )
  }
  
  return( return.object )
}

#------------------------------------------------------------------------
# 
#------------------------------------------------------------------------
multiview = function( block,
                      lib               = NULL,
                      pred              = NULL,
                      norm              = 2,
                      E                 = 1,
                      tau               = -1,
                      tp                = 1,
                      max_lag           = 3, 
                      num_neighbors     = "e+1",
                      k                 = "sqrt",
                      na.rm             = FALSE,
                      target_column     = 1, 
                      stats_only        = TRUE,
                      save_lagged_block = FALSE,
                      first_column_time = FALSE, 
                      exclusion_radius  = NULL,
                      silent            = FALSE )
{
  verbose = ! silent

  if ( norm != 2 ) {
    stop( "multiview(): L2-norm is the only metric currently available." )
  }

  #-----------------------------------------------------------------
  # block : either a vector to be used as the time series, or
  #         data.frame or matrix where each column is a time series
  #-----------------------------------------------------------------
  if ( is.null( dim( block ) ) ) {
    # Not a data.frame or matrix, not much sense since E must be 1.
    # Create a time vector and data.frame, play God and make names...
    dataFrame = data.frame( Index = seq( 1:length(block) ),
                            Data  = block )
    columns = "Data"
    target  = "Data"
  }
  else if ( ncol( block ) >= 2 ) {
    if ( first_column_time ) {
      dataFrame = block
    }
    else {
      # Create a time vector and data.frame
      Index     = seq( 1:nrow(block) )
      dataFrame = data.frame( Index = Index, cbind( block ) )
      first_column_time = TRUE # Argh... Bogus
    }
    
    columns = names( dataFrame )[ 2 : ncol( dataFrame ) ] # ignore first col time
    columns = paste( columns, collapse = " " )
    
    if ( is.numeric( target_column ) ) {
      target_column = target_column + 1
      target        = names( dataFrame )[ target_column ]
    }
    else {
      target = target_column
    }
  }
  #-----------------------------------------------------------------

  if ( verbose ) {
    print( paste( 'columns:', columns ) )
    print( paste( 'target:',  target  ) )
  }
  
  if ( is.null( lib ) ) {
    lib = c( 1, floor(NROW(block)/2) )
  }
  if ( is.null( pred ) ) {
    pred = c( floor(NROW(block)/2) + 1, NROW(block) )
  }
  
  if ( "e+1"   %in% num_neighbors || "E+1"   %in% num_neighbors ||
       "e + 1" %in% num_neighbors || "E + 1" %in% num_neighbors ) {
    knn = 0 # cppEDM Simplex() will set knn = E + 1
  }
  else {
    knn = num_neighbors
  }
  
  if ( is.null( exclusion_radius ) ) {
    exclusionRadius = 0
  }
  else {
    exclusionRadius = exclusion_radius
  }

  if ( 'sqrt' %in% k || k == 0 ) {
    multiview = 0
  }
  else {
    multiview = k
  }
  
  mv = Multiview( pathIn          = "./",
                  dataFile        = "",
                  dataFrame       = dataFrame,
                  pathOut         = "./",
                  predictFile     = "",
                  lib             = lib,
                  pred            = pred,
                  D               = E,       # Nice: E is not embedding dim
                  E               = max_lag, # Nice: E is not embedding dim
                  Tp              = tp,
                  knn             = knn,
                  tau             = tau,
                  columns         = columns,
                  target          = target,
                  multiview       = multiview,
                  exclusionRadius = exclusionRadius,
                  verbose         = FALSE,
                  numThreads      = 4,
                  showPlot        = FALSE )

  return( mv )
}

#------------------------------------------------------------------------
# 
#------------------------------------------------------------------------
make_block = function( block,
                       columns         = NULL,
                       t               = NULL,
                       max_lag         = 3,
                       tau             = -1,
                       lib             = NULL,
                       restrict_to_lib = TRUE ) {

  if ( ! is.null( lib ) ) {
    stop( "make_block(): lib selection not available." )
  }
  if ( ! is.null( t ) ) {
    stop( "make_block(): t not used." )
  }

  if ( is.null( columns ) ) {
    # Embed assumes column 1 is Time and is ignored...
    columns = names( block )[ 2:ncol(block) ]

    print( "make_block(): Ignoring first column of block, assumed to be time." )
  }

  print( columns )
  
  embed.block = Embed( path      = "./",
                       dataFile  = "",
                       dataFrame = block,
                       E         = max_lag,
                       tau       = tau, 
                       columns   = columns,
                       verbose   = FALSE )

  return( embed.block )
}

#------------------------------------------------------------------------
# 
#------------------------------------------------------------------------
compute_stats = function( observed, predicted ) {
  err = ComputeError( observed, predicted )

  num_pred = length( which( ! is.na( predicted ) ) )
  rho      = err $ rho
  mae      = err $ MAE
  rmse     = err $ RMSE
  perc     = sum( abs( sign(observed) + sign(predicted) ), na.rm = TRUE ) /
             2 / length( which( ! is.na(observed) ) )
  p_val    = max(1E-10, pnorm( atanh(rho), 0.0, 1 / sqrt(num_pred),FALSE,FALSE))

  stats = data.frame( num_pred = num_pred, rho = rho, mae = mae,
                      rmse = rmse, perc = perc, p_val = p_val )

  return( stats )
}

#------------------------------------------------------------------------
# 
#------------------------------------------------------------------------
DataFrameFromTimeseries = function( time_series ) {
  # time_series
  # either a vector to be used as the time series, or a data.frame or matrix
  # with at least 2 columns (in which case the first column will be used as
  # the time index, and the second column as the time series)
  if ( is.null( dim( time_series ) ) ) {
    # Not a data.frame or matrix
    # Create a time vector and data.frame, play God and make names...
    dataFrame = data.frame( Index = seq( 1:length(time_series) ),
                            Data  = time_series )
    columns = "Data"
    target  = "Data"
  }
  else if ( ncol( time_series ) >= 2 ) {
    dataFrame = time_series
    columns   = names( dataFrame )[2]
    target    = names( dataFrame )[2]
  }
  return( list( dataFrame = dataFrame, columns = columns, target = target ) )
}

#------------------------------------------------------------------------
# Setup stats_only = TRUE data.frame
#------------------------------------------------------------------------
ComputeStats = function( PredictList, E, tau, tp, knn.E, theta ) {
  #----------------------------------------------------------------------
  # rEDM 0.7 simplex stats_only = TRUE : data.frame E rows x 16 columns
  # rEDM 0.7 s_map   stats_only = TRUE : data.frame theta rows x 17 cols
  #----------------------------------------------------------------------
  #  "E"                "tau"                 "tp"
  #  "nn"   ("theta")   "num_pred"            "rho"
  #  "mae"              "rmse"                "perc"
  #  "p_val"            "const_pred_num_pred" "const_pred_rho"
  #  "const_pred_mae"   "const_pred_rmse"     "const_pred_perc"
  #  "const_p_val"
  #---------------------------------------------------------------------
  N = length( PredictList )
  
  # Here's the redundant part...
  stats = data.frame( E = E, tau = rep( tau, N ),
                      tp = rep( tp, N ), nn = knn.E )

  if ( ! is.null( theta ) ) {
    stats $ theta = theta
  }

  numPred      = sapply( PredictList, PredictN,             simplify = TRUE )
  errors       = sapply( PredictList, PredictError,         simplify = TRUE )
  constErrors  = sapply( PredictList, PredictConstError,    simplify = TRUE )
  percent      = sapply( PredictList, PercentSameSign,      simplify = TRUE ) 
  numConstPred = sapply( PredictList, PredictConstN,        simplify = TRUE )
  constPercent = sapply( PredictList, PercentConstSameSign, simplify = TRUE )
  pvals        = sapply( PredictList, PValue,               simplify = TRUE )
  constpvals   = sapply( PredictList, PValueConst,          simplify = TRUE )

  errors      = as.data.frame( t( errors ) )
  constErrors = as.data.frame( t( constErrors ) )

  stats $ num_pred            = numPred
  stats $ rho                 = errors $ rho
  stats $ mae                 = errors $ MAE
  stats $ rmse                = errors $ RMSE
  stats $ perc                = percent
  stats $ p_val               = pvals
  stats $ const_pred_num_pred = numConstPred
  stats $ const_pred_rho      = constErrors $ rho
  stats $ const_pred_mae      = constErrors $ MAE
  stats $ const_pred_rmse     = constErrors $ RMSE
  stats $ const_pred_perc     = constPercent
  stats $ const_p_val         = constpvals

  return( stats )
}
#------------------------------------------------------------------------
# sapply functions for ComputeStats
#------------------------------------------------------------------------
PValue = function( df ) {
  N.pred = length( which( !is.na( df $ Predictions ) ) )
  rho    = ComputeError( df $ Observations, df $ Predictions ) $ rho
  # pnorm( q, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE )
  max( 1E-10, pnorm( atanh(rho), 0.0, 1 / sqrt(N.pred), FALSE, FALSE ) )
}
PValueConst = function( df ) {
  N.pred = length( which( !is.na( df $ Const_Predictions ) ) )
  rho    = ComputeError( df $ Observations, df $ Const_Predictions ) $ rho
  max( 1E-10, pnorm( atanh(rho), 0.0, 1 / sqrt(N.pred), FALSE, FALSE ) )
}
PercentSameSign = function( df ) {
  # Ratio of observations to predictions with same sign
  o = df $ Observations
  p = df $ Predictions
  N = length( which( ! is.na(o) ) )
  sum( abs( sign(o) + sign(p) ), na.rm = TRUE ) / 2 / N
}
PercentConstSameSign = function( df ) {
  o = df $ Observations
  p = df $ Const_Predictions
  N = length( which( ! is.na(o) ) )
  sum( abs( sign(o) + sign(p) ), na.rm = TRUE ) / 2 / N
}
PredictError = function( df ) {
  ComputeError( df $ Observations, df $ Predictions )
}
PredictConstError = function( df ) {
  ComputeError( df $ Observations, df $ Const_Predictions )
}
PredictN = function( df ) {
  length( which( ! is.na( df $ Predictions ) ) )
}
PredictConstN = function( df ) {
  length( which( ! is.na( df $ Const_Predictions ) ) )
}

#--------------------------------------------------------------------------
# make_surrogate_data()  
#--------------------------------------------------------------------------
make_surrogate_data = function( ts,
                                method   = c("random_shuffle", "ebisuzaki",
                                             "seasonal" ),
                                num_surr = 100,
                                T_period = 1,
                                alpha    = 0 )
{
  data = SurrogateData( ts,
                        method   = method,
                        num_surr = num_surr,
                        T_period = T_period,
                        alpha    = alpha )
  return( data )
}

#--------------------------------------------------------------------------
# ccm_means()  
#--------------------------------------------------------------------------
ccm_means = function( df, FUN = mean, ... ) {
  stop( "ccm_means() deprecated. ccm() returns lib_sizes means." )
}

#--------------------------------------------------------------------------
# tde_gp()  
#--------------------------------------------------------------------------
tde_gp = function( time_series,
                   lib = c(1, NROW(time_series)),
                   pred = lib, 
                   E = 1:10,
                   tau = 1,
                   tp = 1,
                   phi = 0,
                   v_e = 0,
                   eta = 0,
                   fit_params = TRUE, 
                   stats_only = TRUE,
                   save_covariance_matrix = FALSE,
                   silent = FALSE, 
                  ... ) {
  stop( "tde_gp() not implemented." )
}

#--------------------------------------------------------------------------
# block_gp()  
#--------------------------------------------------------------------------
block_gp = function( block,
                     lib = c(1, NROW(block)),
                     pred = lib,
                     tp = 1, 
                     phi = 0,
                     v_e = 0,
                     eta = 0,
                     fit_params = TRUE,
                     columns = NULL, 
                     target_column = 1,
                     stats_only = TRUE,
                     save_covariance_matrix = FALSE, 
                     first_column_time = FALSE,
                     silent = FALSE, ...) {
  stop( "block_gp() not implemented." )
}

#--------------------------------------------------------------------------
# test_nonlinearity()  
#--------------------------------------------------------------------------
test_nonlinearity = function( ts,
                              method = "ebisuzaki",
                              num_surr = 200,
                              T_period = 1, 
                              E = 1,
                              ... ) {
  stop( "test_nonlinearity() not implemented." )
}
