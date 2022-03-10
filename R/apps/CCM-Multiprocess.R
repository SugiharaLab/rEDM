
library( rEDM )

library( foreach )
library( doParallel )

#-------------------------------------------------------------------
# CCM for all dataFrame columns against target using foreach %dopar%
# Presumes first column is time/index, not processed
#-------------------------------------------------------------------
CCM_MP_Columns = function(
  dataFrame = NULL,
  target    = 'V5',
  libSizes  = '20 920 100',
  sample    = 10,
  E         = 5,
  Tp        = 0,
  cores     = 4  # CCM uses 2 cores, max is detectCores()/2 - 2
) {

  if ( is.null( dataFrame ) ) { dataFrame = Lorenz5D }

  registerDoParallel( cores = cores )

  dataCols = names( dataFrame )[ 2 : ncol( dataFrame ) ] # Skip first column

  # Parallel process columns using foreach ... %dopar%
  L = foreach ( col = iter( dataCols ) ) %dopar% {

    CCM( dataFrame = dataFrame,
         E         = E,
         Tp        = Tp,
         columns   = col,
         target    = target,
         libSizes  = libSizes,
         sample    = sample )
  }

  # Get names for the returned list L from the CCM data.frame
  keys = c()
  for ( cmap in L ) {
    keys = c( keys, names( cmap )[3] )
  }
  names( L ) = keys

  invisible( L )
}

#---------------------------------------------------------------------
# CCM for single columns : target over a set of libSizes using foreach %dopar%
# libSizeList is a list of partioned libSizes
# libSizesList elements can be any libSizes format used by CCM
#---------------------------------------------------------------------
CCM_MP_LibSizes = function(
  dataFrame    = NULL,
  columns      = 'V1',
  target       = 'V5',
  libSizesList = c( '20 50 70 100', '150 200 250 300', '400 500 600 700 900' ),
  sample       = 10,
  E            = 5,
  Tp           = 0,
  cores        = 4  # CCM uses 2 cores, max is detectCores()/2 - 2
 ) {

  if ( is.null( dataFrame ) ) { dataFrame = Lorenz5D }

  registerDoParallel( cores = cores )

  # Parallel process libSizesList using foreach ... %dopar%
  L = foreach ( libSize = iter( libSizesList ) ) %dopar% {

    CCM( dataFrame = dataFrame,
         E         = E,
         Tp        = Tp,
         columns   = columns,
         target    = target,
         libSizes  = libSize,
         sample    = sample )
  }

  # Set names
  names( L ) = libSizesList

  invisible( L )
}

#---------------------------------------------------------------------
# CCM for single columns : target over a set of libSizes using clusterApply
# libSizeList is a list of partioned libSizes
# libSizesList elements can be any libSizes format used by CCM
#
# DO NOT USE mclapply
#    From ?mclapply:
#    It is _strongly discouraged_ to use these functions with
#    multi-threaded libraries or packages (see ‘mcfork’ for more
#    details).  If in doubt, it is safer to use a non-FORK cluster
#    (see ‘makeCluster’, ‘clusterApply’).
#---------------------------------------------------------------------
CCM_MP_LibSizes_cluster = function(
  dataFrame    = NULL,
  columns      = 'V1',
  target       = 'V5',
  libSizesList = c( '20 30 40 50 60 70 80 90 100',
                    '120 150 200 250 300',
                    '400 500 600 700 900' ),
  sample       = 20,
  E            = 5,
  Tp           = 0,
  cores        = 4

) {

  if ( is.null( dataFrame ) ) { dataFrame = Lorenz5D }

  cl = makeCluster( cores )

  clusterExport( cl, list("CCM") )

  cmap = clusterApply( cl = cl, x = libSizesList, fun = CrossMapFunc,
                       dataFrame, E, Tp, columns, target, sample )

  stopCluster( cl )

  invisible( cmap )
}

#---------------------------------------------------------------------
# Call rEDM CCM on behalf of CCM_MP_LibSizes_cluster() clusterApply()
#---------------------------------------------------------------------
CrossMapFunc = function(
  libSizes, # First argument : from clusterApply( x = libSizesList )
  dataFrame, E, Tp, columns, target, sample
) {

  cm = CCM( dataFrame = dataFrame,
            E         = E,
            columns   = columns,
            target    = target,
            libSizes  = libSizes,
            sample    = sample )
}
