
source("R/EDM_AuxFuncs.R")

#------------------------------------------------------------------------
# Takens time-delay embedding on columnNames in dataFrame.
# Truncates the timeseries by tau * (E-1) rows.
#------------------------------------------------------------------------
MakeBlock = function( dataFrame,
                      E             = 0, 
                      tau           = -1,
                      columns       = c(),  # vector of strings
                      deletePartial = FALSE) {

  if ( ! isValidDataFrame( dataFrame ) ) {
    stop( "MakeBlock(): dataFrame argument is not valid data.frame." )
  }

  # Mapped to MakeBlock_rcpp() (Embed.cpp) in RcppEDMCommon.cpp
  block = RtoCpp_MakeBlock( dataFrame, E, tau, columns, deletePartial )

  return ( block )
}

#------------------------------------------------------------------------
# Takens time-delay embedding on path/file.
# Embed DataFrame columns (subset) in E dimensions.
# Calls MakeBlock() after validation and column subset selection.
#------------------------------------------------------------------------
Embed = function( path      = "./",
                  dataFile  = "",
                  dataFrame = NULL,
                  E         = 0, 
                  tau       = -1,
                  columns   = "",
                  verbose   = FALSE ) {

  if ( ! is.null( dataFrame ) ) {
    if ( ! isValidDataFrame( dataFrame ) ) {
      stop( "Embed(): dataFrame argument is not valid data.frame." )
    }
  }

  # If columns are vectors/list, convert to string for cppEDM
  # NOTE: columns joined on ',' to enable names with whitespace in cppEDM
  if ( ! is.character( columns ) || length( columns ) > 1 ) {
    columns = FlattenToString( columns, "," )
  }

  # Mapped to Embed_rcpp() (Embed.cpp) in RcppEDMCommon.cpp
  df = RtoCpp_Embed( path,
                     dataFile,
                     dataFrame,
                     E, 
                     tau,
                     columns,
                     verbose )

  return ( df )
}

#------------------------------------------------------------------------
#
#------------------------------------------------------------------------
Simplex = function( pathIn          = "./",
                    dataFile        = "",
                    dataFrame       = NULL,
                    pathOut         = "./",
                    predictFile     = "",
                    lib             = "",
                    pred            = "",
                    E               = 0, 
                    Tp              = 1,
                    knn             = 0,
                    tau             = -1,
                    exclusionRadius = 0,
                    columns         = "",
                    target          = "", 
                    embedded        = FALSE,
                    verbose         = FALSE,
                    validLib        = vector(),
                    generateSteps   = 0,
                    parameterList   = FALSE,
                    showPlot        = FALSE,
                    noTime          = FALSE ) {

  if ( noTime ) {
    if ( nchar( dataFile ) ) {
      # Read the dataFile into a data.frame, insert index
      # Override dataFile argument to "", assign dataFrame with index column
      #   Disable check.names that calls make.names() on column names.
      dataFrame = read.csv( paste( pathIn, dataFile, sep = '/' ),
                            as.is = TRUE, check.names = FALSE )
      dataFile  = ""
    }
    # Insert an index column as first column for cppEDM
    index     = seq( 1, nrow( dataFrame ) )
    dataFrame = cbind( index, dataFrame )
  }

  if ( ! ValidateDataFrame( pathIn, dataFile, dataFrame,
                            columns, target, noTime ) ) {
    stop( "Simplex(): dataFrame validation failed." )
  }

  # If lib, pred, columns are vectors/list, convert to string for cppEDM
  # NOTE: columns joined on ',' to enable names with whitespace in cppEDM
  if ( ! is.character( lib ) || length( lib ) > 1 ) {
    lib = FlattenToString( lib )
  }
  if ( ! is.character( pred ) || length( pred ) > 1 ) {
    pred = FlattenToString( pred )
  }
  if ( ! is.character( columns ) || length( columns ) > 1 ) {
    columns = FlattenToString( columns, "," )
  }
  if ( length( strsplit( target, ' ' )[[1]] ) > 1 ) {
    target = paste0( target, ',' ) # space in target: add , for cppEDM
  }

  # NOTE: Rcpp has a 20 argument limit!
  # Mapped to Simplex_rcpp() (Simplex.cpp) in RcppEDMCommon.cpp
  smplx = RtoCpp_Simplex( pathIn,
                          dataFile,
                          dataFrame,
                          pathOut,
                          predictFile,
                          lib,
                          pred,
                          E,
                          Tp,
                          knn,
                          tau,
                          exclusionRadius,
                          columns,
                          target,
                          embedded,
                          verbose,
                          validLib,
                          generateSteps,
                          parameterList )

  if ( showPlot ) {
    PlotObsPred( smplx[['predictions']], dataFile, E, Tp )
  }

  if ( parameterList ) {
    return( smplx )
  }
  else {
    # Return the exposed dataFrame
    return( smplx[['predictions']] )
  }
}

#------------------------------------------------------------------------
#
#------------------------------------------------------------------------
SMap = function( pathIn          = "./",
                 dataFile        = "",
                 dataFrame       = NULL,
                 # pathOut       = "./", # Rcpp 20 arg limit
                 # predictFile   = "",   # Rcpp 20 arg limit
                 lib             = "",
                 pred            = "",
                 E               = 0, 
                 Tp              = 1,
                 knn             = 0,
                 tau             = -1,
                 theta           = 0,
                 exclusionRadius = 0,
                 columns         = "",
                 target          = "",
                 # smapCoefFile  = "",   # Rcpp 20 arg limit
                 # smapSVFile    = "",   # Rcpp 20 arg limit
                 embedded        = FALSE,
                 # const_pred    = FALSE,# Rcpp 20 arg limit
                 verbose         = FALSE,
                 validLib        = vector(),
                 ignoreNan       = TRUE,
                 generateSteps   = 0,
                 parameterList   = FALSE,
                 showPlot        = FALSE,
                 noTime          = FALSE ) {

  if ( noTime ) {
    if ( nchar( dataFile ) ) {
      # Read the dataFile into a data.frame, insert index
      # Override dataFile argument to "", assign dataFrame with index column
      dataFrame = read.csv( paste( pathIn, dataFile, sep = '/' ),
                            as.is = TRUE, check.names = FALSE )
      dataFile  = ""
    }
    # Insert an index column as first column for cppEDM
    index     = seq( 1, nrow( dataFrame ) )
    dataFrame = cbind( index, dataFrame )
  }

  if ( ! ValidateDataFrame( pathIn, dataFile, dataFrame,
                            columns, target, noTime ) ) {
    stop( "SMap(): dataFrame validation failed." )
  }
  
  # If lib, pred, columns are vectors/list, convert to string for cppEDM
  # NOTE: columns joined on ',' to enable names with whitespace in cppEDM
  if ( ! is.character( lib ) || length( lib ) > 1 ) {
    lib = FlattenToString( lib )
  }
  if ( ! is.character( pred ) || length( pred ) > 1 ) {
    pred = FlattenToString( pred )
  }
  if ( ! is.character( columns ) || length( columns ) > 1 ) {
    columns = FlattenToString( columns, "," )
  }
  if ( length( strsplit( target, ' ' )[[1]] ) > 1 ) {
    target = paste0( target, ',' ) # space in target: add , for cppEDM
  }

  # NOTE: Rcpp has a 20 argument limit!
  # Mapped to SMap_rcpp() (SMap.cpp) in RcppEDMCommon.cpp
  # smapList has data.frames of: predictions, coefficients, singularValues
  smapList = RtoCpp_SMap( pathIn,
                          dataFile,
                          dataFrame,
                          # pathOut,     # Rcpp 20 arg limit
                          # predictFile, # Rcpp 20 arg limit
                          lib,
                          pred,
                          E,
                          Tp,
                          knn,
                          tau,
                          theta,
                          exclusionRadius,
                          columns,
                          target,
                          # smapCoefFile, # Rcpp 20 arg limit
                          # smapSVFile,   # Rcpp 20 arg limit
                          embedded,
                          # const_pred,   # Rcpp 20 arg limit
                          verbose,
                          validLib,
                          ignoreNan,
                          generateSteps,
                          parameterList )

  if ( showPlot ) {
    PlotSmap( smapList, dataFile, E, Tp )
  }

  return( smapList )
}

#------------------------------------------------------------------------
#
#------------------------------------------------------------------------
Multiview = function( pathIn          = "./",
                      dataFile        = "",
                      dataFrame       = NULL,
                      # pathOut       = "./", # Rcpp 20 arg limit
                      # predictFile   = "",   # Rcpp 20 arg limit
                      lib             = "",
                      pred            = "",
                      D               = 0,
                      E               = 1,
                      Tp              = 1,
                      knn             = 0,
                      tau             = -1,
                      columns         = "",
                      target          = "",
                      multiview       = 0,
                      exclusionRadius = 0,
                      trainLib        = TRUE,
                      excludeTarget   = FALSE,
                      parameterList   = FALSE,
                      verbose         = FALSE,
                      numThreads      = 4,
                      showPlot        = FALSE,
                      noTime          = FALSE ) {

  if ( noTime ) {
    if ( nchar( dataFile ) ) {
      # Read the dataFile into a data.frame, insert index
      # Override dataFile argument to "", assign dataFrame with index column
      dataFrame = read.csv( paste( pathIn, dataFile, sep = '/' ),
                            as.is = TRUE, check.names = FALSE )
      dataFile  = ""
    }
    # Insert an index column as first column for cppEDM
    index     = seq( 1, nrow( dataFrame ) )
    dataFrame = cbind( index, dataFrame )
  }

  if ( ! ValidateDataFrame( pathIn, dataFile, dataFrame,
                            columns, target, noTime ) ) {
    stop( "Multiview(): dataFrame validation failed." )
  }

  # If lib, pred, columns are vectors/list, convert to string for cppEDM
  # NOTE: columns joined on ',' to enable names with whitespace in cppEDM
  if ( ! is.character( lib ) || length( lib ) > 1 ) {
    lib = FlattenToString( lib )
  }
  if ( ! is.character( pred ) || length( pred ) > 1 ) {
    pred = FlattenToString( pred )
  }
  if ( ! is.character( columns ) || length( columns ) > 1 ) {
    columns = FlattenToString( columns, "," )
  }
  if ( length( strsplit( target, ' ' )[[1]] ) > 1 ) {
    target = paste0( target, ',' ) # space in target: add , for cppEDM
  }

  # NOTE: Rcpp has a 20 argument limit!
  # Mapped to Multiview_rcpp() (Multiview.cpp) in RcppEDMCommon.cpp
  # mvList has data.frames "Combo_rho" and  "Predictions" 
  mvList = RtoCpp_Multiview( pathIn,
                             dataFile,
                             dataFrame,
                             # pathOut,     # Rcpp 20 arg limit
                             # predictFile, # Rcpp 20 arg limit
                             lib,
                             pred,
                             D,
                             E, 
                             Tp,
                             knn,
                             tau,
                             columns,
                             target,
                             multiview,
                             exclusionRadius,
                             trainLib,
                             excludeTarget,
                             parameterList,
                             verbose,
                             numThreads )

  if ( showPlot ) {
    PlotObsPred( mvList $ Predictions, dataFile, E, Tp )
  }

  # mvList: [[ "ComboRho" : data.frame, "Predictions" : data.frame,
  #            "ColumnNames" : List, "parameters" : List ]]
  # Append ColumnNames to ComboRho to create View
  colNames = names( mvList $ ColumnNames )  # dict keys
  for ( i in 1 : length( mvList $ ColumnNames ) ) {
    mvList $ ComboRho[[ colNames[i] ]] = mvList $ ColumnNames[[ i ]]
  }

  MV = list( "View"        = mvList $ ComboRho,
             "Predictions" = mvList $ Predictions )

  if ( parameterList ) {
    MV[['parameters']] = mvList $ parameters
  }

  return( MV )
}

#------------------------------------------------------------------------
#
#------------------------------------------------------------------------
CCM = function( pathIn          = "./",
                dataFile        = "",
                dataFrame       = NULL,
                # pathOut         = "./", # Rcpp 20 param limit
                # predictFile     = "",   # Rcpp 20 param limit
                E               = 0, 
                Tp              = 0,
                knn             = 0,
                tau             = -1,
                exclusionRadius = 0,
                columns         = "",
                target          = "",
                libSizes        = "",
                sample          = 0,
                random          = TRUE,
                # replacement   = FALSE,  # Rcpp 20 param limit
                seed            = 0,
                embedded        = FALSE,
                includeData     = FALSE,
                parameterList   = FALSE,
                verbose         = FALSE,
                showPlot        = FALSE,
                noTime          = FALSE ) {

  if ( noTime ) {
    if ( nchar( dataFile ) ) {
      # Read the dataFile into a data.frame, insert index
      # Override dataFile argument to "", assign dataFrame with index column
      dataFrame = read.csv( paste( pathIn, dataFile, sep = '/' ),
                            as.is = TRUE, check.names = FALSE )
      dataFile  = ""
    }
    # Insert an index column as first column for cppEDM
    index     = seq( 1, nrow( dataFrame ) )
    dataFrame = cbind( index, dataFrame )
  }

  if ( ! ValidateDataFrame( pathIn, dataFile, dataFrame,
                            columns, target, noTime ) ) {
    stop( "CCM(): dataFrame validation failed." )
  }

  # If libSizes, columns are vectors/list, convert to string for cppEDM
  # NOTE: columns joined on ',' to enable names with whitespace in cppEDM
  # NOTE: CCM can have multiple target
  if ( ! is.character( libSizes ) || length( libSizes ) > 1 ) {
    libSizes = FlattenToString( libSizes )
  }
  if ( ! is.character( columns ) || length( columns ) > 1 ) {
    columns = FlattenToString( columns, "," )
  }
  if ( ! is.character( target ) || length( target ) > 1 ) {
    columns = FlattenToString( target, "," )
  }
  else {
    if ( length( strsplit( target, ' ' )[[1]] ) > 1 ) {
      target = paste0( target, ',' ) # space in target: add , for cppEDM
    }
  }

  # NOTE: Rcpp has a 20 argument limit!
  # Mapped to CCM_rcpp() (CCM.cpp) in RcppEDMCommon.cpp
  # CCMList has "LibSize" and columns:target target:columns
  CCMList = RtoCpp_CCM( pathIn,
                        dataFile,
                        dataFrame,
                        # pathOut,     # Rcpp 20 param limit
                        # predictFile, # Rcpp 20 param limit
                        E, 
                        Tp,
                        knn,
                        tau,
                        exclusionRadius,
                        columns,
                        target,
                        libSizes,
                        sample,
                        random,
                        # replacement, # Rcpp 20 param limit
                        seed,
                        embedded,
                        includeData,
                        parameterList,
                        verbose )

  if ( showPlot ) {
    ccm.df = CCMList[[ 'LibMeans' ]]
    libSize = ccm.df $ LibSize
    V1      = names( ccm.df )[2]
    V2      = names( ccm.df )[3]
                                    
    title = paste( V1, " : ", V2, "\nE=" , E )
    
    plot( libSize, ccm.df[ , V1 ],
          ylim = range( ccm.df[ , V1 ], ccm.df[ , V2 ] ),
          main = title, col = "blue", type = "l", lwd = 3,
          xlab = 'Library Size',
          ylab = expression( "Prediction Skill (" * rho * ")" ) )
    lines( libSize, ccm.df[ , V2 ], col = "red", lwd = 3 )
    abline( h = 0 )
    legend( 'topright', c( V1, V2 ), 
            fill = c( 'blue', 'red' ), bty = 'n', cex = 1.2 )
  }

  if ( includeData ) {
    output = CCMList  # return list with all data
  }
  else {
    output = CCMList $ LibMeans # return data.frame LibMeans
  }

  return( output )
}

#------------------------------------------------------------------------
#
#------------------------------------------------------------------------
EmbedDimension = function ( pathIn          = "./",
                            dataFile        = "",
                            dataFrame       = NULL,
                            pathOut         = "",
                            predictFile     = "",
                            lib             = "",
                            pred            = "",
                            maxE            = 10,
                            Tp              = 1,
                            tau             = -1,
                            exclusionRadius = 0,
                            columns         = "",
                            target          = "",
                            embedded        = FALSE,
                            verbose         = FALSE,
                            validLib        = vector(),
                            numThreads      = 4,
                            showPlot        = TRUE,
                            noTime          = FALSE ) {

  if ( noTime ) {
    if ( nchar( dataFile ) ) {
      # Read the dataFile into a data.frame, insert index
      # Override dataFile argument to "", assign dataFrame with index column
      dataFrame = read.csv( paste( pathIn, dataFile, sep = '/' ),
                            as.is = TRUE, check.names = FALSE )
      dataFile  = ""
    }
    # Insert an index column as first column for cppEDM
    index     = seq( 1, nrow( dataFrame ) )
    dataFrame = cbind( index, dataFrame )
  }

  if ( ! ValidateDataFrame( pathIn, dataFile, dataFrame,
                            columns, target, noTime ) ) {
    stop( "EmbedDimension(): dataFrame validation failed." )
  }

  # If lib, pred, columns are vectors/list, convert to string for cppEDM
  # NOTE: columns joined on ',' to enable names with whitespace in cppEDM
  if ( ! is.character( lib ) || length( lib ) > 1 ) {
    lib = FlattenToString( lib )
  }
  if ( ! is.character( pred ) || length( pred ) > 1 ) {
    pred = FlattenToString( pred )
  }
  if ( ! is.character( columns ) || length( columns ) > 1 ) {
    columns = FlattenToString( columns, "," )
  }
  if ( length( strsplit( target, ' ' )[[1]] ) > 1 ) {
    target = paste0( target, ',' ) # space in target: add , for cppEDM
  }

  # Mapped to EmbedDimension_rcpp() (EmbedDim.cpp) in RcppEDMCommon.cpp
  df = RtoCpp_EmbedDimension( pathIn,
                              dataFile,
                              dataFrame,
                              pathOut,
                              predictFile,
                              lib,
                              pred, 
                              maxE,
                              Tp,
                              tau,
                              exclusionRadius,
                              columns,
                              target,
                              embedded,
                              verbose,
                              validLib,
                              numThreads )

  if ( showPlot ) {
    title = paste(dataFile , "\nTp=" , Tp )
    plot( df $ E, df $ rho, main = title, xlab = "Embedding Dimension",
          ylab = expression( "Prediction Skill (" * rho * ")" ),
          type = "l", lwd = 3 )
  }

  return ( df )
}

#------------------------------------------------------------------------
#
#------------------------------------------------------------------------
PredictInterval = function( pathIn          = "./",
                            dataFile        = "",
                            dataFrame       = NULL,
                            pathOut         = "./",
                            predictFile     = "",
                            lib             = "",
                            pred            = "",
                            maxTp           = 10,
                            E               = 1,
                            tau             = -1,
                            exclusionRadius = 0,
                            columns         = "",
                            target          = "",
                            embedded        = FALSE,
                            verbose         = FALSE,
                            validLib        = vector(),
                            numThreads      = 4,
                            showPlot        = TRUE,
                            noTime          = FALSE ) {

  if ( noTime ) {
    if ( nchar( dataFile ) ) {
      # Read the dataFile into a data.frame, insert index
      # Override dataFile argument to "", assign dataFrame with index column
      dataFrame = read.csv( paste( pathIn, dataFile, sep = '/' ),
                            as.is = TRUE, check.names = FALSE )
      dataFile  = ""
    }
    # Insert an index column as first column for cppEDM
    index     = seq( 1, nrow( dataFrame ) )
    dataFrame = cbind( index, dataFrame )
  }

  if ( ! ValidateDataFrame( pathIn, dataFile, dataFrame,
                            columns, target, noTime ) ) {
    stop( "PredictInterval(): dataFrame validation failed." )
  }

  # If lib, pred, columns are vectors/list, convert to string for cppEDM
  # NOTE: columns joined on ',' to enable names with whitespace in cppEDM
  if ( ! is.character( lib ) || length( lib ) > 1 ) {
    lib = FlattenToString( lib )
  }
  if ( ! is.character( pred ) || length( pred ) > 1 ) {
    pred = FlattenToString( pred )
  }
  if ( ! is.character( columns ) || length( columns ) > 1 ) {
    columns = FlattenToString( columns, "," )
  }
  if ( length( strsplit( target, ' ' )[[1]] ) > 1 ) {
    target = paste0( target, ',' ) # space in target: add , for cppEDM
  }

  # Mapped to PredictInterval_rcpp() (PredictInterval.cpp) in RcppEDMCommon.cpp
  df = RtoCpp_PredictInterval( pathIn,
                               dataFile,
                               dataFrame,
                               pathOut,
                               predictFile,
                               lib,
                               pred, 
                               maxTp,
                               E,
                               tau,
                               exclusionRadius,
                               columns,
                               target,
                               embedded,
                               verbose,
                               validLib,
                               numThreads )

  if ( showPlot ) {
    title = paste( dataFile , "\nE=" , E )
    plot( df $ Tp, df $ rho, main = title, xlab = "Forecast Interval",
          ylab = expression( "Prediction Skill (" * rho * ")" ),
          type = "l", lwd = 3 )
  }

  return( df )
}

#------------------------------------------------------------------------
#
#------------------------------------------------------------------------
PredictNonlinear = function( pathIn          = "./",
                             dataFile        = "",
                             dataFrame       = NULL,
                             pathOut         = "./",
                             predictFile     = "",
                             lib             = "",
                             pred            = "",
                             theta           = "",
                             E               = 1,
                             Tp              = 1,
                             knn             = 0,
                             tau             = -1,
                             exclusionRadius = 0,
                             columns         = "",
                             target          = "",
                             embedded        = FALSE,
                             verbose         = FALSE,
                             validLib        = vector(),
                             ignoreNan       = TRUE,
                             numThreads      = 4,
                             showPlot        = TRUE,
                             noTime          = FALSE ) {

  if ( noTime ) {
    if ( nchar( dataFile ) ) {
      # Read the dataFile into a data.frame, insert index
      # Override dataFile argument to "", assign dataFrame with index column
      dataFrame = read.csv( paste( pathIn, dataFile, sep = '/' ),
                            as.is = TRUE, check.names = FALSE )
      dataFile  = ""
    }
    # Insert an index column as first column for cppEDM
    index     = seq( 1, nrow( dataFrame ) )
    dataFrame = cbind( index, dataFrame )
  }

  if ( ! ValidateDataFrame( pathIn, dataFile, dataFrame,
                            columns, target, noTime ) ) {
    stop( "PredictNonlinear(): dataFrame validation failed." )
  }

  # If lib, pred, theta, columns are vectors/list, convert to string for cppEDM
  # NOTE: columns joined on ',' to enable names with whitespace in cppEDM
  if ( ! is.character( lib ) || length( lib ) > 1 ) {
    lib = FlattenToString( lib )
  }
  if ( ! is.character( pred ) || length( pred ) > 1 ) {
    pred = FlattenToString( pred )
  }
  if ( ! is.character( theta ) || length( theta ) > 1 ) {
    theta = FlattenToString( theta )
  }
  if ( ! is.character( columns ) || length( columns ) > 1 ) {
    columns = FlattenToString( columns, "," )
  }
  if ( length( strsplit( target, ' ' )[[1]] ) > 1 ) {
    target = paste0( target, ',' ) # space in target: add , for cppEDM
  }

  # Mapped to PredictNonlinear_rcpp() (PredictNL.cpp) in RcppEDMCommon.cpp
  df = RtoCpp_PredictNonlinear( pathIn,
                                dataFile,
                                dataFrame,
                                pathOut,
                                predictFile,
                                lib,
                                pred, 
                                theta,
                                E,
                                Tp,
                                knn,
                                tau,
                                exclusionRadius,
                                columns,
                                target,
                                embedded,
                                verbose,
                                validLib,
                                ignoreNan,
                                numThreads )

  if ( showPlot ) {
    title = paste(dataFile , "\nE=", E )
    plot( df $ Theta, df $ rho, main=title, 
          xlab = "S-map Localisation",
          ylab = expression( "Prediction Skill (" * rho * ")" ),
          type = "l", lwd = 3 )
  }

  return( df )
}
