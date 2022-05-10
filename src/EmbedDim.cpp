
#include "RcppEDMCommon.h"

//---------------------------------------------------------------
// 
//---------------------------------------------------------------
r::DataFrame EmbedDimension_rcpp( std::string  pathIn,
                                  std::string  dataFile,
                                  r::DataFrame dataFrame,
                                  std::string  pathOut,
                                  std::string  predictFile,
                                  std::string  lib,
                                  std::string  pred,
                                  int          maxE,
                                  int          Tp,
                                  int          tau,
                                  int          exclusionRadius,
                                  std::string  columns,
                                  std::string  target,
                                  bool         embedded,
                                  bool         verbose,
                                  std::vector<bool> validLib,
                                  unsigned     numThreads ) {

    DataFrame< double > EmbedDimDF;

    if ( dataFile.size() ) {
        // dataFile specified, dispatch overloaded EmbedDimension,
        // ignore dataFrame
        EmbedDimDF = EmbedDimension( pathIn,
                                     dataFile,
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
                                     numThreads );
    }
    else if ( dataFrame.size() ) {
        DataFrame< double > dataFrame_ = DFToDataFrame( dataFrame );

        EmbedDimDF = EmbedDimension( dataFrame_,
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
                                     numThreads );
    }
    else {
        Rcpp::warning( "EmbedDimension_rcpp(): Invalid input.\n" );
    }

    return DataFrameToDF( EmbedDimDF );
}
