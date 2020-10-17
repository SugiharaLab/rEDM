
#include "RcppEDMCommon.h"

//---------------------------------------------------------------
// Input data path and file
//---------------------------------------------------------------
r::DataFrame PredictInterval_rcpp( std::string  pathIn,
                                   std::string  dataFile,
                                   r::DataFrame dataFrame,
                                   std::string  pathOut,
                                   std::string  predictFile,
                                   std::string  lib,
                                   std::string  pred,
                                   int          maxTp,
                                   int          E,
                                   int          tau,
                                   std::string  columns,
                                   std::string  target,
                                   bool         embedded,
                                   bool         verbose,
                                   unsigned     numThreads ) {

    DataFrame< double > PredictDF;

    if ( dataFile.size() ) {
        // dataFile specified, dispatch overloaded PredictInterval,
        // ignore dataFrame
        PredictDF = PredictInterval( pathIn,
                                     dataFile,
                                     pathOut,
                                     predictFile,
                                     lib,
                                     pred,
                                     maxTp,
                                     E,
                                     tau,
                                     columns,
                                     target,
                                     embedded,
                                     verbose,
                                     numThreads );
    }
    else if ( dataFrame.size() ) {
        DataFrame< double > dataFrame_ = DFToDataFrame( dataFrame );

        PredictDF = PredictInterval( dataFrame_,
                                     pathOut,
                                     predictFile,
                                     lib,
                                     pred,
                                     maxTp,
                                     E,
                                     tau,
                                     columns,
                                     target,
                                     embedded,
                                     verbose,
                                     numThreads );
    }
    else {
        Rcpp::warning("PredictInterval_rcpp(): Invalid input.\n");
    }
    
    return DataFrameToDF( PredictDF );
}
