
#include "RcppEDMCommon.h"

//---------------------------------------------------------------
// Input data path and file
//---------------------------------------------------------------
r::DataFrame PredictInterval_rcpp( std::string  pathIn,
                                   std::string  dataFile,
                                   r::DataFrame dataList,
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
        // ignore dataList
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
    else if ( dataList.size() ) {
        DataFrame< double > dataFrame = DFToDataFrame( dataList );

        PredictDF = PredictInterval( dataFrame,
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
