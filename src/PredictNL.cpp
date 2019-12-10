
#include "RcppEDMCommon.h"

//---------------------------------------------------------------
// Input data path and file
//---------------------------------------------------------------
r::DataFrame PredictNonlinear_rcpp( std::string  pathIn,
                                    std::string  dataFile,
                                    r::DataFrame dataList,
                                    std::string  pathOut,
                                    std::string  predictFile,
                                    std::string  lib,
                                    std::string  pred,
                                    std::string  theta,
                                    int          E,
                                    int          Tp,
                                    int          knn,
                                    int          tau,
                                    std::string  columns,
                                    std::string  target,
                                    bool         embedded,
                                    bool         verbose,
                                    unsigned     numThreads ) {

    DataFrame< double > PredictDF;

    if ( dataFile.size() ) {
        // dataFile specified, dispatch overloaded PredictNonlinear,
        // ignore dataList
        PredictDF  = PredictNonlinear( pathIn,
                                       dataFile,
                                       pathOut,
                                       predictFile,
                                       lib,
                                       pred,
                                       theta,
                                       E,
                                       Tp,
                                       knn,
                                       tau,
                                       columns,
                                       target,
                                       embedded,
                                       verbose,
                                       numThreads );
    }
    else if ( dataList.size() ) {
        DataFrame< double > dataFrame = DFToDataFrame( dataList );
        
        PredictDF  = PredictNonlinear( dataFrame,
                                       pathOut,
                                       predictFile,
                                       lib,
                                       pred,
                                       theta,
                                       E,
                                       Tp,
                                       knn,
                                       tau,
                                       columns,
                                       target,
                                       embedded,
                                       verbose,
                                       numThreads );
    }
    else {
        Rcpp::warning("PredictNonlinear_rcpp(): Invalid input.\n");
    }

    return DataFrameToDF( PredictDF );
}
