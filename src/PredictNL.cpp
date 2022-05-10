
#include "RcppEDMCommon.h"

//---------------------------------------------------------------
// Input data path and file
//---------------------------------------------------------------
r::DataFrame PredictNonlinear_rcpp( std::string  pathIn,
                                    std::string  dataFile,
                                    r::DataFrame dataFrame,
                                    std::string  pathOut,
                                    std::string  predictFile,
                                    std::string  lib,
                                    std::string  pred,
                                    std::string  theta,
                                    int          E,
                                    int          Tp,
                                    int          knn,
                                    int          tau,
                                    int          exclusionRadius,
                                    std::string  columns,
                                    std::string  target,
                                    bool         embedded,
                                    bool         verbose,
                                    std::vector<bool> validLib,
                                    unsigned     numThreads ) {

    DataFrame< double > PredictDF;

    if ( dataFile.size() ) {
        // dataFile specified, dispatch overloaded PredictNonlinear,
        // ignore dataFrame
        PredictDF = PredictNonlinear( pathIn,
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

        PredictDF = PredictNonlinear( dataFrame_,
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
                                      numThreads );
    }
    else {
        Rcpp::warning("PredictNonlinear_rcpp(): Invalid input.\n");
    }

    return DataFrameToDF( PredictDF );
}
