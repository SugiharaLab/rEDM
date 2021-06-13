
#include "RcppEDMCommon.h"

//-------------------------------------------------------------
// 
//-------------------------------------------------------------
r::DataFrame Simplex_rcpp( std::string       pathIn,
                           std::string       dataFile,
                           r::DataFrame      dataFrame,
                           std::string       pathOut,
                           std::string       predictFile,
                           std::string       lib,
                           std::string       pred,
                           int               E,
                           int               Tp,
                           int               knn,
                           int               tau,
                           int               exclusionRadius,
                           std::string       columns,
                           std::string       target,
                           bool              embedded,
                           bool              const_predict,
                           bool              verbose,
                           std::vector<bool> validLib ) {

    DataFrame< double > S;

    if ( dataFile.size() ) {
        // dataFile specified, dispatch overloaded Simplex, ignore dataFrame
        S = Simplex( pathIn,
                     dataFile,
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
                     const_predict,
                     verbose,
                     validLib );
    }
    else if ( dataFrame.size() ) {
        DataFrame< double > dataFrame_ = DFToDataFrame( dataFrame );

        S = Simplex( dataFrame_,
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
                     const_predict,
                     verbose,
                     validLib );
    }
    else {
        Rcpp::warning( "Simplex_rcpp(): Invalid input.\n" );
    }

    return DataFrameToDF( S );
}
