
#include "RcppEDMCommon.h"

//-------------------------------------------------------------
// 
//-------------------------------------------------------------
r::DataFrame Simplex_rcpp( std::string  pathIn,
                           std::string  dataFile,
                           r::DataFrame dataList,
                           std::string  pathOut,
                           std::string  predictFile,
                           std::string  lib,
                           std::string  pred, 
                           int          E,
                           int          Tp,
                           int          knn,
                           int          tau, 
                           int          exclusionRadius, 
                           std::string  columns,
                           std::string  target,
                           bool         embedded,
                           bool         const_predict,
                           bool         verbose ) {

    DataFrame< double > S;
    
    if ( dataFile.size() ) {
        // dataFile specified, dispatch overloaded Simplex, ignore dataList
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
                     verbose );
    }
    else if ( dataList.size() ) {
        DataFrame< double > dataFrame = DFToDataFrame( dataList );
        
        S = Simplex( dataFrame,
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
                     verbose );
    }
    else {
        Rcpp::warning( "Simplex_rcpp(): Invalid input.\n" );
    }
    
    return DataFrameToDF( S );
}
