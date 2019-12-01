#include "RcppEDMCommon.h"

//----------------------------------------------------------
// 
//----------------------------------------------------------
r::List SMap_rcpp( std::string  pathIn, 
                   std::string  dataFile,
                   r::DataFrame dataFrameIn,
                   std::string  pathOut,
                   std::string  predictFile,
                   std::string  lib,
                   std::string  pred, 
                   int          E,
                   int          Tp,
                   int          knn,
                   int          tau,
                   double       theta,
                   int          exlusionRadius,
                   std::string  columns,
                   std::string  target,
                   std::string  smapFile,
                   std::string  jacobians,
                   bool         embedded,
                   bool         const_predict,
                   bool         verbose ) {
    
    SMapValues SM;
    
    if ( dataFile.size() ) {
        // dataFile specified, dispatch overloaded SMap, ignore dataFrameIn
        
        SM = SMap( pathIn,
                   dataFile,
                   pathOut,
                   predictFile,
                   lib,
                   pred,
                   E, 
                   Tp,
                   knn,
                   tau,
                   theta,
                   exlusionRadius,
                   columns, 
                   target,
                   smapFile,
                   jacobians,
                   embedded,
                   const_predict,
                   verbose );
    }
    else if ( dataFrameIn.size() ) {
        DataFrame< double > dataFrame = DFToDataFrame( dataFrameIn );
        
        SM = SMap( dataFrame,
                   pathOut,
                   predictFile,
                   lib,
                   pred,
                   E, 
                   Tp,
                   knn,
                   tau,
                   theta,
                   exlusionRadius,
                   columns, 
                   target,
                   smapFile,
                   jacobians,
                   embedded,
                   const_predict,
                   verbose );
    }
    else {
        Rcpp::warning( "SMap_rcpp(): Invalid input.\n" );
    }
    
    r::DataFrame df_pred = DataFrameToDF( SM.predictions  );
    r::DataFrame df_coef = DataFrameToDF( SM.coefficients );
    r::List output = r::List::create( r::Named("predictions")  = df_pred,
                                      r::Named("coefficients") = df_coef );

    return output;
}
