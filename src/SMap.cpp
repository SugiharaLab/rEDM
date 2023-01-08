#include "RcppEDMCommon.h"

//----------------------------------------------------------
// 
//----------------------------------------------------------
r::List SMap_rcpp( std::string       pathIn,
                   std::string       dataFile,
                   r::DataFrame      dataFrame,
                   //std::string       pathOut,     // Rcpp has 20 arg limit
                   //std::string       predictFile, // Rcpp has 20 arg limit
                   std::string       lib,
                   std::string       pred,
                   int               E,
                   int               Tp,
                   int               knn,
                   int               tau,
                   double            theta,
                   int               exlusionRadius,
                   std::string       columns,
                   std::string       target,
                   std::string       smapFile,
                   // std::string    jacobians, // Rcpp has 20 arg limit
                   bool              embedded,
                   bool              const_predict,
                   bool              verbose,
                   std::vector<bool> validLib,
                   int               generateSteps,
                   //bool            generateLibrary, // Rcpp has 20 arg limit
                   bool              parameterList ) {

    SMapValues SM;
    
    std::string pathOut("./");    // Rcpp has 20 arg limit
    std::string predictFile("");  // Rcpp has 20 arg limit
    std::string jacobians("");    // Rcpp has 20 arg limit
    bool generateLibrary = false; // Rcpp has 20 arg limit

    if ( dataFile.size() ) {
        // dataFile specified, dispatch overloaded SMap, ignore dataFrame
        
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
                   verbose,
                   validLib,
                   generateSteps,
                   generateLibrary,
                   parameterList );
    }
    else if ( dataFrame.size() ) {
        DataFrame< double > dataFrame_ = DFToDataFrame( dataFrame );

        SM = SMap( dataFrame_,
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
                   verbose,
                   validLib,
                   generateSteps,
                   generateLibrary,
                   parameterList );
    }
    else {
        Rcpp::warning( "SMap_rcpp(): Invalid input.\n" );
    }

    r::DataFrame df_pred = DataFrameToDF( SM.predictions  );
    r::DataFrame df_coef = DataFrameToDF( SM.coefficients );
    r::List output = r::List::create( r::Named("predictions")  = df_pred,
                                      r::Named("coefficients") = df_coef );

    if ( parameterList ) {
        r::List paramList = ParamMaptoList( SM.parameterMap );
        output["parameters"] = paramList;
    }

    return output;
}
