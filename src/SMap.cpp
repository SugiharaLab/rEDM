#include "RcppEDMCommon.h"

//----------------------------------------------------------
// 
//----------------------------------------------------------
r::List SMap_rcpp( std::string       pathIn,
                   std::string       dataFile,
                   r::DataFrame      dataFrame,
                   //std::string     pathOut,       // Rcpp 20 arg limit
                   //std::string     predictFile,   // Rcpp 20 arg limit
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
                   //std::string     smapCoefFile,  // Rcpp 20 arg limit
                   //std::string     smapSVFile,    // Rcpp 20 arg limit
                   bool              embedded,
                   //bool            const_predict, // Rcpp 20 arg limit
                   bool              verbose,
                   std::vector<bool> validLib,
                   bool              ignoreNan,
                   int               generateSteps,
                   //bool            generateLibrary, // Rcpp 20 arg limit
                   bool              parameterList ) {

    SMapValues SM;

    std::string pathOut("./");    // Rcpp 20 arg limit
    std::string predictFile("");  // Rcpp 20 arg limit
    std::string smapCoefFile(""); // Rcpp 20 arg limit
    std::string smapSVFile("");   // Rcpp 20 arg limit
    bool generateLibrary = false; // Rcpp 20 arg limit
    bool const_predict   = false; // Rcpp 20 arg limit

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
                   smapCoefFile,
                   smapSVFile,
                   embedded,
                   const_predict,
                   verbose,
                   validLib,
                   ignoreNan,
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
                   smapCoefFile,
                   smapSVFile,
                   embedded,
                   const_predict,
                   verbose,
                   validLib,
                   ignoreNan,
                   generateSteps,
                   generateLibrary,
                   parameterList );
    }
    else {
        Rcpp::warning( "SMap_rcpp(): Invalid input.\n" );
    }

    r::DataFrame df_pred = DataFrameToDF( SM.predictions    );
    r::DataFrame df_coef = DataFrameToDF( SM.coefficients   );
    r::DataFrame df_SV   = DataFrameToDF( SM.singularValues );
    r::List output = r::List::create( r::Named("predictions")  = df_pred,
                                      r::Named("coefficients") = df_coef,
                                      r::Named("singularValues") = df_SV );

    if ( parameterList ) {
        r::List paramList = ParamMaptoList( SM.parameterMap );
        output["parameters"] = paramList;
    }

    return output;
}
