
#include "RcppEDMCommon.h"

//--------------------------------------------------------------
// 
//--------------------------------------------------------------
r::List Multiview_rcpp ( std::string  pathIn,
                         std::string  dataFile,
                         r::DataFrame dataFrame,
                         //std::string  pathOut,     // Rcpp 20 arg limit
                         //std::string  predictFile, // Rcpp 20 arg limit
                         std::string  lib,
                         std::string  pred,
                         int          D,
                         int          E,
                         int          Tp,
                         int          knn,
                         int          tau, 
                         std::string  columns,
                         std::string  target,
                         int          multiview,
                         int          exclusionRadius,
                         bool         trainLib,
                         bool         excludeTarget,
                         bool         parameterList,
                         bool         verbose,
                         unsigned int numThreads ) {

    MultiviewValues MV;

    std::string pathOut("./");   // Rcpp has 20 arg limit
    std::string predictFile(""); // Rcpp has 20 arg limit

    if ( dataFile.size() ) {
        // dataFile specified, dispatch overloaded Multiview, ignore dataFrame

        MV = Multiview( pathIn,
                        dataFile,
                        pathOut,
                        predictFile,
                        lib,
                        pred,
                        D,
                        E,
                        Tp,
                        knn,
                        tau,
                        columns,
                        target,
                        multiview,
                        exclusionRadius,
                        trainLib,
                        excludeTarget,
                        parameterList,
                        verbose,
                        numThreads );
    }
    else if ( dataFrame.size() ) {
        DataFrame< double > dataFrame_ = DFToDataFrame( dataFrame );

        MV = Multiview( dataFrame_,
                        pathOut,
                        predictFile,
                        lib,
                        pred,
                        D,
                        E,
                        Tp,
                        knn,
                        tau,
                        columns,
                        target,
                        multiview,
                        exclusionRadius,
                        trainLib,
                        excludeTarget,
                        parameterList,
                        verbose,
                        numThreads );
    }
    else {
        Rcpp::warning( "Multiview_rcpp(): Invalid input.\n" );
    }

    // Copy ComboRhoTable into a Rcpp::StringVector
    r::StringVector comboLines( MV.ComboRhoTable.size() );
    for ( size_t row = 0; row < MV.ComboRhoTable.size(); row++ ) {
        comboLines[ row ] = MV.ComboRhoTable[ row ];
    }

    r::DataFrame predictions = DataFrameToDF( MV.Predictions );

    r::List output = r::List::create(
        r::Named("Views")       = comboLines,
        r::Named("Predictions") = predictions );

    if ( parameterList ) {
        r::List paramList = ParamMaptoList( MV.parameterMap );
        output["parameters"] = paramList;
    }

    // Multiview.R in EDM.R will convert comboLines into an R data.frame
    return output;
}
