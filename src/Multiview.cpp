
#include "RcppEDMCommon.h"

//--------------------------------------------------------------
// 
//--------------------------------------------------------------
r::List Multiview_rcpp ( std::string  pathIn,
                         std::string  dataFile,
                         r::DataFrame dataFrame,
                         std::string  pathOut,
                         std::string  predictFile,
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
                         bool         verbose,
                         unsigned int numThreads ) {

    MultiviewValues MV;

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

    // Multiview.R in EDM.R will convert comboLines into an R data.frame
    return output;
}
