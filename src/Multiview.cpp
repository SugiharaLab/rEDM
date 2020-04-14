
#include "RcppEDMCommon.h"

//--------------------------------------------------------------
// 
//--------------------------------------------------------------
r::List Multiview_rcpp ( std::string  pathIn,
                         std::string  dataFile,
                         r::DataFrame dataList,
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
                         bool         verbose,
                         unsigned int numThreads ) {

    MultiviewValues MV;

    if ( dataFile.size() ) {
        // dataFile specified, dispatch overloaded Multiview, ignore dataList
        
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
                        verbose,
                        numThreads );
    }
    else if ( dataList.size() ) {
        DataFrame< double > dataFrame = DFToDataFrame( dataList );
        
        MV = Multiview( dataFrame,
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
                        verbose,
                        numThreads );
    }
    else {
        Rcpp::warning( "Multiview_rcpp(): Invalid input.\n" );
    }
    
    // Copy Combo_rho_table into a Rcpp::StringVector
    r::StringVector comboLines( MV.Combo_rho_table.size() );
    for ( size_t row = 0; row < MV.Combo_rho_table.size(); row++ ) {
        comboLines[ row ] = MV.Combo_rho_table[ row ];
    }

    r::DataFrame predictions = DataFrameToDF( MV.Predictions );
    
    r::List output = r::List::create(
        r::Named("Views")       = comboLines,
        r::Named("Predictions") = predictions );

    // Multiview.R in EDM.R will convert comboLines into an R data.frame
    return output;
}
