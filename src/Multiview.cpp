
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

    r::DataFrame comboRho    = DataFrameToDF( MV.ComboRho    );
    r::DataFrame predictions = DataFrameToDF( MV.Predictions );

    // ColumnNames are: map< string, vector<string> >, convert to List
    r::List columnNames;
    for ( auto cni  = MV.ColumnNames.begin();
               cni != MV.ColumnNames.end(); cni++ ) {
        r::StringVector strVec;
        std::vector< std::string > names = cni->second;
        for ( auto ni = names.begin(); ni != names.end(); ni++ ) {
            strVec.push_back( *ni );
        }
        columnNames[ cni->first ] = strVec;
    }

    r::List output = r::List::create(
        r::Named("ComboRho")    = comboRho,
        r::Named("ColumnNames") = columnNames,
        r::Named("Predictions") = predictions );

    if ( parameterList ) {
        r::List paramList = ParamMaptoList( MV.parameterMap );
        output["parameters"] = paramList;
    }

    // Multiview.R in EDM.R will convert comboLines into an R data.frame
    return output;
}
