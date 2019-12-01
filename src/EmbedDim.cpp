
#include "RcppEDMCommon.h"

//---------------------------------------------------------------
// 
//---------------------------------------------------------------
r::DataFrame EmbedDimension_rcpp(   std::string  pathIn,
                                    std::string  dataFile,
                                    r::DataFrame dataList,
                                    std::string  pathOut,
                                    std::string  predictFile,
                                    std::string  lib,
                                    std::string  pred,
                                    int          maxE,
                                    int          Tp,
                                    int          tau,
                                    std::string  columns,
                                    std::string  target,
                                    bool         embedded,
                                    bool         verbose,
                                    unsigned     numThreads ) {
    
    DataFrame< double > EmbedDimDF;

    if ( dataFile.size() ) {
        // dataFile specified, dispatch overloaded EmbedDimension,
        // ignore dataList
        EmbedDimDF = EmbedDimension( pathIn,
                                     dataFile,
                                     pathOut,
                                     predictFile,
                                     lib,
                                     pred,
                                     maxE,
                                     Tp,
                                     tau,
                                     columns,
                                     target,
                                     embedded,
                                     verbose,
                                     numThreads );
    }
    else if ( dataList.size() ) {
        DataFrame< double > dataFrame = DFToDataFrame( dataList );
        
        EmbedDimDF = EmbedDimension( dataFrame,
                                     pathOut,
                                     predictFile,
                                     lib,
                                     pred,
                                     maxE,
                                     Tp,
                                     tau,
                                     columns,
                                     target,
                                     embedded,
                                     verbose,
                                     numThreads );
    }
    else {
        Rcpp::warning( "EmbedDimension_rcpp(): Invalid input.\n" );
    }
    
    return DataFrameToDF( EmbedDimDF );
}
