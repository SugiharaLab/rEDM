
#include "RcppEDMCommon.h"

//-----------------------------------------------------------
// 
//-----------------------------------------------------------
r::DataFrame CCM_rcpp( std::string  pathIn, 
                       std::string  dataFile,
                       r::DataFrame dataList,
                       std::string  pathOut,
                       std::string  predictFile,
                       int          E,
                       int          Tp,
                       int          knn,
                       int          tau, 
                       std::string  columns,
                       std::string  target,
                       std::string  libSizes,
                       int          sample,
                       bool         random,
                       bool         replacement,
                       unsigned     seed, 
                       bool         verbose ) {
    
    DataFrame< double > ccmOutput;

    if ( dataFile.size() ) {
        // dataFile specified, dispatch overloaded CCM, ignore dataList
        ccmOutput = CCM( pathIn,
                         dataFile,
                         pathOut,
                         predictFile,
                         E, 
                         Tp,
                         knn,
                         tau,
                         columns,
                         target, 
                         libSizes,
                         sample,
                         random,
                         replacement,
                         seed,
                         verbose );
    }
    else if ( dataList.size() ) {
        DataFrame< double > dataFrame = DFToDataFrame( dataList );

        ccmOutput = CCM( dataFrame,
                         pathOut,
                         predictFile,
                         E, 
                         Tp,
                         knn,
                         tau,
                         columns,
                         target, 
                         libSizes,
                         sample,
                         random,
                         replacement,
                         seed,
                         verbose );
    }
    else {
        Rcpp::warning( "CCM_rcpp(): No dataFile or dataFrame.\n" );
    }
    
    return DataFrameToDF( ccmOutput );
}
