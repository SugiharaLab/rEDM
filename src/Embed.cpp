
#include "RcppEDMCommon.h"
#include "Embed.h"

//---------------------------------------------------------------
// 
//---------------------------------------------------------------
r::DataFrame Embed_rcpp( std::string  path,
                         std::string  dataFile,
                         r::DataFrame df,
                         int          E,
                         int          tau,
                         std::string  columns,
                         bool         verbose ) {

    DataFrame< double > embedded;
    
    if ( dataFile.size() ) {
        // dataFile specified, dispatch overloaded Embed, ignore df
        embedded = Embed( path,
                          dataFile,
                          E,
                          tau,
                          columns,
                          verbose );
    }
    else if ( df.ncol() ) {
        DataFrame< double > dataFrame = DFToDataFrame( df );
        
        embedded = Embed( dataFrame,
                          E,
                          tau,
                          columns,
                          verbose );
    }
    else {
        Rcpp::warning( "Embed_rcpp(): Invalid input.\n" );
    }

    return DataFrameToDF( embedded );
}

//---------------------------------------------------------------
// 
//---------------------------------------------------------------
r::DataFrame MakeBlock_rcpp( r::DataFrame             dataList,
                             int                      E,
                             int                      tau,
                             std::vector<std::string> columnNames,
                             bool                     verbose ) {
    
    DataFrame< double > dataFrame = DFToDataFrame( dataList );
    
    DataFrame< double > block = MakeBlock( dataFrame,
                                           E,
                                           tau,
                                           columnNames,
                                           verbose );
    
    return DataFrameToDF( block );
}
