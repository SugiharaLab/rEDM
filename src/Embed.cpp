
#include "RcppEDMCommon.h"
#include "API.h"

//---------------------------------------------------------------
// 
//---------------------------------------------------------------
r::DataFrame Embed_rcpp( std::string  path,
                         std::string  dataFile,
                         r::DataFrame dataFrame,
                         int          E,
                         int          tau,
                         std::string  columns,
                         bool         verbose ) {

    DataFrame< double > embedded;

    if ( dataFile.size() ) {
        // dataFile specified, dispatch overloaded Embed, ignore dataFrame
        embedded = Embed( path,
                          dataFile,
                          E,
                          tau,
                          columns,
                          verbose );
    }
    else if ( dataFrame.ncol() ) {
        DataFrame< double > dataFrame_ = DFToDataFrame( dataFrame );

        embedded = Embed( dataFrame_,
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
r::DataFrame MakeBlock_rcpp( r::DataFrame             dataFrame,
                             int                      E,
                             int                      tau,
                             std::vector<std::string> columnNames,
                             bool                     deletePartial ) {

    DataFrame< double > dataFrame_ = DFToDataFrame( dataFrame );

    DataFrame< double > block = MakeBlock( dataFrame_,
                                           E,
                                           tau,
                                           columnNames,
                                           deletePartial );

    return DataFrameToDF( block );
}
