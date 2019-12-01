
#include "RcppEDMCommon.h"

//----------------------------------------------------------------
// Compute Error Wrapper method
// @param vec1      : the first vector to compare
// @param vec2      : the second vector to compare
// @return          : map/dictionary with the rho, mae, rmse
//----------------------------------------------------------------
r::List ComputeError_rcpp ( std::vector<double> vec1, 
                            std::vector<double> vec2 ) {
    
    std::valarray<double> val1 ( vec1.data(), vec1.size() );  
    std::valarray<double> val2 ( vec2.data(), vec2.size() );  
    
    VectorError vecErr = ComputeError( val1, val2 );
    
    // Setup as map instead of vecErr struct
    return r::List::create( r::Named( "MAE"  ) = vecErr.MAE,
                            r::Named( "rho"  ) = vecErr.rho,
                            r::Named( "RMSE" ) = vecErr.RMSE );
}
