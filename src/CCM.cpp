
#include "RcppEDMCommon.h"

//-----------------------------------------------------------
// 
//-----------------------------------------------------------
Rcpp::List CCM_rcpp( std::string  pathIn, 
                     std::string  dataFile,
                     r::DataFrame dataFrame,
                     std::string  pathOut,
                     std::string  predictFile,
                     int          E,
                     int          Tp,
                     int          knn,
                     int          tau,
                     int          exclusionRadius,
                     std::string  columns,
                     std::string  target,
                     std::string  libSizes,
                     int          sample,
                     bool         random,
                     bool         replacement,
                     unsigned     seed,
                     bool         includeData,
                     bool         verbose ) {
    
    CCMValues ccmValues;

    if ( dataFile.size() ) {
        // dataFile specified, dispatch overloaded CCM, ignore dataFrame
        ccmValues = CCM( pathIn,
                         dataFile,
                         pathOut,
                         predictFile,
                         E, 
                         Tp,
                         knn,
                         tau,
                         exclusionRadius,
                         columns,
                         target, 
                         libSizes,
                         sample,
                         random,
                         replacement,
                         seed,
                         includeData,
                         verbose );
    }
    else if ( dataFrame.size() ) {
        DataFrame< double > dataFrame_ = DFToDataFrame( dataFrame );

        ccmValues = CCM( dataFrame_,
                         pathOut,
                         predictFile,
                         E, 
                         Tp,
                         knn,
                         tau,
                         exclusionRadius,
                         columns,
                         target, 
                         libSizes,
                         sample,
                         random,
                         replacement,
                         seed,
                         includeData,
                         verbose );
    }
    else {
        Rcpp::warning( "CCM_rcpp(): No dataFile or dataFrame.\n" );
    }

    // Ouput Rcpp DataFrames
    r::DataFrame allLibStat = DataFrameToDF( ccmValues.AllLibStats );

    r::List output;
    if ( includeData ) {
        // Have to unroll and convert CCMValues.Predictions forward_list
        // to Rcpp::DataFrame for output.
        r::List PredictionsList1;
        for ( auto pi =  ccmValues.CrossMap1.Predictions.begin();
              pi != ccmValues.CrossMap1.Predictions.end(); ++pi ) {
            PredictionsList1.push_back( DataFrameToDF( *pi ) );
        }
        r::List PredictionsList2;
        for ( auto pi =  ccmValues.CrossMap2.Predictions.begin();
              pi != ccmValues.CrossMap2.Predictions.end(); ++pi ) {
            PredictionsList2.push_back( DataFrameToDF( *pi ) );
        }
        
        r::DataFrame cm1_PredStat =
            DataFrameToDF( ccmValues.CrossMap1.PredictStats );
        r::DataFrame cm2_PredStat =
            DataFrameToDF( ccmValues.CrossMap2.PredictStats );
        
        output =
            r::List::create(r::Named( "LibMeans"         ) = allLibStat,
                            r::Named( "CCM1_PredictStat" ) = cm1_PredStat,
                            r::Named( "CCM1_Predictions" ) = PredictionsList1,
                            r::Named( "CCM2_PredictStat" ) = cm2_PredStat,
                            r::Named( "CCM2_Predictions" ) = PredictionsList2);
    }
    else {
        output = r::List::create( r::Named( "LibMeans" ) = allLibStat);
    }
    return output;
}
