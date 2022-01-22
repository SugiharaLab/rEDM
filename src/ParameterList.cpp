
#include "RcppEDMCommon.h"

//----------------------------------------------------------------
//
//----------------------------------------------------------------
r::List ParamMaptoList( std::map< std::string, std::string > m ) {

    r::List L;

    for ( auto pi = m.begin(); pi != m.end(); pi++ ) {
        // string types
        if ( pi->first == "version"  or
             pi->first == "method"   or pi->first == "columns" or
             pi->first == "target"   or pi->first == "pathIn"  or
             pi->first == "dataFile" or pi->first == "pathOut" or
             pi->first == "predictOutputFile" or
             pi->first == "SmapOutputFile"    or
             pi->first == "blockOutputFile" ) {

            L[ pi->first ] = pi->second;
        }
        // int types
        else if ( pi->first == "E"   or pi->first == "Tp"  or
                  pi->first == "knn" or pi->first == "tau" or
                  pi->first == "exclusionRadius"   or
                  pi->first == "seed"              or
                  pi->first == "subSamples"        or
                  pi->first == "multiviewEnsemble" or
                  pi->first == "multiviewD"        or
                  pi->first == "generateSteps" ) {

            L[ pi->first ] = std::stoi( pi->second );
        }
        // boolean types
        else if ( pi->first == "randomLib"   or
                  pi->first == "replacement" or
                  pi->first == "includeData" or
                  pi->first == "multiviewTrainLib"      or
                  pi->first == "multiviewExcludeTarget" or
                  pi->first == "embedded"      or
                  pi->first == "const_predict" or
                  pi->first == "parameterList" or
                  pi->first == "verbose" ) {

            if ( pi->second == "0" ) {
                L[ pi->first ] = false;
            }
            if ( pi->second == "1" ) {
                L[ pi->first ] = true;
            }
        }
        // vector of int
        else if ( pi->first == "lib"      or pi->first == "pred" or
                  pi->first == "libSizes" or pi->first == "validLib" ) {
            std::stringstream iss( pi->second );
            std::vector< int > intVector;
            int                value;

            while ( iss >> value ) {
                intVector.push_back( value );
            }

            L[ pi->first ] = intVector;
        }
        // float type
        else if ( pi->first == "theta" ) {
            L[ pi->first ] = std::stof( pi->second );
        }
    }

    return L;
}
