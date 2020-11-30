
#include "Parameter.h"

//----------------------------------------------------------------
// Constructor
// Default values set in Parameter.h declaration
//----------------------------------------------------------------
Parameters::Parameters(
    Method      method,
    std::string pathIn,
    std::string dataFile,
    std::string pathOut,
    std::string predictOutputFile,

    std::string lib_str,
    std::string pred_str,

    int         E,
    int         Tp,
    int         knn,
    int         tau,
    double      theta,
    int         exclusionRadius,

    std::string columns_str,
    std::string target_str,

    bool        embedded,
    bool        const_predict,
    bool        verbose,

    std::string SmapOutputFile,
    std::string blockOutputFile,

    int         multiviewEnsemble,
    int         multiviewD,
    bool        multiviewTrainLib,
    bool        multiviewExcludeTarget,

    std::string libSizes_str,
    int         subSamples,
    bool        randomLib,
    bool        replacement,
    unsigned    seed,
    bool        includeData
    ) :
    // Variable initialization from Parameters arguments
    method           ( method ),
    pathIn           ( pathIn ),
    dataFile         ( dataFile ),
    pathOut          ( pathOut ),
    predictOutputFile( predictOutputFile ),

    lib_str          ( lib_str ),
    pred_str         ( pred_str ),

    E                ( E ),
    Tp               ( Tp ),
    knn              ( knn ),
    tau              ( tau ),
    theta            ( theta ),
    exclusionRadius  ( exclusionRadius ),

    columns_str      ( columns_str ),
    target_str       ( target_str ),
    targetIndex      ( 0 ),

    embedded         ( embedded ),
    const_predict    ( const_predict ),
    verbose          ( verbose ),

    SmapOutputFile   ( SmapOutputFile ),
    blockOutputFile  ( blockOutputFile ),

    multiviewEnsemble     ( multiviewEnsemble ),
    multiviewD            ( multiviewD ),
    multiviewTrainLib     ( multiviewTrainLib ),
    multiviewExcludeTarget( multiviewExcludeTarget ),

    libSizes_str     ( libSizes_str ),
    subSamples       ( subSamples ),
    randomLib        ( randomLib ),
    replacement      ( replacement ),
    seed             ( seed ),
    includeData      ( includeData ),

    // Set validated flag and instantiate Version
    validated        ( false ),
    version          ( 1, 7, 1, "2020-11-29" )
{
    // Constructor code
    if ( method != Method::None ) {

        Validate();

        if ( verbose ) {
            version.ShowVersion();
        }
    }
}

//----------------------------------------------------------------
// Destructor
//----------------------------------------------------------------
Parameters::~Parameters() {}

//----------------------------------------------------------------
// Index offsets, generate library and prediction indices,
// and parameter validation
//----------------------------------------------------------------
void Parameters::Validate() {

    validated       = true;
    disjointLibrary = false;

    if ( not embedded and tau == 0 ) {
        std::string errMsg( "Parameters::Validate(): "
                            "tau must be non-zero.\n" );
        throw std::runtime_error( errMsg );
    }

    //--------------------------------------------------------------
    // Generate library indices: Apply zero-offset
    //--------------------------------------------------------------
    if ( lib_str.size() ) {
        // Parse lib_str into vector of strings
        std::vector<std::string> lib_vec = SplitString( lib_str, " \t," );
        if ( lib_vec.size() % 2 != 0 ) {
            std::string errMsg( "Parameters::Validate(): "
                                "library must be even number of integers.\n" );
            throw std::runtime_error( errMsg );
        }

        // Generate libPairs vector of start, stop index pairs
        std::vector< std::pair< size_t, size_t > > libPairs; // numeric bounds
        for ( size_t i = 0; i < lib_vec.size(); i = i + 2 ) {
            libPairs.emplace_back( std::make_pair( std::stoi( lib_vec[i] ),
                                                   std::stoi( lib_vec[i+1] ) ) );
        }
        if ( libPairs.size() > 1 ) {
            disjointLibrary = true;
        }

        size_t nLib = 0; // Count of lib items

        // Get number of lib indices, validate end > start
        for ( auto thisPair : libPairs ) {
            size_t lib_start = thisPair.first;
            size_t lib_end   = thisPair.second;

            nLib += lib_end - lib_start + 1;

            // Validate end > stop indices
            if ( method == Method::Simplex or method == Method::SMap ) {
                // Don't check if method == None, Embed or CCM since default
                // of "1 1" is used.
                if ( lib_start >= lib_end ) {
                    std::stringstream errMsg;
                    errMsg << "Parameters::Validate(): library start "
                           << lib_start << " exceeds end " << lib_end << ".\n";
                    throw std::runtime_error( errMsg.str() );
                }
            }
            // Disallow indices < 0, the user may have specified 0 start
            if ( lib_start < 1 or lib_end < 1 ) {
                std::stringstream errMsg;
                errMsg << "Parameters::Validate(): Library indices less than "
                       << "1 not allowed.\n";
                throw std::runtime_error( errMsg.str() );
            }
        }

        // Create library vector of indices
        library = std::vector< size_t >( nLib );
        size_t i = 0;
        for ( auto thisPair : libPairs ) {
            for ( size_t li = thisPair.first; li <= thisPair.second; li++ ) {
                library[ i ] = li - 1; // apply zero-offset
                i++;
            }
        }

        // If disjointLibrary create vector of disallowed library rows
        if ( disjointLibrary ) {
            int maxLibrary_i = (int) library.back();
            int minLibrary_i = (int) library.front();
            int allLibSize   = maxLibrary_i - minLibrary_i + 1;

            // A "full" library with no missing values, for set difference
            std::vector< size_t > allLibrary( allLibSize );

            // Fill with indices 0, 1, ... allLibSize - 1
            i = 0;
            for ( int lib_i = minLibrary_i; lib_i <= maxLibrary_i; lib_i++ ) {
                allLibrary[ i ] = lib_i;
                i++;
            }

            // Find library rows that are "gaps", not in library
            disjointLibraryRows = std::vector< size_t > ( allLibSize );
            std::vector< size_t >::iterator li;
            li = std::set_difference ( allLibrary.begin(), allLibrary.end(),
                                       library.begin(),    library.end(),
                                       disjointLibraryRows.begin());

            disjointLibraryRows.resize( li - disjointLibraryRows.begin() );

            // Add rows that span gaps by virture of embedding
            int shift = tau * ( E - 1 ); // default tau = -1
            std::vector< size_t > disjointLibraryRowsCopy( disjointLibraryRows );
            
            for ( auto disjointLibRow : disjointLibraryRowsCopy ) {
                int thisRow = (int) disjointLibRow;
                int newRow  = embedded ? thisRow : thisRow - shift;

                if ( newRow > 0 ) {
                    li = std::find( disjointLibraryRowsCopy.begin(),
                                    disjointLibraryRowsCopy.end(), newRow );
                    if ( li == disjointLibraryRowsCopy.end() ) {
                        // newRow not in disjointLibraryRows
                        disjointLibraryRows.push_back( (size_t) newRow );
                    }
                }
            }

            // Remove from library
            std::vector< size_t >::iterator dli;
            for ( dli  = disjointLibraryRows.begin();
                  dli != disjointLibraryRows.end(); dli++ ) {
                
                li = std::find( library.begin(), library.end(), *dli );
                if ( li != library.end() ) {
                    library.erase( li );
                }   
            }
        } // disjointLibrary

        if ( method == Method::Simplex or
             method == Method::SMap or method == Method::CCM ) {

            // Validate lib: E, tau, Tp combination
            int vectorStart  = std::max( (E - 1) * tau, 0 );
            vectorStart      = std::max( vectorStart, Tp );
            int vectorEnd    = std::min( (E - 1) * tau, Tp );
            vectorEnd        = std::min( vectorEnd, 0 );
            int vectorLength = std::abs( vectorStart - vectorEnd ) + 1;

            int maxLibrarySegment = 0;
            for ( auto thisPair : libPairs ) {
                int libPairSpan = thisPair.second - thisPair.first + 1;
                maxLibrarySegment = std::max( maxLibrarySegment, libPairSpan );
            }

            if ( vectorLength > maxLibrarySegment ) {
                std::stringstream errMsg;
                errMsg << "Parameters::Validate(): Combination of E = "
                       << E << " Tp = " << Tp << " tau = " << tau
                       << " is invalid.\n";
                throw std::runtime_error( errMsg.str() );
            }
        }
    }

    //--------------------------------------------------------------
    // Generate prediction indices: Apply zero-offset
    //--------------------------------------------------------------
    if ( pred_str.size() ) {
        // Parse pred_str into vector of strings
        std::vector<std::string> pred_vec = SplitString( pred_str, " \t," );
        if ( pred_vec.size() % 2 != 0 ) {
            std::string errMsg( "Parameters::Validate(): "
                                "prediction must be even number of integers.\n");
            throw std::runtime_error( errMsg );
        }

        // Generate vector of start, stop index pairs
        std::vector< std::pair< size_t, size_t > > predPairs;
        for ( size_t i = 0; i < pred_vec.size(); i = i + 2 ) {
            predPairs.emplace_back( std::make_pair( std::stoi( pred_vec[i] ),
                                                    std::stoi( pred_vec[i+1])));
        }

        size_t nPred = 0; // Count of pred items

        // Get number of pred indices, validate end > start
        for ( auto thisPair : predPairs ) {
            size_t pred_start = thisPair.first;
            size_t pred_end   = thisPair.second;

            nPred += pred_end - pred_start + 1;

            // Validate end > stop indices
            if ( method == Method::Simplex or method == Method::SMap ) {
                // Don't check if method == None, Embed or CCM since default
                // of "1 1" is used.
                if ( pred_start >= pred_end ) {
                    std::stringstream errMsg;
                    errMsg << "Parameters::Validate(): prediction start "
                           << pred_start << " exceeds end " << pred_end << ".\n";
                    throw std::runtime_error( errMsg.str() );
                }
            }
            // Disallow indices < 0, the user may have specified 0 start
            if ( pred_start < 1 or pred_end < 1 ) {
                std::stringstream errMsg;
                errMsg << "Parameters::Validate(): Prediction indices less than "
                       << "1 not allowed.\n";
                throw std::runtime_error( errMsg.str() );
            }
        }

        // Create prediction vector of indices
        prediction = std::vector< size_t >( nPred );
        size_t i = 0;
        for ( auto thisPair : predPairs ) {
            for ( size_t li = thisPair.first; li <= thisPair.second; li++ ) {
                prediction[ i ] = li - 1; // apply zero-offset
                i++;
            }
        }
    }

    if ( method == Method::Simplex or method == Method::SMap ) {
        if ( not library.size() ) {
            std::string errMsg( "Parameters::Validate(): "
                                "library indices not found.\n" );
            throw std::runtime_error( errMsg );
        }
        if ( not prediction.size() ) {
            std::string errMsg( "Parameters::Validate(): "
                                "prediction indices not found.\n" );
            throw std::runtime_error( errMsg );
        }
    }
    else {
        // Defaults if Method is None, Embed or CCM
        if ( not library.size() ) {
            library = std::vector<size_t>( 1, 0 );
        }
        if ( not prediction.size() ) {
            prediction = std::vector<size_t>( 1, 0 );
        }
    }

#ifdef DEBUG_ALL
    PrintIndices( library, prediction );
#endif

    //--------------------------------------------------------------
    // Convert multi argument parameters from string to vectors
    //--------------------------------------------------------------

    //--------------------------------------------------------------
    // Columns
    // If columns are purely integer, set vector<size_t> columnIndex
    // Otherwise fill in vector<string> columnNames
    //--------------------------------------------------------------
    if ( columns_str.size() ) {

        std::vector<std::string> columns_vec = SplitString( columns_str,
                                                            " \t,\n" );

        bool onlyDigits = false;

        for ( auto ci = columns_vec.begin(); ci != columns_vec.end(); ++ci ) {
            onlyDigits = OnlyDigits( *ci, true );
            
            if ( not onlyDigits ) { break; }
        }

        if ( onlyDigits ) {
            for ( auto ci =  columns_vec.begin();
                       ci != columns_vec.end(); ++ci ) {
                columnIndex.push_back( std::stoi( *ci ) );
            }
        }
        else {
            columnNames = columns_vec;
        }
    }

    if ( not columnIndex.size() and not columnNames.size() ) {
        std::stringstream errMsg;
        errMsg << "Parameters::Validate(): Simplex/CCM: "
               << " No valid columns found." << std::endl;
        throw std::runtime_error( errMsg.str() );
    }

    //--------------------------------------------------------------
    // target
    //--------------------------------------------------------------
    if ( target_str.size() ) {
        bool onlyDigits = OnlyDigits( target_str, true );
        if ( onlyDigits ) {
            targetIndex = std::stoi( target_str );
        }
        else {
            targetName = target_str;
        }
    }

    //--------------------------------------------------------------------
    // CCM sample not 0 if random is true
    //--------------------------------------------------------------
    if ( method == Method::CCM ) {
        if ( randomLib ) {
            if ( subSamples < 1 ) {
                std::string errMsg( "Parameters::Validate(): "
                                    "CCM samples must be > 0.\n" );
                throw std::runtime_error( errMsg );
            }
        }
    }

    //--------------------------------------------------------------------
    // CCM librarySizes
    //--------------------------------------------------------------------
    if ( libSizes_str.size() > 0 ) {
        std::vector<std::string> libsize_vec = SplitString(libSizes_str," \t,");

        bool   libSizeSequence = false;
        size_t start;
        size_t stop;
        size_t increment;
        
        if ( libsize_vec.size() == 3 ) {
            // Presume ( start, stop, increment ) sequence arguments
            start     = std::stoi( libsize_vec[0] );
            stop      = std::stoi( libsize_vec[1] );
            increment = std::stoi( libsize_vec[2] );

            // However, it might just be 3 library sizes...
            // If increment < stop, then presume 3 sequence arguments
            if ( increment < stop ) {
                libSizeSequence = true;
            }
        }

        if ( libSizeSequence ) {
            if ( increment < 1 ) {
                std::stringstream errMsg;
                errMsg << "Parameters::Validate(): "
                       << "CCM librarySizes increment " << increment
                       << " is invalid.\n";
                throw std::runtime_error( errMsg.str() );
            }

            if ( start > stop ) {
                std::stringstream errMsg;
                errMsg << "Parameters::Validate(): "
                       << "CCM librarySizes start " << start
                       << " stop " << stop  << " are invalid.\n";
                throw std::runtime_error( errMsg.str() );
            }

            size_t N_lib = std::floor((stop-start)/increment + 1/increment) + 1;

            if ( (int) start < E ) {
                std::stringstream errMsg;
                errMsg << "Parameters::Validate(): "
                       << "CCM librarySizes start < E = " << E << "\n";
                throw std::runtime_error( errMsg.str() );
            }
            else if ( (int) start < 3 ) {
                std::string errMsg( "Parameters::Validate(): "
                                    "CCM librarySizes start < 3.\n" );
                throw std::runtime_error( errMsg );
            }

            // Allocate the librarySizes vector
            librarySizes = std::vector< size_t >( N_lib, 0 );

            // Fill in the sizes
            size_t libSize = start;
            for ( size_t i = 0; i < librarySizes.size(); i++ ) {
                librarySizes[i] = libSize;
                libSize = libSize + increment;
            }
        }
        else {
            // Presume a list of lib sizes
            librarySizes = std::vector< size_t >( libsize_vec.size(), 0 );
            for ( size_t i = 0; i < librarySizes.size(); i++ ) {
                size_t libSize;
                std::stringstream libStringStream( libsize_vec[i] );
                libStringStream >> libSize;
                librarySizes[i] = libSize;
            }
        }
    }

    //--------------------------------------------------------------------
    // Simplex: knn not specified: knn set to E+1
    //    embedded true : E set to number columns
    //--------------------------------------------------------------------
    if ( method == Method::Simplex or method == Method::CCM ) {

        // embedded = true: Set E to number of columns if not already set
        if ( embedded ) {
            if ( columnIndex.size() ) {
                if ( E != (int) columnIndex.size() ) {
                    E = columnIndex.size();
                }
            }
            else if ( columnNames.size() ) {
                if ( E != (int) columnNames.size() ) {
                    E = columnNames.size();
                }
            }
        }

        if ( knn < 1 ) {
            knn = E + 1;

            if ( verbose ) {
                std::stringstream msg;
                msg << "Parameters::Validate(): Set knn = " << knn
                    << " (E+1) for Simplex. " << std::endl;
                std::cout << msg.str();
            }
        }

        if ( knn < E + 1 ) {
            std::stringstream errMsg;
            errMsg << "Parameters::Validate(): Simplex knn of " << knn
                   << " is less than E+1 = " << E + 1 << std::endl;
            throw std::runtime_error( errMsg.str() );
        }
    }
    //--------------------------------------------------------------------
    // SMap and knn not specified: knn set to library.size()
    // embedded true : E set to size( columns ) for output processing
    //--------------------------------------------------------------------
    else if ( method == Method::SMap ) {
        if ( embedded and columnNames.size() > 1 ) {
            if ( columnIndex.size() ) {
                E = columnIndex.size();
            }
            else if ( columnNames.size() ) {
                E = columnNames.size();
            }
        }

        if ( knn == 0 ) { // default knn = 0, set knn value
            knn = library.size() - E;

            if ( verbose ) {
                std::stringstream msg;
                msg << "Parameters::Validate(): Set knn = " << knn
                    << " for SMap. " << std::endl;
                std::cout << msg.str();
            }
        }
        if ( knn < 2 ) {
            std::stringstream errMsg;
            errMsg << "Parameters::Validate() S-Map knn must be at least 2.\n";
            throw std::runtime_error( errMsg.str() );
        }

        if ( not embedded and columnNames.size() > 1 ) {
            std::stringstream errMsg( "Parameters::Validate(): "
                                      "Multivariable S-Map must use "
                                      "embedded = true to ensure "
                                      "data/dimension correspondance.\n" );
            throw std::runtime_error( errMsg.str() );
        }
    }
    else if ( method == Method::Embed ) {
        // no-op
    }
    else {
        throw std::runtime_error( "Parameters::Validate() "
                                  "Prediction method error.\n" );
    }

    if ( method == Method::Simplex or method == Method::SMap ) {
        if ( E < 1 ) {
            std::stringstream errMsg;
            errMsg << "Parameters::Validate() E = " << E
                   << " is invalid with embedded = true.\n" ;
            throw std::runtime_error( errMsg.str() );
        }
    }
}

//------------------------------------------------------------
// Adjust lib/pred concordant with Embed() removal of tau(E-1)
// rows, and DeletePartialDataRow()
//------------------------------------------------------------
void Parameters::DeleteLibPred() {

    size_t shift          = abs( tau ) * ( E - 1 );
    size_t library_len    = library.size();
    size_t prediction_len = prediction.size();

    // If [0, 1, ... shift]  (negative tau) or
    // [N-shift, ... N-1, N] (positive tau) are in library or prediction
    // those rows were deleted, delete these index elements.
    // First, create vectors of indices to delete.
    std::vector< size_t > deleted_pred_elements( shift, 0 );
    std::vector< size_t > deleted_lib_elements ( shift, 0 );

    if ( tau < 0 ) {
        std::iota(deleted_pred_elements.begin(), deleted_pred_elements.end(),0);
        std::iota(deleted_lib_elements.begin(),  deleted_lib_elements.end(), 0);
    }
    else {
        std::iota( deleted_pred_elements.begin(),
                   deleted_pred_elements.end(), prediction_len - shift );
        std::iota( deleted_lib_elements.begin(),
                   deleted_lib_elements.end(), library_len - shift );
    }

    // Now that we have the indices that could have been deleted,
    // check to see if any are in lib and pred
    bool deleteLibIndex = false;
    std::vector< size_t >::iterator it;
    for ( auto element  = deleted_lib_elements.begin();
               element != deleted_lib_elements.end(); element++ ) {
        it = std::find( library.begin(), library.end(), *element );
        if ( it != library.end() ) {
            deleteLibIndex = true;
            break;
        }
    }

    bool deleteDisjointLibIndex = false;
    for ( auto element  = deleted_lib_elements.begin();
               element != deleted_lib_elements.end(); element++ ) {
        it = std::find( disjointLibraryRows.begin(),
                        disjointLibraryRows.end(), *element );
        if ( it != prediction.end() ) {
            deleteDisjointLibIndex = true;
            break;
        }
    }

    bool deletePredIndex = false;
    for ( auto element  = deleted_pred_elements.begin();
               element != deleted_pred_elements.end(); element++ ) {
        it = std::find( prediction.begin(), prediction.end(), *element );
        if ( it != prediction.end() ) {
            deletePredIndex = true;
            break;
        }
    }

#ifdef DEBUG_ALL
    std::cout << "DeleteLibPred(): shift: " << shift
              << " deleteLibIndex: " << deleteLibIndex
              << " deletePredIndex: " << deletePredIndex << std::endl;
    std::cout << " deleted_lib_elements: ";
    for ( auto element  = deleted_lib_elements.begin();
               element != deleted_lib_elements.end(); element++ ) {
        std::cout << *element << ", ";
    } std::cout << std::endl << " deleted_pred_elements: ";
    for ( auto element  = deleted_pred_elements.begin();
               element != deleted_pred_elements.end(); element++ ) {
        std::cout << *element << ", ";
    } std::cout << std::endl;
#endif

    // Erase elements of row indices that were deleted
    if ( deleteLibIndex ) {
        for ( auto element  = deleted_lib_elements.begin();
                   element != deleted_lib_elements.end(); element++ ) {

            std::vector< size_t >::iterator it;
            it = std::find( library.begin(), library.end(), *element );
            
            if ( it != library.end() ) {
                library.erase( it );
            }
        }
    }

    if ( deleteDisjointLibIndex ) {
        for ( auto element  = deleted_lib_elements.begin();
                   element != deleted_lib_elements.end(); element++ ) {

            std::vector< size_t >::iterator it;
            it = std::find( disjointLibraryRows.begin(),
                            disjointLibraryRows.end(), *element );
            
            if ( it != disjointLibraryRows.end() ) {
                disjointLibraryRows.erase( it );
            }
        }
    }

    if ( deletePredIndex ) {
        for ( auto element  = deleted_pred_elements.begin();
              element != deleted_pred_elements.end(); element++ ) {

            std::vector< size_t >::iterator it;
            it = std::find( prediction.begin(),
                            prediction.end(), *element );

            if ( it != prediction.end() ) {
                prediction.erase( it );
            }
        }
    }

    // Now offset all values by shift so that vectors indices
    // in library and prediction refer to the same data rows
    // before the deletion/shift.
    if ( tau < 0 ) {
        for ( auto li = library.begin(); li != library.end(); li++ ) {
            *li = *li - shift;
        }
        for ( auto li  = disjointLibraryRows.begin();
                   li != disjointLibraryRows.end(); li++ ) {
            *li = *li - shift;
        }
        for ( auto pi = prediction.begin(); pi != prediction.end(); pi++ ) {
            *pi = *pi - shift;
        }
    }
    // tau > 0  : Forward shifting: no adjustment needed from origin

    return;
}

//------------------------------------------------------------------
// Overload << to output to ostream
//------------------------------------------------------------------
std::ostream& operator<< ( std::ostream &os, Parameters &p ) {

    // print info about the dataframe
    os << "Parameters: -------------------------------------------\n";

    std::string method("Unknown");
    if      ( p.method == Method::Simplex ) { method = "Simplex"; }
    else if ( p.method == Method::SMap    ) { method = "SMap";    }
    else if ( p.method == Method::CCM     ) { method = "CCM";     }
    else if ( p.method == Method::None    ) { method = "None";    }
    else if ( p.method == Method::Embed   ) { method = "Embed";   }

    os << "Method: " << method
       << " E=" << p.E << " Tp=" << p.Tp
       << " knn=" << p.knn << " tau=" << p.tau << " theta=" << p.theta
       << std::endl;

    if ( p.columnNames.size() ) {
        os << "Column Names : [ ";
        for ( auto ci = p.columnNames.begin();
              ci != p.columnNames.end(); ++ci ) {
            os << *ci << " ";
        } os << "]" << std::endl;
    }

    if ( p.targetName.size() ) {
        os << "Target: " << p.targetName << std::endl;
    }

    os << "Library: [" << p.library[0] << " : "
       << p.library[ p.library.size() - 1 ] << "]  "
       << "Prediction: [" << p.prediction[0] << " : "
       << p.prediction[ p.prediction.size() - 1 ]
       << "] " << std::endl;

    os << "-------------------------------------------------------\n";

    return os;
}

#ifdef DEBUG_ALL
//------------------------------------------------------------------
// 
//------------------------------------------------------------------
void Parameters::PrintIndices( std::vector<size_t> library,
                               std::vector<size_t> prediction )
{
    std::cout << "Parameters(): library: ";
    for ( auto li = library.begin(); li != library.end(); ++li ) {
        std::cout << *li << " ";
    } std::cout << std::endl;
    std::cout << "Parameters(): prediction: ";
    for ( auto pi = prediction.begin(); pi != prediction.end(); ++pi ) {
        std::cout << *pi << " ";
    } std::cout << std::endl;
}
#endif
