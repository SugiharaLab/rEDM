#include "RcppEDMCommon.h"

//-----------------------------------------------------------------------
// Convert R DataFrame to cppEDM DataFrame<double>
//-----------------------------------------------------------------------
DataFrame< double > DFToDataFrame ( Rcpp::DataFrame df ) {

    // Get number of valarray rows from first pair
    size_t numRows = df.nrow();

    // ensure that we have > 1 column for reading
    if ( df.ncol() == 1 ) {
        std::string err = "DFToDataFrame(): Input must have > 1 column, "
                          "first column is interpreted as a time vector.\n";
        throw Rcpp::exception( err.c_str() );
    }

    // Get column names
    // JP: Are df.names() ensured to be in order as accessed by index?
    //     If not, this will give incorrect results.
    std::vector< std::string > colNames;
    r::CharacterVector tmp_colNames = df.names();

    for ( int idx = 1; idx < tmp_colNames.size(); idx++ ) {
        colNames.push_back( r::as<std::string>( tmp_colNames[idx] ) );
    }

    // Create cpp DataFrame
    DataFrame< double > dataFrame ( numRows, df.ncol()-1, colNames ); 

    // Setup time column and time name for dataframe
    // It is assumed that the first column is a time vector !!!
    r::CharacterVector tmp = r::as<r::CharacterVector>( df[0] );
    dataFrame.Time()       = r::as< std::vector<std::string> >( tmp );
    dataFrame.TimeName()   = r::as<std::string>( 
                             ((r::CharacterVector)df.names())[0] );  

    // read in data columns to the cppEDM DF
    // JP: Are df.names() ensured to be in order as accessed by index?
    //     If not, this will give incorrect results.
    for ( int idx = 1; idx < df.ncol(); idx++ ) {
        // unfortunately we can't convert numeric vec to valarray
        std::vector<double> tmp = r::as<std::vector<double>>( df[idx] );
        std::valarray<double> col ( tmp.data(), tmp.size() );
        dataFrame.WriteColumn( idx-1, col ); 
    }

    return dataFrame;
}

//---------------------------------------------------------------
// Convert cppEDM DataFrame<double> to R DataFrame
//---------------------------------------------------------------
r::DataFrame DataFrameToDF ( DataFrame< double > dataFrame ) {

    r::List columnList; // List of columns to create new R data.frame

    // NOTE: cppEDM DataFrame columnNames are data only, not time
    std::vector<std::string> columnNamesIn = dataFrame.ColumnNames();

#ifdef _WIN32 // Rcpp UTF8 encoding fails on Windows
    r::StringVector columnNames; 
#else
    std::vector<std::string> columnNames;
#endif

    bool hasTime = false;

    // If dataFrame has time vector and timeName, add to columnList
    if ( dataFrame.Time().size() ) {
        hasTime = true; // Skip time column in dataFrame.VectorColumnName()

#ifdef _WIN32 // Rcpp UTF8 encoding fails on Windows
        r::String timeName( dataFrame.TimeName() );
        columnNames.push_back( timeName );
#else
        columnNames.push_back( dataFrame.TimeName() );
#endif

        // Probe dataFrame.Time() to see if we can convert it to
        // a numeric, Date, or Datetime...
        std::string firstTime = dataFrame.Time()[0];

        // Is firstTime purely numeric characters (not Date or DateTime)?
        // We presume time is not negative, or exponential 
        bool numericTime = strspn( firstTime.c_str(),
                                   ".0123456789" ) == firstTime.size();

        // Does firstTime have two hyphens as in "%Y-%m-%d" Date format?
        size_t nHyphen  = std::count(firstTime.begin(), firstTime.end(), '-');
        bool   dateTime = nHyphen == 2 ? true : false;

        // Does firstTime have two '-' and two ':' as in DateTime format?
        size_t nColon     = std::count(firstTime.begin(), firstTime.end(), ':');
        bool dateTimeTime = dateTime and nColon == 2 ? true : false;
        if ( dateTimeTime ) { dateTime = false; }

        if ( numericTime and not dateTime and not dateTimeTime ) {
            // Convert the dataFrame.Time() vector to numeric/double
            r::NumericVector timeVec( dataFrame.Time().size() );

            char *pEnd;
            for ( size_t i = 0; i < dataFrame.Time().size(); i++ ) {
                timeVec[ i ] = strtod( dataFrame.Time().at( i ).c_str(), &pEnd );
                // JP: check pEnd?
            }
            columnList.push_back( timeVec );
        }

        if ( not numericTime and dateTime and not dateTimeTime ) {
            // Convert to Date
            r::DateVector dateVec( dataFrame.Time().size() );
            
            for ( size_t i = 0; i < dataFrame.Time().size(); i++ ) {
                dateVec[ i ] = r::Date( dataFrame.Time().at( i ),
                                        "%Y-%m-%d" );
            }
            columnList.push_back( dateVec );
        }

        if ( not numericTime and not dateTime and dateTimeTime )  {
            // Convert to Datetime
            r::DatetimeVector datetimeVec( dataFrame.Time().size() );

            for ( size_t i = 0; i < dataFrame.Time().size(); i++ ) {
                datetimeVec[ i ] = r::Datetime( dataFrame.Time().at( i ),
                                                "%Y-%m-%d %H:%M:%OS" );
            }
            columnList.push_back( datetimeVec );
        }

        if ( not numericTime and not dateTime and not dateTimeTime )  {
            // Couldn't convert dataFrame.Time(), just push it as-is
            // R will see it as a vector of factors... why not strings?
            columnList.push_back( dataFrame.Time() );
        }
    } // if ( dataFrame.Time().size() ) 

    // Copy data: NOTE in cppEDM data and time vector are separate
    //            data are in a row-major valarray with NColumns().
    for ( auto ci = columnNamesIn.begin(); ci != columnNamesIn.end(); ci++ ) {
        if ( hasTime and (*ci).compare( dataFrame.TimeName() ) == 0 ) {
            continue;  // skip time. It's a vector< std::string > 
        }

        // Unfortunately we have to copy to vector first
        std::valarray<double> col_val = dataFrame.VectorColumnName( *ci );
        std::vector<double> col_vec(std::begin(col_val), std::end(col_val));
        columnList.push_back( col_vec );

#ifdef _WIN32 // Rcpp UTF8 encoding fails on Windows
        r::String colName( *ci );
        colName.replace_all( "∂", "d" ); // ∂ set in cppEDM/src/SMap.cc
        // colName.set_encoding( cetype_t::CE_UTF8 );
        columnNames.push_back( colName );
#else
        columnNames.push_back( *ci );
#endif
    }

    r::DataFrame df ( columnList );
    df.attr("names") = columnNames;

    return df;
}

//---------------------------------------------------------------
// Load path/file into cppEDM DataFrame, convert to Python
// dict{ column : array }
//---------------------------------------------------------------
r::DataFrame ReadDataFrame ( std::string path, std::string file ) {
    return DataFrameToDF( DataFrame< double >( path, file ) );
}
