
library( rEDM )

#-------------------------------------------------------------------
# EDM Embed wrapper
# Create time-delay embedding with time column for EDM. 
# Useful to create mixed multivariate embeddings for SMap and
# embeddings with time-advanced vectors.  
# Rename V(t-0), V(t+0) to V. Add Time column.
# If columns is NULL, embedd all except the first (time) column.
# If plusminus create time-advanced & time-delayed columns.
#-------------------------------------------------------------------
Embedding = function(
  dataFrame = NULL,
  dataFile  = NULL,        
  outFile   = NULL,
  plusminus = FALSE,
  columns   = NULL,
  E         = 2,
  tau       = -1,
  verbose   = FALSE
) {

  if ( is.null( dataFrame ) & is.null( dataFile ) ) {
    stop( 'dataFrame and dataFile are empty, specify one.' )
  }
  if ( tau > 0 & plusminus ) {
    # Convert to negative
    tau = -tau
  }

  if ( is.null( dataFrame ) ) {
    # Load from dataFile
    data = read.csv( dataFile )
  }
  else {
    data = dataFrame
  }

  # Presume time is first column
  timeName   = colnames( data )[1]
  timeSeries = data[ , timeName ]

  # If no columns specified, use all except first
  if ( is.null( columns ) ) {
    columns = colnames( data )[ 2 : ncol( data ) ]
  }

  if ( verbose ) {
    print( paste( "Time column: ", timeName ) )
    print( "Embed columns: " ); print( columns )
  }

  # Create embeddings of columns
  # There will be redundancies vis V1(t-0), V1(t+0)
  if ( plusminus ) {
    embed_minus = Embed( dataFrame = data, E = E, tau = tau, columns = columns )
    embed_plus  = Embed( dataFrame = data, E = E, tau = abs( tau ),
                         columns = columns )
    embed = cbind( timeSeries, embed_minus, embed_plus, stringsAsFactors=FALSE )

    # TRUE / FALSE vector
    cols_tplus0 = grepl( '(t+0)', colnames( embed ), fixed = TRUE )
    # Remove *(t+0) : redunant with *(t-0)
    embed = embed[ , !cols_tplus0 ]
  }
  else {
    embed_ = Embed( dataFrame = data, E = E, tau = tau, columns = columns )
    embed  = cbind( timeSeries, embed_, stringsAsFactors = FALSE )
  }

  # Rename *(t-0) to original column names
  columnNames = colnames( embed )
  for ( i in 1:length( columnNames ) ) {
    if ( grepl( '(t-0)', columnNames[i], fixed = TRUE ) ) {
      columnNames[i] = sub( '(t-0)', '', columnNames[i], fixed = TRUE )
    }
  }

  # Rename *(t+0) to original column names
  for ( i in 1:length( columnNames ) ) {
    if ( grepl( '(t+0)', columnNames[i], fixed = TRUE ) ) {
      columnNames[i] = sub( '(t+0)', '', columnNames[i], fixed = TRUE )
    }
  }

  # Rename first column to original time column name
  columnNames[ 1 ]  = timeName
  colnames( embed ) = columnNames

  if ( verbose ) {
    print( head( embed, 4 ) )
    print( tail( embed, 4 ) )
  }

  if ( ! is.null( outFile ) ) {
    write.csv( embed, file = outFile, row.names = FALSE )
  }

  return( embed )
}
