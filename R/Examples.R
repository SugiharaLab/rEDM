
#------------------------------------------------------------------------
# 
#------------------------------------------------------------------------
Examples = function() {

  library( rEDM )
  
  # make sure data is loaded
  tryCatch(
    expr = {
      data( TentMap,             envir = environment() )
      data( TentMapNoise,        envir = environment() )
      data( block_3sp,           envir = environment() )
      data( circle,              envir = environment() )
      data( sardine_anchovy_sst, envir = environment() )
    },
    error = function( err ) {
      print( err )
      stop("Examples(): Failed to load package data.")
    }
  )

  if ( is.null( dev.list() ) ) {
    newPlot( mfrow = c( 4, 2 ) )
  }
  else {
    par( mfrow = c( 4, 2 ) )
  }
  
  # EmbedDimension()
  cmd = paste0('EmbedDimension( dataFrame=TentMap,',
               ' lib="1 100", pred="201 500",',
               ' columns="TentMap", target="TentMap")' )
  df = eval( parse( text = cmd ) )
  
  # PredictInterval()
  cmd = paste0('PredictInterval( dataFrame=TentMap,',
               ' lib="1 100", pred="201 500", E = 2,',
               ' columns="TentMap", target="TentMap") ')
  df = eval( parse( text = cmd ) )
  
  # PredictNonlinear()
  cmd = paste0('PredictNonlinear( dataFrame=TentMapNoise,',
               ' E=2,lib="1 100", pred="201 500", ',
               ' columns="TentMap",target="TentMap") ')
  df = eval( parse( text = cmd ) )

  # Simplex() 
  # Tent map : specify multivariable columns embedded = TRUE
  cmd = paste0('Simplex( dataFrame=block_3sp,',
               ' lib="1 99", pred="100 195", ',
               ' E=3, embedded=TRUE, showPlot=TRUE, const_pred=TRUE,',
               ' columns="x_t y_t z_t", target="x_t") ')
  df = eval( parse( text = cmd ) )

  # Simplex() 
  # Tent map : Embed column x_t to E=3, embedded = False
  cmd = paste0('Simplex( dataFrame=block_3sp,',
               ' lib="1 99", pred="105 190", ',
               ' E=3, showPlot=TRUE, const_pred=TRUE,',
               ' columns="x_t", target="x_t") ')
  df = eval( parse( text = cmd ) )

  # Multiview()
  cmd = paste0('Multiview( dataFrame=block_3sp,',
               ' lib="1 99", pred="105 190", ',
               ' E=3, columns="x_t y_t z_t", target="x_t",',
               ' showPlot=TRUE) ')
  df = eval( parse( text = cmd ) )

  # CCM demo
  cmd = paste0('CCM( dataFrame=sardine_anchovy_sst,',
               ' E=3, Tp=0, columns="anchovy", target="np_sst",',
               ' libSizes="10 70 10", sample=100, verbose=TRUE, ',
               ' showPlot=TRUE) ')
  df = eval( parse( text = cmd ) )

  par( ask = TRUE )
  # SMap circle : specify multivariable columns embedded = TRUE
  cmd = paste0('SMap( dataFrame=circle,',
               ' lib="1 100", pred="110 190", theta=4, E=2,',
               ' verbose=TRUE, showPlot=TRUE, embedded=TRUE,',
               ' columns="x y", target="x") ')
  df = eval( parse( text = cmd ) )
  par( ask = FALSE )
}

#------------------------------------------------------------------------
# 
#------------------------------------------------------------------------
newPlot = function(
  mar      = c( 4, 4, 1, 1 ),
  mgp      = c( 2.3, 0.8, 0 ),
  cex      = 1.5, 
  cex.axis = 1.3,
  cex.lab  = 1.3,
  mfrow    = c( 1, 1 )
) {
    dev.new()
    par( mar = mar )
    par( mgp = mgp )
    par( cex = cex )
    par( cex.axis = cex.axis )
    par( cex.lab = cex.lab )
    par( mfrow = mfrow )
    invisible( dev.cur() )
}
