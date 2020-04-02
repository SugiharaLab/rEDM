#------------------------------------------------------------------------
# 
#------------------------------------------------------------------------
ComputeError <- function( obs, pred ) {
  # Pearson rho, RMSE, MAE.
  return ( RtoCpp_ComputeError( obs, pred ) )
}

#------------------------------------------------------------------------
# Validate columns are present
#------------------------------------------------------------------------
ColumnsInDataFrame <- function( pathIn, dataFile, dataFrame, columns, target ) {

  if ( nchar( dataFile ) ) {
    # Shame to have to read the data just for this...
    # Perhaps pass the df back ?  Gets messy.
    #
    # R data.frame names are presumed to be valid variable names.
    #   Numeric fist character or hyphen/dash/minus "-" are not valid.
    #   So we disable check.names that calls make.names() on column names.
    df = read.csv( paste( pathIn, dataFile, sep = '/' ),
                   as.is = TRUE, check.names = FALSE )
  }
  else {
    df = dataFrame
  }

  if ( ! isValidDF( df ) ) {
    print( "Error: ColumnsInDataFrame(): dataFrame is not valid." )
    return( FALSE )
  }

  columnNames = names( df )
  columnVec   = strsplit( columns, "\\s+" )[[1]] # split on whitespace
  
  if ( ! (target %in% columnNames) ) {
    print( paste( "Error: ColumnsInDataFrame(): Target",
                  target, "not found." ) )
    return( FALSE )
  }
  
  for ( column in columnVec ) {
    if ( length( df[,column] ) == 0 ) {
      print( paste("Error: ColumnsInDataFrame(): Column", column, "is empty."))
      return( FALSE )
    }

    if ( ! (column %in% columnNames) ) {
      print( paste( "Error: ColumnsInDataFrame(): Column",
                    column, "not found." ) )
      return( FALSE )
    }
  }

  return( TRUE )
}

#------------------------------------------------------------------------
# Is dataFrame a non-empty data.frame?  TRUE : FALSE
#------------------------------------------------------------------------
isValidDF <- function( dataFrame ) {
  if ( inherits( dataFrame, "data.frame" ) ) {
    if ( nrow( dataFrame ) == 0 || ncol( dataFrame ) == 0 ) { 
      print( paste( "isValidDF(): dataFrame is empty." ) )
      return( FALSE )
    }
    return( TRUE )
  }
  else {
    return( FALSE )
  }
}

#------------------------------------------------------------------------
# Plot data.frame with "time" "Observations" "Predictions"
#------------------------------------------------------------------------
PlotObsPred <- function( df,
                         dataFile = NULL,
                         E        = NULL,
                         Tp       = NULL ) {

  if ( ncol( df ) < 3 ) {
    print( "PlotObsPred: at least 3 columns are expected." )
    return( 0 )
  }
  if ( ! "Observations" %in% names( df ) ) {
    print( "PlotObsPred: unable to find Observations." )
    return( 0 )
  }
  if ( ! "Predictions" %in% names( df ) ) {
    print( "PlotObsPred: unable to find Predictions." )
    return( 0 )
  }

  # Try to convert first column to Date or POSIXlt or numeric
  time = NULL
  if ( is.numeric( df[,1] ) ) {
    time = df[,1]
  }
  else {
    time = try( as.Date( df[,1] ), silent = TRUE )
    if ( "try-error" %in% class( time ) ) {
      time = try( as.POSIXlt( df[,1] ), silent = TRUE )
      if ( "try-error" %in% class( time ) ) {
        time = try(as.numeric(levels(df[,1]))[df[,1]], silent = TRUE)
      }
    }
  }
  if ( "try-error" %in% class( time ) ) {
    # Create a bogus time vector
    time = seq( 1, nrow( df ) )
  }

  # stats: {'MAE': 0., 'RMSE': 0., 'rho': 0. }
  stats = ComputeError( df $ Observations,
                        df $ Predictions )
  
  title = paste( "\nE=", E, " Tp=", Tp,
                 " rho=",  round( stats[['rho']],  2 ),    
                 " RMSE=", round( stats[['RMSE']], 2 ) )
  
  plot( time, df $ Observations, main = title,
        xlab = names(df)[1], ylab = "",
        type = "l", col = "blue", lwd = 3,
        cex.axis = 1.3, cex.lab = 1.3 )
  
  lines( time, df $ Predictions, col = "red", lwd = 3 )
  
  legend( 'topright', c( "Predictions", "Observations" ), 
          fill = c('red', 'blue' ), bty = 'n', cex = 1.2 )
}

#------------------------------------------------------------------------
# Plot S-Map coefficients
#------------------------------------------------------------------------
PlotSmap <- function( SmapList,
                      dataFile = NULL,
                      E        = NULL,
                      Tp       = NULL ) {

  if ( ! "predictions" %in% names( SmapList ) ) {
    print( "PlotSmap: unable to find predictions." )
    return( 0 )
  }
  if ( ! "coefficients" %in% names( SmapList ) ) {
    print( "PlotSmap: unable to find coefficients." )
    return( 0 )
  }

  p = SmapList[[ "predictions"  ]]
  c = SmapList[[ "coefficients" ]]

  if ( ncol( p ) < 3 ) {
    print( "PlotSmap: expected at least 3 columns in predictions." )
    return( 0 )
  }
  
  # Try to convert first column to Date or POSIXlt or numeric
  time = NULL
  if ( is.numeric( p[,1] ) ) {
    time = p[,1]
  }
  else {
    time = try( as.Date( p[,1] ), silent = TRUE )
    if ( "try-error" %in% class( time ) ) {
      time = try( as.POSIXlt( p[,1] ), silent = TRUE )
      if ( "try-error" %in% class( time ) ) {
        time = try(as.numeric(levels(p[,1]))[p[,1]], silent = TRUE)
      }
    }
  }
  if ( "try-error" %in% class( time ) ) {
    # Create a bogus time vector
    time = seq( 1, nrow( p ) )
  }

  numCoeff = ncol( c ) - 1 
  
  old.par = par( no.readonly = TRUE )
  
  par( mfrow = c( numCoeff + 1, 1 ), mar = c( 3.5, 4, 0.5, 1 ),
       mgp = c( 1.5, 0.5, 0 ), cex.axis = 1.3, cex.lab = 1.3 )
  
  # Observations & Predictions
  plot( time, p $ Observations,
        xlab = names(p)[1], ylab = "",
        type = "l", col = "blue", lwd = 3,
        cex.axis = 1.3, cex.lab = 1.3 )
  
  lines( time, p $ Predictions, col = "red", lwd = 3 )
  legend( 'topright', c( "Predictions", "Observations" ), 
          fill = c('red', 'blue' ), bty = 'n', cex = 1.5 )

  # Coefficients
  title = paste( dataFile , 'S-Map Coefficients', '\nE=' , E, ' Tp=', Tp )

  coefName = names( c )
  for ( coef in 2:ncol(c) ) {
    plot( time, c[,coef], xlab = "Time", ylab = coefName[coef],
          col = "blue", type = "l", lwd = 3 )
    mtext( title, side = 3, line = -1.5, cex = 1.2 )
  }

  par( old.par )
}

#------------------------------------------------------------------------
# 
# Generate surrogate data for permutation/randomization tests
# 
# Method "random_shuffle" creates surrogates by randomly permuting the values 
# of the original time series.
# 
# Method "Ebisuzaki" creates surrogates by randomizing the phases of a Fourier 
# transform, preserving the power spectra of the null surrogates.
# 
# Method "seasonal" creates surrogates by computing a mean seasonal trend of 
# the specified period and shuffling the residuals.
# 
#------------------------------------------------------------------------
SurrogateData <- function(
  ts,
  method = c("random_shuffle", "ebisuzaki", "seasonal"), 
  num_surr = 100, T_period = 1, alpha = 0 )
{
  
  method = match.arg(method)
  if( method == "random_shuffle" ) {
    return( sapply( 1:num_surr, function(i) {
      sample(ts, size = length(ts))
    } ) )
  }
  else if( method == "ebisuzaki" ) {
    if( any( ! is.finite(ts) ) ) {
      stop("SurrogateData(): input time series contained invalid values")
    }
    
    n  = length(ts)
    n2 = floor(n/2)
    
    mu    = mean(ts)
    sigma = sd(ts)
    a     = fft(ts)
    amplitudes    = abs(a)
    amplitudes[1] = 0
    
    return( sapply(1:num_surr, function(i) {
      if(n %% 2 == 0) # even length
      {
        thetas   = 2*pi*runif(n2-1)
        angles   = c(0, thetas, 0, -rev(thetas))
        recf     = amplitudes * exp(complex(imaginary = angles))
        recf[n2] = complex(real = sqrt(2) * amplitudes[n2] * cos(runif(1)*2*pi))
      }
      else # odd length
      {
        thetas = 2*pi*runif(n2)
        angles = c(0, thetas, -rev(thetas))
        recf   = amplitudes * exp(complex(imaginary = angles))
      }
      temp = Re( fft(recf, inverse = T) / n )
      
      # adjust variance of the surrogate time series to match the original
      return(temp / sd(temp) * sigma)
    }))
  }
  else {
    if( any(!is.finite(ts)) ) {
      stop("SurrogateData(): input time series contained invalid values")
    }
    
    n = length(ts)
    I_season = suppressWarnings( matrix( 1:T_period, nrow = n, ncol = 1 ) )
    
    # Calculate seasonal cycle using smooth.spline
    seasonal_F = smooth.spline(
      c(I_season - T_period, I_season, I_season + T_period), c(ts, ts, ts) )
    
    seasonal_cyc   = predict( seasonal_F, I_season ) $ y
    seasonal_resid = ts - seasonal_cyc
    
    return(sapply(1:num_surr, function(i) {
      seasonal_cyc + sample(seasonal_resid, n) + rnorm(n, 0, alpha)
    }))
  }
}
