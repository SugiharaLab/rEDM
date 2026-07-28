#------------------------------------------------------------------------
# Plotting helpers for showPlot = TRUE (base graphics).
#------------------------------------------------------------------------

#------------------------------------------------------------------------
#' Observed vs predicted time series.
#' @param df data.frame with Time, Observations, Predictions.
#' @param dataName optional label used in the title.
#' @param E,Tp embedding dimension and forecast interval, shown in the title.
#' @return Invisibly \code{NULL}.
#' @keywords internal
#------------------------------------------------------------------------
PlotObsPred <- function(df, dataName = "", E = NULL, Tp = NULL) {
  err   <- ComputeError(df$Observations, df$Predictions)
  title <- sprintf("%s E=%d Tp=%d  rho=%.2f  RMSE=%.2f",
                   dataName, E, Tp,
                   ifelse(is.na(err$rho),  NA_real_, err$rho ),
                   ifelse(is.na(err$RMSE), NA_real_, err$RMSE))
  time <- PlotTimeAxis(df[[1]], nrow(df))
  yl <- range(c(df$Observations, df$Predictions), na.rm = TRUE)

  plot(time, df$Observations, type = "l", lwd = 3, ylim = yl,
       col = "blue", xlab = names(df)[1], ylab = "", main = title,
       cex.axis = 1.3, cex.lab = 1.3)

  graphics::lines(time, df$Predictions, col = "red", lwd = 3)
  graphics::legend("topright", legend = c("Observations", "Predictions"),
                   fill = c("blue", "red"), bty = "n", cex = 1.2)
  invisible(NULL)
}

#------------------------------------------------------------------------
#' S-map observed/predicted plus one panel per coefficient.
#' @param SmapList list returned by \code{SMap} (predictions, coefficients).
#' @param dataName optional label used in the title.
#' @param E,Tp embedding dimension and forecast interval, shown in the title.
#' @return Invisibly \code{NULL}.
#' @keywords internal
#------------------------------------------------------------------------
PlotSmap <- function(SmapList, dataName = "", E = NULL, Tp = NULL) {
  if (!all(c("predictions", "coefficients") %in% names(SmapList))) {
    message("PlotSmap(): need predictions and coefficients.")
    return(invisible(NULL))
  }
  p        <- SmapList$predictions
  co       <- SmapList$coefficients
  time     <- PlotTimeAxis(p[[1]], nrow(p))
  numCoeff <- ncol(co) - 1L
  oldPar   <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldPar), add = TRUE)
  
  graphics::par(mfrow = c(numCoeff + 1L, 1L), mar = c(3.5, 4, 0.5, 1),
                mgp = c( 1.5, 0.5, 0 ), cex.axis = 1.3, cex.lab = 1.3)
  yl <- range(c(p$Observations, p$Predictions), na.rm = TRUE)
  
  # Observations & Predictions
  title = sprintf("S-Map %s E=%d Tp=%d", dataName, E, Tp)
  
  plot(time, p$Observations, type = "l", ylim = yl, xlab = names(p)[1],
       ylab = "", col = "blue", lwd = 3, cex.axis = 1.3, cex.lab = 1.3,
       main = title )
  
  graphics::lines(time, p$Predictions, col = "red", lwd = 3)

  graphics::legend("topright", legend = c("Observations", "Predictions"),
                   fill = c("blue", "red"), bty = "n", cex = 1.5)

  # Coefficients
  coefName <- names(co)
  title = paste( dataName, 'S-Map Coefficients', '\nE=' , E, ' Tp=', Tp )
  for (j in 2:ncol(co)) {
    plot(time, co[,j], type = "l", lwd = 3, col = "blue", xlab = names(p)[1],
         ylab = gsub("\u2202", "d", coefName[j]))  # ASCII-safe label (any locale)
    graphics::mtext( title, side = 3, line = -1.5, cex = 1.2 )
  }

  invisible(NULL)
}

#------------------------------------------------------------------------
#' rho versus a swept parameter (EmbedDimension / PredictInterval / ...).
#' @keywords internal
#------------------------------------------------------------------------
PlotSweep <- function(df, xlab, title = NULL) {
  plot(df[[1]], df$rho, type = "l", lwd = 3, main = title,
       xlab = xlab, ylab = expression("Prediction Skill (" * rho * ")"))

  invisible(NULL)
}

#------------------------------------------------------------------------
#' CCM cross-map skill versus library size, both directions.
#' @param df CCM result data.frame (LibSize plus two rho columns).
#' @param E embedding dimension, shown in the title.
#' @return Invisibly \code{NULL}.
#' @keywords internal
#------------------------------------------------------------------------
PlotCCM <- function(df, E = NULL) {
  yl <- range(c(df[[2]], df[[3]]), na.rm = TRUE)

  title = if (!is.null(E)) paste0(names(df)[2], " : ", names(df)[3], "  E=", E) else ""
  plot(df$LibSize, df[[2]], type = "l", col = "blue", ylim = yl, lwd = 3, 
       main = title,
       xlab = "Library Size",
       ylab = expression("Cross Map Skill (" * rho * ")"))
  
  graphics::lines(df$LibSize, df[[3]], col = "red", lwd = 3)
  graphics::abline(h = 0)
  graphics::legend("topright", legend = names(df)[2:3],
                   fill = c("blue", "red"), bty = "n", cex = 1.2)
  invisible(NULL)
}

#------------------------------------------------------------------------
#' Coerce a first column into a plottable time axis (numeric, Date, or index).
#' @keywords internal
#' @noRd
#------------------------------------------------------------------------
PlotTimeAxis <- function(col, n) {
  if (is.numeric(col)) return(col)
  t <- try(as.Date(col), silent = TRUE)
  if (inherits(t, "try-error") || all(is.na(t)))
    t <- try(as.POSIXct(col), silent = TRUE)
  if (inherits(t, "try-error") || all(is.na(t))) t <- seq_len(n)
  t
}

