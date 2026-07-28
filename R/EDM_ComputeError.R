#------------------------------------------------------------------------
# ComputeError : Pearson rho, MAE, CAE, RMSE.
#
# Faithful port of pyEDM AuxFunc.ComputeError, including its quirks:
#   - MAE is round(max(obs - pred)) : the SIGNED maximum error (pyEDM labels
#     it MAE but does not take the mean or absolute value). Preserved for
#     bit-compatibility with pyEDM outputs.
#   - CAE is the summed absolute error.
#   - Fewer than 5 finite pairs -> all NaN.
#   - NaN removed from pred first, then obs (Pearson needs finite pairs).
#
# Naming: functions/constants uppercase-first, variables lowercase-first.
#------------------------------------------------------------------------

#------------------------------------------------------------------------
#' Prediction skill statistics.
#'
#' @param obs numeric vector of observations.
#' @param pred numeric vector of predictions (aligned with obs).
#' @param digits rounding applied to each statistic (pyEDM default 6).
#' @return named list(rho, MAE, CAE, RMSE).
#' @keywords internal
#------------------------------------------------------------------------
ComputeError <- function(obs, pred, digits = 6) {
  notNan <- is.finite(pred)
  if (any(!notNan)) { pred <- pred[notNan]; obs <- obs[notNan] }
  notNan <- is.finite(obs)
  if (any(!notNan)) { pred <- pred[notNan]; obs <- obs[notNan] }

  if (length(pred) < 5)
    return(list(rho = NA_real_, MAE = NA_real_, CAE = NA_real_, RMSE = NA_real_))

  err  <- obs - pred
  rho  <- round(suppressWarnings(stats::cor(obs, pred)), digits)   # Pearson
  MAE  <- round(max(err), digits)                                  # signed max (pyEDM)
  CAE  <- round(sum(abs(err)), digits)
  RMSE <- round(sqrt(mean(err^2)), digits)

  list(rho = rho, MAE = MAE, CAE = CAE, RMSE = RMSE)
}
