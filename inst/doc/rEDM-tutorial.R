## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, 
                      tidy = TRUE, cache = TRUE, 
                      fig.width = 5, fig.height = 3.5)
knitr::opts_knit$set(global.par = TRUE, progress = TRUE)
options(digits = 2)
par(mar = c(4, 4, 1, 1), oma = c(0, 0, 0, 0), mgp = c(2.5, 1, 0))

## ----CRAN installation instructions, eval = FALSE------------------------
#  install.packages("rEDM")

## ----GitHub installation instructions, eval = FALSE----------------------
#  devtools::install_github("ha0ye/rEDM")

## ----fig_time_series_projection, echo = FALSE, fig.cap = "Time Series Projection from the Lorenz Attractor"----
knitr::include_graphics(here::here("vignettes", "vignette_figs", "figure_1.svg"))

## ----fig_attractor_reconstruction, echo = FALSE, fig.cap = "Attractor Reconstruction from Lagged Coordinates"----
knitr::include_graphics(here::here("vignettes", "vignette_figs", "figure_2.svg"))

## ----load package--------------------------------------------------------
library(rEDM)

## ----load tentmap data---------------------------------------------------
data(tentmap_del)
str(tentmap_del)

## ----lib and pred for tentmap--------------------------------------------
ts <- tentmap_del
lib <- c(1, 100)
pred <- c(201, 500)

## ----simplex on tentmap--------------------------------------------------
simplex_output <- simplex(ts, lib, pred)
str(simplex_output)

## ----rho vs. E for tentmap-----------------------------------------------
par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0)) # set up margins for plotting
plot(simplex_output$E, simplex_output$rho, type = "l",  
     xlab = "Embedding Dimension (E)", ylab = "Forecast Skill (rho)")

## ----simplex varying tp for tentmap--------------------------------------
simplex_output <- simplex(ts, lib, pred, E = 2, tp = 1:10)

## ----rho vs. tp for tentmap----------------------------------------------
plot(simplex_output$tp, simplex_output$rho, type = "l", 
     xlab = "Time to Prediction (tp)", ylab = "Forecast Skill (rho)")

## ----smap for tentmap----------------------------------------------------
smap_output <- s_map(ts, lib, pred, E = 2)

## ----rho vs. theta for tentmap-------------------------------------------
plot(smap_output$theta, smap_output$rho, type = "l", 
     xlab = "Nonlinearity (theta)", ylab = "Forecast Skill (rho)")

## ----rho vs. theta with noise--------------------------------------------
ts_err <- ts + rnorm(length(ts), sd = sd(ts) * 0.2)
smap_output_err <- s_map(ts_err, lib, pred, E = 2)
plot(smap_output_err$theta, smap_output_err$rho, type = "l", 
     xlab = "Nonlinearity (theta)", ylab = "Forecast Skill (rho)")

## ----load block_3sp data-------------------------------------------------
data(block_3sp)
str(block_3sp)

## ----block_lnlp for block_3sp, warning = FALSE---------------------------
lib <- c(1, NROW(block_3sp))
pred <- c(1, NROW(block_3sp))

cols <- c(1, 2, 4)
target <- 1

block_lnlp_output <- block_lnlp(block_3sp, lib = lib, pred = pred, 
                                columns = cols, target_column = target, 
                                stats_only = FALSE, first_column_time = TRUE)

## ------------------------------------------------------------------------
block_lnlp_output_2 <- block_lnlp(block_3sp, lib = lib, pred = pred, 
                                columns = c("x_t", "x_t-1", "y_t"), target_column = "x_t", 
                                stats_only = FALSE, first_column_time = TRUE)

# test for equality
stopifnot(identical(block_lnlp_output, block_lnlp_output_2))

## ------------------------------------------------------------------------
str(block_lnlp_output)

## ------------------------------------------------------------------------
list_of_model_predictions <- block_lnlp_output$model_output
first_data_frame_of_predictions <- list_of_model_predictions[[1]]

observed <- first_data_frame_of_predictions$obs
predicted <- first_data_frame_of_predictions$pred

## ---- echo = FALSE-------------------------------------------------------
par(pty = "s")

## ----observed vs predicted for block_lnlp, fig.width = 4, fig.height = 4----
plot_range <- range(c(observed, predicted), na.rm = TRUE)
plot(observed, predicted, xlim = plot_range, ylim = plot_range, 
     xlab = "Observed", ylab = "Predicted", asp = 1)
abline(a = 0, b = 1, lty = 2, col = "blue")

## ---- echo = FALSE-------------------------------------------------------
par(pty = "m")

## ----fig_cross_mapping, echo = FALSE, fig.cap = "Cross Mapping Between Reconstructions of the Lorenz Attractor"----
knitr::include_graphics(here::here("vignettes", "vignette_figs", "figure_3.svg"))

## ----sardine anchovy ccm, warning = FALSE--------------------------------
data(sardine_anchovy_sst)
anchovy_xmap_sst <- ccm(sardine_anchovy_sst, E = 3, 
                        lib_column = "anchovy", target_column = "np_sst", 
                        lib_sizes = seq(10, 80, by = 10), num_samples = 100, 
                        random_libs = TRUE, replace = TRUE)
sst_xmap_anchovy <- ccm(sardine_anchovy_sst, E = 3, 
                        lib_column = "np_sst", target_column = "anchovy", 
                        lib_sizes = seq(10, 80, by = 10), num_samples = 100, 
                        random_libs = TRUE, replace = TRUE)
str(anchovy_xmap_sst)

## ----sardine anchovy ccm plot--------------------------------------------
a_xmap_t_means <- ccm_means(anchovy_xmap_sst)
t_xmap_a_means <- ccm_means(sst_xmap_anchovy)

plot(a_xmap_t_means$lib_size, pmax(0, a_xmap_t_means$rho), type = "l", col = "red", 
     xlab = "Library Size", ylab = "Cross Map Skill (rho)", ylim = c(0, 0.25))
lines(t_xmap_a_means$lib_size, pmax(0, t_xmap_a_means$rho), col = "blue")
legend(x = "topleft", legend = c("anchovy xmap SST", "SST xmap anchovy"), 
       col = c("red", "blue"), lwd = 1, bty = "n", inset = 0.02, cex = 0.8)

## ------------------------------------------------------------------------
data(paramecium_didinium)
str(paramecium_didinium)

## ------------------------------------------------------------------------
vars <- names(paramecium_didinium)[2:3] # c("paramecium", "didinium")

# generate all combinations of lib_column, target_column, tp
params <- expand.grid(lib_column = vars, 
                      target_column = vars, 
                      tp = -10:10)

# throw out cases where lib == target
params <- params[params$lib_column != params$target_column, ]

## ------------------------------------------------------------------------
params$E <- 3

## ---- warning = FALSE----------------------------------------------------
output <- do.call(rbind, lapply(seq_len(NROW(params)), function(i) {
    ccm(paramecium_didinium, E = params$E[i], 
        lib_sizes = NROW(paramecium_didinium), random_libs = FALSE, 
        lib_column = params$lib_column[i], 
        target_column = params$target_column[i], 
        tp = params$tp[i])
}))

## ------------------------------------------------------------------------
output$direction <- paste(output$lib_column, "xmap to\n", output$target_column)

library(ggplot2)
time_delay_ccm_fig <- ggplot(output, aes(x = tp, y = rho, color = direction)) + 
    geom_line() + theme_bw()

## ------------------------------------------------------------------------
print(time_delay_ccm_fig)

## ----load e120 data------------------------------------------------------
data(e120_biodiversity)
str(e120_biodiversity)

## ----setup e120 data-----------------------------------------------------
normalize <- function(x, ...) {(x - mean(x, ...))/sd(x, ...)}

# separate time column from data
vars <- c("AbvBioAnnProd", "noh020tot", "invrichness", "SummerPrecip.mm.")
composite_ts <- e120_biodiversity[, vars]

# normalize each time series within a plot
data_by_plot <- split(composite_ts, e120_biodiversity$Plot)
normalized_data <- lapply(data_by_plot, function(df) sapply(df, normalize))
composite_ts <- cbind(Year = e120_biodiversity$Year, 
                      data.frame(do.call(rbind, normalized_data)))

## ----make composite library----------------------------------------------
segments_end <- cumsum(sapply(data_by_plot, NROW))
segments_begin <- c(1, segments_end[-length(segments_end)] + 1)
segments <- cbind(segments_begin, segments_end)

# Choose random segments for prediction
set.seed(2312)
rndlib <- sample(1:NROW(segments), floor(NROW(segments) * 0.75))
composite_lib <- segments[rndlib, ]
composite_pred <- segments[-rndlib, ]

## ------------------------------------------------------------------------
precip_ts <- unique(e120_biodiversity[, c("Year", "SummerPrecip.mm.")])
precip_ts <- precip_ts[order(precip_ts$Year), ]

## ----simplex on e120, warning = FALSE, fig.width = 6---------------------
simplex_out <- lapply(c("AbvBioAnnProd", "noh020tot", "invrichness"), 
                      function(var) {
                          simplex(composite_ts[, c("Year", var)], E = 2:4, 
                                  lib = composite_lib, pred = composite_pred)
                      })
simplex_out[[4]] <- simplex(precip_ts, E = 2:5)
names(simplex_out) <- names(composite_ts)[-1]

par(mfrow = c(2, 2))
for(var in names(simplex_out))
{
    plot(simplex_out[[var]]$E, simplex_out[[var]]$rho, type = "l", 
         xlab = "Embedding Dimension (E)", ylab = "Forecast Skill (rho)", 
         main = var)
}

## ----best E for e120-----------------------------------------------------
best_E <- sapply(simplex_out, function(df) {df$E[which.max(df$rho)]})
best_E

## ----smap on e120, warning = FALSE, fig.width = 6, results = "hide"------
smap_out <- lapply(c("AbvBioAnnProd", "noh020tot", "invrichness"), 
                   function(var) {
                       s_map(composite_ts[, c("Year", var)], E = best_E[var], 
                             lib = composite_lib, pred = composite_pred)
                   })
smap_out[[4]] <- s_map(precip_ts, E = best_E[4])
names(smap_out) <- names(simplex_out)

par(mfrow = c(2, 2))
for(var in names(simplex_out))
{    plot(smap_out[[var]]$theta, smap_out[[var]]$rho, type = "l", 
          xlab = "Nonlinearity (theta)", ylab = "Forecast Skill (rho)", 
          main = var)
}

## ----make block for e120-------------------------------------------------
block_data <- make_block(composite_ts[, 2:5], t = composite_ts$Year, 
                         max_lag = 4, lib = segments)
str(block_data)

## ----block_lnlp for e120, warning = FALSE--------------------------------
AB_columns <- c("AbvBioAnnProd", "AbvBioAnnProd_1", "AbvBioAnnProd_2")
AB_output <- block_lnlp(block_data, lib = composite_lib, pred = composite_pred, 
                        columns = AB_columns, target_column = "AbvBioAnnProd", 
                        stats_only = FALSE, first_column_time = TRUE)

ABNO_columns <- c(AB_columns, "noh020tot", "noh020tot_1", "noh020tot_2")
ABNO_output <- block_lnlp(block_data, lib = composite_lib, pred = composite_pred, 
                          columns = ABNO_columns, target_column = "AbvBioAnnProd", 
                          stats_only = FALSE, first_column_time = TRUE)

## ---- echo = FALSE-------------------------------------------------------
par(mfrow = c(1, 1), pty = "s")

## ----block_lnlp on e120, warning = FALSE, fig.width = 4, fig.height = 4----
observed_AB <- AB_output$model_output[[1]]$obs
predicted_AB <- AB_output$model_output[[1]]$pred

observed_ABNO <- ABNO_output$model_output[[1]]$obs
predicted_ABNO <- ABNO_output$model_output[[1]]$pred

plot_range <- range(c(observed_AB, predicted_AB), na.rm = TRUE)
plot(observed_AB, predicted_AB, xlim = plot_range, ylim = plot_range, 
     xlab = "Observed", ylab = "Predicted")
abline(a = 0, b = 1, lty = 2, col = "darkgrey", lwd = 2)
abline(lm(predicted_AB ~ observed_AB), col = "black", lty = 3, lwd = 2)

points(observed_ABNO, predicted_ABNO, pch = 2, col = "red")
abline(lm(predicted_ABNO ~ observed_ABNO), col = "red", lty = 3, lwd = 2)

legend("topleft", legend = c(paste("(biomass alone) rho =", round(AB_output$rho, 2)), 
                             paste("(biomass and nitrate) rho =", round(ABNO_output$rho, 2))), 
       lty = 3, lwd = 2, col = c("black", "red"), bty = "n")

## ---- echo = FALSE-------------------------------------------------------
par(mfrow = c(1, 1), pty = "m")

## ----ccm on e120, warning = FALSE----------------------------------------
# A. repens:
lib_sizes <- c(seq(5, 55, by = 5), seq(55, 400, by = 50))

no_xmap_inv <- ccm(composite_ts, lib = segments, pred = segments, 
                   lib_column = "noh020tot", target_column = "invrichness", 
                   E = best_E["noh020tot"], lib_sizes = lib_sizes)
inv_xmap_no <- ccm(composite_ts, lib = segments, pred = segments, 
                   lib_column = "invrichness", target_column = "noh020tot", 
                   E = best_E["invrichness"], lib_sizes = lib_sizes)

n_xmap_i_means <- ccm_means(no_xmap_inv)
i_xmap_n_means <- ccm_means(inv_xmap_no)

plot(n_xmap_i_means$lib_size, pmax(0, n_xmap_i_means$rho), type = "l", 
     xlab = "Library Size", ylab = "Cross Map Skill (rho)", 
     col = "red", ylim = c(0, 0.4), lwd = 2)
lines(i_xmap_n_means$lib_size, pmax(0, i_xmap_n_means$rho), 
      col = "blue", lwd = 2)
legend(x = "topleft", col = c("red", "blue"), lwd = 2, 
       legend = c("Nitrate xmap Inv. Richness", "Inv. Richness xmap Nitrate"), 
       inset = 0.02, bty = "n", cex = 0.8)

## ----ccm on e120 with biological productivity, warning = FALSE-----------
ab_xmap_inv <- ccm(composite_ts, lib = segments, pred = segments,
                   lib_column = "AbvBioAnnProd", target_column = "invrichness", 
                   E = best_E["AbvBioAnnProd"], lib_sizes = lib_sizes, silent = TRUE)
inv_xmap_ab <- ccm(composite_ts, lib = segments, pred = segments, 
                   lib_column = "invrichness", target_column = "AbvBioAnnProd", 
                   E = best_E["invrichness"], lib_sizes = lib_sizes, silent = TRUE)

a_xmap_i_means <- ccm_means(ab_xmap_inv)
i_xmap_a_means <- ccm_means(inv_xmap_ab)

plot(a_xmap_i_means$lib_size, pmax(0, a_xmap_i_means$rho), type = "l", 
     xlab = "Library Size", ylab = "Cross Map Skill (rho)", 
     col = "orange", ylim = c(0, 0.4), lwd = 2)
lines(i_xmap_a_means$lib_size, pmax(0, i_xmap_a_means$rho), 
      col = "blue", lwd = 2)
legend(x = "topleft", col = c("orange", "blue"), 
       legend = c("Abv. Biomass xmap Inv. Richness", 
                  "Inv. Richness xmap Abv. Biomass"), 
       lwd = 2, inset = 0.02, bty = "n", cex = 0.8)

## ----thrips data---------------------------------------------------------
data(thrips_block)
str(thrips_block)

## ----thrips plot, echo = FALSE, fig.width = 6, fig.height = 7------------
par(mfrow = c(4, 1), mar = c(2, 4, 1, 1), oma = c(2, 0, 0, 0), 
    mgp = c(2.5, 1, 0))

time_dec <- thrips_block$Year + (thrips_block$Month) / 12
plot(time_dec, thrips_block$Thrips_imaginis, type = "l", col = "green", ylab = "Thrips", xlab = "")
plot(time_dec, thrips_block$maxT_degC, type = "l", col = "red", ylab = "maxT (oC)", xlab = "")
plot(time_dec, thrips_block$Rain_mm, type = "l", col = "blue", ylab = "Rain (mm)", xlab = "")
plot(time_dec, thrips_block$Season, type = "l", col = "magenta", ylab = "Season", xlab = "")
mtext("Year", side = 1, outer = TRUE, line = 1)

par(mfrow = c(1, 1))

## ----univariate thrips, warning = FALSE----------------------------------
ts <- thrips_block$Thrips_imaginis
lib <- c(1, length(ts))
pred <- c(1, length(ts))
simplex_output <- simplex(ts, lib, pred)

## ----rho vs. e for thrips, echo = FALSE----------------------------------
plot(simplex_output$E, simplex_output$rho, type = "l", 
     xlab = "Embedding Dimension (E)", ylab = "Forecast Skill (rho)")

## ----smap for thrips, warning = FALSE------------------------------------
smap_output <- list(s_map(ts, lib, pred, E = 4), 
                    s_map(ts, lib, pred, E = 8))

## ---- echo = FALSE-------------------------------------------------------
par(mfrow = c(1, 2))

## ----rho vs. theta for thrips, echo = FALSE------------------------------
plot(smap_output[[1]]$theta, smap_output[[1]]$rho, type = "l", xlim = c(0, 4), 
     xlab = "Nonlinearity (theta)", ylab = "Forecast Skill (rho)")
plot(smap_output[[2]]$theta, smap_output[[2]]$rho, type = "l", xlim = c(0, 4), 
     xlab = "Nonlinearity (theta)", ylab = "Forecast Skill (rho)")

## ---- echo = FALSE-------------------------------------------------------
par(mfrow = c(1, 1))

## ----compute ccm matrix for thrips, results = "hold", warning = FALSE----
vars <- colnames(thrips_block[3:6])
n <- NROW(thrips_block)
ccm_matrix <- matrix(NA, nrow = length(vars), ncol = length(vars), 
                     dimnames = list(vars, vars))

for(ccm_from in vars)
{
    for(ccm_to in vars[vars != ccm_from]) # ignore self-interactions
    {
        out_temp <- ccm(thrips_block, E = 8, 
                        lib_column = ccm_from, target_column = ccm_to,
                        lib_sizes = n, replace = FALSE)
        ccm_matrix[ccm_from, ccm_to] <- out_temp$rho
    }
}

## ----compute corr matrix for thrips--------------------------------------
corr_matrix <- array(NA, dim = c(length(vars), length(vars)), 
                     dimnames = list(vars, vars))

for(ccm_from in vars)
{
    for(ccm_to in vars[vars != ccm_from])
    {
        cf_temp <- ccf(thrips_block[, ccm_from], thrips_block[, ccm_to], 
                       type = "correlation", lag.max = 6, plot = FALSE)$acf
        corr_matrix[ccm_from, ccm_to] <- max(abs(cf_temp))
    }
}

## ----xmap vs. corr matrix for thrips-------------------------------------
head(ccm_matrix)
head(corr_matrix)

## ----ccm on thrips, results = "hide", warning = FALSE--------------------
thrips_xmap_maxT <- ccm(thrips_block, E = 8, random_libs = TRUE, 
                        lib_column = "Thrips_imaginis", target_column = "maxT_degC", 
                        lib_sizes = seq(10, 75, by = 5), num_samples = 300)
maxT_xmap_thrips <- ccm(thrips_block, E = 8, random_libs = TRUE, 
                        lib_column = "maxT_degC", target_column = "Thrips_imaginis", 
                        lib_sizes = seq(10, 75, by = 5), num_samples = 300)

ccm_out <- list(ccm_means(thrips_xmap_maxT), ccm_means(maxT_xmap_thrips))

## ----ccm plot, echo = FALSE----------------------------------------------
plot(ccm_out[[1]]$lib_size, pmax(0, ccm_out[[1]]$rho), type = "l", col = "red",  
     xlab = "Library Size", ylab = "Cross Map Skill (rho)", ylim = c(0, 1))
lines(ccm_out[[2]]$lib_size, pmax(0, ccm_out[[2]]$rho), col = "blue")
abline(h = corr_matrix['Thrips_imaginis', 'maxT_degC'], col = "black", lty = 2)
legend(x = "bottomright", legend = c("Thrips xmap maxT", "maxT xmap Thrips"), 
       col = c("red", "blue"), lwd = 1, inset = 0.02)

## ----ccm on thrips and rainfall, results = "hide", warning = FALSE-------
thrips_xmap_Rain <- ccm(thrips_block, E = 8, random_libs = TRUE, 
                        lib_column = "Thrips_imaginis", target_column = "Rain_mm", 
                        lib_sizes = seq(10, 75, by = 5), num_samples = 300)
Rain_xmap_thrips <- ccm(thrips_block, E = 8, random_libs = TRUE, 
                        lib_column = "Rain_mm", target_column = "Thrips_imaginis", 
                        lib_sizes = seq(10, 75, by = 5), num_samples = 300)

ccm_out <- list(ccm_means(thrips_xmap_Rain), ccm_means(Rain_xmap_thrips))

## ----rainfall and thrips ccm plot, echo = FALSE--------------------------
plot(ccm_out[[1]]$lib_size, pmax(0, ccm_out[[1]]$rho), type = "l", col = "red", 
     xlab = "Library Size", ylab = "Cross Map Skill (rho)", ylim = c(0, 1))
lines(ccm_out[[2]]$lib_size, pmax(0, ccm_out[[2]]$rho), col = "blue")
abline(h = corr_matrix['Thrips_imaginis', 'Rain_mm'], col = 'black', lty = 2)
legend(x = "topleft", legend = c("Thrips xmap Rain", "Rain xmap Thrips"),
       col = c("red", "blue"), lwd = 1, inset = 0.02)

## ----ccm on thrips and season, results = "hide", warning = FALSE---------
thrips_xmap_Season <- ccm(thrips_block, E = 8, random_libs = TRUE, 
                          lib_column = "Thrips_imaginis", target_column = "Season", 
                          lib_sizes = seq(10, 75, by = 5), num_samples = 300)
Season_xmap_thrips <- ccm(thrips_block, E = 8, random_libs = TRUE, 
                          lib_column = "Season", target_column = "Thrips_imaginis", 
                          lib_sizes = seq(10, 75, by = 5), num_samples = 300)

ccm_out <- list(ccm_means(thrips_xmap_Season), ccm_means(Season_xmap_thrips))

## ----season and thrips ccm plot, echo = FALSE----------------------------
par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0)) # set up margins for plotting
plot(ccm_out[[1]]$lib_size, pmax(0, ccm_out[[1]]$rho), type = "l", col = "red", 
     xlab = "Library Size", ylab = "Cross Map Skill (rho)", ylim = c(0, 1))
lines(ccm_out[[2]]$lib_size, pmax(0, ccm_out[[2]]$rho), col = "blue")
abline(h = corr_matrix['Thrips_imaginis', 'Season'], col = 'black', lty = 2)
legend(x = "bottomright", legend = c("Thrips xmap Season", "Season xmap Thrips"), 
       col = c("red", "blue"), lwd = 1, inset = 0.02)

## ----seasonal surrogates for thrips, warning = FALSE---------------------
num_surr <- 1000
surr_maxT <- make_surrogate_data(thrips_block$maxT_degC, method = "seasonal", 
                                 T_period = 12, num_surr = num_surr)
surr_Rain <- make_surrogate_data(thrips_block$Rain_mm, method = "seasonal", 
                                 T_period = 12, num_surr = num_surr)

rho_surr <- data.frame(maxT = numeric(num_surr), Rain = numeric(num_surr))

for (i in 1:num_surr) {
    rho_surr$maxT[i] <- ccm(cbind(thrips_block$Thrips_imaginis, surr_maxT[,i]), 
                            E = 8, lib_column = 1, target_column = 2, 
                            lib_sizes = NROW(thrips_block), replace = FALSE)$rho

    rho_surr$Rain[i] <- ccm(cbind(thrips_block$Thrips_imaginis, surr_Rain[,i]), 
                            E = 8, lib_column = 1, target_column = 2, 
                            lib_sizes = NROW(thrips_block), replace = FALSE)$rho
    }

## ----significance of randomization test, tidy = FALSE--------------------
(sum(ccm_matrix['Thrips_imaginis', 'Rain_mm'] < rho_surr$Rain) + 1) / 
    (length(rho_surr$Rain) + 1)
(sum(ccm_matrix['Thrips_imaginis', 'maxT_degC'] < rho_surr$maxT) + 1) / 
    (length(rho_surr$maxT) + 1)

