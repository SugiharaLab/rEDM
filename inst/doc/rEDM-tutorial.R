## ----CRAN installation instructions, eval = FALSE------------------------
#  install.packages("rEDM")

## ----GitHub installation instructions, eval = FALSE----------------------
#  devtools::install_github("ha0ye/rEDM")

## ----load tentmap data---------------------------------------------------
library(rEDM)
data(tentmap_del)
str(tentmap_del)

## ----lib and pred for tentmap--------------------------------------------
ts <- tentmap_del
lib <- c(1, 100)
pred <- c(201, 500)

## ----simplex on tentmap--------------------------------------------------
simplex_output <- simplex(ts, lib, pred)
str(simplex_output)

## ----rho vs. E for tentmap, tidy = TRUE, fig.width = 5, fig.height = 3.5----
par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0)) # set up margins for plotting
plot(simplex_output$E, simplex_output$rho, type = "l", xlab = "Embedding Dimension (E)", ylab = "Forecast Skill (rho)")

## ----simplex varying tp for tentmap--------------------------------------
simplex_output <- simplex(ts, lib, pred, E = 2, tp = 1:10)

## ----rho vs. tp for tentmap, tidy = TRUE, fig.width = 5, fig.height = 3.5----
par(mar = c(4,4,1,1))
plot(simplex_output$tp, simplex_output$rho, type = "l", 
     xlab = "Time to Prediction (tp)", ylab = "Forecast Skill (rho)")

## ----smap for tentmap----------------------------------------------------
smap_output <- s_map(ts, lib, pred, E = 2)

## ----rho vs. theta for tentmap, tidy = TRUE, fig.width = 5, fig.height = 3.5----
par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0))
plot(smap_output$theta, smap_output$rho, type = "l", 
     xlab = "Nonlinearity (theta)", ylab = "Forecast Skill (rho)")

## ----rho vs. theta with noise, tidy = TRUE, fig.width = 5, fig.height = 3.5----
ts <- ts + rnorm(length(ts), sd = sd(ts) * 0.2)
smap_output <- s_map(ts, lib, pred, E = 2)
par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0))
plot(smap_output$theta, smap_output$rho, type = "l", 
     xlab = "Nonlinearity (theta)", ylab = "Forecast Skill (rho)")

## ----load block_3sp data-------------------------------------------------
data(block_3sp)
str(block_3sp)

## ----block_lnlp for block_3sp, tidy = TRUE, warning = FALSE--------------
lib <- c(1, NROW(block_3sp))
pred <- c(1, NROW(block_3sp))

cols <- c(1, 2, 4) # c("x_t", "x_t-1", "y_t")
target <- 1 # "x_t"

block_lnlp_output <- block_lnlp(block_3sp, lib = lib, pred = pred, 
                                columns = cols, target_column = target, 
                                stats_only = FALSE, first_column_time = TRUE)

## ------------------------------------------------------------------------
str(block_lnlp_output)

## ----observed vs predicted for block_lnlp, tidy = TRUE, fig.width = 4, fig.height = 4----
observed <- block_lnlp_output[[1]]$model_output$obs
predicted <- block_lnlp_output[[1]]$model_output$pred

par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0), pty = "s")
plot_range <- range(c(observed, predicted), na.rm = TRUE)
plot(observed, predicted, xlim = plot_range, ylim = plot_range, 
     xlab = "Observed", ylab = "Predicted")
abline(a = 0, b = 1, lty = 2, col = "blue")

## ----sardine anchovy ccm plot, tidy = TRUE, fig.width = 5, fig.height = 3.5----
a_xmap_t_means <- ccm_means(anchovy_xmap_sst)
t_xmap_a_means <- ccm_means(sst_xmap_anchovy)

par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0)) # set up margins for plotting
y1 <- pmax(0, a_xmap_t_means$rho)
y2 <- pmax(0, t_xmap_a_means$rho)

plot(a_xmap_t_means$lib_size, y1, type = "l", col = "red", 
     xlab = "Library Size", ylab = "Cross Map Skill (rho)", ylim = c(0, 0.25))
lines(t_xmap_a_means$lib_size, y2, col = "blue")
legend(x = "topleft", legend = c("anchovy xmap SST", "SST xmap anchovy"), 
       col = c("red", "blue"), lwd = 1, bty = "n", inset = 0.02, cex = 0.8)

## ----load e120 data, tidy = TRUE-----------------------------------------
data(e120_biodiversity)
head(e120_biodiversity)

# separate time column from data
composite_ts <- e120_biodiversity[,c(7:9,12)]

# normalize each time series
n <- NCOL(composite_ts)
blocks <- e120_biodiversity$Plot
blocks_index <- sort(unique(blocks))
for (j in 1:n) {
    for (i in 1:length(blocks_index)) {
        subs <- which(blocks == blocks_index[i])
        composite_ts[subs,j] <- (composite_ts[subs,j] - mean(composite_ts[subs,j])) / sd(composite_ts[subs,j])
        }
    }

composite_ts <- cbind(year = e120_biodiversity$Year, composite_ts)

## ----make composite ts for e120------------------------------------------
# make composite library
segments <- NULL
startpos <- 1
for (i in 2:nrow(composite_ts)) {
  if (composite_ts$year[i] < composite_ts$year[i - 1]) {
    segments <- rbind(segments, c(startpos, i))
    startpos <- i + 1
  }
}
segments <- rbind(segments, c(max(segments) + 1, nrow(composite_ts)))

# choose random segments for prediction
set.seed(2312)
rndlib <- sort(sample(1:nrow(segments), round(nrow(segments)/2, 0), replace = FALSE))
composite_lib <- segments[rndlib, ]
composite_pred <- segments[-rndlib, ]

## ----precip for e120-----------------------------------------------------
precip_ts <- unique(e120_biodiversity[, c("Year", "SummerPrecip.mm.")])
precip_ts <- precip_ts[order(precip_ts$Year), ]

## ----best E for e120-----------------------------------------------------
bestE <- sapply(simplex_output_list, function(simplex_output) {
    simplex_output$E[which.max(simplex_output$rho)]
    })
bestE

## ----make block for e120, tidy = TRUE------------------------------------
n <- NROW(composite_ts)

# make lags
block_data <- data.frame(year = composite_ts$year)
block_data$AB_tm <- composite_ts$AbvBioAnnProd
block_data$AB_tm1 <- c(NA, block_data$AB_tm[-n])
block_data$AB_tm2 <- c(NA, block_data$AB_tm1[-n])
block_data$AB_tm3 <- c(NA, block_data$AB_tm2[-n])

block_data$NO_tm <- composite_ts$noh020tot
block_data$NO_tm1 <- c(NA, block_data$NO_tm[-n])
block_data$NO_tm2 <- c(NA, block_data$NO_tm1[-n])
block_data$NO_tm3 <- c(NA, block_data$NO_tm2[-n])

block_data$IV_tm <- composite_ts$invrichness
block_data$IV_tm1 <- c(NA, block_data$IV_tm[-n])
block_data$IV_tm2 <- c(NA, block_data$IV_tm1[-n])
block_data$IV_tm3 <- c(NA, block_data$IV_tm2[-n])

block_data$PR_tm <- composite_ts$SummerPrecip.mm
block_data$PR_tm1 <- c(NA, block_data$PR_tm[-n])
block_data$PR_tm2 <- c(NA, block_data$PR_tm1[-n])
block_data$PR_tm3 <- c(NA, block_data$PR_tm2[-n])

# remove overlaps from other plots
startyear <- 1996
for (i in 2:nrow(block_data)) {
    if (block_data$year[i] < block_data$year[i - 1]) {
        startyear <- block_data$year[i]
        }
    if (block_data$year[i] == startyear) {
        block_data[i,c("AB_tm1", "NO_tm1", "IV_tm1", "PR_tm1")] <- NA
        block_data[i,c("AB_tm2", "NO_tm2", "IV_tm2", "PR_tm2")] <- NA
        block_data[i,c("AB_tm3", "NO_tm3", "IV_tm3", "PR_tm3")] <- NA
        }
    if (block_data$year[i] == (startyear + 1)) {
        block_data[i,c("AB_tm2", "NO_tm2", "IV_tm2", "PR_tm2")] <- NA
        block_data[i,c("AB_tm3", "NO_tm3", "IV_tm3", "PR_tm3")] <- NA
        }
    if (block_data$year[i] == (startyear + 2)) {
        block_data[i,c("AB_tm3", "NO_tm3", "IV_tm3", "PR_tm3")] <- NA
        }
    }
head(block_data[,1:5],20)

## ----block_lnlp on e120, tidy = TRUE, warning = FALSE, fig.width = 4, fig.height = 4----
observed_AB <- block_lnlp_output_AB[[1]]$model_output$obs
predicted_AB <- block_lnlp_output_AB[[1]]$model_output$pred

observed_ABPR <- block_lnlp_output_ABPR[[1]]$model_output$obs
predicted_ABPR <- block_lnlp_output_ABPR[[1]]$model_output$pred

par(mar = c(4,4,1,1), pty = "s", mgp = c(2.5, 1, 0))
plot_range <- range(c(observed_AB, predicted_AB), na.rm = TRUE)
plot(observed_AB, predicted_AB, xlim = plot_range, ylim = plot_range, xlab = "Observed", ylab = "Predicted")
abline(a = 0, b = 1, lty = 2, col = "darkgrey", lwd = 2)
abline(lm(predicted_AB~observed_AB), col = "black", lty = 3, lwd = 2)

points(observed_ABPR, predicted_ABPR, pch = 2, col = "red")
abline(lm(predicted_ABPR~observed_ABPR), col = "red", lty = 3, lwd = 2)

## ----thrips data---------------------------------------------------------
data(thrips_block)
colnames(thrips_block)

## ----thrips plot, echo = FALSE, fig.width = 6, fig.height = 7------------
par(mar = c(4,4,1,1), mfrow = c(4,1), mgp = c(2.5,1,0))
time_dec <- thrips_block$Year + (thrips_block$Month)/12
plot(time_dec, thrips_block$Thrips_imaginis, type = 'l', col = 'green', ylab = 'Thrips')
plot(time_dec, thrips_block$maxT_degC, type = 'l', col = 'red', ylab = 'maxT (oC)')
plot(time_dec, thrips_block$Rain_mm, type = 'l', col = 'blue', ylab = 'Rain (mm)')
plot(time_dec, thrips_block$Season, type = 'l', col = 'magenta', ylab = 'Season')

## ----univariate thrips, warning = FALSE----------------------------------
ts <- thrips_block$Thrips_imaginis
lib <- c(1, length(ts))
pred <- c(1, length(ts))
simplex_output <- simplex(ts, lib, pred, tau = 1)

## ----rho vs. e for thrips, echo=FALSE, fig.width = 5, fig.height = 3.5, tidy = TRUE----
par(mar = c(4,4,1,1))
plot(simplex_output$E, simplex_output$rho, type = "l", xlab = "Embedding Dimension (E)", ylab = "Forecast Skill (rho)")

## ----smap for thrips, warning = FALSE------------------------------------
smap_output <- list()
smap_output[[1]] <- s_map(ts, lib, pred, E = 4)
smap_output[[2]] <- s_map(ts, lib, pred, E = 8)

## ----rho vs. theta for thrips, echo=FALSE, tidy = TRUE, fig.width = 6, fig.height = 3.5----
par(mar = c(4,4,1,1), mfrow = c(1,2))
plot(smap_output[[1]]$theta, smap_output[[1]]$rho, type = "l", xlim = c(0,4), xlab = "Nonlinearity (theta)", ylab = "Forecast Skill (rho)")
plot(smap_output[[2]]$theta, smap_output[[2]]$rho, type = "l", xlim = c(0,4), xlab = "Nonlinearity (theta)", ylab = "Forecast Skill (rho)")

## ----compute corr matrix for thrips, tidy=TRUE---------------------------
M_corr <- array(NA, dim = c(ncol, ncol), dimnames = list(colnames(thrips_block[3:6]), colnames(thrips_block[3:6])))

for (i in 1:ncol) {
    for (j in 1:ncol) {
        if (i != j) {
            cf_temp <- ccf(x = thrips_block[, 2 + i], y = thrips_block[, 2 + j],
                           type = "correlation", lag.max = 6, plot = FALSE)$acf
            M_corr[i, j] <- max(abs(cf_temp))
            }
        }
}

## ----xmap matrix for thrips, echo=FALSE----------------------------------
head(M_rho)

## ----corr matrix for thrips, echo=FALSE----------------------------------
head(M_corr)

## ----ccm plot, echo=FALSE, tidy = TRUE, fig.width = 5, fig.height = 3.5----
par(mar = c(4,4,1,1), mgp = c(2.5, 1, 0))
plot(xmap_means[[1]]$lib_size, pmax(0, xmap_means[[1]]$rho), type = "l", col = "red",  xlab = "Library Size", ylab = "Cross Map Skill (rho)", ylim = c(0, 1))
lines(xmap_means[[2]]$lib_size, pmax(0, xmap_means[[2]]$rho), col = "blue")
abline(h = M_corr['Thrips_imaginis','maxT_degC'], col = "black", lty = 2)
legend(x = "bottomright", legend = c("Thrips xmap maxT", "maxT xmap Thrips"), col = c("red", "blue"), lwd = 1, inset = 0.02)

## ----rainfall and thrips ccm plot, echo=FALSE, tidy = TRUE, fig.width = 5, fig.height = 3.5----
par(mar = c(4,4,1,1), mgp = c(2.5, 1, 0))
plot(xmap_means[[1]]$lib_size, pmax(0, xmap_means[[1]]$rho), type = "l", col = "red", xlab = "Library Size", ylab = "Cross Map Skill (rho)", ylim = c(0, 1))
lines(xmap_means[[2]]$lib_size, pmax(0, xmap_means[[2]]$rho), col = "blue")
abline(h = M_corr['Thrips_imaginis','Rain_mm'], col = 'black', lty = 2)
legend(x = "topleft", legend = c("Thrips xmap Rain", "Rain xmap Thrips"),
 col = c("red", "blue"), lwd = 1, inset = 0.02)

## ----season and thrips ccm plot, echo=FALSE, fig.width = 5, fig.height = 3.5, tidy = TRUE----
par(mar = c(4,4,1,1), mgp = c(2.5, 1, 0))
plot(xmap_means[[1]]$lib_size, pmax(0, xmap_means[[1]]$rho), type = "l", col = "red", xlab = "Library Size", ylab = "Cross Map Skill (rho)", ylim = c(0, 1))
lines(xmap_means[[2]]$lib_size, pmax(0, xmap_means[[2]]$rho), col = "blue")
abline(h = M_corr['Thrips_imaginis','Season'], col = 'black', lty = 2)
legend(x = "bottomright", legend = c("Thrips xmap Season", "Season xmap Thrips"), col = c("red", "blue"), lwd = 1, inset = 0.02)

## ----significance of randomization test----------------------------------
(sum(M_rho['Thrips_imaginis','Rain_mm'] < rho_surr$Rain) + 1) / (length(rho_surr$Rain) + 1)
(sum(M_rho['Thrips_imaginis','maxT_degC'] < rho_surr$maxT) + 1) / (length(rho_surr$maxT) + 1)

