#!/usr/bin/Rscript
library(rEDM)

hold <- function(k=0){
    message("Press Return To Continue")
    invisible(readLines("stdin", n=k))
}

block <- read.csv( "lotka.csv", nrows = 2*T )


## Create the sum of the three species
block$sum <- block$y + block$w  

x11()
## plot(block_3sp$time,
##      block_3sp$total,
##      xlab = "time", 
##      ylab = "abundance",
##      type = "l",
##      col = "black" )

## lines(block_3sp$time,
##       block_3sp$x_t,
##       type = "l",
##       col = "red" )


T = 1000
lib  <- c( 1  , T   )
pred <- c( T+1, 2*T )
simplex_sum <- simplex(block$sum, lib, pred)
simplex_y   <- simplex(block$y,   lib, pred)
simplex_w   <- simplex(block$w,   lib, pred)

plot(simplex_sum$E,
     simplex_sum$rho,
     type = "l",
     col  = "black",
     xlab = "Embedding Dimension (E)", 
     ylab = "Forecast Skill",
     ylim = c(0.85,1) )

lines(simplex_y$E,
      simplex_y$rho,
      col = "red" )


lines(simplex_w$E,
      simplex_w$rho,
      col = "green" )


hold(1)

block_lnlp_output <- block_lnlp(block,
                                lib = lib,
                                pred = pred,
                                columns = c("x_t","x_t-1","y_t"),
                                target_column = "x_t",
                                stats_only = FALSE,
                                first_column_time = TRUE)

observed <- block_lnlp_output[[1]]$model_output$obs
predicted <- block_lnlp_output[[1]]$model_output$pred


par(mar = c(4, 4, 1, 1), pty = "s")

plot_range <- range(c(observed, predicted),
                    na.rm = TRUE)

x11()
plot(observed,
     predicted,
     xlim = plot_range,
     ylim = plot_range,
     xlab = "Observed", 
     ylab = "Predicted")

abline(a = 0, b = 1, lty = 2, col = "blue")

hold()


