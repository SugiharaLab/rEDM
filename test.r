#!/usr/bin/Rscript
hold <- function(k=1){
    message("Press Return To Continue")
    invisible(readLines("stdin", n=k))
}

library(rEDM)
## data(tentmap_del)
## lib <- c(1, 100)
## pred <- c(201, 500)
## ts <- tentmap_del
## glm    <- s_map(ts, lib, pred, E = 2, glm = TRUE )
## no_glm <- s_map(ts, lib, pred, E = 2 )
## par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0))
## x11()
## plot(glm$theta,
##      glm$rho,
##      type = "l",
##      col = "red",
##      xlab = "Nonlinearity (theta)", 
##      ylab = "Forecast Skill (rho)",
##      ylim = c(0.8,1) )
## lines(no_glm$theta,
##       no_glm$rho,
##       type = "l",
##       col = "blue" )

## hold(0)





data(thrips_block)
ts <- thrips_block$Thrips_imaginis
## data(sardine_anchovy_sst)
## names( sardine_anchovy_sst )
## ts <- sardine_anchovy_sst$sardine## anchovy
output <- simplex( ts, E = c(1:15) )

x11()
plot(output$E,
     output$rmse,
     type = "l",
     col = "red",
     ylab = "Forecast Skill (rho)" )
hold(1)

## E <- 4
## theta <- seq(0.0, 8, 0.05)
## glm    <- s_map(ts, E = E, theta = theta, glm = TRUE  )
## no_glm <- s_map(ts, E = E, theta = theta, glm = FALSE )
## x11()
## plot(glm$theta,
##      glm$rho,
##      type = "l",
##      col = "red",
##      ylab = "Forecast Skill (rho)",
##      ylim = c(0,1) )

## lines(no_glm$theta,
##       no_glm$rho,
##       type = "l",
##       col = "blue" )
## hold(1)
warnings()
## hold(1)
## data(sardine_anchovy_sst)
## anchovy_xmap_sst <- ccm(sardine_anchovy_sst, E = 3, lib_column = "anchovy", 
##     target_column = "np_sst", lib_sizes = seq(10, 80, by = 10), random_libs = FALSE, glm = TRUE )
## sst_xmap_anchovy <- ccm(sardine_anchovy_sst, E = 3, lib_column = "np_sst", target_column = "anchovy", 
##     lib_sizes = seq(10, 80, by = 10), random_libs = FALSE, glm = TRUE )
