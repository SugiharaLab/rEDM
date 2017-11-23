## ------------------------------------------------------------------------
library(rEDM)
data(block_3sp)

x <- block_3sp$x_t
y <- block_3sp$y_t

## ------------------------------------------------------------------------
concatenated_xy <- c(x, y)
lib_x <- c(1, length(x))
lib_y <- length(x) + c(1, length(y))

## ------------------------------------------------------------------------
simplex_out_x <- simplex(concatenated_xy, lib = lib_x, pred = lib_x)
best_E_x <- simplex_out_x$E[which.max(simplex_out_x$rho)]

copred_x_to_y <- simplex(concatenated_xy, lib = lib_x, pred = lib_y, E = best_E_x)

## ------------------------------------------------------------------------
simplex_out_y <- simplex(concatenated_xy, lib = lib_y, pred = lib_y)
best_E_y <- simplex_out_y$E[which.max(simplex_out_y$rho)]

copred_y_to_x <- simplex(concatenated_xy, lib = lib_y, pred = lib_x, E = best_E_y)

## ------------------------------------------------------------------------
to_plot <- data.frame(label = c("prediction of x (from x)", 
                                "coprediction of x (from y)", 
                                "prediction of y (from y)", 
                                "coprediction of y (from x)"), 
                      rbind(simplex_out_x[which.max(simplex_out_x$rho), ], 
                            copred_y_to_x, 
                            simplex_out_y[which.max(simplex_out_y$rho), ], 
                            copred_x_to_y)
                      )

## ---- fig.height = 5-----------------------------------------------------
library(ggplot2)
ggplot(to_plot, aes(x = label, y = rho)) + 
    geom_col() + theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

