## ------------------------------------------------------------------------
library(rEDM)
data(paramecium_didinium)

## ------------------------------------------------------------------------
vars <- names(paramecium_didinium)[2:3] # c("paramecium", "didinium")

# generate all combinations of lib_column, target_column, tp
params <- expand.grid(lib_column = vars, 
                      target_column = vars, 
                      tp = -10:10)

# throw out cases where lib == target
params <- params[params$lib_column != params$target_column, ]

# E = 3 is optimal or very close to optimal for both vars
# In other circumstances, we should use the best univariate E for each lib_column
E <- 3

## ---- warning = FALSE----------------------------------------------------
output <- do.call(rbind, lapply(seq_len(NROW(params)), function(i) {
    ccm(paramecium_didinium, E = 3, 
        lib_sizes = NROW(paramecium_didinium), random_libs = FALSE, 
        lib_column = params$lib_column[i], 
        target_column = params$target_column[i], 
        tp = params$tp[i], silent = TRUE)
}))

## ------------------------------------------------------------------------
output$direction <- paste(output$lib_column, "xmap to", output$target_column)

## ---- fig.width = 6------------------------------------------------------
library(ggplot2)
ggplot(output, aes(x = tp, y = rho, color = direction)) + 
    geom_line() + theme_bw()

