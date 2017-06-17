#!/usr/bin/Rscript
% This script finds the best predictors and their number
for prediction of the abundancy of Dinoflagellate A, Dinoflagellate B
and Cyanobacteria, based on the model described in the data
entered from the file "monod_2box_2bloom.species.csv". The idea
is to first learn a model that predicts algal bloom for each animal
and then see that unless some serios changes are made, even these models
cannot predict algal blooms when combined (???).
%

library(rEDM)

predict <- function(target, dataFile)
{
    
    dat <- read.csv( dataFile )
    
    ## This is the name of the variable we are trying to predict
    ## Hence, we "train" our model to predict algal blooms for
    ## this (and only this) variable
    ## target = commandArgs(trailingOnly = TRUE)
   
    
    ## Size of library and pred...
    lib <- c(1, NROW(dat))
    pred <- c(1, NROW(dat))

    ## Removing variables. Remove day variable because we don't need it
    ## explicitly.  Remove the all animals except the one we focus on
    ## since here we're trying to predict only dino A blooms.
    drops <- list( "day", "DinoA", "DinoB", "CynB")
    ind <- match( target, drops )
    drops <- drops[ -ind ]
    dat[, (names(dat) %in% drops) ] <- NULL
    
    ## Generally, a log transform is not a good idea but you should try it
    ## and see for yourself!!!
    ## dat$DinoA <- log(dat$ DinoA )
    
    ## Normalize data to heave mean = 0, sd = 1.
    dat <- data.frame(lapply( dat, scale ))
    
    ## Find the best time lag. This is the "embedding dimension" for the
    ## system.
    simplex_output <- simplex(dat[, target], lib, pred)
    E <- which.max( simplex_output$rho )
    
    ## Create an extended data frame with a bunch of lagged
    ## coordinates. Then we choose only the ones that are actually useful
    ## for prediciton:
    
    ## For every column in the data...
    for( col in names( dat ) ){
        
        ## ...place the column in the "current" variable ...
        current <- dat[,col]
        
        ## if( identical(col,target) ) 
        ##     n = E-1
        ## else 
        ##     n = E-1
    
        ## ...and for every lag between 1 and 10;
        for( lag in seq( 1, E-1 ) ) {

            ## Lag the temp variable by one time step ...
            tmp <- c(NA,current[1:nrow(dat)-1])
            
            ## ... then add it to the data fram with its appropriate name.
            dat[, paste0(col,"_",lag)] <- tmp
            
            ## Then, place the new lagged variable in the "current" and repeat. 
            current <- tmp
        }
        
    }
    
    ## Show the new data frame.
    ## dat <- dat[ , order(names(dat))]
    ## print( names( dat ) )
    
    ## Intially, the best rho is zero and the best column is empty
    bestRho = 0
    bestCols <- NULL
    
    ## Generate an empty file to hold all our results
    filename <- paste0("best_",target,"_predictors.txt")
    if (file.exists(filename))
        file.remove(filename)
    
    ## This list contains all possible combinations of column names
    ## (without replacement and disregarding ordering).
    colList <- combn( names(dat), E, simplify = FALSE )

    ## For every possible combination of variables...
    for( cols in colList) {

        ## ... use these variables to make predicitions...
        output <- block_lnlp(dat,
                             lib = lib,
                             pred = pred,
                             columns = cols,
                             target_column = target,
                             stats_only = TRUE,
                             first_column_time = FALSE)

        ## ... and record the performance of the predictor.
        currRho <- output$rho

        ## Update best rho and columns and save data (only if we need to,
        ## of course).
        if (currRho > bestRho){
            bestCols <- cols
            bestRho  <- currRho
            
            ## Save (and print?) the top performing columns.
            bestCols[ E + 1 ] <- bestRho
            newline <- paste0( paste( bestCols, collapse = ' , ' ) )
            ## print( newline )
            write.table(newline,
                        file = filename,
                        append = TRUE,
                        row.names = FALSE,
                        col.names = FALSE )        
            
        }
    }
}

dataFile <- "monod_2box_2bloom.species.csv"
predict( "DinoA", dataFile )
predict( "DinoB", dataFile )
predict( "CynB" , dataFile )
## obs  <- predictor[[1]]$model_output$obs
## pred <- predictor[[1]]$model_output$pred
    
## par(mar = c(4, 4, 1, 1), pty = "s")
## plot_range <- range(c(obs, pred), na.rm = TRUE)

## x11()
## plot(obs,
##      pred,
##      xlim = plot_range,
##      ylim = plot_range,
##      xlab = "Observed dino A", 
##      ylab = "Predicted dino A")
## abline(a = 0, b = 1, lty = 2, col = "blue")



## message("Press Return To Continue")
## invisible(readLines("stdin", n=1))
# print( dat[1:4,] )


