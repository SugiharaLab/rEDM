#!/usr/bin/Rscript
library(rEDM)
library(data.table)


## Estimate log likelihood of a time series using SSE.
logLikelihood <- function(ts,
                          E = c(1:15),
                          lib = c(1,nrows(ts)),
                          pred = lib )
{
    output <- simplex(ts,
                      E = E,
                      lib = lib,
                      pred = pred )
       
    ## How much data we actually have (subtract since we use time
    ## lags)
    n <- length(ts) - E
    
    ## assume the variance is 1 but we may need to estimate it.
    sig2 <- 1
    
    ## Since sum-squared-errors is hopefully approximating the part
    ## inside the exponent in a Gaussian likelihood, we use that (and
    ## enter sigma^2 later).
    sse <- n * output$rmse ** 2 
    
    ## Log likelihood. We take the likelihood to be a product of
    ## n Gaussians with noise sig2 and sum-squared-errors
    ## as it is usually defined.
    return( n / 2 * log(2*pi*sig2)  - sse / (2*sig2) )
}


ic <- function(ts,
               E = c(1:15),
               lib = c(1,length(ts)),
               pred = lib,
               criterion = "AIC",
               showPlot = FALSE,
               title = "")
    
    ## This funciton estimates the Akaike / Bayesian Information
    ## Criterion. Likelihood is estimated in the function logLikelihood
    ## above. Number of model parameters is taken to equal E the embedding
    ## dimension.
    ##     Parameters used are the same as the ones used elsewhere in the
    ## rEDM library or just self explanatory.
{
    
    ## Get the log likelihood. Calculating this is:
    LL <- logLikelihood(ts,
                        E,
                        lib = lib,
                        pred = pred )

    n <- length(ts) - E
    
    ## The parameters are lagged coordinates (we have E of those).
    k <- E

    ## We don't need to include other variables that are constants
    ## independent of E and theta.
    c = substr(criterion,1,1)
    if( c == "A" || c == "a" ){
        ic <- 2*k - 2*LL + 2*k*(k+1)/(n-k-1) 
    } else if ( c == "B" || c == "b") {
        ic <- log(n)*k - 2*LL 
    }

    if ( showPlot || title != "")
    {
        x11()
        plot(E,
             ic,
             type = "l", 
             xlab = "Embedding dimension",
             ylab = criterion )
        title( main = title )
        message("Press Return To Continue")
        invisible(readLines("stdin", n=k))
    }

    ## Return ic but do not print if it is not caught
    ic <- ic
}



## Tests \ examples
###################
data(thrips_block)
ts <- thrips_block$Thrips_imaginis
ic( ts, title = "Thrips imaginis" ) 

###################
data(tentmap_del)
ts <- tentmap_del
ic(ts,
   lib = c(1,100),
   pred = c(201,500),
   title = "Tent Map" )
 
###################
data(block_3sp)
ts <- block_3sp$x_t
ic( ts, title = "Block Three Species - x" )

## x11()
## plot(E,
##      simplex( ts, E = E )$rho,
##      type = "l", 
##      xlab = "Embedding dimension",
##      ylab = "rho" )
## title( main= "Block Three Speceis - x (correlation)" )
## hold(1)

################
data(sardine_anchovy_sst)
ts <- sardine_anchovy_sst$anchovy
ic( ts, title = "Sardine, Anchovy, SST - Anchovy" )

ts <- sardine_anchovy_sst$sardine
ic( ts, title = "Sardine, Anchovy, SST - Sardine" )

ts <- sardine_anchovy_sst$sio_sst
ic( ts, title = "Sardine, Anchovy, SST - SIO SST" )

ts <- sardine_anchovy_sst$np_sst
ic( ts, title = "Sardine, Anchovy, SST - NP(??) SST" )

    
