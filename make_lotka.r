#!/usr/bin/Rscript
make_lotka_volterra <- function(num_points = 1000,
                                alpha,
                                r )
                                    
{    
    n <- length(r)
    names <- c("x", "y", "z", "u", "v", "w" )
    names <- names[1:n]
    
    data <- matrix(0, nrow = num_points, ncol = n)
    data[1,] <- rep.int(0.2, times = n)
    
    for(i in 1:(num_points-1))
    {
        data[i+1,] <- r * data[i,] * (1 - alpha %*% data[i,]) 
    }
    data <- data.frame(data)
    names(data) <- names
    data$time <- c(1:length(data$x))
    return(data)
}

## A generalized Lotka Volterra matrix with
## DECOUPLED coordinates - first three do not
## interact with the last two.
a <- 0.2
b <- 0.34
r <- c(3.6, 3, 3, 4, 2.5, 2.5)   
alpha <- matrix(c(1 ,  a ,  a ,      0 , 0 ,  0 , 
                  a ,  1 , -a ,      0 , 0 ,  0 , 
                  a , -a ,  1 ,      0 , 0 ,  0 ,
                  ##
                  ##
                  0 ,  0 ,  0 ,      1 , b ,  b ,
                  0 ,  0 ,  0 ,      b , 1 , -b ,
                  0 ,  0 ,  0 ,      b ,-b ,  1),
                nrow = length(r),
                byrow = TRUE)

T = 100000
df <- make_lotka_volterra( T, alpha, r )
write.csv(df, file='lotka.csv', row.names=FALSE, quote=FALSE)
## print( df[ T-2, ] )
## print( df[ T-1, ] )
## print( df[ T  , ] )

