simplex <- function(data, lib = c(1, NROW(data)), pred = c(1, NROW(data)), 
                    norm_type = c("L2 norm", "L1 norm"), exclusion_radius = NULL, 
                    E = 1:10, tau = 1, tp = 1, num_neighbors = "e+1", epsilon = NULL, 
                    stats_only = FALSE)
{
    # check inputs?
    
    # make new model object
    my_lnlp <- new(LNLP)
    
    # setup data
    if(is.vector(data))
    {
        time <- seq_along(data)
        ts <- data
    } else {
        time <- data[,1]
        ts <- data[,2]
    }
    my_lnlp$set_time(time)
    my_lnlp$set_time_series(ts)
           
    # setup norm and pred types
    my_lnlp$set_norm_type(switch(match.arg(norm_type), "L2 norm" = 2, "L1 norm" = 1))
    my_lnlp$set_pred_type(2) # 2 = simplex
    
    # setup lib and pred ranges
    if(is.vector(lib))
        lib <- matrix(lib, ncol = 2, byrow = TRUE)
    if(is.vector(pred))
        pred <- matrix(pred, ncol = 2, byrow = TRUE)
    my_lnlp$set_lib(lib)
    my_lnlp$set_pred(pred)
    
    # TODO: handle exclusion radius
    # TODO: handle epsilon
    
    # setup other params in data.frame
    params = expand.grid(E, tau, tp, num_neighbors)
    names(params) = c("E", "tau", "tp", "nn")
    e_plus_1_index = match(num_neighbors, c("e+1", "E+1", "e + 1", "E + 1"))
    if(any(e_plus_1_index, na.rm = TRUE))
        params$nn <- params$E+1
        
    # apply model prediction function to params
    if(stats_only)
    {
        stats <- lapply(1:NROW(params), function(i) {
            my_lnlp$set_params(params$E[i], params$tau[i], params$tp[i], params$nn[i])
            my_lnlp$run()
            return(my_lnlp$get_stats())
        })
        return(cbind(params, do.call(rbind, stats)))
    }
    
    # else
    output <- lapply(1:NROW(params), function(i) {
        my_lnlp$set_params(params$E[i], params$tau[i], params$tp[i], params$nn[i])
        my_lnlp$run()
        return(list(params = params[i,], 
                    model_output = my_lnlp$get_output(), 
                    stats = my_lnlp$get_stats()))
    })
    return(output)
}