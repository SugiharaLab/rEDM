simplex <- function(data, lib = c(1, NROW(data)), pred = c(1, NROW(data)), 
                    norm_type = c("L2 norm", "L1 norm"), exclusion_radius = NULL, 
                    E = 1:10, tau = 1, tp = 1, num_neighbors = "e+1", epsilon = NULL)
{
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
           
    # global params
    my_lnlp$set_norm_type(switch(match.arg(norm_type), "L2 norm" = 2, "L1 norm" = 1))
    
    # setup params to run
    if(is.vector(lib))
        lib <- matrix(lib, ncol = 2, byrow = TRUE)
    if(is.vector(pred))
        pred <- matrix(pred, ncol = 2, byrow = TRUE)
    my_lnlp$set_lib(lib)
    my_lnlp$set_pred(pred)
    
    # apply model prediction function to params
    model_output <- lapply(seq_along(params), function(i) {
        my_lnlp$set_params(params$E[i], params$tau[i], params$tp[i], params$nn[i])
        return(my_lnlp$run())
    })
    
    # compute stats
    stats <- cbind(params, do.call(rbind, lapply(model_output, compute_stats)))
    
    return(list(model = my_lnlp, model_output = model_output, stats = stats))
}