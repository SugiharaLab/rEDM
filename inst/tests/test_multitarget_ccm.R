library(rEDM)
library(purrr)
data("sardine_anchovy_sst")

n <- NROW(sardine_anchovy_sst)

ccm_results <- rbind(
    map_df(1:10, function(E) {
        ccm(sardine_anchovy_sst, E = E, tp = -1, 
            lib_column = "anchovy", target_column = c("sardine", "np_sst", "sio_sst"), 
            lib_sizes = n, num_samples = 1, random_libs = FALSE, 
            silent = TRUE)
    }), 
    map_df(1:10, function(E) {
        ccm(sardine_anchovy_sst, E = E, tp = 0, 
            lib_column = "anchovy", target_column = c("sardine", "np_sst", "sio_sst"),  
            lib_sizes = n, num_samples = 1, random_libs = FALSE, 
            silent = TRUE)
    })
)

#save(ccm_results, file = "~/Desktop/ccm_out_dev.Rdata")
