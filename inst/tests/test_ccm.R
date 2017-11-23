library(rEDM)
data("sardine_anchovy_sst")

anchovy_xmap_sst <- ccm(sardine_anchovy_sst, E = 3, 
                        lib_column = "anchovy", target_column = "np_sst", 
                        lib_sizes = seq(10, 80, by = 10), num_samples = 100, 
                        RNGseed = 42)

a_xmap_t_means <- ccm_means(anchovy_xmap_sst)
plot(a_xmap_t_means$lib_size, a_xmap_t_means$rho, type = "l", 
     xlab = "Library Size", ylab = "Cross Map rho")
