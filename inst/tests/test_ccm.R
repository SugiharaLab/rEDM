data("sardine_anchovy_sst")

anchovy_xmap_sst <- ccm(sardine_anchovy_sst, E = 3, 
                            lib_column = 2, target_column = 5, 
                            lib_sizes = seq(10, 80, by = 10), 
                            random_libs = FALSE)

a_xmap_t_means <- ccm_means(anchovy_xmap_sst)
