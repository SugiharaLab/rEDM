#!/usr/bin/Rscript
library(rEDM)
data(tentmap_del)
lib <- c(1, 10)
pred <- c(201, 204)
ts <- tentmap_del
s_map(ts, lib, pred, E = 2)
## warnings()
