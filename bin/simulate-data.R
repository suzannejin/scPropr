library(scDesign2)
library(copula)   

# usage: simulate-data.R <model> <size factor> <output> --------------------------
args = commandArgs(trailingOnly = T)
model = readRDS(args[1])
size_factor = as.numeric(args[2])  # define ratio = (size simulated data) / (size original data)
out = args[3]

# simulate data ------------------------------------------------------------------
message("simulate data")
simulated = simulate_count_scDesign2(model, 
                                     model$ncell, 
                                     total_count_new= model$n_read * size_factor,   # set sequencing depth. Note that there is some stochasticity, so simulated n_read won't match exactly the original n_read
                                     sim_method = 'copula')
simulated = t(simulated)
fwrite(simulated, out, row.names=F, col.names=F)