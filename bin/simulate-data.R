library(scDesign2)
library(copula)   
library(data.table)

# usage: simulate-data.R <model> <size factor> <output> --------------------------
args = commandArgs(trailingOnly = T)
model = readRDS(args[1])
size_factor = as.numeric(args[2])  # define ratio = (size simulated data) / (size original data)
comb_factor = as.numeric(args[3])  # define ratio = (ncell original data) / (ncell simulated data)
output = args[4]

# simulate data ------------------------------------------------------------------
message("simulate data")
new_read = as.numeric(model$n_read) * size_factor / comb_factor   # here it is also divided by comb factor, because otherwise, with less number of cells, each cell will receive more reads
new_cell = floor( as.numeric(model$n_cell) / comb_factor )
simulated = simulate_count_scDesign2(list(model),   # this function reads model as a list of models
                                     model$n_cell, 
                                     total_count_new = new_read,   # set sequencing depth. Note that there is some stochasticity, so simulated n_read won't match exactly the original n_read
                                     n_cell_new = new_cell,
                                     sim_method = 'copula')
simulated = t(simulated)
fwrite(simulated, file=output, quote=F, sep=",", row.names=F, col.names=F, compress="gzip")