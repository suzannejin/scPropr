#!/usr/bin/env Rscript

library(scDesign2)
library(argparse)
library(copula)   
library(data.table)

# parse arguments -----------------------------------------------------------------------

parser <- ArgumentParser(description='Simulate data using scDesign2')
parser$add_argument('model', type='character', help="Model")
parser$add_argument('output', type='character', help="Output simulated count data")
parser$add_argument('--slope', type='integer', help="Slope of sequencing depth")
parser$add_argument('--ndata', type='integer', help="Number of simulations for a given slope")
parser$add_argument('--cell_factor', type='double', default=1, help="define ratio = (ncell simulated data) / (ncell original data)")
parser$add_argument('--depth_factor', type='double', default=1, help="define ratio = (seq depth simulated data) / (seq depth original data)")
parser = parser$parse_args()

# functions -----------------------------------------------------------------------------

get_seqdepth <- function(slope, ndata){
    step = (slope - 1) / (ndata - 1)
    seqdepth = c(1)
    v = 1
    for (i in 2:ndata){
        v = v + step
        seqdepth[i] = v
    }
    seqdepth = seqdepth / sum(seqdepth)
    return(seqdepth)
}

# calculate nreads and ncells -----------------------------------------------------------

# read model
model = readRDS(parser$model)
message("loading model with [nread:", model$n_read, "][ncell:", model$n_cell, "]")

# calculate
seqdepth = get_seqdepth(parser$slope, parser$ndata)
nreads   = model$n_read * seqdepth * parser$depth_factor
ncells   = rep( model$n_cell, length(seqdepth) ) / parser$ndata
ncells   = ncells * parser$cell_factor   # TODO think better how we should deal with cell_factor change. Should I also modify the seqdepth, when less cells or so? 

# simulate data -------------------------------------------------------------------------

simulated = list()
for (i in 1:length(seqdepth)){
    message("simulating data with [seqdepth:", seqdepth[i], "][cell_factor:", parser$cell_factor, "][nread:", nreads[i], "][ncell:", ncells[i], "]")
    simulated[[i]] = simulate_count_scDesign2(list(model),   # this function reads model as a list of models
                                     model$n_cell, 
                                     total_count_new = nreads[i],   
                                     n_cell_new = ncells[i],
                                     sim_method = 'copula')
    simulated[[i]] = t(simulated[[i]])   # since scDesign creates matrix of row=gene, col=cell
}

# merge simulated data
merged = simulated[[1]]
if (length(simulated) > 1){
    for (i in 2:length(simulated)){
        merged = rbind(merged, simulated[[i]])
    }
}

# write output
message("writing simulated data [", nrow(merged), "][", ncol(merged), "]")
fwrite(merged, file=parser$output, quote=F, sep=",", row.names=F, col.names=F, compress="gzip")