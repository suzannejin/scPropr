#!/usr/bin/env Rscript

library(scDesign2)
library(argparse)
library(copula)   
library(data.table)

# parse arguments -----------------------------------------------------------------------

parser <- ArgumentParser(description='Simulate data using scDesign2')
parser$add_argument('model', type='character', help="Model")
parser$add_argument('output', type='character', help="Output simulated count data")
parser$add_argument('--size_factor', type='integer', help="define ratio = (size simulated data) / (size original data)")
parser$add_argument('--slope', type='integer', help="Slope of sequencing depth")
parser$add_argument('--ndata', type='integer', help="Number of simulations for a given slope")
parser$add_argument('--bymax', action='store_true', help="The datasets with smallest and highest seq depth will have specifically <slope> times difference")
parser$add_argument('--bystep', action='store_true', help="One dataset and the one with next higher seq depth will have specifically <slope> times difference. So the datasets with smallest and highest seq depth will have < (ndata - 1) * slope > times difference")
parser = parser$parse_args()

# functions -----------------------------------------------------------------------------

get_seqdepth_bymax <- function(slope, ndata){
    # data = c(1:ndata)
    # seqdepth = c(1, data[1:ndata-1] * slope ) 
    # seqdepth = 1 + (seqdepth - min(seqdepth)) * (slope - 1) / (max(seqdepth) - min(seqdepth))
    step = (slope - 1) / (ndata - 1)
    seqdepth = c(1)
    v = 1
    for (i in 2:ndata){
        v = v + step
        seqdepth[i] = v
    }
    return(seqdepth)
}
get_seqdepth_bystep <- function(slope, ndata){
    data = c(1:ndata)
    seqdepth = c(1, data[1:ndata-1] * slope )
    return(seqdepth)
}
get_nreads_ncells <- function(model, seqdepth){
    ndata = length(seqdepth)
    nreads   = c()
    ncells   = c()
    for (i in 1:length(seqdepth)){
        nreads[i] = round( as.numeric(model$n_read) * seqdepth[i] / ndata )
        ncells[i] = as.numeric(model$n_cell) / ndata
    }
    return(list(nreads, ncells))
}

# calculate nreads and ncells -----------------------------------------------------------

# read model
model = readRDS(parser$model)
message("load model with [nread:", model$n_read, "][ncell:", model$n_cell, "]")

# calculate
if (!is.null(parser$size_factor)){
    seqdepth = c(parser$size_factor)
    nreads = c(as.numeric(model$n_read) * parser$size_factor)
    ncells = c(as.numeric(model$n_cell))
}else if(!is.null(parser$slope)){
    if (parser$bymax){
        seqdepth = get_seqdepth_bymax(parser$slope, parser$ndata)
    }else if(parser$bystep){
        seqdepth = get_seqdepth_bystep(parser$slope, parser$ndata)
    }else{
        stop("You should specify --bystep or --bymax")
    }
    l = get_nreads_ncells(model, seqdepth)
    nreads = l[[1]]
    ncells = l[[2]]
}

# simulate data -------------------------------------------------------------------------

simulated = list()
for (i in 1:length(seqdepth)){
    message("simulate data with [seqdepth:", seqdepth[i], "][nread:", nreads[i], "][ncell:", ncells[i], "]")
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
print(dim(merged))
fwrite(merged, file=parser$output, quote=F, sep=",", row.names=F, col.names=F, compress="gzip")