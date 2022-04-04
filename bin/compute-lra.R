#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(easyCODA)

# parse arguments
parser = ArgumentParser(description='Compute procrustes analysis')
parser$add_argument('input', type='character', help="Input count data (with zero already imputed)")
parser$add_argument('output', type='character', help="Output stat filename")
parser = parser$parse_args()

# read data
data = fread(parser$input)
if(sum(data[1,]!=1)) data <- data/rowSums(data)

# compute lra
message("computing the exact logratio geometrcy")
data.lra = LRA(data, weight=F)
saveRDS(data.lra, file=parser$output)
message("finished")