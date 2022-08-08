#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(reshape2)

# parse arguments -----------------------------------------------------------------------
parser = ArgumentParser(description='Get one table with columns: gene i, gene j, cor abs, cor rel.')
parser$add_argument('--abs', type='character', help="Absolute count data")
parser$add_argument('--rel', type='character', help="Relative count data")
parser$add_argument('--out', type='character', help="Output table")
parser$add_argument('--ref', type='character', help="Reference gene")
parser$add_argument('--features', type='character', help="Features list")
parser = parser$parse_args()

# read count data
message('reading input data')
abs = as.matrix(fread(parser$abs))
rel = as.matrix(fread(parser$rel))
features = fread(parser$features, header=F)$V1
if (!is.null(parser$ref)){
    nref = which( toupper(features) == toupper(parser$refgene) )
    abs  = abs[-nref, -nref]
    features = features[-nref]
}
if (nrow(abs) != nrow(rel) || ncol(abs) != ncol(rel) ) stop(dim(abs), ' != ', dim(rel))
colnames(abs) = features; rownames(abs) = features
colnames(rel) = features; rownames(rel) = features
print(dim(abs))

# transforming multdimensional array into table
abs = reshape2::melt(abs)
rel = reshape2::melt(rel)


