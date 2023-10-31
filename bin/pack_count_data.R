#!/usr/bin/env Rscript

library(argparse)
library(data.table)

# parse arguments
parser = ArgumentParser(description='Add features and barcodes to count matrix, and transpose.')
parser$add_argument('-c', '--count', type='character', help="Input count data")
parser$add_argument('-f', '--features', type='character', help="Features list")
parser$add_argument('-b', '--barcodes', type='character', help="Barcodes list")
parser$add_argument('-o', '--output', type='character', help="Output data")
parser = parser$parse_args()

# read input data
count    = as.matrix( fread(parser$count, header=F) )
features = fread(parser$features, header=F)$V1
barcodes = fread(parser$barcodes, header=F)$V1

# add features and barcodes
count = t(count)
count = cbind(features, count)
colnames(count) = c('GeneID', barcodes)

# write output
write.csv(count, parser$output, row.names=F, quote=F)