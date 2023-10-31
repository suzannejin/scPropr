#!/usr/bin/env Rscript

library(argparse)

# parse arguments
parser = ArgumentParser(description='Get the count matrix, transpose, and compress')
parser$add_argument('-i', '--input', type='character', help="Input data, tab separated, with row and column names")
parser$add_argument('-o', '--output', type='character', help="Transposed and compressed output data, comma separated, without row and column names")
parser = parser$parse_args()

# organize data
mat = read.table(parser$input, sep='\t', header=T, row.names=1)
mat = t(mat)

# write output
write.table(mat, gzfile(parser$output), sep=',', row.names=F, col.names=F, quote=F)