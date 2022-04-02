#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(stringr)
library(propr)

# parse arguments -----------------------------------------------------------------------
parser = ArgumentParser(description='Compute correlation coefficients given count data.')
parser$add_argument('-i', '--input', type='character', help="Input count data")
parser$add_argument('-o', '--output', type='character', help="Output count data")
parser$add_argument('-o2', '--output2', type='character', help="Standard output")
parser$add_argument('-m', '--method', type='character', help="Correlation method")
parser = parser$parse_args()

# read count data
count = fread(parser$input)
print(dim(count))

# compute association coefficients ------------------------------------------------------
message("calculating coefficients with method ", parser$method)
stdout = vector('character')  # connect standard output
con = textConnection('stdout', 'wr', local=T)
sink(con)
corr   = propr(count, metric=parser$method, ivar=NA, p=0)@matrix  
# corr = propr(dat, metric=coef, ivar=NA, p=20)
# corr = updateCutoffs(corr, seq(min(corr@matrix), max(corr@matrix), length=102)[2:101])
sink()
close(con)
print(dim(corr))

# write output
fwrite(corr, file=parser$output, quote=F, sep=",", row.names=F, col.names=F, compress="gzip")
writeLines(stdout, parser$output2)
message('finished')