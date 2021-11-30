library(data.table)
library(stringr)
library(propr)
# library(propr, lib.loc=paste0(.libPaths()[1], "/propr_sjin"))

# usage: get-corr.R <count> <method> <output> -------------------------------------------
args = commandArgs(trailingOnly = T)
countfile = args[1]
coef = args[2]
output = args[3]

# read count data
count = fread(countfile)

# compute association coefficients ------------------------------------------------------
corr = propr(count, metric=coef, ivar=NA, p=0)@matrix  
# corr = propr(dat, metric=coef, ivar=NA, p=20)
# corr = updateCutoffs(corr, seq(min(corr@matrix), max(corr@matrix), length=102)[2:101])
fwrite(corr, file=output, quote=F, sep=",", row.names=F, col.names=F, compress="gzip")