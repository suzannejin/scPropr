library(data.table)
library(stringr)
library(propr)
# library(propr, lib.loc=paste0(.libPaths()[1], "/propr_sjin"))

# usage: get-corr.R <count> <method> <output> -------------------------------------------
args = commandArgs(trailingOnly = T)
countfile = args[1]
coef = args[2]
output = args[3]
refgene = args[4]

# read count data
count = fread(countfile)
ngenes = dim(count)[1]

# compute association coefficients ------------------------------------------------------
corr = propr(count, metric=coef, ivar=NA, p=0)@matrix  
# corr = propr(dat, metric=coef, ivar=NA, p=20)
# corr = updateCutoffs(corr, seq(min(corr@matrix), max(corr@matrix), length=102)[2:101])

# keep the same matrix dimensions to number of genes, so that we can compare different matrices
### since proprALR returns a matrix with n-1 genes, therefore propr will also generate a correlation matrix of n-1 
### we need to add the lost column and row to keep the same matrix dimensions than other cases.
if (refgene != "NA"){
    corr = as.data.frame(corr)
    refgene = as.numeric(refgene)
    corr = tibble::add_column(corr, .before=refgene, ref=NA)
    corr = tibble::add_row(corr, .before=refgene, ref=NA)
    corr = as.matrix(corr)
}

# write output
fwrite(corr, file=output, quote=F, sep=",", row.names=F, col.names=F, compress="gzip")