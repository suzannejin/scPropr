#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(ggplot2)

# parse arguments
parser = ArgumentParser(description='Plot mean expression vs variance')
parser$add_argument('input', type='character', help="Input count data")
parser$add_argument('name', type='character', help="Dataset name - it will be used as plot title")
parser$add_argument('outdir', type='character', help="Output directory")
parser$add_argument('--bygene', action='store_true', help="Plot gene mean expression vs variance")
parser$add_argument('--bycell', action='store_true', help="Plot cell mean expression vs variance")
parser = parser$parse_args()

count = fread(parser$input)
if (parser$bygene){
    by = 2
    filename = paste0(parser$outdir,'/mean-vs-variance-bygene.png')
}else{
    by = 1
    filename = paste0(parser$outdir,'/mean-vs-variance-bycell.png')
}

# calculate mean expression, and variance
expr_mean = apply(count, by, function(x) mean(x))
expr_var  = apply(count, by, function(x)  var(x))

# plot
png(filename=filename, width=12, height=4, units="in", res=100)
par(mfrow=c(1,3))
plot(expr_mean, expr_var, xlab="expression mean",  ylab="variance")
plot(expr_mean, expr_var, xlab="expression mean",  ylab="variance", log='y')
plot(expr_mean, expr_var, xlab="expression mean",  ylab="variance", log='xy')
mtext(parser$name, side = 3, line = - 2, outer = TRUE)