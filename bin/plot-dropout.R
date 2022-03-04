#!/usr/bin/env Rscript

library(data.table)

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
    stop("Usage: plot-dropout.R <count data> <dataset name> <output directory>", call.=FALSE)
}
count   = args[1]   # cells x genes
dataset = args[2]
outdir  = args[3]

# read count file
count = fread(count)

# total dropout 
dropout = mean(count==0, na.rm=T)
write(dropout, file=paste0(outdir,'/dropout'))

# dropout histogram
r = rowMeans(count==0, na.rm=T)  # dropout per cell
c = colMeans(count==0, na.rm=T)  # dropout per gene
png(filename=paste0(outdir,'/dropout.png'), width=500, height=280)
par(mfrow=c(1,2), oma=c(0,0,2,0))
hist(r, xlab="Dropout per cell", main=NA)
hist(c, xlab="Dropout per gene", main=NA)
mtext(dataset, line=0, side=3, outer=TRUE, cex=1.2, font=2)

# no zero genes/cells
nozero_r = length(which(r==0))
nozero_c = length(which(c==0))
write(nozero_r, file=paste0(outdir,'/nozero_cells'))
write(nozero_c, file=paste0(outdir,'/nozero_genes'))