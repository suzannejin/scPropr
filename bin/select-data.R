#!/usr/bin/env Rscript

library(argparse)
library(data.table)

# parse arguments
parser = ArgumentParser(description='Select genes/cells from count data')
parser$add_argument('-i', '--input', type='character', help="Input count data")
parser$add_argument('-o', '--output', type='character', help="Output count data")
parser$add_argument('--filter_gene', type='double', help="Remove the genes with dropout > threshold")
parser$add_argument('--filter_cell', type='double', help="Remove the cells with dropout > threshold")
parser$add_argument('--select_gene', type='character', help="File specifying the genes to be used")
parser$add_argument('--select_cell', type='character', help="File specifying the cells to be used")
parser = parser$parse_args()

# read input count data
count = fread(parser$input)

# filter genes 
if (is.numeric(parser$filter_gene)){
    message("filter genes with dropout > ", parser$filter_gene)
    dropout = colMeans(count==0)
    pos     = which(dropout <= parser$filter_gene)
    count   = count[,..pos]
}
# filter cells
if (is.numeric(parser$filter_cell)){
    message("filter cells with dropout > ", parser$filter_cell)
    dropout = rowMeans(count==0)
    pos     = which(dropout <= parser$filter_cell)
    count   = count[pos,]
}
# select genes
if (!is.null(parser$select_gene)){
    message("select genes from ", parser$select_gene)
    pos   = fread(parser$select_gene, header=F)$V1
    count = count[,..pos]
}
# select cells
if (!is.null(parser$select_cell)){
    message("select cells from ", parser$select_cell)
    pos   = fread(parser$select_cell, header=F)$V1
    count = count[pos,]
}

print(dim(count))
print(mean(count==0))

# write output file
fwrite(count, file=parser$output, quote=F, sep=",", row.names=F, col.names=F, compress="gzip")