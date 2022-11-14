#!/usr/bin/env Rscript

library(argparse)
library(data.table)

# parse arguments
parser = ArgumentParser(description='Select genes/cells from count data')
parser$add_argument('-i', '--input', type='character', help="Input count data")
parser$add_argument('-f', '--features', type='character', help="Features (col) name")
parser$add_argument('-b', '--barcodes', type='character', help="Barcodes (row) id")
parser$add_argument('-o', '--output', type='character', help="Output count data")
parser$add_argument('-o2', '--output2', type='character', help="Updated features name")
parser$add_argument('-o3', '--output3', type='character', help="Updated barcodes id")
parser$add_argument('--filter_gene', type='double', help="Remove the genes with dropout > threshold")
parser$add_argument('--filter_cell', type='double', help="Remove the cells with dropout > threshold")
parser$add_argument('--select_gene', type='character', help="File specifying the genes to be used")
parser$add_argument('--select_cell', type='character', help="File specifying the cells to be used")
parser$add_argument('--ncell', type='character', help="Number of cells wanted")
parser = parser$parse_args()


# read input count data
count    = fread(parser$input, header=F)
features = fread(parser$features, header=F)$V1
barcodes = fread(parser$barcodes, header=F)$V1

# filter genes 
if (is.numeric(parser$filter_gene)){
    parser$filter_gene = parser$filter_gene / 100
    message("filter genes with dropout > ", parser$filter_gene)
    dropout  = colMeans(count==0)
    pos      = which(dropout <= parser$filter_gene)
    count    = count[,..pos]
    features = features[pos]
}
# filter cells
if (is.numeric(parser$filter_cell)){
    parser$filter_cell = parser$filter_cell / 100
    message("filter cells with dropout > ", parser$filter_cell)
    dropout  = rowMeans(count==0)
    pos      = which(dropout <= parser$filter_cell)
    count    = count[pos,]
    barcodes = barcodes[pos]
}
# select genes
if (!is.null(parser$select_gene)){
    message("select genes from ", parser$select_gene)
    tmp      = fread(parser$select_gene, header=F)$V1
    pos      = which(features %in% tmp)
    count    = count[,..pos]
    features = features[pos]
}
# select cells
if (!is.null(parser$select_cell)){
    message("select cells from ", parser$select_cell)
    tmp      = fread(parser$select_cell, header=F)$V1
    pos      = which(barcodes %in% tmp)
    count    = count[pos,]
    barcodes = barcodes[pos]
}
# specific number of cells
if (!is.null(parser$ncell)){
    message("get only ", parser$ncell, " cells")
    set.seed(0); pos = sample(1:nrow(count), parser$ncell)
    count    = count[pos,]
    barcodes = barcodes[pos]
}

print(dim(count))
print(mean(count==0))

# write output files
fwrite(count, file=parser$output, quote=F, sep=",", row.names=F, col.names=F, compress="gzip")
if (!is.null(parser$output2)) fwrite(list(features), file=parser$output2, quote=F, sep=",", row.names=F, col.names=F)
if (!is.null(parser$output3)) fwrite(list(barcodes), file=parser$output3, quote=F, sep=",", row.names=F, col.names=F)