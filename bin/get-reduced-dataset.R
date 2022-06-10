#!/usr/bin/env Rscript

library(argparse)
library(data.table)

# parse arguments
parser = ArgumentParser(description='Get a reduced set of data with no zero values')
parser$add_argument('-i', '--input', type='character', help="Input count data")
parser$add_argument('-f', '--features', type='character', help="Features (col) name")
parser$add_argument('-b', '--barcodes', type='character', help="Barcodes (row) id")
parser$add_argument('-o', '--output', type='character', help="Output count data")
parser$add_argument('-o2', '--output2', type='character', help="Updated features name")
parser$add_argument('-o3', '--output3', type='character', help="Updated barcodes id")
parser$add_argument('--ncell', type='character', help="Number of cells in the output matrix wanted")
parser$add_argument('--ngene', type='character', help="Number of genes in the output matrix wanted")
parser = parser$parse_args()


# read input files ----------------------------------------------------------------------

count    = as.matrix(fread(parser$input))
features = fread(parser$features, header=F)$V1
barcodes = fread(parser$barcodes, header=F)$V1


# find the nonzero set of samples for each gene individually ----------------------------

message('Getting nonzero set of samples for each gene individually')
set_nonzero_samples_per_gene = list()
n_nonzero_samples_per_gene   = c()
for (j in 1:ncol(count)){
    pos = which(count[,j] > 0)
    set_nonzero_samples_per_gene[[j]] = pos
    n_nonzero_samples_per_gene[j] = length(pos)
}
idx = order(n_nonzero_samples_per_gene, decreasing=T)
n_nonzero_samples_per_gene_ordered = n_nonzero_samples_per_gene[idx]


# find a set of genes that has *common* nonzero samples --------------------------------
# as a result one gets a matrix of genes and samples with no zero values

message('Finding set of genes with common nonzero samples')

# start from the genes with most nonzero samples, 
# recursively get the intersection of *common* nonzero samples and new nonzero samples from the gene i
first = idx[[1]]
common_nonzero_sample_set = set_nonzero_samples_per_gene[[first]]
n_common_nonzero_sample_set = c()
n_reduced_gene_set = c()
history_common_nonzero_sample_set = list()
for (c in 1:length(idx)){
    j = idx[c]
    current_nonzero_sample_set =  set_nonzero_samples_per_gene[[j]]
    common_nonzero_sample_set = intersect(common_nonzero_sample_set, current_nonzero_sample_set)
    n_common_nonzero_sample_set[c] = length(common_nonzero_sample_set)
    history_common_nonzero_sample_set[[c]] = common_nonzero_sample_set
}


# get the set with similar number of genes and samples that have no zeros ---------------

df = data.frame(
    ii=n_common_nonzero_sample_set, 
    jj=1:length(n_common_nonzero_sample_set)
    )
df = df[which(df$ii > 500),]
df$dist = abs(df$ii - df$jj)
df = df[df$dist <= 5,]
if (nrow(df) == 0) stop("Error: conditions not met")
pos = which(df$dist == min(df$dist))
df  = df[pos,]
if (length(pos) > 1){
    pos = which(df[,1] == max(df[,1]))
    df  = df[pos,]
}
p = df[,2]


# get matrix ----------------------------------------------------------------------------

ipos = history_common_nonzero_sample_set[[p]]
jpos = idx[1:p]
count    = count[ipos,jpos]
features = features[jpos]
barcodes = barcodes[ipos]
dim(count) 

if (any(count == 0)) stop("Error: obtained data with zero values")


# restrinct size, if wanted
if (!is.null(parser$ngene)){
    message("getting only ", parser$ngene, " genes")
    set.seed(0); pos = sample(1:ncol(count), parser$ngene)
    count    = count[,pos]
    features = features[pos]
}
if (!is.null(parser$ncell)){
    message("getting only ", parser$ncell, " cells")
    set.seed(0); pos = sample(1:nrow(count), parser$ncell)
    count    = count[pos,]
    barcodes = barcodes[pos]
}
if (!is.null(parser$ngene) && !is.null(parser$ncell)){
    if (nrow(count) != parser$ncell && ncol(count) != parser$ngene) stop("Error with parsing the correct [ncell,ngene] sizes.")
}

# write output files
message("writing output with size [", nrow(count), ",", ncol(count), "]")
fwrite(count, file=parser$output, quote=F, sep=",", row.names=F, col.names=F, compress="gzip")
if (!is.null(parser$output2)) fwrite(list(features), file=parser$output2, quote=F, sep=",", row.names=F, col.names=F)
if (!is.null(parser$output3)) fwrite(list(barcodes), file=parser$output3, quote=F, sep=",", row.names=F, col.names=F)