#!/usr/bin/env Rscript
# script adapted from get-reduced-dataset-ionas.R

library(data.table)


# read input count data -----------------------------------------------------------------
message('Reading input count data')
count    = '/users/cn/sjin/projects/proportionality/run_scPropr/data/mouseStemCells/mouseStemCells_experimental_full_absolute.csv.gz'
features = '/users/cn/sjin/projects/proportionality/run_scPropr/data/mouseStemCells/mouseStemCells_experimental_full_features.csv'
barcodes = '/users/cn/sjin/projects/proportionality/run_scPropr/data/mouseStemCells/mouseStemCells_experimental_full_barcodes.csv'
count    = as.matrix(fread(count))
features = fread(features, header=F)$V1
barcodes = fread(barcodes, header=F)$V1
colnames(count) = features


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

# plot number of nonzero samples vs number of genes
plot(
    n_nonzero_samples_per_gene_ordered,
    xlab = 'number of genes',
    ylab = 'number of nonzero samples',
    main = 'reverse cumulative of nonzero samples'
)


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

# plot number of common nonzero samples vs number of genes considered
plot(
    n_common_nonzero_sample_set,
    xlab = 'number of genes considered',
    ylab = 'number of common nonzero samples',
    xlim = c(0, 2000)
)
p = 1200
abline(h=n_common_nonzero_sample_set[p], v=p, col='red')


# get the reduced matrix ----------------------------------------------------------------
ipos = history_common_nonzero_sample_set[[p]]
jpos = idx[1:p]
count2 = count[ipos,jpos]
features2 = features[jpos]
barcodes2 = barcodes[ipos]
dim(count2)  
