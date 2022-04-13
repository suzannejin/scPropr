#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(ggplot2)
library(ggpubr)


# parse arguments
parser = ArgumentParser(description='Compare coefficients computed on relative vs absolute data.')
parser$add_argument('--ori', type='character', help="Original (raw) absolute count data")
parser$add_argument('--abs', type='character', help="Coefficients computed on absolute data")
parser$add_argument('--rel', type='character', help="Coefficients computed on relative data")
parser$add_argument('--features', type='character', help="List of gene names")
parser$add_argument('--refgene', type='character', help="Reference gene name")
parser$add_argument('--output', type='character', help="Output figure filename")
parser$add_argument('--npoint', type='integer', default=1e4, help="Number of points to plot")
parser = parser$parse_args()

# read input files
message('reading input data')
ori      = as.matrix(fread(parser$ori))
abs      = as.matrix(fread(parser$abs))
rel      = as.matrix(fread(parser$rel))
features = fread(parser$features, header=F)$V1
if ( !is.null(parser$refgene) ) {
    pos = which(features == parser$refgene)
    ori = ori[, -pos]
    abs = abs[-pos, -pos]
    features = features[-pos]
}
colnames(ori) = features
colnames(abs) = features
rownames(abs) = features
colnames(rel) = features
rownames(rel) = features
if ( (nrow(abs) != nrow(rel)) || (ncol(abs) != ncol(rel)) ) stop("abs and rel matrices have different dimensiones") 
if ( (ncol(abs) != ncol(ori)) ) stop("abs and ori are based on different number of genes (ncol)") 
if ( length(features) != ncol(abs) ) stop("length(features) != ncol(abs)")

# get dropout and variance 
message('computing dropout and variance for each gene')
tmp           = ori / rowSums(ori) * mean(rowSums(ori))
dropout       = colMeans(ori == 0)
var_genes_abs = apply(ori, 2, var)
var_genes_rel = apply(tmp, 2, var)

# organize data frame
message('organizing data frame and get random sample')
ind  = which( upper.tri(abs, diag=F) , arr.ind = TRUE )
comb = data.frame(
    i         = features[ind[,1]],
    j         = features[ind[,2]],
    abs       = abs[ind],
    rel       = rel[ind],
    dropout_i = dropout[ind[,1]],
    dropout_j = dropout[ind[,2]],
    var_abs_i = var_genes_abs[ind[,1]],
    var_abs_j = var_genes_abs[ind[,2]],
    var_rel_i = var_genes_rel[ind[,1]],
    var_rel_j = var_genes_rel[ind[,2]]
)
set.seed(0); pos = sample(c(1:nrow(comb)), size=parser$npoint)
comb = comb[pos,]
comb[,'dropout'] = rowMeans( data.frame(comb[,'dropout_i'], comb[,'dropout_j']) )
comb[,'var_abs'] = rowMeans( data.frame(comb[,'var_abs_i'], comb[,'var_abs_j']) )
comb[,'var_rel'] = rowMeans( data.frame(comb[,'var_rel_i'], comb[,'var_rel_j']) )
dim(comb)
head(comb)

# plot abs vs relative
message('plotting figure')
g1 = ggplot(comb, aes_string(x='abs', y='rel', color='dropout')) +
     geom_point(alpha=.1, size=.5) + 
     geom_abline(intercept = 0, slope = 1, linetype = 2, color="gray46") +
     xlab("On absolute data") +
     ylab("On relative data") +
     scale_colour_viridis_c() +
     theme(strip.text = element_text(size = 12)) 
g2 = ggplot(comb, aes_string(x='abs', y='rel', color='var_abs')) +
     geom_point(alpha=.1, size=.5) + 
     geom_abline(intercept = 0, slope = 1, linetype = 2, color="gray46") +
     xlab("On absolute data") +
     ylab("On relative data") +
     scale_colour_viridis_c() +
     theme(strip.text = element_text(size = 12))
g3 = ggplot(comb, aes_string(x='abs', y='rel', color='var_rel')) +
     geom_point(alpha=.1, size=.5) + 
     geom_abline(intercept = 0, slope = 1, linetype = 2, color="gray46") +
     xlab("On absolute data") +
     ylab("On relative data") +
     scale_colour_viridis_c() +
     theme(strip.text = element_text(size = 12))
g = ggarrange(g1, g2, g3, ncol=3)
ggsave(g, filename=parser$output, width = 15, height = 4)
message('finished')