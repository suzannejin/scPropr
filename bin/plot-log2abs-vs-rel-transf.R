#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(ggplot2)
library(ggpubr)


# parse arguments
parser = ArgumentParser(description='Compare transformed count on relative vs absolute data.')
parser$add_argument('--abs', type='character', help="log2 absolute data")
parser$add_argument('--rel', type='character', nargs='+', help="Transformed count on relative data")
parser$add_argument('--method_zero', type='character', nargs='+', help="Zero handling method names")
parser$add_argument('--method_transf', type='character', nargs='+', help="Transformation method names")
parser$add_argument('--refgene', type='character', nargs='+', help="Reference gene names")
parser$add_argument('--features', type='character', help="List of gene names")
parser$add_argument('--output', type='character', help="Output figure filename")
parser$add_argument('--npoint', type='integer', default=1e4, help="Number of points to plot")
parser = parser$parse_args()


get_mat <- function(filename, nrefgene){
    count = fread(filename)
    count = as.matrix(count)
    if (!is.na(nrefgene)) {count = count[,-nrefgene]}   # alr transformed matrices have one less column. So in order to match the positions, we should remove the corresponding column from the original matrix
    return(count)
}
get_df <- function(abs, rel, method, nrefgene, npoint){
    abs = get_mat(abs, nrefgene)
    rel = get_mat(rel, nrefgene=NA)
    if ( (nrow(abs) != nrow(rel)) || (ncol(abs) != ncol(rel)) ) { stop("abs and rel matrices have different dimensiones") }
    abs = as.vector(abs)
    rel = as.vector(rel)
    if (length(abs) > npoint){
        set.seed(0); pos = sample(c(1:length(abs)), size=npoint)
        abs = abs[pos]
        rel = rel[pos]
    }
    df  = data.frame(
        abs    = abs,
        rel    = rel,
        method = method
    )
    return(df)
}

# organize data
features = fread(parser$features, header=F)$V1
df = data.frame()
for (i in 1:length(parser$rel)){
    if ( is.na(parser$refgene[i]) || parser$refgene[i] == 'NA' ) {
        refgene  = ''
        nrefgene = NA
    } else {
        refgene  = parser$refgene[i]
        nrefgene = which(features == refgene)
    }
    method = paste0(parser$method_zero[i], ' - ', parser$method_transf[i], ' - ', parser$refgene[i])
    tmp = get_df(parser$abs, parser$rel[i], method, nrefgene, parser$npoint)
    df  = rbind(df, tmp)
}

# set figure size
nmethod = length(unique(df$method))
if (nmethod == 1){
    width = height = 4
    nrow = ncol = 1
}else{
    if (nmethod == 2){
        nrow   = 1
        ncol   = 2
        width  = 8
        height = 4
    }else if(nmethod == 3){
        nrow   = 1
        ncol   = 3
        width  = 12
        height = 4
    }else if(nmethod >= 4){
        nrow   = ceiling(nmethod / 4)
        ncol   = 4
        width  = 12
        height = 3 * nrow
    }
}

# plot abs vs relative
g = ggplot(df, aes(x=abs, y=rel)) +
    facet_wrap(~method, scales="free", nrow=nrow, ncol=ncol) +
    geom_point(alpha=.1, size=.5) + 
    geom_abline(intercept = 0, slope = 1, linetype = 2, color="gray46") +
    xlab("On absolute data") +
    ylab("On relative data") +
    theme(strip.text = element_text(size = 12))
ggsave(
    parser$output, 
    plot   = g, 
    width  = width,
    height = height
)