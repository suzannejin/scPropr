#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(ggplot2)
library(ggpubr)


# parse arguments
parser = ArgumentParser(description='Compare transformed count on relative vs absolute data.')
parser$add_argument('--abs', type='character', nargs='+', help="log2 absolute data")
parser$add_argument('--rel', type='character', nargs='+', help="Transformed count on relative data")
parser$add_argument('--features', type='character', help="List of gene names")
parser$add_argument('--output', type='character', help="Output figure filename")
parser$add_argument('--npoint', type='integer', default=1e4, help="Number of points to plot")
parser = parser$parse_args()

# functions
parse_files <- function(files){
    l = list(
        name   = c(),
        transf = c(),
        ref    = c(),
        cor    = c(),
        file   = c()
    )
    for (i in 1:length(files)){
        file    = files[i]
        base    = sub('\\.csv.gz', '', file)
        base    = strsplit(base, split='_')[[1]]
        transf  = base[6]
        refgene = base[7]
        cor     = base[8]
        if (refgene == "NA"){ transf = transf; refgene = NA } else { transf = paste0(transf, ' - ', refgene) }
        l$name[i]   = paste0(transf, ' - ', cor)
        l$transf[i] = transf
        l$ref[i]    = refgene
        l$cor[i]    = cor
        l$file[i]   = file
    }
    return(l)
}
get_mat <- function(filename, nrefgene){
    count = fread(filename)
    count = as.matrix(count)
    if (!is.na(nrefgene)) { count = count[-nrefgene,-nrefgene] }   # alr transformed matrices have one less column. So in order to match the positions, we should remove the corresponding column from the original matrix
    return(count)
}

# parse files
abs_l = parse_files(parser$abs)
rel_l = parse_files(parser$rel)

# organize data
features = fread(parser$features, header=F)$V1
df = data.frame()
for (i in 1:length(parser$rel)){
    if ( !is.na(rel_l$ref[i]) ){ nrefgene = which(features == rel_l$ref[i]) } else { nrefgene=NA }
    xpos = which( abs_l$name == paste0('log2 - ', rel_l$cor[i]) )
    x    = get_mat(abs_l$file[xpos], nrefgene)
    y    = get_mat(rel_l$file[i], NA)
    if ( (nrow(x) != nrow(y)) || (ncol(x) != ncol(y)) ) { stop("abs and rel matrices have different dimensiones") }
    x = x[lower.tri(x)]
    y = y[lower.tri(y)]
    set.seed(0); pos = sample(c(1:length(x)), size=parser$npoint)
    x = x[pos]
    y = y[pos]
    df2  = data.frame(
        abs    = x,
        rel    = y,
        transf = rel_l$transf[i],
        cor    = rel_l$cor[i]
    )
    df  = rbind(df, df2)
}

# set figure size
ncol   = length( unique(df$transf) )
nrow   = length( unique(df$cor) )
width  = 3 * ncol
height = 3 * nrow

# plot abs vs relative
g = ggplot(df, aes(x=abs, y=rel)) +
    facet_wrap(~cor + transf, scales="free", nrow=nrow, ncol=ncol) +
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