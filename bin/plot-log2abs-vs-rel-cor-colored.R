#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(ggplot2)
library(ggpubr)


# parse arguments
parser = ArgumentParser(description='Compare coefficients computed on relative vs absolute data.')
parser$add_argument('--ori', type='character', help="Original (raw) absolute count data")
parser$add_argument('--abs', type='character', nargs='+', help="Coefficients computed on absolute data")
parser$add_argument('--rel', type='character', nargs='+', help="Coefficients computed on relative data")
parser$add_argument('--features', type='character', help="List of gene names")
parser$add_argument('--color', type='character', default='dropout', help="Color by dropout or var. Default=dropout")
parser$add_argument('--filter', type='double', default=0.2, help="Keep only the pairs whose dropout <= filtering threshold. Default = 0.2")
parser$add_argument('--outdir', type='character', help="Output directory")
parser$add_argument('--npoint', type='integer', default=1e4, help="Number of points to plot")
parser$add_argument('--seed', type='integer', default=0, help="Random seed. Default=0")
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
message("parsing files")
abs_l = parse_files(parser$abs)
rel_l = parse_files(parser$rel)

# get dropout and variance 
message('computing dropout and variance for each gene')
features      = fread(parser$features, header=F)$V1
ori           = fread(parser$ori)
tmp           = ori / rowSums(ori) * mean(rowSums(ori))
dropout       = colMeans(ori == 0)
# var_genes_abs = apply(ori, 2, var)
# var_genes_rel = apply(tmp, 2, var)
rm(ori)
rm(tmp)

# organize data
message("organizing data")
df = data.frame()
df_filtered = data.frame()
for (i in 1:length(parser$rel)){
    if ( !is.na(rel_l$ref[i]) ){ 
        nrefgene  = which(features == rel_l$ref[i]) 
        features2 = features[-nrefgene]
        dropout2  = dropout[-nrefgene]
        # var_genes_abs2 = var_genes_abs[-nrefgene]
        # var_genes_rel2 = var_genes_rel[-nrefgene]
    } else { 
        nrefgene  = NA 
        features2 = features
        dropout2  = dropout
        # var_genes_abs2 = var_genes_abs
        # var_genes_rel2 = var_genes_rel
    }
    xpos = which( abs_l$name == paste0('log2 - ', rel_l$cor[i]) )
    abs    = get_mat(abs_l$file[xpos], nrefgene)
    rel    = get_mat(rel_l$file[i], NA)
    if ( (nrow(abs) != nrow(rel)) || (ncol(abs) != ncol(rel)) ) stop("abs and rel matrices have different dimensiones") 
    if ( length(features2) != ncol(abs) ) stop("length(features) != ncol(abs)")
    ind  = which( upper.tri(abs, diag=F) , arr.ind = TRUE )
    df2  = data.frame(
        transf    = rel_l$transf[i],
        cor       = rel_l$cor[i],
        i         = features2[ind[,1]],
        j         = features2[ind[,2]],
        abs       = abs[ind],
        rel       = rel[ind],
        dropout_i = dropout2[ind[,1]],
        dropout_j = dropout2[ind[,2]]
        # var_abs_i = var_genes_abs2[ind[,1]],
        # var_abs_j = var_genes_abs2[ind[,2]],
        # var_rel_i = var_genes_rel2[ind[,1]],
        # var_rel_j = var_genes_rel2[ind[,2]]
    )
    
    pos = which( df2[,'dropout_i'] <= parser$filter & df2[,'dropout_j'] <= parser$filter )
    if (length(pos) == 0) stop('please make the filtering less stringent')
    df2_filtered = df2[pos,]
    set.seed(parser$seed); pos = sample(c(1:nrow(df2)), size=parser$npoint)
    df2 = df2[pos,]
    df  = rbind(df, df2)
    set.seed(parser$seed); pos = sample(c(1:nrow(df2_filtered)), size=parser$npoint)
    df2_filtered = df2_filtered[pos,]
    df_filtered  = rbind(df_filtered, df2_filtered)
}
df[,'dropout'] = rowMeans( data.frame(df[,'dropout_i'], df[,'dropout_j']) )
df_filtered[,'dropout'] = rowMeans( data.frame(df_filtered[,'dropout_i'], df_filtered[,'dropout_j']) )
# df[,'var_abs'] = rowMeans( data.frame(df[,'var_abs_i'], df[,'var_abs_j']) )
# df[,'var_rel'] = rowMeans( data.frame(df[,'var_rel_i'], df[,'var_rel_j']) )
dim(df)
head(df)
dim(df_filtered)
head(df_filtered)

# release memory
rm(abs_l)
rm(rel_l)
rm(abs)
rm(rel)
rm(features)
rm(dropout)
# rm(var_genes_abs)
# rm(var_genes_rel)

# set figure size
ncol    = length( unique(df$transf) )
nrow    = length( unique(df$cor) )
width   = 3 * ncol
height  = 3 * nrow
ncol2   = length( unique(df_filtered$transf) )
nrow2   = length( unique(df_filtered$cor) )
width2  = 3 * ncol2
height2 = 3 * nrow2

# plot abs vs relative
message('plotting figure')
g1 = ggplot(df, aes_string(x='abs', y='rel', color=parser$color)) +
    facet_wrap(~cor + transf, scales="free", nrow=nrow, ncol=ncol) +
    geom_point(alpha=.1, size=.5) + 
    geom_abline(intercept = 0, slope = 1, linetype = 2, color="gray46") +
    scale_colour_viridis_c() +
    xlab("On absolute data") +
    ylab("On relative data") +
    theme(strip.text = element_text(size = 12))
g2 = ggplot(df, aes_string(x='abs', y='rel')) +
    facet_wrap(~cor + transf, scales="free", nrow=nrow, ncol=ncol) +
    geom_point(alpha=.1, size=.5) + 
    geom_abline(intercept = 0, slope = 1, linetype = 2, color="gray46") +
    xlab("On absolute data") +
    ylab("On relative data") +
    theme(strip.text = element_text(size = 12))
g3 = ggplot(df_filtered, aes_string(x='abs', y='rel', color=parser$color)) +
    facet_wrap(~cor + transf, scales="free", nrow=nrow2, ncol=ncol2) +
    geom_point(alpha=.1, size=.5) + 
    geom_abline(intercept = 0, slope = 1, linetype = 2, color="gray46") +
    scale_colour_viridis_c() +
    xlab("On absolute data") +
    ylab("On relative data") +
    theme(strip.text = element_text(size = 12))
g4 = ggplot(df_filtered, aes_string(x='abs', y='rel')) +
    facet_wrap(~cor + transf, scales="free", nrow=nrow2, ncol=ncol2) +
    geom_point(alpha=.1, size=.5) + 
    geom_abline(intercept = 0, slope = 1, linetype = 2, color="gray46") +
    xlab("On absolute data") +
    ylab("On relative data") +
    theme(strip.text = element_text(size = 12))
message('saving colored figure')
ggsave(
    paste0(parser$outdir, '/log2abs-vs-rel-cor-colored.png'), 
    plot   = g1, 
    width  = width,
    height = height
)
message('saving black figure')
ggsave(
    paste0(parser$outdir, '/log2abs-vs-rel-cor-black.png'), 
    plot   = g2, 
    width  = width,
    height = height
)
message('saving colored and filtered figure')
ggsave(
    paste0(parser$outdir, '/log2abs-vs-rel-cor-colored-filtered.png'), 
    plot   = g3, 
    width  = width2,
    height = height2
)
message('saving black and filtered figure')
ggsave(
    paste0(parser$outdir, '/log2abs-vs-rel-cor-black-filtered.png'), 
    plot   = g4, 
    width  = width2,
    height = height2
)
message('finished')