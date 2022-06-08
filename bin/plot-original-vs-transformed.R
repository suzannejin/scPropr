#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(ggplot2)
library(ggpubr)

# parse arguments
parser = ArgumentParser(description='Compare original vs transformed counts.')
parser$add_argument('--original', type='character', help="Original count data")
parser$add_argument('--transformed', type='character', nargs='+', help="Transformed count data(s)")
parser$add_argument('--method_zero', type='character', nargs='+', help="Zero handling method(s) used")
parser$add_argument('--method_transf', type='character', nargs='+', help="Transformation method(s) used")
parser$add_argument('--refgene', type='character', nargs='+', help="Reference gene(s) used")
parser$add_argument('--features', type='character', help="List of gene names")
parser$add_argument('--outdir', type='character', help="Output directory")
parser$add_argument('--prefix', type='character', help="Output filename prefix")
parser$add_argument('--npoint', type='integer', default=1e4, help="Number of points to plot")
parser$add_argument('--xlab', type='character', default='original')
parser$add_argument('--ylab', type='character', default='transformed')
parser = parser$parse_args()

if ( length(parser$transformed) != length(parser$method_zero) || length(parser$transformed) != length(parser$method_transf) || length(parser$method_transf) != length(parser$refgene) ) stop("please make sure you give the correct inputs")

get_mat <- function(filename, nrefgene){
    count = fread(filename)
    count = as.matrix(count)
    if (!is.na(nrefgene)) {count = count[,-nrefgene]}   # alr transformed matrices have one less column. So in order to match the positions, we should remove the corresponding column from the original matrix
    return(count)
}
get_df <- function(original, transformed, method, nrefgene, npoint){
    original    = get_mat(original, nrefgene)
    transformed = get_mat(transformed, nrefgene=NA)
    if ( (nrow(original) != nrow(transformed)) || (ncol(original) != ncol(transformed)) ) { stop("original and transformed matrices have different dimensiones") }
    original    = as.vector(original)
    transformed = as.vector(transformed)
    if (length(original) > npoint){
        set.seed(0); pos = sample(c(1:length(original)), size=npoint)
        original    = original[pos]
        transformed = transformed[pos]
    }
    df = data.frame(
        original    = original,
        transformed = transformed,
        method      = method
    )
    return(df)
}

# organize data
message('--organize data')
features = fread(parser$features, header=F)$V1
df = data.frame()
for (i in 1:length(parser$transformed)){
    if (is.na(parser$refgene[i]) || parser$refgene[i] == 'NA') {
        refgene  = ''
        nrefgene = NA
        method   = paste0( parser$method_zero[i], '_', parser$method_transf[i], '_NA' )
    } else {
        refgene  = parser$refgene[i]
        nrefgene = which(features == refgene)
        method   = paste0( parser$method_zero[i], '_', parser$method_transf[i], '_', refgene )
    }
    tmp = get_df(parser$original, parser$transformed[i], method, nrefgene, parser$npoint)
    df  = rbind(df, tmp)
}

# set figure size
nmethod = length(parser$transformed)
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

# plot figure
message('--plot figure')
g = ggplot(df, aes(x=original, y=transformed)) +
        facet_wrap(~method, nrow=nrow, ncol=ncol, scales='free') +
        geom_point(size=1, alpha=.4) +
        geom_abline(intercept=1, slope=1, linetype='dashed') +
        theme(strip.text.x = element_text(size = 14)) + 
        xlab(parser$xlab) +
        ylab(parser$ylab)
ggsave(
    paste0(parser$outdir, '/', parser$prefix, '.png'), 
    plot   = g, 
    width  = width,
    height = height
)
# logx
message('--plot figure with log10(x)')
g = ggplot(df, aes(x=original, y=transformed)) +
        facet_wrap(~method, nrow=nrow, ncol=ncol, scales='free') +
        geom_point(size=1, alpha=.4) +
        geom_abline(intercept=1, slope=1, linetype='dashed') +
        scale_x_continuous(trans='log10') +
        theme(strip.text.x = element_text(size = 14)) +
        xlab(parser$xlab) +
        ylab(parser$ylab)
ggsave(
    paste0(parser$outdir, '/', parser$prefix, '-logx.png'), 
    plot   = g, 
    width  = width,
    height = height
)
# logy
message('--plot figure with log10(y)')
g = ggplot(df, aes(x=original, y=transformed)) +
        facet_wrap(~method, nrow=nrow, ncol=ncol, scales='free') +
        geom_point(size=1, alpha=.4) +
        geom_abline(intercept=1, slope=1, linetype='dashed') +
        scale_y_continuous(trans='log10') +
        theme(strip.text.x = element_text(size = 14)) +
        xlab(parser$xlab) +
        ylab(parser$ylab)
ggsave(
    paste0(parser$outdir, '/', parser$prefix, '-logy.png'), 
    plot   = g, 
    width  = width,
    height = height
)
# logxy
message('--plot figure with log10(x) and log10(y)')
g = ggplot(df, aes(x=original, y=transformed)) +
        facet_wrap(~method, nrow=nrow, ncol=ncol, scales='free') +
        geom_point(size=1, alpha=.4) +
        geom_abline(intercept=1, slope=1, linetype='dashed') +
        scale_x_continuous(trans='log10') +
        scale_y_continuous(trans='log10') +
        theme(strip.text.x = element_text(size = 14)) +
        xlab(parser$xlab) +
        ylab(parser$ylab)
ggsave(
    paste0(parser$outdir, '/', parser$prefix, '-logxy.png'), 
    plot   = g, 
    width  = width,
    height = height
)