#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(ggplot2)
library(ggpubr)


# parse arguments
parser = ArgumentParser(description='Compare transformed count on relative vs absolute data.')
parser$add_argument('--abs', type='character', nargs='+', help="Transformed count on absolute data")
parser$add_argument('--rel', type='character', nargs='+', help="Transformed count on relative data")
parser$add_argument('--method_zero', type='character', nargs='+', help="Zero handling method names")
parser$add_argument('--method_transf', type='character', nargs='+', help="Transformation method names")
parser$add_argument('--refgene', type='character', nargs='+', help="Reference gene names")
parser$add_argument('--method_cor', type='character', help="Correlation method name")
parser$add_argument('--output', type='character', help="Output figure filename")
parser$add_argument('--npoint', type='integer', default=1e4, help="Number of points to plot")
parser = parser$parse_args()


# Organize data frame
df = data.frame(matrix(ncol=4, nrow=0))
colnames(df) = c("abs", "rel", "transf", "cor")
for (i in 1:length(parser$abs)){
    # read x
    x = parser$abs[i]
    x = as.matrix(fread(x))
    x = x[lower.tri(x)]
    
    # read y
    y = parser$rel[i]
    y = as.matrix(fread(y))
    y = y[lower.tri(y)]

    if( length(x) != length(y) ) stop("length abs != length rel")

    # random sample
    set.seed(0); pos = sample(c(1:length(x)), parser$npoint)
    x = x[pos]
    y = y[pos]
    
    # update data frame
    transf = paste0(parser$method_zero[i], " - ", parser$method_transf[i], " - ", parser$refgene[i])
    df2 = data.frame(
        "abs"    = x, 
        "rel"    = y, 
        "transf" = rep(transf, length(x)),
        "cor"    = rep(parser$method_cor, length(x))
    )
    df = rbind(df, df2)
}

# set figure size
ncol   = length( unique(df$transf) )
nrow   = 1
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