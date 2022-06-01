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
parser$add_argument('--output', type='character', help="Output figure filename")
parser$add_argument('--npoint', type='integer', default=1e4, help="Number of points to plot")
parser = parser$parse_args()


# Organize data frame
df = data.frame(matrix(ncol=3, nrow=0))
colnames(df) = c("x", "y", "z")
methods = c()
for (i in 1:length(parser$method_transf)){
    # read x
    x = parser$abs[i]
    x = as.matrix(fread(x))
    x = as.vector(x)
    
    # read y
    y = parser$rel[i]
    y = as.matrix(fread(y))
    y = as.vector(y)

    if( length(x) != length(y) ) stop("length abs != length rel")

    # random sample
    set.seed(0); pos = sample(c(1:length(x)), parser$npoint)
    x = x[pos]
    y = y[pos]
    
    # update data frame
    method  = paste0(parser$method_zero[i], ' - ', parser$method_transf[i], ' - ', parser$refgene[i])
    methods = c(methods, method)
    df2 = data.frame("x"=x, "y"=y, "z"=rep(method, length(x)))
    df = rbind(df, df2)
}
df$z = factor(df$z, levels=methods)

# set figure size
nmethod = length(methods)
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
g = ggplot(df, aes(x=x, y=y)) +
    facet_wrap(~z, scales="free", nrow=nrow, ncol=ncol) +
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