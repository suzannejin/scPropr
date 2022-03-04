#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(ggplot2)
library(ggpubr)

# parse arguments
parser <- ArgumentParser(description='Plot library size (geometric mean) histogram')
parser$add_argument('count', type="character", help="Input count file (cells x genes)")
parser$add_argument('name' , type="character", help="Dataset name")
parser$add_argument('outdir',type="character", help="Output directory")
parser$add_argument('--libsize', action='store_true', help="Compute lib size histogram")
parser$add_argument('--geomean', action='store_true', help="Compute geometric mean histogram")
parser$add_argument('--mean', action='store_true', help="Compute mean")
parser$add_argument('--var', action='store_true', help="Compute variance")
parser$add_argument('--bygene', action='store_true', help="Compute per gene. Default compute per cell")
parser = parser$parse_args()

# functions
calculate_lib_size <- function(count){
    libsize = rowSums(count)
    return(libsize)
}
calculate_geom_mean <-  function(count){
    gm = apply( count, 1, function(x) exp(mean(log(x[x!=0]))) )
    return(gm)
}
calculate_mean <- function(count, by){
    val = apply( count, by, function(x) mean(x) )
    return(val)
}
calculate_var <- function(count, by){
    val = apply( count, by, function(x) var(x) )
    return(val)
}


# compute per cell (default) or per gene
by = 1
bypref = 'cell'
if (parser$bygene){
    by = 2
    bypref = 'gene'
}

# read count data and compute
count = fread(parser$count)
if (parser$libsize){
    x    = calculate_lib_size(count)
    xlab = 'library size'
    val  = 'libsize'
}else if(parser$geomean){
    x    = calculate_geom_mean(count)
    xlab = 'geometric mean'
    val  = 'geomean'
}else if(parser$mean){
    x    = calculate_mean(count, by)
    xlab = paste0(bypref, ' expression mean')
    val  = paste0(bypref, 'mean')
}else if(parser$var){
    x    = calculate_var(count, by)
    xlab = paste0(bypref, ' expression variance')
    val  = paste0(bypref, 'var')
}

# plot histogram
# png(filename=paste0(parser$outdir,'/histogram-', val, '.png'), width=5, height=5, units="in", res=100)
# hist(x, main=parser$name, xlab=xlab)
# png(filename=paste0(parser$outdir,'/histogram-', val, '-log.png'), width=5, height=5, units="in", res=100)
# hist(x, main=parser$name, xlab=xlab, log='x')
df = data.frame(x)
g1 = ggplot(df, aes(x)) + 
        geom_histogram(colour="black", fill="white") +
        xlab(xlab) 
g2 = ggplot(df, aes(x)) +               # Histogram with log10 axis
        geom_histogram(colour="black", fill="white") +
        scale_x_continuous(
            trans='log10', 
            breaks = scales::trans_breaks("log10", function(x) 10^x),
            labels=scales::trans_format("log10", scales::math_format(10^.x))) +
        xlab(xlab) 
g = ggarrange(g1, g2, ncol=2, nrow=1)
ggsave(g, filename=paste0(parser$outdir, '/histogram-', val, '.png'), width = 9, height = 4)