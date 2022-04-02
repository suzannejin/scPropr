#!/usr/bin/env Rscript

library(argparse)
library(data.table)

# parse arguments
parser <- ArgumentParser(description='Plot scatter values vs cell/gene index ')
parser$add_argument('count', type="character", help="Input count file (cells x genes)")
parser$add_argument('name' , type="character", help="Dataset name")
parser$add_argument('outdir',type="character", help="Output directory")
parser$add_argument('--libsize', action='store_true', help="Compute lib size vs index")
parser$add_argument('--geomean', action='store_true', help="Compute geometric mean vs index")
parser$add_argument('--mean', action='store_true', help="Compute arithmetic mean vs index")
parser$add_argument('--var', action='store_true', help="Compute variance vs index")
parser$add_argument('--bygene', action='store_true', help="Use gene index")
parser = parser$parse_args()

# functions
calculate_lib_size <- function(count){
    libsize = rowSums(count)
    return(libsize)
}
calculate_geom_mean <-  function(count, by){
    gm = apply( count, by, function(x) exp(mean(log(x[x!=0]))) )
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

# read count data
count = fread(parser$count)

# compute per cell (default) or per gene
by = 1
bypref = 'cell'
x = c(1:nrow(count))
if (parser$bygene){
    by = 2
    bypref = 'gene'
    x = c(1:ncol(count))
}

# compute values
if (parser$libsize){
    y    = calculate_lib_size(count)
    ylab = 'library size'
    metric = 'libsize'
}else if(parser$geomean){
    y    = calculate_geom_mean(count, by)
    ylab = 'geometric mean'
    metric = 'geomean'
}else if(parser$mean){
    y    = calculate_mean(count, by)
    ylab = 'mean'
    metric  = 'mean'
}else if(parser$var){
    y    = calculate_var(count, by)
    ylab = 'variance'
    metric  = 'var'
}

# plot 
png(filename=paste0(parser$outdir,'/scatter-', metric, '-', bypref, 'index.png'), width=12, height=5, units="in", res=100)
plot(x,y, main=parser$name, xlab=paste0(bypref, ' index'), ylab=ylab)

# plot sorted figure
y = sort(y)
png(filename=paste0(parser$outdir,'/scatter-', metric, '-', bypref, 'index-sorted.png'), width=12, height=5, units="in", res=100)
plot(x,y, main=parser$name, xlab=paste0(bypref, ' sorted by ', ylab), ylab=ylab)