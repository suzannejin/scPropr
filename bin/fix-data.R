#!/usr/bin/env Rscript

library(data.table)
library(tibble)

# parse arguments -----------------------------------------------------------------------
parser = ArgumentParser(description='Organize data dimensions. Eg when a reference gene is used in ALR computation, the corresponding column is removed. Here we fill the column with NA to recover the dimensions.')
parser$add_argument('--count', type='character', help="Input count data")
parser$add_argument('--coeff', type='character', help="Correlation coefficient matrix")
parser$add_argument('-f', '--features', type='character', help="File containing the feature (gene) names")
parser$add_argument('-r', '--refgene', type='character', nargs='+', help="Reference gene(s)")
parser$add_argument('-o', '--output', type='character', help="Output filename")
parser = parser$parse_args()


# read input data -----------------------------------------------------------------------
if (!is.null(parser$count)){
    input = fread(parser$count)
    type  = "count"
}else if (!is.null(parser$coeff)){
    input = fread(parser$coeff)
    type  = "coeff"
}else{
    stop("plase provide the input data")
}
features = fread(parser$features, header=F)$V1


# fill data dimension, if needed --------------------------------------------------------
if ( length(features) == nrow(input) ){
    out = input
}else if ( length(features) > nrow(input) ){
    pos = which(features %in% parser$refgene)
    out = tibble::add_column(input, .before=pos, ref=NA)
    if (type == "coeff"){
        out = tibble::add_row(out, .before=refgene, ref=NA)
    }
}else { 
    stop("feature list provided has smaller length")
}


# output --------------------------------------------------------------------------------

# check new dimensions
if ( length(features) != nrow(out) ) stop("something went wrong. Still wrong dimensions.")

# save output
fwrite(out, file=parser$output, quote=F, sep=",", row.names=F, col.names=F, compress="gzip")