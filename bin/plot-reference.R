#!/usr/bin/env Rscript

library(argparse)
library(data.table)

# parse arguments
parser <- ArgumentParser(description='Plot true reference vs norm factor')
# parser$add_argument('--abs_raw', type="character", help="Absolute raw count data")
parser$add_argument('--abs', type="character", help="Absolute count data")
parser$add_argument('--rel', type="character", help="Relative count data")
parser$add_argument('--features', type="character", help="File with list of col names (features)")
parser$add_argument('--ref', type="character", help="Reference gene id. Or normalisation name, if normalisation factors provided")
parser$add_argument('--absnorm', type="character", help="File with normalisation factors on absolute data")
parser$add_argument('--relnorm', type="character", help="File with normalisation factors on relative data")
parser$add_argument('--out', type="character", help="Output figure filename")
parser$add_argument('--libsize', action="store_true", help="If true, just use as reference the libsize of the absolute data. Otherwise, the real factor, 1/libsize * mean(libsize)")
parser = parser$parse_args()

# read input
message('reading input')
abs = as.matrix( fread(parser$abs) )
rel = as.matrix( fread(parser$rel) )
features = fread(parser$features, header=F)$V1
# if ( !is.null(parser$abs_raw) ){  # if raw is provided, it means that the other abs and rel count data were already processed, no need of zero replacing
#     abs_raw = fread(parser$abs_raw)
#     abs_raw = log2(abs_raw + 1)
#     abs = log2(abs)
#     rel = log2(rel)
# }else{
#     abs = log2(abs + 1)
#     rel = log2(rel + 1)
#     abs_raw = abs
# }
# ref = 1 / rowSums(abs) * mean(rowSums(abs))
abs = log2(abs + 1)
rel = log2(rel + 1)
abs_raw = abs

# true reference
if (parser$libsize){
    xlab = 'Abs libsize'
    ref  = rowSums(abs_raw)
}else{
    xlab = 'True reference'
    ref = 1 / rowSums(abs_raw) * mean(rowSums(abs_raw)) 
}

# get norm factor
message('reading norm factor')
if ( !is.null(parser$absnorm) ){
    name  = parser$ref
    fac_abs = fread(parser$absnorm, header=F)$V1
    fac_rel = fread(parser$relnorm, header=F)$V1
} else {
    name = parser$ref
    n    = which(features == parser$ref)
    fac_abs = abs[,n]
    fac_rel = rel[,n]
}

# plot figure
message('plotting')
png(filename=parser$out, width=8, height=4, units="in", res=100)
par(mfrow=c(1,2))
plot(
    x=ref, 
    y=fac_abs,
    xlab=xlab,
    ylab=paste0('abs_', name)
)
plot(
    x=ref, 
    y=fac_rel,
    xlab=xlab,
    ylab=paste0('rel_', name)
)
message('finished')