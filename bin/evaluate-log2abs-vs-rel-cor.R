#!/usr/bin/env Rscript

library(argparse)
library(data.table)

# parse arguments
parser = ArgumentParser(description='Calculate a metric evaluating the coefficients on absolute vs relative data.')
parser$add_argument('--abs', type='character', help="Coefficients on log2 absolute data")
parser$add_argument('--rel', type='character', help="Coefficients on transformed relative data")
parser$add_argument('--features', type='character', help="List of gene names")
parser$add_argument('--out', type='character', help="Output filename containing the metric")
parser$add_argument('--method', type='character', help="Metric method. Choices: [pearson, spearman, kendall]")
parser = parser$parse_args()

parse_file <- function(filename){
    l = c()
    base    = sub('\\.csv.gz', '', filename)
    base    = strsplit(base, split='_')[[1]]
    l['transf']  = base[6]
    l['refgene'] = base[7]
    l['cor']     = base[8]
    l['file']    = filename
    return(l)
}

get_mat <- function(filename, nrefgene){
    count = fread(filename)
    count = as.matrix(count)
    if (!is.na(nrefgene)) { count = count[-nrefgene,-nrefgene] }   # alr transformed matrices have one less column. So in order to match the positions, we should remove the corresponding column from the original matrix
    return(count)
}

compute.cor <- function(abs, rel, method){
    res = cor.test(abs, rel, method=method)
    out = data.frame(method = method, coef = as.numeric(res$estimate), pvalue = res$p.value, stat = as.numeric(res$statistic), test = names(res$statistic))
    return(list(res, out))
}

# parse file
abs_l    = parse_file(parser$abs)
rel_l    = parse_file(parser$rel)
features = fread(parser$features, header=F)$V1
refgene  = rel_l['refgene']
if(is.na(refgene) || refgene == "NA"){ 
    nrefgene=NA 
} else { 
    nrefgene = which( toupper(features) == toupper(refgene) )
    if (length(nrefgene) == 0) stop("wrong reference gene")
}

# read data
message("reading input data")
abs = get_mat(abs_l['file'], nrefgene)
rel = get_mat(rel_l['file'], NA)
if ( (nrow(abs) != nrow(rel)) || (ncol(abs) != ncol(rel)) ) { stop("abs and rel matrices have different dimensiones") }
abs = abs[lower.tri(abs)]
rel = rel[lower.tri(rel)]

# calculate cor
if (parser$method %in% c('pearson', 'spearman', 'kendall')){
    message("calculating pairwise [", parser$method, "] correlation as evaluation method")
    res = compute.cor(abs, rel, parser$method)
}else{
    stop("method not supported")
}
write.table(res[[2]], file=parser$out, quote=F, sep = ',', row.names=F)
message('finished')