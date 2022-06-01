#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(propr)

# parse arguments
parser = ArgumentParser(description='Calculate a metric evaluating the coefficients on absolute vs relative data.')
parser$add_argument('--ori', type='character', help="Original (raw) absolute count data")
parser$add_argument('--abs', type='character', help="Coefficients on log2 absolute data")
parser$add_argument('--rel', type='character', help="Coefficients on transformed relative data")
parser$add_argument('--features', type='character', help="List of gene names")
parser$add_argument('--out', type='character', help="Output filename containing the metric")
parser$add_argument('--filter', type='integer', default=20, help="Keep only the pairs whose dropout <= filtering threshold (in percentage). Default = 20")
parser$add_argument('--method', type='character', help="Metric method. Choices: [pearson, spearman, kendall, rho]")
parser = parser$parse_args()

filter = parser$filter / 100

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
    if (method %in% c('pearson','spearman','kendall')){
        res = cor.test(abs, rel, method=method)
        out = data.frame(method = method, coef = as.numeric(res$estimate), pvalue = res$p.value, stat = as.numeric(res$statistic), test = names(res$statistic))
    }else if (method == 'rho'){
        df  = data.frame(abs=abs, rel=rel)
        pro = propr(df, metric=method, ivar=NA, p=0)
        out = data.frame(method = method, coef = pro@matrix[1,2], pvalue=NA, stat=NA, test=NA)
    }
    return(out)
}

# parse file
message('parsing data')
abs_l    = parse_file(parser$abs)
rel_l    = parse_file(parser$rel)
features = fread(parser$features, header=F)$V1
refgene  = rel_l['refgene']
if(is.na(refgene) || refgene == "NA"){ 
    nrefgene=NA 
} else { 
    nrefgene = which( toupper(features) == toupper(refgene) )
    if (length(nrefgene) == 0) stop("wrong reference gene")
    features = features[-nrefgene]
}

# read data
message("reading input data")
ori = get_mat(parser$ori, nrefgene)
abs = get_mat(abs_l['file'], nrefgene)
rel = get_mat(rel_l['file'], NA)
if ( (nrow(abs) != nrow(rel)) || (ncol(abs) != ncol(rel)) ) stop("abs and rel matrices have different dimensiones") 
if ( ncol(abs) != ncol(ori) ) stop("abs and ori matrices have different number of genes") 
if ( ncol(abs) != length(features) ) stop("ncol(abs) != length(features)")

# get dropout
message('filtering genes that have dropout > ', filter)
dropout = colMeans(ori == 0)
pos = which(dropout <= filter)
abs = abs[pos, pos]
rel = rel[pos, pos]
abs = abs[lower.tri(abs)]
rel = rel[lower.tri(rel)]
features = features[pos]

# calculate cor
if (parser$method %in% c('pearson', 'spearman', 'kendall', 'rho')){
    message("calculating pairwise [", parser$method, "] correlation as evaluation method")
    out = compute.cor(abs, rel, parser$method)
}else{
    stop("method not supported")
}
write.table(out, file=parser$out, quote=F, sep = ',', row.names=F)
message('finished')