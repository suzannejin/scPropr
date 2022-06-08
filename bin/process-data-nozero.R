#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(propr)

# parse arguments
parser = ArgumentParser(description='Transform and normalize count data. Also save the normalization factors.')
parser$add_argument('-i', '--input', type='character', help="Input count data")
parser$add_argument('-o', '--output', type='character', help="Output count data, if required")
parser$add_argument('-o2', '--output2', type='character', help="Output size factor file")
parser$add_argument('--features', type='character')
parser$add_argument('--method_transf', type='character', help="Transformation or normalization method. Choices = [log2, clr, alr, tmm, scran]")
parser$add_argument('--refgene', type='character', help="Reference gene index. This is required to compute alr.")
parser = parser$parse_args()

# check arguments
if (!is.null(parser$method_transf) && !( parser$method_transf %in% c("log2", "clr", "alr", "tmm", "scran") )){
    stop("wrong transformation or normalization method - ", parser$method_transf)
}

# read and process data --------------------------------------------------------

# read input data
count    = as.matrix( fread(parser$input, header=F) )
features = fread(parser$features, header=F)$V1
if ( length(features) != ncol(count) ) stop("length(features) != ncol(count)")
if (!is.null(parser$refgene)){
    nrefgene  = which( toupper(features) == toupper(parser$refgene) )
    if (length(nrefgene) == 0) stop("wrong reference gene")
}

# count dropout
dropout = round(mean(count==0, na.rm=T), 3)
if (any(count == 0)){
    print(paste0("dropout: ", dropout))
}
dim(count)
class(count)


# functions -------------------------------------------------------------------

calculate_geom_mean <-  function(count){
    gm = apply( count, 1, function(x) exp(mean(log(x[x!=0]))) )
    return(gm)
}


# transform/normalize data --------------------------------------------------------------

if (parser$method_transf == 'log2'){
    message("log2-transforming data")
    count = log2(count)

}else if (parser$method_transf == 'alr'){
    message("alr-transforming data")
    size.factor = count[,nrefgene]
    count  = count[,-nrefgene] / size.factor 
    count  = log(count)

}else if (parser$method_transf == 'clr'){
    message("clr-transforming data")
    size.factor = calculate_geom_mean(count)
    count  = count / size.factor
    count  = log(count)

}else if (parser$method_transf == 'tmm'){
    message("normalizing data with method tmm")
    require(edgeR)
    count       = t(count) # it requires row=genes, col=cells
    y           = DGEList(counts=count)   
    y           = calcNormFactors(y, method="TMM")
    size.factor = y$samples$norm.factors
    eff.libsize = y$samples$lib.size * size.factor
    count       = (t(count)) / eff.libsize * 1e6  
    count       = log2(count)

}else if (parser$method_transf == 'scran'){
    message(" normalizing data with method scran")
    require(scran)
    require(scuttle)
    require(SingleCellExperiment)
    count       = t(count)   # it requires row=genes, col=cells
    sce         = SingleCellExperiment(list(counts=count))
    clusters    = quickCluster(sce)
    sce         = computeSumFactors(sce, clusters=clusters)  
    size.factor = sce@colData@listData$sizeFactor
    count       = t(count) / size.factor
    count       = log2(count)

}else{
    stop("make sure you introduced the correct processing methods")
}


# save output --------------------------------------------------------------
message("saving output files")

# store processed count data
fwrite(count, file=parser$output, quote=F, sep=",", row.names=F, col.names=F, compress="gzip")

# store size factor
if (parser$method_transf != 'log2') {
    fwrite(data.frame(size.factor), file=parser$output2, quote=F, sep=",", row.names=F, col.names=F)
} else {
    file.create(parser$output2)
}