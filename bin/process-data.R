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
parser$add_argument('--method_zero', type='character', help="Zero handling method. Choices = [zcompositions, one, min]")
parser$add_argument('--method_transf', type='character', help="Transformation or normalization method. Choices = [log2, clr, alr, tmm, scran]")
parser$add_argument('--refgene', type='character', help="Reference gene index. This is required to compute alr.")
parser = parser$parse_args()

# check arguments
if (!is.null(parser$method_zero) && !( parser$method_zero %in% c("zcompositions", "one", "min", "pseudocount") )){
    stop("wrong zero replacement method - ", parser$method_zero)
}
if (!is.null(parser$method_transf) && !( parser$method_transf %in% c("log2", "clr", "alr", "clr2", "alr2", "tmmpre","tmmpos", "tmmcpm", "tmmscaled", "scran") )){
    stop("wrong transformation or normalization method - ", parser$method_transf)
}
pseudo.count = 1

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
print(paste0("dropout: ", dropout))
dim(count)
class(count)

# functions -------------------------------------------------------------------

replace_zero <- function(count, method=c("zcompositions", "one", "min", "pseudocount"), pseudo.count=1){

    method = match.arg(method)
    if (!any(count==0, na.rm=T)) return(count)
    message("replacing zeros with method ", method)

    # replace zero
    if (method == "zcompositions"){
        require(zCompositions)
        pars = as.character(formals(cmultRepl)$output)
        if ("p-counts" %in% pars){
            par = "p-counts"
        }else if ("counts" %in% pars){
            par = "counts"
        }else{
            stop("wrong output parameter for zCompositions")
        }
        count = cmultRepl(count, method="CZM", label=0, output=par)

    }else if (method == "one"){
        count[count == 0] = 1

    }else if (method == "min"){
        count = count
        zeros = count == 0
        count[zeros] = min(count[!zeros])
    }else if (method == "pseudocount"){
        count = count + pseudo.count
    }

    # TODO add scImpute, etc zero handling methods for single cell

    # check zeros and negative values
    if(any(count<=0, na.rm=T)) stop("counts matrix contain zero and/or negative values")

    return(count)
}

calculate_geom_mean <-  function(count){
    gm = apply( count, 1, function(x) exp(mean(log(x[x!=0]))) )
    return(gm)
}


# transform/normalize data --------------------------------------------------------------

if (parser$method_transf == 'log2'){
    message("log2-transforming data")
    count = replace_zero(count, parser$method_zero, pseudo.count=pseudo.count)
    count = log2(count)

}else if (parser$method_transf == 'alr'){
    message("alr-transforming data")
    count = replace_zero(count, parser$method_zero, pseudo.count=pseudo.count)
    size.factor = count[,nrefgene]
    count = count[,-nrefgene] / size.factor 
    count = log(count)

}else if (parser$method_transf == 'alr2'){
    message("alr2-transforming data")
    size.factor = replace_zero(count, parser$method_zero, pseudo.count=pseudo.count)[,nrefgene]
    count = count[,-nrefgene] / size.factor 
    count = replace_zero(count, parser$method_zero, pseudo.count=pseudo.count)
    count = log(count)

}else if (parser$method_transf == 'clr'){
    message("clr-transforming data")
    count = replace_zero(count, parser$method_zero, pseudo.count=pseudo.count)
    size.factor = calculate_geom_mean(count)
    count = count / size.factor
    count = log(count)

}else if (parser$method_transf == 'clr2'){
    message("clr2-transforming data")
    size.factor = calculate_geom_mean( replace_zero(count, parser$method_zero, pseudo.count=pseudo.count) )
    count = count / size.factor
    count = replace_zero(count, parser$method_zero, pseudo.count=pseudo.count)
    count = log(count)

}else if (parser$method_transf == 'tmmcpm' && parser$method_zero == 'pseudocount'){
    message("normalizing data with method tmmcpm")
    require(edgeR)
    count       = t(count) # it requires row=genes, col=cells
    y           = DGEList(counts=count)   
    y           = calcNormFactors(y, method="TMM")
    size.factor = y$samples$norm.factors
    eff.libsize = y$samples$lib.size * size.factor
    adj.prior   = pseudo.count * eff.libsize / mean(eff.libsize)
    adj.libsize = eff.libsize + 2 * adj.prior
    count       = (t(count) + adj.prior) / adj.libsize * 1e6  ## note that adj.prior has length libsize, therefore you should be careful when transposing and making matrix calculations
    count       = log2(count)

}else if (parser$method_transf == 'tmmpre'){
    message("normalizing data with method tmmpre")
    require(edgeR)
    count       = t(count) # it requires row=genes, col=cells
    y           = DGEList(counts=count)   
    y           = calcNormFactors(y, method="TMM")
    size.factor = y$samples$norm.factors
    effective.libsize = y$samples$lib.size * size.factor
    count       = replace_zero(t(count), parser$method_zero, pseudo.count=pseudo.count)   
    count       = count / effective.libsize * 1e6
    count       = log2(count)  

}else if (parser$method_transf == 'tmmpos'){
    message("normalizing data with method tmmpos")
    require(edgeR)
    count       = t(count) # it requires row=genes, col=cells
    y           = DGEList(counts=count)   
    y           = calcNormFactors(y, method="TMM")
    size.factor = y$samples$norm.factors
    effective.libsize = y$samples$lib.size * size.factor   
    count       = t(count) / effective.libsize * 1e6
    count       = replace_zero(count, parser$method_zero, pseudo.count=pseudo.count)
    count       = log2(count)  

}else if (parser$method_transf == 'tmmscaled'){
    message("normalizing data with method tmmscaled")
    require(edgeR)
    count       = t(count) # it requires row=genes, col=cells
    y           = DGEList(counts=count)   
    y           = calcNormFactors(y, method="TMM")
    size.factor = y$samples$norm.factors
    count       = t(count) / size.factor
    count       = replace_zero(count, parser$method_zero, pseudo.count=pseudo.count)
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
    count       = replace_zero(count, parser$method_zero, pseudo.count=pseudo.count)
    count       = log2(count)

}else{
    stop("make sure you introduced the correct processing methods")
}


# save output --------------------------------------------------------------
message("saving output files")
fwrite(count, file=parser$output, quote=F, sep=",", row.names=F, col.names=F, compress="gzip")
if (parser$method_transf != 'log2') {
    fwrite(data.frame(size.factor), file=parser$output2, quote=F, sep=",", row.names=F, col.names=F)
} else {
    file.create(parser$output2)
}