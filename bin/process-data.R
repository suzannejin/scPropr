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
    stop("wrong zero replacement method - ", method_zero)
}
if (!is.null(parser$method_transf) && !( parser$method_transf %in% c("log2", "clr", "alr", "alr2", "clr2", "tmm", "scran","tmm2","tmmcpm") )){
    stop("wrong transformation or normalization method - ", method_transf)
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

normalize_data <- function(count, method=c("tmm", "scran","tmm2","tmmcpm")){
    method = match.arg(method)
    message("normalizing data with method ", method)
    if (method == "tmm"){
        require(edgeR)
        # Note that normalization in edgeR is model-based, and the original read counts are not themselves transformed. This means that users should not transform the read counts in any way
        # before inputing them to edgeR. For example, users should not enter RPKM or FPKM values to edgeR in place of read counts. Such quantities will prevent edgeR from correctly
        # estimating the mean-variance relationship in the data, which is a crucial to the statistical
        # strategies underlying edgeR. Similarly, users should not add artificial values to the counts
        # before inputing them to edgeR

        # For further information check: https://www.biostars.org/p/84087/ and https://www.biostars.org/p/317701/
        count       = t(count) # it requires row=genes, col=cells
        y           = DGEList(counts=count)   
        y           = calcNormFactors(y, method="TMM")
        size.factor = y$samples$norm.factors
        effective.libsize = y$samples$lib.size * size.factor   
        count       = t(count) / effective.libsize * 1e6
    }else if (method == "tmm2"){
        require(edgeR)
        count       = t(count) # it requires row=genes, col=cells
        y           = DGEList(counts=count)   
        y           = calcNormFactors(y, method="TMM")
        size.factor = y$samples$norm.factors
        effective.libsize = y$samples$lib.size * size.factor   
        size.factor = effective.libsize
        count       = t(count) 
    }else if (method == "tmmcpm"){
        require(edgeR)
        count       = t(count)
        y           = DGEList(counts=count)   
        y           = calcNormFactors(y, method="TMM")
        norm        = cpm(y, log=T, prior.count=1)
        count       = t(norm)
        size.factor = y$samples$norm.factors
        effective.libsize = y$samples$lib.size * size.factor   
        size.factor = effective.libsize
    }else if (method == "scran"){
        require(scran)
        require(scuttle)
        require(SingleCellExperiment)
        count       = t(count)   # it requires row=genes, col=cells
        sce         = SingleCellExperiment(list(counts=count))
        clusters    = quickCluster(sce)
        sce         = computeSumFactors(sce, clusters=clusters)  
        size.factor = sce@colData@listData$sizeFactor
        count       = t(count) / size.factor
        ## TODO count = log1p( t( t(count) / centered.size.factor ) ) / log(2)
    }

    # TODO add other popular normalization procedures
    norm = list(count = count, size.factor = size.factor)
    return(norm)
}

calculate_geom_mean <-  function(count){
    gm = apply( count, 1, function(x) exp(mean(log(x[x!=0]))) )
    return(gm)
}


# transform/normalize data --------------------------------------------------------------

if (parser$method_transf %in% c('log2', 'clr', 'alr')){
    if(!is.null(parser$method_zero)) count = replace_zero(count, parser$method_zero, pseudo.coount=1)

    if (parser$method_transf == "log2"){
        message("log-transforming data")
        count = log2(count)

    }else if(parser$method_transf == "clr"){
        message("clr-transforming data")
        # count = propr:::proprCLR(count)
        size.factor = calculate_geom_mean(count)
        count = log(count / size.factor)

    }else if(parser$method_transf == "alr"){
        if (!is.numeric(nrefgene)) stop("Make sure you are giving the correct reference gene index for ALR-transformation")
        message("alr-transforming data")
        size.factor = count[,nrefgene]
        count = propr:::proprALR(count, nrefgene)
    }

}else if(parser$method_transf == 'alr2'){
    count2 = replace_zero(count, parser$method_zero, pseudo.count=1)
    size.factor = count2[,nrefgene]
    count = count[,-nrefgene] / size.factor
    count = replace_zero(count, parser$method_zero, pseudo.count=1)
    count = log(count)

}else if (parser$method_transf == 'clr2'){
    count2 = replace_zero(count, parser$method_zero, pseudo.count=1)
    size.factor = calculate_geom_mean(count2)
    count = count / size.factor
    count = replace_zero(count, parser$method_zero, pseudo.count=1)
    count = log(count)

}else if (parser$method_transf %in% c('tmm','scran')){
    norm  = normalize_data(count, parser$method_transf)
    count = norm$count
    if(!is.null(parser$method_zero)) count = replace_zero(count, parser$method_zero, pseudo.count=1)
    count = log2(count)
    size.factor = norm$size.factor

}else if (parser$method_transf == 'tmm2'){  ## corrected version,  transcript per milion
    norm  = normalize_data(count, 'tmm2')
    size.factor = norm$size.factor
    if(!is.null(parser$method_zero)) count = replace_zero(count, parser$method_zero, pseudo.count=1)
    count = count / size.factor * 1e6
    count = log2(count)

}else if (parser$method_transf == 'tmmcpm' & parser$method_zero == 'pseudocount'){
    norm  = normalize_data(count, 'tmmcpm')
    count = norm$count
    size.factor = norm$size.factor

}else if (parser$method_transf == 'alr2' & parser$method_zero == 'pseudocount'){
    if (!is.numeric(nrefgene)) stop("Make sure you are giving the correct reference gene index for ALR-transformation")
    message("alr-transforming data")
    size.factor = count[,nrefgene]
    count = count[,-nrefgene]
    count = count / size.factor
    count = replace_zero(count, parser$method_zero, pseudo.count=1)
    count = log2(count)
}

# save output
message("saving output files")
fwrite(count, file=parser$output, quote=F, sep=",", row.names=F, col.names=F, compress="gzip")
if (parser$method_transf != 'log2') {
    fwrite(data.frame(size.factor), file=parser$output2, quote=F, sep=",", row.names=F, col.names=F)
} else {
    file.create(parser$output2)
}