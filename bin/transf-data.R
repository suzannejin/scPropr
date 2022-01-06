library(data.table)
library(propr)
# library(propr, lib.loc=paste0(.libPaths()[1], "/propr_sjin"))

# get arguments ---------------------------------------------------------------

# get arguments
args = commandArgs(trailingOnly = T)
countfile = args[1] # observations x features matrix
output = args[2]
method_zero = args[3]    # zcompositions, one, min, none
method_transf = args[4]  # tmm, log2, clr
genesfile = args[5]   # list with the positions of the genes of interest
refgene = args[6]    # reference gene index

# check arguments
if (! method_zero %in% c("zcompositions", "one", "min", "none")){
    stop("wrong zero replacement method - ", method_zero)
}else if (! method_transf %in% c("log2", "clr", "alr", "tmm", "scran")){
    stop("wrong transformation or normalization method - ", method_transf)
}

# functions -------------------------------------------------------------------

replace_zero <- function(count, method=c("zcompositions", "one", "min")){

    method = match.arg(method)

    # check if zero
    if (!any(count==0, na.rm=T)) return(count)
    message("replacing zeros with method ", method)

    # replace zero
    if (method == "zcompositions"){
        require(zCompositions)
        count = cmultRepl(count, method="CZM", label=0, output="counts")  

    }else if (method == "one"){
        count[count == 0] = 1

    }else{
        count = as.matrix(count) 
        zeros = count == 0
        count[zeros] = min(count[!zeros])
    }

    # check negative values
    if(any(count<0, na.rm=T)) stop("counts matrix contain negative values")

    return(count)
}

normalize_data <- function(count, method=c("tmm", "scran"), pseudo.count=0){
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
        count = t(count) # it requires row=genes, col=cells
        y = DGEList(counts=count)   
        y = calcNormFactors(y, method="TMM")
        count = cpm(y, log=T, prior.count=pseudo.count)   # it calculates the normalized count per million. So normalized counts = pseudocount / effective library size * 10E6, then the normalized count is log2-transformed
        count = t(count)
    }else if (method == "scran"){
        require(scran)
        require(scuttle)
        require(SingleCellExperiment)
        count = t(count)   # it requires row=genes, col=cells
        sce = SingleCellExperiment(list(counts=count))
        clusters = quickCluster(sce)
        sce = computeSumFactors(sce, clusters=clusters)  #compute the size factors
        sce = logNormCounts(sce, pseudo.count=pseudo.count)   # compute log2(pseudocount / size factor)
        count = assays(sce)[["logcounts"]]
        count = t(count)
    }
    # TODO add other popular normalization procedures
    return(count)
}

# process and transform data --------------------------------------------------

# read count data
count = fread(countfile)

# filter genes
if (file.exists(genesfile)){
    g = fread(genesfile, header=F)$V1
    count = count[,..g]
}

# count dropout
dropout = round(mean(count==0, na.rm=T), 3)
print(paste0("dropout: ", dropout))
dim(count)
class(count)

# handle zero
if (dropout > 0){
    handle_zero = TRUE
    pseudo.count = 1
}else{
    handle_zero = FALSE
}

# clr-transform or normalize
if (method_transf == "log2"){
    if (handle_zero) count = replace_zero(count, method_zero)
    message("log-transforming data")
    count = log2(count)
}else if(method_transf == "clr"){
    if (handle_zero) count = replace_zero(count, method_zero)
    message("clr-transforming data")
    count = propr:::proprCLR(count)
}else if(method_transf == "alr"){
    if (refgene == "NA"){
        stop("Make sure you are giving the correct reference gene index for ALR-transformation")
    }else{
        refgene = as.numeric(refgene)
    }
    if (handle_zero) count = replace_zero(count, method_zero)
    message("alr-transforming data")
    count = propr:::proprALR(count, refgene)
}else if (method_transf %in% c("tmm", "scran")){
    count = normalize_data(count, method_transf, pseudo.count=1)   # THINK. Do we need to use zcompositions instead of a pseudocount of 1, to unify comparison with other methods? And control for zero imputation effect on the analysis?
}

# save output
fwrite(count, file=output, quote=F, sep=",", row.names=F, col.names=F, compress="gzip")
