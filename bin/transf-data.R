library(data.table)
library(propr)
# library(propr, lib.loc=paste0(.libPaths()[1], "/propr_sjin"))

# usage: clr-data.R <count> <method> <output> ---------------

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
}else if (! method_transf %in% c("log2", "clr", "alr", "tmm")){
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

normalize_data <- function(count, method=c("tmm")){
    method = match.arg(method)
    if (method == "tmm"){
        require(edgeR)
        message("normalizing data with method ", method)
        # Note that normalization in edgeR is model-based, and the original read counts are not themselves transformed. This means that users should not transform the read counts in any way
        # before inputing them to edgeR. For example, users should not enter RPKM or FPKM values to edgeR in place of read counts. Such quantities will prevent edgeR from correctly
        # estimating the mean-variance relationship in the data, which is a crucial to the statistical
        # strategies underlying edgeR. Similarly, users should not add artificial values to the counts
        # before inputing them to edgeR

        # For further information check: https://www.biostars.org/p/84087/ and https://www.biostars.org/p/317701/
        y = DGEList(counts=t(count))   # it requires a matrix of features x samples
        y = calcNormFactors(y, method="TMM")
        m = cpm(y)   # it calculates the normalized count per million. So normalized counts = counts / effective library size * 10E6
        count = t(m)
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
print(paste0("dropout: ", round(mean(count==0, na.rm=T), 3)))
dim(count)
class(count)

# clr-transform or normalize
if (method_transf == "log2"){
    count = replace_zero(count, method_zero)
    message("log-transforming data")
    count = log2(count)
}else if(method_transf == "clr"){
    count = replace_zero(count, method_zero)
    message("clr-transforming data")
    count = propr:::proprCLR(count)
}else if(method_transf == "alr"){
    if (refgene == "NA"){
        stop("Make sure you are giving the correct reference gene index for ALR-transformation")
    }else{
        refgene = as.numeric(refgene)
    }
    count = replace_zero(count, method_zero)
    message("alr-transforming data")
    count = propr:::proprALR(count, refgene)
}else if (method_transf == "tmm"){
    count = normalize_data(count, method_transf)
}

# save output
fwrite(count, file=output, quote=F, sep=",", row.names=F, col.names=F, compress="gzip")
