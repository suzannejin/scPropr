library(data.table)
library(zCompositions)
library(propr, lib.loc=paste0(.libPaths()[1], "/propr_sjin"))

# usage: clr-data.R <count> <method> <output> ---------------
args = commandArgs(trailingOnly = T)
countfile = args[1] # observations x features matrix
method = args[2]    # zcompositions, one, min, tmm
output = args[3]
clr_methods = c("zcompositions","one","min")
norm_methods = c("tmm")
valid_methods = c(clr_methods, norm_methods)
if (!method %in% valid_methods){
    stop("Please give a proper zero handling method {", toString(valid_methods), "}")
}

# functions -------------------------------------------------------------------

clr_data <- function(count, method=c("zcompositions", "one", "min")){

    # replace zeros
    message("replace zeros")
    count = replace_zero(count, method)

    # check negative values
    if(any(count<0, na.rm=T)) stop("counts matrix contain negative values")

    # clr-transform data
    message("CLR-transform data")
    clr = propr:::proprCLR(count)

    return(clr)
}

replace_zero <- function(count, method=c("zcompositions", "one", "min")){
    if (method == "zcompositions"){
        count = cmultRepl(count, method="CZM", label=0, output="counts")  
    }else if (method == "one"){
        count[count == 0] = 1
    }else{
        count = as.matrix(count) 
        zeros = count == 0
        count[zeros] = min(count[!zeros])
    }
    return(count)
}

normalize_data <- function(count, method=c("tmm")){
    if (method == "tmm"){
        require(edgeR)
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
dim(count)
class(count)

# count dropout
print(paste0("dropout: ", round(mean(count==0, na.rm=T), 3)))

# clr-transform or normalize
if (method %in% clr_methods){
    count = clr_data(count, method)
}else{
    count = normalize_data(count, method)
}

# save output
fwrite(count, file=output, quote=F, sep=",", row.names=F, col.names=F, compress="gzip")