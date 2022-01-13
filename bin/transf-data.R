library(data.table)
library(propr)
# library(propr, lib.loc=paste0(.libPaths()[1], "/propr_sjin"))

# get arguments ---------------------------------------------------------------

# get arguments
args = commandArgs(trailingOnly = T)
countfile = args[1] # observations x features matrix
output = args[2]
method_zero = args[3]    # zcompositions, one, min, none, pseudocount
method_transf = args[4]  # tmm, log2, clr
genesfile = args[5]   # list with the positions of the genes of interest
refgene = args[6]    # reference gene index

# check arguments
if (! method_zero %in% c("zcompositions", "one", "min", "none")){
    stop("wrong zero replacement method - ", method_zero)
}else if (! method_transf %in% c("log2", "clr", "alr", "tmm", "scran")){
    stop("wrong transformation or normalization method - ", method_transf)
}

# read and process data --------------------------------------------------------

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

# functions -------------------------------------------------------------------

replace_zero <- function(count, method=c("zcompositions", "one", "min", "none")){

    method = match.arg(method)
    if (!any(count==0, na.rm=T)) return(count)
    if (method == "none") return(count)
    message("replacing zeros with method ", method)

    # replace zero
    if (method == "zcompositions"){
        require(zCompositions)
        count = cmultRepl(count, method="CZM", label=0, output="counts")  

    }else if (method == "one"){
        count[count == 0] = 1

    }else if (method == "min"){
        count = as.matrix(count) 
        zeros = count == 0
        count[zeros] = min(count[!zeros])
    }

    # check zeros and negative values
    if(any(count<0, na.rm=T)) stop("counts matrix contain negative values")

    return(count)
}

normalize_data <- function(count, method=c("tmm", "scran")){
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
        effective.libsize = y$samples$lib.size * y$samples$norm.factors
        count = t(count) / effective.libsize   # normalized count per million
    }else if (method == "scran"){
        require(scran)
        require(scuttle)
        require(SingleCellExperiment)
        count = t(count)   # it requires row=genes, col=cells
        sce = SingleCellExperiment(list(counts=count))
        clusters = quickCluster(sce)
        sce = computeSumFactors(sce, clusters=clusters)  
        size.factor = sce@colData@listData$sizeFactor
        count = t(count) / size.factor
    }
    # TODO add other popular normalization procedures
    return(count)
}

# transform/normalize data --------------------------------------------------------------

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

}else if (method_transf %in% c("tmm", "scran")){
    count = normalize_data(count, method_transf)
    count = replace_zero(count, method_zero)
    count = log2(count)
}

# save output
fwrite(count, file=output, quote=F, sep=",", row.names=F, col.names=F, compress="gzip")