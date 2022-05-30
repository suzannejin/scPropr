#!/usr/bin/env Rscript

library(argparse)
library(data.table)

# parse arguments -----------------------------------------------------------------------
parser = ArgumentParser(description='Replace zeros')
parser$add_argument('-i', '--input', type='character', help="Input count data")
parser$add_argument('-o', '--output', type='character', help="Output count data")
parser$add_argument('--method', type='character', help="Zero handling method. Choices = [zcompositions, one, min]")
parser = parser$parse_args()


# functions -----------------------------------------------------------------------------
replace_zero <- function(count, method=c("zcompositions", "one", "min")){

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
        count = as.matrix(count) 
        zeros = count == 0
        count[zeros] = min(count[!zeros])
    }

    # TODO add scImpute, etc zero handling methods for single cell

    # check zeros and negative values
    if(any(count<=0, na.rm=T)) stop("counts matrix contain zero and/or negative values")

    return(count)
}

# read and process data -----------------------------------------------------------------
count = fread(parser$input)
count = replace_zero(count, method=parser$method)
fwrite(count, file=parser$output, quote=F, sep=",", row.names=F, col.names=F, compress="gzip")
message('finished')