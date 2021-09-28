library(zCompositions)
library(propr, lib.loc=paste0(.libPaths()[1], "/propr_sjin"))

# usage: clr-data.R <count> <zero handling method> <output clr> <output zero replaced> ---------------
args = commandArgs(trailingOnly = T)
countfile = args[1]
zerohandling = args[2]
output = args[3]
if (length(args) >= 4){output2 = args[4]}

# read count data
count = as.matrix(read.csv(gzfile(countfile)))
dim(count)
class(count)

# count dropout
print(paste0("dropout: ", round(mean(count==0, na.rm=T), 3)))

# quit if clr is not required
if (zerohandling == "none"){
    quit(status=0)
}
if (!zerohandling %in% c("one", "min", "zcompositions")){
    stop("Error: please give a proper zero handling method {one, min, zcompositions}")
}

# process and transform data --------------------------------------------------

# replace zeros 
message("replace zeros")
if (zerohandling == "zcompositions"){
    count = cmultRepl(count, method="CZM", label=0, output="counts")  
}else if (zerohandling == "one"){
    count[count == 0] = 1
}else{
    zeros = count == 0
    count[zeros] = min(count[!zeros])
}

# check negative values
if(any(count<0, na.rm=T)) stop("counts matrix contain negative values")

# clr-transform data
message("CLR-transform data")
clr = propr:::proprCLR(count)

# save output
write.table(clr, file=gzfile(output), quote=F, sep=",")
if (exists("output2")){
    write.table(count, file=gzfile(output2), quote=F, sep=",")
}