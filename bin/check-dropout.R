library(data.table)

# usage: check-dropout.R <simulated data> <features pos> <output> --------------------------
args = commandArgs(trailingOnly = T)
count = fread(args[1])
features_pos = args[2]
output = args[3]

if (file.exists(features_pos)){
    g = fread(features_pos, header=F)$V1
    count = count[,..g]
}

# check dropout ----------------------------------------------------------------------------
dropout = mean(count==0, na.rm=T)
cat(dropout, file=output, sep="\n")