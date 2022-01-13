library(data.table)

# usage: get-relative.R <count data> <output> ---------------------------------
args = commandArgs(trailingOnly = T)
count = fread(args[1])
output = args[2]

# get relative data -----------------------------------------------------------
relative = count / rowSums(count) * mean(rowSums(count))   # "Reclose" the data
                                                           # The closed data is multiplied by the average read/cell to keep the same seq depth across cells
fwrite(relative, file=output, quote=F, sep=",", row.names=F, col.names=F, compress="gzip")