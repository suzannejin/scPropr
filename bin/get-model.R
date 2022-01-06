library(scDesign2)
library(copula)   
library(data.table)

# usage: get-model.R <count> <feature positions> <features> <cell type>  <output file> --
args = commandArgs(trailingOnly = T)
count = fread(args[1])   
pos = fread(args[2])$V1
features = fread(args[3])$V1[pos]
cell_type = args[4]
out = args[5]

# rename and reorganize count data
# sdDesign2 requires: row=gene, column=cell (labeled with cell type)
count = t(count)   
rownames(count) = features
colnames(count) = rep(cell_type, ncol(count))

# fit model and save data --------------------------------------------------------------
message("fit model")
RNGkind("L'Ecuyer-CMRG")  # set the random generator
set.seed(1)
model <- fit_model_scDesign2(count, cell_type, sim_method = 'copula')[[1]]
saveRDS(model, out)