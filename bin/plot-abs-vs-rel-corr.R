library(data.table)
library(ggplot2)
library(ggpubr)

# usage: plot-abs-vs-rel.R 
args = commandArgs(trailingOnly = T)
cell.type = args[1]
data.type = args[2]
zero.method = args[3]
alr.ref = args[4]
transf.list = strsplit(args[5], ",")[[1]]
corr.list = strsplit(args[6], ",")[[1]]
out = args[7]
random = as.numeric(args[8])

# match association method names
dic = list("cor"   = "pearson",
           "rho"   = "proportionality",
           "vlr"   = "vlr",
           "pcor"  = "pcor",
           "pcor.shrink" = "pcorshrink")

# parse methods and filenames
methods = c(); abs = c(); rel = c()
n = 1
for (i in corr.list){
    for (j in transf.list){
        if (j == "alr") { k = alr.ref } else { k = "NA" }
        name = paste(n, j, dic[[i]], sep=".")
        name2_abs = paste(data.type, zero.method, "log2_NA", i, sep="_")
        name2_rel = paste(data.type, zero.method, j, k, i, sep="_")
        methods[n] = name
        abs[n] = paste0("absolute_", name2_abs, ".csv.gz")
        rel[n] = paste0("relative_", name2_rel, ".csv.gz")
        n = n + 1
    }
}

# Organize data frame
df = data.frame(matrix(ncol=3, nrow=0))
colnames(df) = c("x", "y", "z")
for (i in 1:length(methods)){
    # read x
    x = abs[i]
    x = as.matrix(fread(x))
    x = x[lower.tri(x)]
    
    # read y
    y = rel[i]
    y = as.matrix(fread(y))
    y = y[lower.tri(y)]

    # random sample
    if (random > 0){
        set.seed(0); x = sample(x, random)
        set.seed(0); y = sample(y, random)
    }
    
    # update data frame
    df2 = data.frame("x"=x, "y"=y, "z"=rep(methods[i], length(x)))
    df = rbind(df, df2)
}
df$z = factor(df$z, levels=methods)

# set figure dimensions
if (length(methods) == 25){
    nrow = 5
    width = height = 15
} else if (length(methods) == 15) {
    nrow = 3
    width = 15
    height = 9
}

# plot abs vs relative
g = ggplot(df, aes(x=x, y=y)) +
    facet_wrap(~z, scales="free", nrow=nrow) +
    geom_point(alpha=.1, size=.5) + 
    geom_abline(intercept = 0, slope = 1, linetype = 2, color="gray46") +
    xlab("On absolute data") +
    ylab("On relative data") +
    ggtitle( paste0(cell.type, " - ", data.type) ) + 
    theme(strip.text = element_text(size = 12))
ggsave(out, width = width, height = height)