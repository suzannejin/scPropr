#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(dplyr)
library(egg)
library(ggplot2)
library(scales)

# parse arguments
parser <- ArgumentParser(description='Plot sequencing depth info: rna per cell, genes per cell, cells per sample, reads per sample')
parser$add_argument('-i', '--input', type="character", nargs='+', help="Input file(s)")
parser$add_argument('-n', '--name', type="character", nargs='+', help="Dataset name(s)")
parser$add_argument('-o', '--output', type="character", help="Output filename")
parser = parser$parse_args()

# functions
compute_values <- function(count, name){
    df = data.frame('dataset'         = name, 
                    'RNA.per.cell'    = rowSums(count), 
                    'Genes.per.cell'  = rowSums(count!=0), 
                    'Cells.per.sample'= nrow(count),
                    'Reads.per.sample'= sum(count))
    return(df)
}
plot_violin <- function(df, xvar){
    if (mean(df[,xvar]) > 1000){ 
        g = ggplot(df, aes_string(y='dataset', x=xvar)) +
            geom_violin() +
            geom_boxplot(width=0.1, color="grey", alpha=0.2, outlier.shape = NA) + 
            scale_x_continuous(position = "top", labels = label_number(suffix = "k", scale = 1e-3)) +
            ylab("")
    }else{
        g = ggplot(df, aes_string(y='dataset', x=xvar)) +
            geom_violin() +
            geom_boxplot(width=0.1, color="grey", alpha=0.2, outlier.shape = NA) + 
            scale_x_continuous(position = "top") +
            ylab("")
    }
    return(g)
}
plot_bar_cells <- function(df){
    g = ggplot(df, aes(y=dataset, x=Cells.per.sample)) + 
            geom_bar(stat = "identity", width=.75) + 
            scale_x_continuous(position = "top") +
            ylab("")
}
plot_bar_reads <- function(df){
    g = ggplot(df, aes(y=dataset, x=Reads.per.sample)) + 
            geom_bar(stat = "identity", width=.75) + 
            scale_x_continuous(position = "top", labels = label_number(suffix = "M", scale = 1e-6)) +
            ylab("")
}

# compute seqdepth per dataset
df = data.frame()
for (i in 1:length(parser$input)){
    count = fread(parser$input[i])
    name  = parser$name[i]
    df2   = compute_values(count, name)
    df    = rbind(df, df2)
}

# plot figure
g1 = plot_violin(df, "RNA.per.cell")
g2 = plot_violin(df, "Genes.per.cell") + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
g3 = plot_bar_cells(distinct(df[,c(1,4)])) + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
g4 = plot_bar_reads(distinct(df[,c(1,5)])) + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
gs = list(g1,g2,g3,g4)
g  = ggarrange(plots=gs, nrow=1)
ggsave(parser$output, plot=g, width = 12, height = length(parser$input))