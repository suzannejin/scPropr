#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(ggplot2)

# parse arguments
parser = ArgumentParser(description='Compare transformed count on relative vs absolute data.')
parser$add_argument('input', type='character', help="Input table")
parser$add_argument('output', type='character', help="Output figure filename")
parser$add_argument('--cor', type='character', nargs='+', help="Association methods to be considered")
parser$add_argument('--transf', type='character', nargs='+', help="Transformation/normalisation methods to be considered")
parser$add_argument('--name', type='character', nargs='+', help="Names for transformation/normalisation methods")
parser = parser$parse_args()


message('reading and organizing input')
df = as.data.frame(fread(parser$input, header=T))
pos= which(df[,"type"] == "experimental")
df = df[pos,]
if ( !is.null(parser$cor) ){
    pos= which(df[,"cor"] %in% parser$cor)
    df = df[pos,]
}
if ( !is.null(parser$transf) ){
    dic = c(parser$name)
    names(dic) = parser$transf
    pos= which(df[,"transf"] %in% parser$transf)
    df = df[pos,]
    df[,"transf"] = as.vector(dic[df[,"transf"]])
    df[,"transf"] = factor(df[,"transf"], levels=parser$name)
}

message('plotting figure')
g  = ggplot(df, aes(x=transf, y=val, color=cor)) +
     geom_point() +
     ylim(-0.05,1) +
     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(
    parser$output, 
    plot   = g, 
    width  = 4,
    height = 3 
)
message('finished')