#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(ggplot2)

# parse arguments
parser = ArgumentParser(description='Compare transformed count on relative vs absolute data.')
parser$add_argument('input', type='character', help="Input table")
parser$add_argument('output', type='character', help="Output figure filename")
parser = parser$parse_args()

dic = c(
    "phaseS+simulate+s2+n5+c3+d3"     = "slope=2",
    "phaseS+simulate+s4+n5+c3+d3"     = "slope=4",
    "phaseS+simulate+s6+n5+c3+d3"     = "slope=6",
    "phaseS+simulate+s8+n5+c3+d3"     = "slope=8"
)

df = as.data.frame(fread(parser$input, header=T))
df[,"type"] = as.vector(dic[df[,"type"]])
g  = ggplot(df, aes(x=type, y=val)) +
     geom_point() +
     facet_wrap( ~ cor + transf, ncol=length(unique(df$transf))) +
     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(
    parser$output, 
    plot   = g, 
    width  = 2 * length(unique(df$transf)),
    height = 2 * length(unique(df$cor))
)