#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(ggplot2)

# parse arguments
parser = ArgumentParser(description='Compare transformed count on relative vs absolute data.')
parser$add_argument('input1', type='character', help="Input table - not filtered")
parser$add_argument('input2', type='character', help="Input table - filtered")
parser$add_argument('output', type='character', help="Output figure filename")
parser = parser$parse_args()

dic = c(
    "phaseS"                       = "phaseS",
    "phaseS+simulate+s2+n5+c1"     = "slope=2",
    "phaseS+simulate+s4+n5+c1"     = "slope=4",
    "phaseS+simulate+s6+n5+c1"     = "slope=6",
    "phaseS+simulate+s8+n5+c1"     = "slope=8",
    "phaseS+rel"                   = "phaseS+rel",
    "phaseS+rel+simulate+s2+n5+c1" = "slope=2",
    "phaseS+rel+simulate+s4+n5+c1" = "slope=4",
    "phaseS+rel+simulate+s6+n5+c1" = "slope=6",
    "phaseS+rel+simulate+s8+n5+c1" = "slope=8"
)

df1 = as.data.frame(fread(parser$input1, header=T))
df2 = as.data.frame(fread(parser$input2, header=T))
df1[,'filt'] = 'not filtered'
df2[,'filt'] = 'filtered'
df = rbind(df1, df2)
df[,"type"] = as.vector(dic[df[,"type"]])
g  = ggplot(df, aes(x=type, y=val, color=filt)) +
     geom_point() +
     scale_colour_manual(values=c('darkorange3','darkslategray')) +
     facet_wrap( ~ cor + transf, ncol=length(unique(df$transf))) +
     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(
    parser$output, 
    plot   = g, 
    width  = 2 * length(unique(df$transf)),
    height = 2 * length(unique(df$cor))
)