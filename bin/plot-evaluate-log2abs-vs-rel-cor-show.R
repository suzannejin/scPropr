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
parser$add_argument('--color', type='character', nargs='+')
parser = parser$parse_args()

# dic = c(
#     "experimental"                    = "experimental",
#     "phaseS+simulate+s1+n5+c1+d1"     = "slope=1",
#     "phaseS+simulate+s2+n5+c1+d1"     = "slope=2",
#     "phaseS+simulate+s4+n5+c1+d1"     = "slope=4",
#     "phaseS+simulate+s6+n5+c1+d1"     = "slope=6",
#     "phaseS+simulate+s8+n5+c1+d1"     = "slope=8"
# )

dic = c(
    "experimental"                    = "experimental",
    "phaseS+simulate+s1+n5+c1+d1"     = "1",
    "phaseS+simulate+s2+n5+c1+d1"     = "2",
    "phaseS+simulate+s4+n5+c1+d1"     = "4",
    "phaseS+simulate+s6+n5+c1+d1"     = "6",
    "phaseS+simulate+s8+n5+c1+d1"     = "8"
)

dic2 = c(
    "NA_log2_NA"                       = "log2",
    "NA_tmm_NA"                        = "TMM",
    "NA_scran_NA"                      = "scran",
    "NA_clr_NA"                        = "clr",
    "NA_alr_MT-RNR2"                   = "alr (MT-RNR2)",
    "NA_alr_RPL8"                      = "alr (RPL8)"
)

dic_cor = c(
    "cor"                              = "Pearson correlation",
    "rho"                              = "Proportionality",
    "pcor"                             = "Partial correlation",
    "pcor.shrink"                      = "Partial correlation (regularized)"
)

message('reading and organizing input')
df = as.data.frame(fread(parser$input, header=T))
df = df[which(df[,"type"] %in% names(dic)), ]
df[,"type"] = as.vector(dic[df[,"type"]])
if ( !is.null(parser$cor) ){
    pos= which(df[,"cor"] %in% parser$cor)
    df = df[pos,]
    df[,"cor"] = as.vector(dic_cor[df[,"cor"]])
    df[,"cor"] = factor(df[,"cor"], level=dic_cor[parser$cor])
}
if ( !is.null(parser$transf) ){
    pos= which(df[,"transf"] %in% parser$transf)
    df = df[pos,]
    df[,"transf"] = as.vector(dic2[df[,"transf"]])
    df$transf = factor(df$transf, levels=dic2[parser$transf])
}

message('plotting figure')
g  = ggplot(df, aes(x=type, y=val, color=cor)) +
     geom_point() +
     facet_wrap( ~ transf, nrow = 1) +
     theme(legend.position="top",
       strip.text.x = element_text(face = "bold")) +
     labs(colour="") #+
    #  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
if (!is.null(parser$color)){
    g = g + scale_color_manual(values=parser$color) 
}
head(df)
ggsave(
    parser$output, 
    plot   = g, 
    width  = 2 * length(unique(df$transf)),
    height = 3 # 2 * length(unique(df$cor))
)
message('finished')