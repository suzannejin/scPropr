#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(ggplot2)
library(ggpubr)

# parse arguments
parser = ArgumentParser(description='Plot log expression vs sorted library size, and coloring the normalization factor (scale). The figure is obtained for absolute and relative data')
parser$add_argument('--abs', type='character', help="Absolute count data")
parser$add_argument('--rel', type='character', help="Relative count data")
parser$add_argument('--abs_size', type='character', help="File with normalization factors for absolute data")
parser$add_argument('--rel_size', type='character', help="File with normalization factors for relative data")
parser$add_argument('--norm', type='character', help="Name of the normalization factor")
parser$add_argument('--outdir', type='character', help="Output directory")
parser$add_argument('--ncell', type='integer', help="Number of cells to plot")
parser = parser$parse_args()

# functions -----------------------------------------------------------------------------
compute_boxplot_stats <- function(x, coef=1.5){
    stats = as.numeric(quantile(x, na.rm = TRUE))
    iqr   = diff(stats[c(2, 4)])
    out   = x < (stats[2L] - coef * iqr) | x > (stats[4L] + coef * iqr)
    stats[c(1, 5)] = range(x[!out], na.rm = TRUE)
    return(stats)
}
get_boxplot_df <- function(count, x, sizefactor, norm_name){
    count = log10(count)  # we want to compute the boxplot on log10 values
    x     = log10(x)   
    stats = apply(count, 1, function(x) compute_boxplot_stats(x))
    df = data.frame(
        group = 1:nrow(count),
        x     = x,
        min   = stats[1,],
        low   = stats[2,],
        mid   = stats[3,],
        top   = stats[4,],
        max   = stats[5,],
        norm  = log10(sizefactor),
        norm_name = norm_name
    )
    if (norm_name %in% c('tmm', 'scran')) df$norm = sizefactor
    return(df)
}
plot_boxplot <- function(df, xlab, ylab, title, norm, pointsize=1, ylimits=c(1, 1e4)){

    if (norm == 'tmm'){
        sc = .5
    }else if (norm == 'scran'){
        sc = 2
    }else{
        sc =1
    }
    
    # plot
    g = ggplot(df, aes(x=x, ymin = min, lower = low, middle = mid, upper = top, ymax = max, group=group)) +
        geom_boxplot(stat="identity") +
        geom_point(aes(x=x, y=norm/sc, color=norm_name), size=.7)
    
    # y axis
    if (norm %in% c('clr', 'alr')){
        g = g +
            scale_y_continuous(
                breaks = seq(0,4,by = 1),
                limits = c(0,4),
                labels = scales::math_format(),
                sec.axis = sec_axis(~.*sc, name = "Normalization factor", labels = scales::math_format())
            )
    }else{
        g = g +
            scale_y_continuous(
                breaks = seq(0,4,by = 1),
                limits = c(0,4),
                labels = scales::math_format(),
                sec.axis = sec_axis(~.*sc, name = "Normalization factor")
            ) 
    }
    # x axis
    g = g +
        scale_x_continuous(
            breaks = seq(2,5,by = 1),
            limits = c(2,5),
            labels = scales::math_format()
        )
    
    # modify theme
    g = g +
        xlab(xlab) +
        ylab(ylab) +
        ggtitle(title) + 
        theme(legend.title = element_blank())
    
    return(g)
}

# read data -----------------------------------------------------------------------------
abs = fread(parser$abs)
rel = fread(parser$rel)
abs_size = fread(parser$abs_size)$V1
rel_size = fread(parser$rel_size)$V1

# limit the number of cells (boxplot) showing
if (is.numeric(parser$ncell)){
    if (nrow(abs) > parser$ncell){
        set.seed(0); 
        pos = sample(1:nrow(abs), parser$ncell)
        abs = abs[pos,]
        rel = rel[pos,]
        abs_size = abs_size[pos]
        rel_size = rel_size[pos]
    }
}

print(dim(abs))
print(mean(abs==0))

# reorder count data
libsize = rowSums(abs)
idx = order(libsize)
abs = abs[idx,]
rel = rel[idx,]
abs_size = abs_size[idx]
rel_size = rel_size[idx]

# parse data
df_abs = get_boxplot_df(abs, rowSums(abs), abs_size, parser$norm)
df_rel = get_boxplot_df(rel, rowSums(abs), rel_size, parser$norm)

# plot ----------------------------------------------------------------------------------
g_abs = plot_boxplot(df_abs, "Library size of absolute data", "log10 expression", "Absolute data", parser$norm)
g_rel = plot_boxplot(df_rel, "Library size of absolute data", "log10 expression", "Relative data", parser$norm)
g = ggarrange(g_abs, g_rel, ncol=1, nrow=2, common.legend=T, legend="bottom")
ggsave(g, filename=paste0(parser$outdir, '/expr-libsize-', parser$norm, '.png'), width = 12, height = 8)