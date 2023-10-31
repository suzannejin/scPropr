#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(ggplot2)

# parse arguments
parser = ArgumentParser(description='Plot histogram: gene pairs vs the %non-zero cells')
parser$add_argument('-i', '--input', type='character', help="Input count data")
parser$add_argument('-o', '--outdir', type='character', help="Output directory")
parser$add_argument('-f', '--filter', type='integer', help="Filter percentage of cells in descending order acording to the dropout rate. So the cells with highest dropout rate will be removed.")
parser$add_argument('-t', '--threshold', type='integer', help="Filter cells with dropout > threshold (in percentage)")
parser$add_argument('--avg', action='store_true')
parser$add_argument('--gm', action='store_true')
parser$add_argument('--both', action='store_true')
parser = parser$parse_args()


# read input count data --------------------------------------------------------

message('reading input count data')
count = as.matrix(fread(parser$input))

# filter cells, if needed
pref = 'null'
if (!is.null(parser$filter)){
    message('Filtering ', parser$filter, '% cells')
    zeros = rowMeans(count == 0)
    ord   = order(zeros, decreasing=F)
    cut   = round( length(ord) * (100 - parser$filter)/100, 0 )
    pos   = ord[1:cut]
    count = count[pos,]   ## NOTE here that the order that original cells appear in the matrix is not perserved
    pref  = paste0('filter', parser$filter)
}else if(!is.null(parser$threshold)){
    message('Filtering cells with dropout > ', parser$threshold, '%')
    pos   = which(rowMeans(count==0) > parser$threshold/100)
    count = count[pos,]
    pref  = paste0('threshold', parser$threshold)
}

# count dropout
dropout = round(mean(count==0, na.rm=T), 3)
print(paste0("dropout: ", dropout))
dim(count)
class(count)


# functions ---------------------------------------------------------------------

calculate_geom_mean <-  function(count){
    gm = apply( count, 1, function(x) exp(mean(log(x[x!=0]))) )
    return(gm)
}

# Get data frame with columns: i, j, nonzero_i, nonzero_j, nonzero_avg, nonzero_gm
# if redundant, do it for all pairs i -> j and j -> i reciprocally
get_df_comb <- function(count, redundant=F){
    message("Compute average percentage of non-zero genes")
    nonzero = colMeans(count != 0)
    idx     = 1:length(nonzero)
    if (redundant){
        idx = expand.grid(idx, idx)
    }else{
        idx = t(combn(idx, 2))
        idx = idx[which(idx[,1] != idx[,2]),]
    }
    df = data.frame(i=idx[,1], j=idx[,2])
    df$nonzero_i = nonzero[idx[,1]]
    df$nonzero_j = nonzero[idx[,2]]
    df$nonzero_avg = rowMeans(df[,3:4])
    df$nonzero_gm = apply(df[,3:4], 1, calculate_geom_mean)
    return(df)
}

# Get data frame with columns: i, j, nonzero
# this function computes the percentage of cells that have no zeros in both genes, instead of averaging like previous function
get_df_ind <- function(count){
    message("For each gene pair, compute the percentage of cells that have nonzero values in both genes")
    ci=c(); cj=c(); cnonzero=c()
    for (i in 1:ncol(count-1)){
        for (j in (i+1):ncol(count)){
            if (i == j) next
            ci       = c(ci, i)
            cj       = c(cj, j)
            both     = count[,c(i,j)]
            both     = rowSums( both != 0 )
            both     = ifelse(both == 2, 1, 0)  # if non-zero cell for both genes
            nonzero  = mean(both == 1)  # average % cells that have non-zero for both genes
            cnonzero = c(cnonzero, nonzero)
        }
    }
    df = data.frame(i=ci, j=cj, nonzero=cnonzero)
    return(df)
}


# plot histogram of number gene pairs vs average % non zero cells
if (parser$avg || parser$gm){
    df1 = get_df_comb(count, redundant=F)
    head(df1)
    message("Plotting histogram")
    if (parser$avg){
        g1  = ggplot(df1, aes(x=nonzero_avg)) + 
            geom_histogram(color="black", fill="white") +
            xlab("Avg % non-zero cells") +
            ylab("Number of gene pairs")
        ggsave(
            g1, 
            filename = paste0(parser$outdir, '/histogram-nonzero-avg-', pref, '.png'), 
            width    = 6, 
            height   = 4
        )
    }
    if (parser$gm){
        g1  = ggplot(df1, aes(x=nonzero_gm)) + 
            geom_histogram(color="black", fill="white") +
            xlab("Gm % non-zero cells") +
            ylab("Number of gene pairs")
        ggsave(
            g1, 
            filename = paste0(parser$outdir, '/histogram-nonzero-gm-', pref, '.png'), 
            width    = 6, 
            height   = 4
        )
    }
}


# plot histogram of number of gene pairs vs average percentage of cells that have nonzero values in both genes
if (parser$both){
    df2 = get_df_ind(count)
    head(df2)
    message("Plotting histogram")
    g2  = ggplot(df2, aes(x=nonzero)) + 
            geom_histogram(color="black", fill="white") +
            xlab("Avg % non-zero cells") +
            ylab("Number of gene pairs")
    ggsave(
        g2, 
        filename = paste0(parser$outdir, '/histogram-nonzero-ind-', pref, '.png'), 
        width    = 6, 
        height   = 4
    )
}