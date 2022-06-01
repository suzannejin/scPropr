#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(ggplot2)
library(ggpubr)


# parse arguments
parser = ArgumentParser(description='Compare coefficients computed on relative vs absolute data.')
parser$add_argument('--ori', type='character', help="Original (raw) absolute count data")
parser$add_argument('--abs', type='character', help="Coefficients computed on absolute data")
parser$add_argument('--rel', type='character', help="Coefficients computed on relative data")
parser$add_argument('--features', type='character', help="List of gene names")
parser$add_argument('--refgene', type='character', help="Reference gene name")
parser$add_argument('--outdir', type='character', help="Output directory")
parser$add_argument('--filter', type='integer', default=20, help="Keep only the pairs whose dropout <= filtering threshold (in percentage). Default = 20")
parser$add_argument('--gm', action='store_true', help="Use geometric mean")
parser$add_argument('--table', action='store_true', help="Save table into csv file")
parser$add_argument('--npoint', type='integer', default=1e4, help="Number of points to plot")
parser$add_argument('--seed', type='integer', default=0, help="Random seed. Default=0")
parser = parser$parse_args()

filter = parser$filter / 100

calculate_geom_mean <-  function(x, y){
     df = data.frame(x,y)
     gm = apply( df, 1, function(x) exp(mean(log(x))) )
     return(gm)
}

# read input files
message('reading input data')
ori      = as.matrix(fread(parser$ori))
abs      = as.matrix(fread(parser$abs))
rel      = as.matrix(fread(parser$rel))
features = fread(parser$features, header=F)$V1
if ( !is.null(parser$refgene) ) {
    pos = which(features == parser$refgene)
    ori = ori[, -pos]
    abs = abs[-pos, -pos]
    features = features[-pos]
}
print(dim(ori))
colnames(ori) = features
colnames(abs) = features
rownames(abs) = features
colnames(rel) = features
rownames(rel) = features
if ( (nrow(abs) != nrow(rel)) || (ncol(abs) != ncol(rel)) ) stop("abs and rel matrices have different dimensiones") 
if ( (ncol(abs) != ncol(ori)) ) stop("abs and ori are based on different number of genes (ncol)") 
if ( length(features) != ncol(abs) ) stop("length(features) != ncol(abs)")

# get dropout and variance 
message('computing dropout and variance for each gene')
tmp           = ori / rowSums(ori) * mean(rowSums(ori))
dropout       = colMeans(ori == 0)
var_genes_abs = apply(ori, 2, var)
var_genes_rel = apply(tmp, 2, var)
rm(ori)
rm(tmp)

# organize data frame
message('organizing data frame and get random sample')
ind  = which( upper.tri(abs, diag=F) , arr.ind = TRUE )
comb = data.frame(
    i         = features[ind[,1]],
    j         = features[ind[,2]],
    abs       = abs[ind],
    rel       = rel[ind],
    dropout_i = dropout[ind[,1]],
    dropout_j = dropout[ind[,2]],
    var_abs_i = var_genes_abs[ind[,1]],
    var_abs_j = var_genes_abs[ind[,2]],
    var_rel_i = var_genes_rel[ind[,1]],
    var_rel_j = var_genes_rel[ind[,2]]
)
pos = which(comb[,'dropout_i'] <= filter & comb[,'dropout_j'] <= filter)
comb_filt = comb[pos,]
if (nrow(comb) > parser$npoint){
     set.seed(parser$seed); pos = sample(c(1:nrow(comb)), size=parser$npoint)
     comb = comb[pos,]
}
if (parser$gm){
     comb[,'dropout'] = calculate_geom_mean( 1-comb[,'dropout_i'], 1-comb[,'dropout_j'] )
     comb[,'var_abs'] = calculate_geom_mean( comb[,'var_abs_i'], comb[,'var_abs_j'] )
     comb[,'var_rel'] = calculate_geom_mean( comb[,'var_rel_i'], comb[,'var_rel_j'] )
     suf = '-gm'
}else{
     comb[,'dropout'] = rowMeans( data.frame(comb[,'dropout_i'], comb[,'dropout_j']) )
     comb[,'var_abs'] = rowMeans( data.frame(comb[,'var_abs_i'], comb[,'var_abs_j']) )
     comb[,'var_rel'] = rowMeans( data.frame(comb[,'var_rel_i'], comb[,'var_rel_j']) )
     suf = ''
}

# comb[,'i_j']     = paste0( comb[,'i'], '_', comb[,'j'] )
if (nrow(comb_filt) > parser$npoint){
     set.seed(parser$seed); pos = sample(c(1:nrow(comb_filt)), size=parser$npoint)
     comb_filt = comb_filt[pos,]
}
if (parser$gm){
     comb_filt[,'dropout'] = calculate_geom_mean( 1-comb_filt[,'dropout_i'], 1-comb_filt[,'dropout_j'] )
     comb_filt[,'var_abs'] = calculate_geom_mean( comb_filt[,'var_abs_i'], comb_filt[,'var_abs_j'] )
     comb_filt[,'var_rel'] = calculate_geom_mean( comb_filt[,'var_rel_i'], comb_filt[,'var_rel_j'] )
}else{
     comb_filt[,'dropout'] = rowMeans( data.frame(comb_filt[,'dropout_i'], comb_filt[,'dropout_j']) )
     comb_filt[,'var_abs'] = rowMeans( data.frame(comb_filt[,'var_abs_i'], comb_filt[,'var_abs_j']) )
     comb_filt[,'var_rel'] = rowMeans( data.frame(comb_filt[,'var_rel_i'], comb_filt[,'var_rel_j']) )
}
dim(comb)
head(comb)
dim(comb_filt)
head(comb_filt)

# save table
if ( parser$table ){
     message('saving table')
     fwrite(
          comb,
          file=paste0(parser$outdir, '/log2abs-vs-rel-cor-table.csv'),
          quote=F,
          sep=',',
          row.names=F,
          col.names=T
     )
}

# plot abs vs relative
message('plotting figure')
g1 = ggplot(comb, aes_string(x='abs', y='rel', color='dropout')) +
     geom_point(alpha=.1, size=.5) + 
     geom_abline(intercept = 0, slope = 1, linetype = 2, color="gray46") +
     xlab("On absolute data") +
     ylab("On relative data") +
     scale_colour_viridis_c( limits=c(0,1) ) +
     theme(strip.text = element_text(size = 12)) 
g2 = ggplot(comb, aes_string(x='abs', y='rel', color='var_abs')) +
     geom_point(alpha=.3, size=.5) + 
     geom_abline(intercept = 0, slope = 1, linetype = 2, color="gray46") +
     xlab("On absolute data") +
     ylab("On relative data") +
     scale_colour_viridis_c(trans = scales::pseudo_log_trans()) +
     # geom_text(data=subset(comb, var_abs > 1500), aes(abs, rel, color=var_abs, label=i_j), size=2) +
     theme(strip.text = element_text(size = 12)) 
g3 = ggplot(comb, aes_string(x='abs', y='rel', color='var_rel')) +
     geom_point(alpha=.3, size=.5) + 
     geom_abline(intercept = 0, slope = 1, linetype = 2, color="gray46") +
     xlab("On absolute data") +
     ylab("On relative data") +
     scale_colour_viridis_c(trans = scales::pseudo_log_trans()) +
     # geom_text(data=subset(comb, var_rel > 1500), aes(abs, rel, color=var_rel, label=i_j), size=2) +
     theme(strip.text = element_text(size = 12)) 
g = ggarrange(g1, g2, g3, ncol=3)
ggsave(
     g, 
     filename=paste0(parser$outdir, '/log2abs-vs-rel-cor-colored', suf, '.png'), 
     width = 15, 
     height = 4
)
g1 = ggplot(comb_filt, aes_string(x='abs', y='rel', color='dropout')) +
     geom_point(alpha=.1, size=.5) + 
     geom_abline(intercept = 0, slope = 1, linetype = 2, color="gray46") +
     xlab("On absolute data") +
     ylab("On relative data") +
     scale_colour_viridis_c() +
     theme(strip.text = element_text(size = 12)) 
g2 = ggplot(comb_filt, aes_string(x='abs', y='rel', color='var_abs')) +
     geom_point(alpha=.3, size=.5) + 
     geom_abline(intercept = 0, slope = 1, linetype = 2, color="gray46") +
     xlab("On absolute data") +
     ylab("On relative data") +
     scale_colour_viridis_c(trans = scales::pseudo_log_trans()) +
     # geom_text(data=subset(comb, var_abs > 1500), aes(abs, rel, color=var_abs, label=i_j), size=2) +
     theme(strip.text = element_text(size = 12)) 
g3 = ggplot(comb_filt, aes_string(x='abs', y='rel', color='var_rel')) +
     geom_point(alpha=.3, size=.5) + 
     geom_abline(intercept = 0, slope = 1, linetype = 2, color="gray46") +
     xlab("On absolute data") +
     ylab("On relative data") +
     scale_colour_viridis_c(trans = scales::pseudo_log_trans()) +
     # geom_text(data=subset(comb, var_rel > 1500), aes(abs, rel, color=var_rel, label=i_j), size=2) +
     theme(strip.text = element_text(size = 12)) 
g = ggarrange(g1, g2, g3, ncol=3)
ggsave(
     g, 
     filename=paste0(parser$outdir, '/log2abs-vs-rel-cor-colored-filter', parser$filter, suf, '.png'), 
     width = 15, 
     height = 4
)
comb = comb[which(comb[,'dropout_i'] <= filter & comb[,'dropout_j'] <= filter),]
g1 = ggplot(comb, aes_string(x='abs', y='rel', color='dropout')) +
     geom_point(alpha=.1, size=.5) + 
     geom_abline(intercept = 0, slope = 1, linetype = 2, color="gray46") +
     xlab("On absolute data") +
     ylab("On relative data") +
     scale_colour_viridis_c() +
     theme(strip.text = element_text(size = 12)) 
g2 = ggplot(comb, aes_string(x='abs', y='rel', color='var_abs')) +
     geom_point(alpha=.3, size=.5) + 
     geom_abline(intercept = 0, slope = 1, linetype = 2, color="gray46") +
     xlab("On absolute data") +
     ylab("On relative data") +
     scale_colour_viridis_c(trans = scales::pseudo_log_trans()) +
     # geom_text(data=subset(comb, var_abs > 1500), aes(abs, rel, color=var_abs, label=i_j), size=2) +
     theme(strip.text = element_text(size = 12)) 
g3 = ggplot(comb, aes_string(x='abs', y='rel', color='var_rel')) +
     geom_point(alpha=.3, size=.5) + 
     geom_abline(intercept = 0, slope = 1, linetype = 2, color="gray46") +
     xlab("On absolute data") +
     ylab("On relative data") +
     scale_colour_viridis_c(trans = scales::pseudo_log_trans()) +
     # geom_text(data=subset(comb, var_rel > 1500), aes(abs, rel, color=var_rel, label=i_j), size=2) +
     theme(strip.text = element_text(size = 12)) 
g = ggarrange(g1, g2, g3, ncol=3)
ggsave(
     g, 
     filename=paste0(parser$outdir, '/log2abs-vs-rel-cor-colored-filter', parser$filter, suf, '-2.png'), 
     width = 15, 
     height = 4
)
message('finished')