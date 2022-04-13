#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(easyCODA)

# parse arguments
parser = ArgumentParser(description='Compute procrustes analysis')
parser$add_argument('--count', type='character', help="Input count data (with zero already imputed)")
parser$add_argument('--lra', type='character', help="Exact log ratio geometric (.rds)")
parser$add_argument('--features', type='character', help="List of gene names")
parser$add_argument('--gene', type='character', help="Gene name to compute procrustes analysis")
parser$add_argument('--output', type='character', help="Output stat filename")
parser = parser$parse_args()

# functions -----------------------------------------------------------------------------

FINDALR <- function(data, data.lra, j, gene, weight = FALSE) { 
  
  ### -------------------------------------------------------------------
  ### function to identify the best reference for a set of ALRs
  ### various statistics are computed for each reference to assist
  ### the choice of best reference
  ### data is a normalized data matrix
  ### equal weighting is default for the logratio geometry here
  ### row (sample) weighting not implemented in this version
  
  if(sum(data[1,]!=1)) data <- data/rowSums(data)
  
  ### get exact logratio geometry
  data.lra.rpc <- data.lra$rowpcoord 
  tot.var <- sum(data.lra$sv^2)
  
  ### compute Procrustes correlation
  message("computing the procrates correlation for gene ", j, " - ", gene)
  dim <- min(nrow(data), ncol(data)) - 1
  # ALR transformation
  alr <- ALR(data, denom=j, weight=weight)
  # ALR geometry using PCA 'by hand' using SVD, without or with weighting
  if(!weight) {
    alr.svd <- svd(sqrt(1/nrow(alr$LR)) * sweep(alr$LR, 2, colMeans(alr$LR)) * sqrt(1/ncol(alr$LR)))
    alr.rpc <- sqrt(nrow(alr$LR)) * alr.svd$u %*% diag(alr.svd$d)
  }
  if(weight) {
    c <- colMeans(data)
    cc <- c*c[j]
    cc <- cc[-j]
    alr.svd <- svd(sqrt(1/nrow(alr$LR)) * sweep(alr$LR, 2, colMeans(alr$LR)) %*%  diag(sqrt(cc)))
    alr.rpc <- sqrt(nrow(alr$LR)) * alr.svd$u %*% diag(alr.svd$d)
  }
  procrust.cor <- protest(alr.rpc[,1:dim],data.lra.rpc, permutations=0)$t0

  message("organizing stats")

  ### the variances of the log-transformed parts
  var.log <- as.numeric( var(log(data[,..j])) ) 

  return(list(
      j=j, 
      id=gene,
      tot.var=tot.var, 
      procrust.cor=procrust.cor, 
      var.log=var.log)
  )
}


# read and process input data -----------------------------------------------------------
data = fread(parser$count)
features = fread(parser$features, header=F)$V1
if (any(data == 0)) stop("please provide a dataset with no zero values (needed for the log ratio analysis)")
if ( ncol(data) != length(features) ) stop("ncol(data) != length(features)")
if (!parser$gene %in% features) stop("incorrect gene")
ngene = which(features == parser$gene)
data.lra = readRDS(parser$lra)

# compute stats -------------------------------------------------------------------------
stat = FINDALR(data, data.lra, ngene, parser$gene, weight=F)
part = data[,..ngene]
dropout = mean(part == 0)
stat$dropout = dropout
stat = data.frame(stat)
fwrite(stat, file=parser$output, quote=F, sep=",", row.names=F, col.names=T)
message("finished")