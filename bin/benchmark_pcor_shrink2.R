#!/usr/bin/env Rscript

library(argparse)
library(corpcor)
library(data.table)
library(ggplot2)
library(ggpubr)

# parse arguments
parser = ArgumentParser(description='Benchmark shrinkage')
parser$add_argument('-i', '--input', type='character', help="Input count data")
parser$add_argument('-o', '--outdir', type='character', help="Output directory")
parser$add_argument('--transf', type='character', help='clr or alr')
parser$add_argument('--ngene', type='integer', help='Number of genes considered')
parser$add_argument('--ncell', type='integer', nargs='+', help='Number of cells considered')
parser$add_argument('--nsamp', type='integer', default=200, help='Number of sampling')
parser = parser$parse_args()


# read input files ----------------------------------------------------------------------

count = fread(parser$input)
count = as.matrix(count)
pos = which(colSums(count == 0) == 0)
count = count[,pos]

if (!file.exists(parser$outdir)){
    dir.create(parser$outdir)
}

# functions -----------------------------------------------------------------------------

bShrink <- function(M, intype, outtype="clr"){
    if (intype=="alr"){
        D=ncol(M)+1
    } else{
        D=ncol(M)
    }
    N=nrow(M)
    
    if (intype=="rco"){
        
        print("Warning: estimating compositions from counts via shrinkage")
        P=matrix(0,N,D)
        for (i in 1:nrow(M)){
            P[i,]=shrink(M[i,])
        }
        B=log(P)
    } else if (intype=="alr"){
        
        P=exp(M)/(1+apply(exp(M),1,sum))
        P=cbind(P,1-apply(P,1,sum))
        B=log(P)
        
    } else if (intype=="clr"){
        f=log(apply(exp(M),1,sum)) #backtransform to log P:
        B=M-f
    }
    
    Cov=cov(B)
    Cov.shrink=cov.shrink(B,verbose=FALSE)
    if (outtype=="alr"){
        F=cbind(diag(rep(1,D-1)),rep(-1,D-1))
        Cov.basis=F%*%Cov.shrink%*%t(F)
    } else if (outtype=="clr"){
        G=diag(rep(1,D))-matrix(1/D,D,D)
        Cov.basis=G%*%Cov.shrink%*%G
    } else{
        die("outtype unsupported")
    }
    Cor=cov2cor(Cov)
    PC=cor2pcor(Cov)
    Cor.shrink=cov2cor(Cov.shrink)
    PC.shrink=cor2pcor(Cov.shrink)
    Cor.basis=cov2cor(Cov.basis)
    PC.basis=cor2pcor(Cov.basis)
    
    
    C=list(Cov,Cov.shrink,Cov.basis,Cor,Cor.shrink,Cor.basis,PC,PC.shrink,PC.basis)
    names(C)=c('cov', 'cov.shrink', 'cov.basis', 'corr', 'corr.shrink', 'corr.basis', 'pcor', 'pcor.shrink', 'pcor.basis')
    return(C)    
}
shrink <- function(n){
        N=sum(n)
        q=n/N       
        D=length(n)
        T=rep(1/D,D)
        
        l=(1-sum(q^2))/((N-1)*sum((T-q)^2))
        if (l<0){
            l=0
        }
        if(l>1){
            l=1
        }
        qs=l*T+(1-l)*q
        
        return(qs)
}
calculate_mse <- function(true, pred){
    l = list()
    for (i in c('cov', 'cov.shrink', 'cov.basis')){
        l[[i]] = mse(true[['cov']], pred[[i]])
    }
    for (i in c('corr', 'corr.shrink', 'corr.basis')){
        l[[i]] = mse(true[['corr']], pred[[i]])
    }
    for (i in c('pcor', 'pcor.shrink', 'pcor.basis')){
        l[[i]] = mse(true[['pcor']], pred[[i]])
    }
    return(l)
}
mse <- function(true, pred){
    true = true[lower.tri(true)]
    pred = pred[lower.tri(pred)]
    # return( mean((true-pred)^2) )
    return(sqrt(sum((true-pred)^2)))
}
compute_mse_boxplot <- function(count, ncell, nsamp, ngene, transf='clr'){
    df = data.frame(mse=numeric(), method=character())
    for (i in 1:nsamp){

        # compute ground truth
        gsamp  = sample(1:ncol(count), ngene)
        count  = count[,gsamp]
        bcor   = bShrink(count, 'rco', transf)

        # compute preds for less cells
        csamp  = sample(1:nrow(count), ncell)
        samp   = count[csamp,]
        bcor_  = bShrink(samp, 'rco', transf)

        # update data frame
        tmp   = calculate_mse(bcor, bcor_)
        for (method in names(tmp)){
            row = data.frame(mse=tmp[[method]], corr=get_cor_name(method), shrink=get_shrink_name(method))
            df = rbind(df, row)
        }
    }
    df$shrink = factor(df$shrink, levels=c('no shrinkage', 'shrink std', 'shrink basis'))
    # covariance
    tmp = df[which(df$corr == 'covariance'),]
    g1  = ggplot(tmp, aes(x=shrink, y=mse)) + geom_boxplot() + xlab('') + ylab('Mean squared error') + ggtitle(paste0('Covariance matrix, D=', ngene, ', N=', ncell))
    # correlation
    tmp = df[which(df$corr == 'correlation'),]
    g2  = ggplot(tmp, aes(x=shrink, y=mse)) + geom_boxplot() + xlab('') + ylab('Mean squared error') + ggtitle(paste0('Correlation matrix, D=', ngene, ', N=', ncell))
    # partial correlation
    tmp = df[which(df$corr == 'partial correlation'),]
    g3  = ggplot(tmp, aes(x=shrink, y=mse)) + geom_boxplot() + xlab('') + ylab('Mean squared error') + ggtitle(paste0('Partial correlation matrix, D=', ngene, ', N=', ncell))
    g   = ggarrange(g1, g2, g3, nrow=1)
    g
}
get_cor_name <- function(name){
    if (name %in% c('cov', 'cov.shrink', 'cov.basis')){
        'covariance'
    }else if (name %in% c('corr', 'corr.shrink', 'corr.basis')){
        'correlation'
    }else if (name %in% c('pcor', 'pcor.shrink', 'pcor.basis')){
        'partial correlation'
    }
}
get_shrink_name <- function(name){
    if (name %in% c('cov', 'corr', 'pcor')){
        'no shrinkage'
    }else if (name %in% c('cov.shrink', 'corr.shrink', 'pcor.shrink')){
        'shrink std'
    }else if (name %in% c('cov.basis', 'corr.basis', 'pcor.basis')){
        'shrink basis'
    }
}

# benchmark -----------------------------------------------------------------------------

for (i in parser$ncell){
    message("benchmarking for n=", i)
    g   = compute_mse_boxplot(count, ncell=i, nsamp=parser$nsamp, ngene=parser$ngene, transf=parser$transf)
    out = paste0(parser$outdir, '/benchmark_pcor_shrink_', i, '_', parser$ngene, '_', parser$transf, '.png')
    ggsave(out, plot=g, width=12, height=3)
}
