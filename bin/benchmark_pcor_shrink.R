#!/usr/bin/env Rscript

library(argparse)
library(corpcor)
library(data.table)
library(ggplot2)

# parse arguments
parser = ArgumentParser(description='Benchmark shrinkage')
parser$add_argument('-i', '--input', type='character', help="Input count data")
parser$add_argument('-o', '--outdir', type='character', help="Output directory")
parser$add_argument('--transf', type='character', help='clr or alr')
parser = parser$parse_args()


# read input files ----------------------------------------------------------------------

count = fread(parser$input)
count = as.matrix(count)

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
        Cov.improved=F%*%Cov.shrink%*%t(F)
    } else if (outtype=="clr"){
        G=diag(rep(1,D))-matrix(1/D,D,D)
        Cov.improved=G%*%Cov.shrink%*%G
    } else{
        die("outtype unsupported")
    }
    Cor=cov2cor(Cov)
    PC=cor2pcor(Cov)
    Cor.shrink=cov2cor(Cov.shrink)
    PC.shrink=cor2pcor(Cov.shrink)
    Cor.improved=cov2cor(Cov.improved)
    PC.improved=cor2pcor(Cov.improved)
    
    
    C=list(Cov,Cov.shrink,Cov.improved,Cor,Cor.shrink,Cor.improved,PC,PC.shrink,PC.improved)
    names(C)=c('cov', 'cov.shrink', 'cov.improved', 'corr', 'corr.shrink', 'corr.improved', 'pcor', 'pcor.shrink', 'pcor.improved')
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
    for (i in c('cov', 'cov.shrink', 'cov.improved')){
        l[[i]] = mse(true[['cov']], pred[[i]])
    }
    for (i in c('corr', 'corr.shrink', 'corr.improved')){
        l[[i]] = mse(true[['corr']], pred[[i]])
    }
    for (i in c('pcor', 'pcor.shrink', 'pcor.improved')){
        l[[i]] = mse(true[['pcor']], pred[[i]])
    }
    return(l)
}
mse <- function(true, pred){
    true = true[lower.tri(true)]
    pred = pred[lower.tri(pred)]
    return( mean((true-pred)^2) )
}
check_mse_boxplot <- function(count, bcor, ncell, nsamp, transf='clr'){
    df = data.frame(mse=numeric(), method=character())
    for (i in 1:nsamp){
        samp  = sample(1:nrow(count), ncell)
        samp  = count[samp,]
        bcor_ = bShrink(samp, 'rco', transf)
        tmp   = calculate_mse(bcor, bcor_)
        for (method in names(tmp)){
            row = data.frame(mse=tmp[[method]], method=method)
            df = rbind(df, row)
        }
    }
    g = ggplot(df, aes(x=method, y=mse)) + geom_boxplot() + ggtitle(paste0('N = ', ncell))
    g
}

# benchmark -----------------------------------------------------------------------------

# calculate ground truth
message("calculating ground truth")
bcor = bShrink(count, 'rco', parser$transf)

# benchmark
for (i in c(2000, 1000, 500, 200, 100, 50)){
    message("benchmarking for n=", i)
    g = check_mse_boxplot(count, bcor, ncell=i, nsamp=200, transf=parser$transf)
    out = paste0(parser$outdir, '/benchmark_pcor_shrink_', i, '_', parser$transf, '.png')
    ggsave(out, plot=g)
}
