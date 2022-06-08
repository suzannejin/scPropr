M=read.csv("../data/mouseStemCells.csv")
#find the nonzero (support) set for each gene individually:
supp=list()
d=rep(0,ncol(M))
for (j in 1:ncol(M)){
    s=which(M[,j]>0)
    d[j]=length(s)
    supp[[j]]=s
}
ind=order(d,decreasing=TRUE)
plot(d[ind],xlab="number of genes",ylab="number of nonzero samples",main="reverse cumulative of nonzero samples")
#find a set of genes that has *common* nonzero samples
#(starting from the genes with most nonzero samples):
first=d[ind[1]]
mysup=supp[[first]]
csup=list()
dc=rep(0,ncol(M))
c=0
for (i in ind){
    c=c+1
    #print(c)
    mysup=intersect(mysup,supp[[i]])
    csup[[c]]=mysup
    dc[c]=length(mysup)
}
plot(dc,xlab="number of genes used",ylab="number of joint nonzero samples",xlim=c(0,2000))
#as an example, choose a cut-off:
abline(v=1000,col="red")
dc[1000]
[1] 1510
#1510 samples left
mygenes=ind[1:1000]
#new nonzero data matrix:
m=M[csup[[1000]],mygenes]
