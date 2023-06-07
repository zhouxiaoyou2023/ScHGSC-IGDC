library(DiffCorr)
#Select expression of core genes
input <- seurat_obj@misc$tutorial$datExpr 
core <- read.table("core_gene.txt",header=T)
name=core[,1]
mut = colnamess(input) %in% name
new=input[,mut]
new1=as.data.frame(t(new))

exp1<- new1[, 1:36]
exp2<- new1[, 37:179]

index<- (apply(exp1,1,sd)!=0)&(apply(exp2,1,sd)!=0)
sum(index)

exp1<- exp1[index,]
exp2<- exp2[index,]

## Clusters on each subset
hc.mol1 <- cluster.molecule(exp1, 'pearson', 'average')  
hc.mol2 <- cluster.molecule(exp2, 'pearson', 'average') 

g1 <- cutree(hc.mol1, k=16)
g2 <- cutree(hc.mol2, k=16)#k设置为1/20的基因数

res1 <- get.eigen.molecule(data = exp1, groups = g1)
res2 <- get.eigen.molecule(data = exp2, groups = g2)

#Visualizing module networks
gg1 <- get.eigen.molecule.graph(res1)
write.modules(g1, res1, outfile="normal_modules.txt")

gg2 <- get.eigen.molecule.graph(res2)
write.modules(g2, res2, outfile="tumor_modules.txt")

pdf('module_network.pdf',width=10,height=6)
par(mfrow=c(1,2))
plot(gg1, layout=layout.fruchterman.reingold(gg1), main='Normal')
plot(gg2, layout=layout.fruchterman.reingold(gg2), main='Tumor')
dev.off()

save(hc.mol1,hc.mol2,g1,g2,res1,res2,file='GSE.RData')

#Export the results (FDR < 0.05)
comp.2.cc.fdr(output.file="res.txt", exp1, exp2, method='pearson', p.adjust.methods='BH', threshold=0.001)
