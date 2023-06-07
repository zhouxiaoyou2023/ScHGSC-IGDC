#导入单细胞表达数据
xx <- read.table("sc_expression.txt")
a <- substr(colnames(xx),1,7)
b <- substr(colnames(xx),1,4)

library(dplyr)
library(Seurat)
library(patchwork)

sce <- CreateSeuratObject(counts = xx)
sce <- NormalizeData(sce)%>%
  FindVariableFeatures()%>%
  ScaleData()%>%
  RunPCA()

library(harmony)
sce <- RunHarmony(sce, "orig.ident",max.iter.harmony = 20)
ElbowPlot(sce,ndims = 50,reduction = "harmony")

sce <- RunUMAP(sce,dims = 1:30,reduction = "harmony") %>%
  FindNeighbors(dims = 1:30,reduction = "harmony") %>%
  FindClusters(resolution =1)

DimPlot(sce,label = T)

meta <-data.frame(celltype=substr(colnames(xx),6,7))
meta <- ifelse(meta$celltype=="FT","Normal","Tumor")
names(meta) <- colnames(sce)
sce <- AddMetaData(sce,metadata = meta,col.name = "Celltype")

DimPlot(sce,label = T,group.by = "Celltype")
save(sce,file = 'sce.Rdata')









