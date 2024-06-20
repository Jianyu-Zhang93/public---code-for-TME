library(Seurat)
library(rmarkdown)
library(tidyverse)
getwd()
setwd('E:/R/glioma-infiltrate')
load('seurat-infiltrate-wanmei.Rdata')

load('infiltrate-qc.Rdata')


load('seurat-infiltrate-metaplus-new.Rdata')

load("hpca.se.RData")

save(a,cluster1.markers,cluster2.markers,cluster3.markers,cluster4.markers,cluster5.markers,cluster6.markers,
     cluster7.markers,cluster8.markers,cluster9.markers,cluster10.markers,cluster11.markers,
     cluster12.markers,cluster13.markers,cluster14.markers,cluster15.markers,
     dat,df,hc,plate,metadata,pbmc,pbmc.data,
     pbmc.markers,pbmc1,pbmc2,pbmc3,plot1,plot2,top10,all.genes,
     test,test1,newmetadata,
     Periphery.markers,Peripheryxx.markers,core.markers,corexx.markers,
     pbmctest,S1.markers,S2.markers,S4.markers,S6.markers,patient_hb.markers,
     S1xx.markers,S2xx.markers,S4xx.markers,S6xx.markers,
     
     pbmc4,pbmc4.markers,pbmc4_top10,
     file = 'seurat-infiltrate-wanmei.Rdata')


save(a,
     pbmc4,pbmc4.markers,
     file = 'infiltrate-qc.Rdata')













pbmc4 <- FindNeighbors(pbmc2, dims = 1:20)
pbmc4 <- FindClusters(pbmc4, resolution = 0.3)
# Look at cluster IDs of the first 5 cells
head(Idents(pbmc4), 5)
#非线性降维
pbmc4 <- RunUMAP(pbmc4, dims = 1:20)
pbmc4<- RunTSNE(pbmc4, dims = 1:20)
## after we run UMAP and TSNE, there are more entries in the reduction slot
str(pbmc4@reductions)
DimPlot(pbmc3, reduction = "umap", label = TRUE)
DimPlot(pbmc3, reduction = "tsne")
DimPlot(pbmc4, reduction = "pca")

cluster1.markers <- FindMarkers(pbmc4, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)
# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc4, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)
# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc4.markers <- FindAllMarkers(pbmc4, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc4.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
# 是log2FC
VlnPlot(pbmc4, features = c("IDH1"))
VlnPlot(pbmc4, features = c("EGFR", "PTPRC"), slot = "counts", log = TRUE)
# MS4A1和CD79A 是B细胞标志
# 可多个基因对比

FeaturePlot(pbmc4,features = c("NOTCH1"))
FeaturePlot(pbmc4, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))
# 基因分布图 + 多基因分布拼图

pbmc4_top10 <- pbmc4.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# top10热图

head(colnames(pbmc.data))


##对比多个基因
cg=c('ABCB5','ALK','BRAF','MITF');cg
DotPlot(cancer,features = cg)
FeaturePlot(cancer,features = cg)
















#自己分组
FeaturePlot(pbmc4, features = c("EGFR", "PTPRC", "MOG", "AGXT2L1", "DCN", "GPR17", "STMN2", "CD14", "CD8A"))



#出图定义亚群
library(ggplot2) 
genes_to_check = c('PTPRC', 'CD3D', 'CD3E','CD4', 'CD8A', ## Tcells
                   'CD68', 'CD163', 'CD14',  'CD86', 'CCL22','S100A4', ## DC(belong to monocyte)
                   'CD68',  'CD163','MRC1','MSR1' ,'S100A8', ## Macrophage (belong to monocyte)
                   'COL1A1','ACTA2','PECAM1', 'VWF' ,'PROX1', ## Fibroblasts,Endothelial
                   'EGFR','PROM1','DCN', ## epi or tumor
                   'MOG','AGXT2L1','GPR17','STMN2'
)
library(stringr)  
p <- DotPlot(pbmc4, features = unique(genes_to_check),
             assay='RNA' ,group.by = 'seurat_clusters' )  + coord_flip()

p
ggsave(plot=p, filename="check_marker_by_seurat_cluster.png")


annotations <- read.csv("annotation.csv")
test=pbmc4.markers
test=data.frame(test)
test=  rownames_to_column(test)
dif=test[,4]-test[,5]
dif=data.frame(dif)
test1=cbind(test,dif)

testwithdisc <- test1 %>% 
        left_join(y = unique(annotations[, c("gene_name", "description")]),
                  by = c("rowname" = "gene_name"))
pbmc4.markers<-testwithdisc










B_cell=c(0)
Plasma_cell=c(4)
T_cell=c(1,2)
NK_cell=c(5)
Unknown=c(13)
Epithelial=c(3, 6, 7, 8, 10, 11, 12)
Endothelial=c(14)
Fibroblasts=c(9)
Doublet=c(15)

current.cluster.ids <- c(B_cell,
                         Plasma_cell,
                         T_cell,
                         NK_cell,
                         Unknown,
                         Epithelial,
                         Endothelial,
                         Fibroblasts,
                         Doublet)

new.cluster.ids <- c(rep("B_cell",length(B_cell)),
                     rep("Plasma_cell",length(Plasma_cell)),
                     rep("T_cell",length(T_cell)),
                     rep("NK_cell",length(NK_cell)),
                     rep("Unknown",length(Unknown)),
                     rep("Epithelial",length(Epithelial)),
                     rep("Endothelial",length(Endothelial)),
                     rep("Fibroblasts",length(Fibroblasts)),
                     rep("Doublet",length(Doublet))
)
















#细胞计数
table(pbmc4@meta.data$seurat_clusters)
Cell_Num=data.frame(table(pbmc3@meta.data$seurat_clusters))
View(Cell_Num)
library(ggplot2)

mytheme <- theme(plot.title=element_text(face="bold", size="14", color="black"), #指定图的标题应该为粗斜体棕色14号
                 axis.title=element_text(face="bold", size=10, color="black"),#轴的标题为粗斜体的棕色10
                 axis.text=element_text(face="bold", size=9, color="black"),#轴标签为粗体的深蓝色9号
                 panel.background=element_rect(fill="white",color="darkblue"),#图片区域有白色的填充和深蓝色的边框
                 panel.grid.major.y=element_line(color="grey",linetype=1),#主水平网格应该是灰色的实线
                 panel.grid.minor.y=element_line(color="grey", linetype=2),#次水平网格应该是灰色的虚线
                 panel.grid.minor.x=element_blank(), #垂直网格不输出
                 legend.position="top",#图例展示在顶部
                 axis.text.x=element_text(angle=50,hjust=0.5, vjust=0.5))

ggplot(data = Cell_Num, mapping = aes(x = Cell_Num$Var1, y = Cell_Num$Freq, fill = Cell_Num$Var1)) + geom_bar(stat = 'identity', position = 'identity')  + xlab('Cell Types')+geom_text(mapping = aes(label = Cell_Num$Freq))+mytheme





#添加meta

library(openxlsx)
newmetadata = read.xlsx('infiltrate-metadata.xlsx')

pbmc4@meta.data <-cbind(pbmc4@meta.data,newmetadata)

view(pbmc4@meta.data)

plate2=str_split(pbmc3@meta.data[,11],'000',simplify = T)[,2] 
plate=str_split(plate2,'[.]',simplify = T)[,1] 
table(plate2)
pbmc3@meta.data <-cbind(pbmc3@meta.data,plate)






#修改ID
# Get cell identity classes
Idents(object = pbmc3)
levels(x = pbmc3)

# Stash cell identity classes 存cluster信息
pbmc3[["old.ident"]] <- Idents(object = pbmc3)
pbmc3 <- StashIdent(object = pbmc3, save.name = "old.ident")

# Set identity classes
Idents(object = pbmc3) <- "CD4 T cells"
Idents(object = pbmc3, cells = 1:10) <- "CD4 T cells"


# Set identity classes to an existing column in meta data  
#Idents(object = pbmc3, cells = 1:10) <- "orig.ident"

#授予分组
Idents(object = pbmc3) <- "characteristics:.cell.type"
Idents(object = pbmc3) <- "old.ident"



#成图
DimPlot(pbmc3, reduction = "umap", label = TRUE)

#找差异
core.markers <- FindMarkers(pbmc3, ident.1 = 'Tumor', ident.2 = 'Periphery',min.pct = 0.25)
Periphery.markers<- FindMarkers(pbmc3, ident.1 = 'Periphery', ident.2 = 'Tumor',min.pct = 0.25)

#重命名
# Rename identity classes
pbmc3 <- RenameIdents(object = pbmc3, `CD4 T cells` = "T Helper cells")





#添加基因注释
annotations <- read.csv("annotation.csv")

test=S6.markers
test=data.frame(test)
test=  rownames_to_column(test)
dif=test[,4]-test[,5]
dif=data.frame(dif)
test1=cbind(test,dif)

testwithdisc <- test1 %>% 
        left_join(y = unique(annotations[, c("gene_name", "description")]),
                  by = c("rowname" = "gene_name"))
S6xx.markers<-testwithdisc






#test 成功分组
pbmctest <- pbmc3
#修改ID
# Get cell identity classes
Idents(object = pbmctest)
levels(x = pbmctest)

# Stash cell identity classes 存cluster信息
pbmctest[["test3.ident"]] <-    Idents(object = pbmctest)

pbmctest <- StashIdent(object = pbmctest, save.name = "test.ident")

# Set identity classes
#Idents(object = pbmctest) <- "CD4 T cells"
#Idents(object = pbmctest, cells = 1:10) <- "CD4 T cells"


# Set identity classes to an existing column in meta data  
#Idents(object = pbmctest, cells = 1:10) <- "orig.ident"

#授予分组
Idents(object = pbmctest) <- "characteristics:.Tissue"
Idents(object = pbmctest) <- "characteristics:.cell.type"
Idents(object = pbmctest) <- "old.ident"
head(Idents(object = pbmctest))

#授予双重分组，重复tissue + cell.type == test
Idents(object = pbmctest) <- "characteristics:.Patient.ID"
pbmctest[["linshi"]]<-    Idents(object = pbmctest)
linshi<- pbmctest[["linshi"]]

Idents(object = pbmctest) <- "characteristics:.cell.type"
pbmctest[["linshi2"]]<-    Idents(object = pbmctest)
linshi2<- pbmctest[["linshi2"]]


linshi3 <- cbind(linshi,linshi2)
linshi3 <- data.frame(linshi3)
linshi4 <- tidyr::unite(linshi3, "pat_tumorcell.ident", sep = "_")

pbmctest[["pat_tumorcell.ident"]] <-   linshi4
Idents(object = pbmctest) <- "pat_tumorcell.ident"
pbmctest <- StashIdent(object = pbmctest, save.name = "pat_tumorcell.ident")

S2hb.markers <- S2.markers
S2hb.markers$patient<-'2'

S6hb.markers <- S6.markers
S6hb.markers$patient<-'6'
#...
S6hb.markers=  rownames_to_column(S6hb.markers)
S2hb.markers=  rownames_to_column(S2hb.markers)
#...
patient_hb.markers<- rbind(S1hb.markers,S6hb.markers,S4hb.markers,S2hb.markers)

heatmap()



##
#成图
DimPlot(pbmctest, reduction = "umap", label = TRUE)

#找差异
S2.markers <- FindMarkers(pbmctest, ident.1 = 'BT_S2_Neoplastic', ident.2 = c('BT_S1_Neoplastic','BT_S4_Neoplastic','BT_S6_Neoplastic'),min.pct = 0.25)
S1.markers <- FindMarkers(pbmctest, ident.1 = 'BT_S1_Neoplastic', ident.2 = c('BT_S2_Neoplastic','BT_S4_Neoplastic','BT_S6_Neoplastic'),min.pct = 0.25)
S4.markers <- FindMarkers(pbmctest, ident.1 = 'BT_S4_Neoplastic', ident.2 = c('BT_S1_Neoplastic','BT_S2_Neoplastic','BT_S6_Neoplastic'),min.pct = 0.25)
S6.markers <- FindMarkers(pbmctest, ident.1 = 'BT_S6_Neoplastic', ident.2 = c('BT_S1_Neoplastic','BT_S4_Neoplastic','BT_S2_Neoplastic'),min.pct = 0.25)




Periphery.markers<- FindMarkers(pbmctest, ident.1 = 'Periphery_Neoplastic', ident.2 = 'Tumor_Neoplastic',min.pct = 0.25)

patient.markers <- FindAllMarkers(pbmctest, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#热图
coretop5 <- patient_hb.markers %>% group_by(patient) %>% top_n(n = 5, wt = avg_log2FC)
DoHeatmap(patient_hb.markers, features = coretop10$gene) + NoLegend()

##
rm(b)












#seurat 常规操作
PBMC <- CreateSeuratObject(pbmc.data, min.cells = 1, min.features = 0, project = '10x_PBMC')
PBMC <- AddMetaData(object = PBMC,   metadata = apply(pbmc.data, 2, sum), col.name = 'nUMI_raw')
PBMC <- AddMetaData(object = PBMC, metadata = timePoints, col.name = 'TimePoints')
PBMC <- ScaleData(object = PBMC, vars.to.regress = c('nUMI_raw'), model.use = 'linear', use.umi = FALSE)
PBMC <- FindVariableFeatures(object = PBMC, mean.function = ExpMean, dispersion.function = LogVMR, mean.cutoff = c(0.0125,3), dispersion.cutoff = c(0.5,Inf))
PBMC <- RunPCA(object = PBMC, pc.genes = PBMC@var.genes)
PBMC <- FindNeighbors(PBMC, reduction = "pca", dims = 1:10, k.param = 35)
PBMC <- FindClusters(object = PBMC, resolution = 0.9, verbose=F) 
PBMC <- RunTSNE(object = PBMC, dims.use = 1:10)
DimPlot(PBMC, cols = c('green4', 'pink', '#FF7F00', 'orchid', '#99c9fb', 'dodgerblue2', 'grey30', 'yellow', 'grey60', 'grey', 'red', '#FB9A99', 'green2','black'))





##singleR注释
library(celldex)
hpca.se <- HumanPrimaryCellAtlasData()
hpca.se
## class: SummarizedExperiment 
## dim: 19363 713 
## metadata(0):
## assays(1): logcounts
## rownames(19363): A1BG A1BG-AS1 ... ZZEF1 ZZZ3
## rowData names(0):
## colnames(713): GSM112490 GSM112491 ... GSM92233 GSM92234
## colData names(3): label.main label.fine label.ont

save(hpca.se,file="hpca.se.RData")



library(SingleR)
library(scater)
library(SummarizedExperiment)
test.count=as.data.frame(pbmc3[["RNA"]]@counts)

load(file="hpca.se.RData")
common_hpca <- intersect(rownames(test.count), rownames(hpca.se))
hpca.se <- hpca.se[common_hpca,]
test.count_forhpca <- test.count[common_hpca,]
test.count_forhpca.se <- SummarizedExperiment(assays=list(counts=test.count_forhpca))
test.count_forhpca.se <- logNormCounts(test.count_forhpca.se)


###main
pred.main.hpca <- SingleR(test = test.count_forhpca.se, ref = hpca.se, labels = hpca.se$label.main)
result_main_hpca <- as.data.frame(pred.main.hpca$labels)
result_main_hpca$CB <- rownames(pred.main.hpca)
colnames(result_main_hpca) <- c('HPCA_Main', 'CB')
write.table(result_main_hpca, file = "HPCA_Main.txt", sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE) #保存下来，方便以后调用

pbmc3b <-pbmc3

pbmc3@meta.data$CB=rownames(pbmc3@meta.data)
pbmc3@meta.data=merge(pbmc3@meta.data,result_main_hpca,by="CB")
rownames(pbmc3@meta.data)=pbmc3@meta.data$CB

p5 <- DimPlot(pbmc3, reduction = "umap", group.by = "HPCA_Main.y", pt.size=0.5)+theme(
        axis.line = element_blank(),
        axis.ticks = element_blank(),axis.text = element_blank()
)
p6 <- DimPlot(pbmc3, reduction = "umap", group.by = "ident",   pt.size=0.5, label = TRUE,repel = TRUE)+theme(
        axis.line = element_blank(),
        axis.ticks = element_blank(),axis.text = element_blank()
)
fig_tsne <- plot_grid(p6, p5, labels = c('ident','HPCA_Main'),rel_widths = c(2,3))
ggsave(filename = "tsne4.pdf", plot = fig_tsne, device = 'pdf', width = 36, height = 12, units = 'cm')



B_cell=c(0)
Plasma_cell=c(4)
T_cell=c(1,2)
NK_cell=c(5)
Unknown=c(13)
Epithelial=c(3, 6, 7, 8, 10, 11, 12)
Endothelial=c(14)
Fibroblasts=c(9)
Doublet=c(15)

current.cluster.ids <- c(B_cell,
                         Plasma_cell,
                         T_cell,
                         NK_cell,
                         Unknown,
                         Epithelial,
                         Endothelial,
                         Fibroblasts,
                         Doublet)

new.cluster.ids <- c(rep("B_cell",length(B_cell)),
                     rep("Plasma_cell",length(Plasma_cell)),
                     rep("T_cell",length(T_cell)),
                     rep("NK_cell",length(NK_cell)),
                     rep("Unknown",length(Unknown)),
                     rep("Epithelial",length(Epithelial)),
                     rep("Endothelial",length(Endothelial)),
                     rep("Fibroblasts",length(Fibroblasts)),
                     rep("Doublet",length(Doublet))
)

pbmc3@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(test.seu@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)

table(pbmc3@meta.data$celltype)

plotCB=as.data.frame(pbmc3@meta.data%>%filter(seurat_clusters!="13" &
                                                         seurat_clusters!="15"))[,"CB"]
DimPlot(pbmc3, reduction = "tsne", group.by = "celltype", pt.size=0.5,cells = plotCB)
ggsave(filename = "tsne5.pdf", width = 15, height = 12, units = 'cm')
saveRDS(pbmc3,file = "pbmc3.rds") #保存test.seu对象，下次可以直接调用，不用重复以上步骤


Idents(pbmc3)="celltype"
markers2 <- FindAllMarkers(pbmc3, logfc.threshold = 0.25, min.pct = 0.1, 
                           only.pos = TRUE)





#分析12组
linshi <-pbmc3[,pbmc3@meta.data$seurat_clusters%in%c(12)]
view(linshi@meta.data)


ng=apply(pbmc.data,2,function(x)sum(x>1))
table(ng)
ng=data.frame(ng)

pbmctest[["ng.ident"]] <-   ng





