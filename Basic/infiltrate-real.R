library(Seurat)
library(rmarkdown) 
library(tidyverse)
getwd()
setwd('E:/R/glioma-infiltrate')

library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)
library(monocle)

load('cancer_ana.Rdata')

load('zs-infil-anno2.Rdata')


load('zs-infil.Rdata')

load('zs-infil-ori.Rdata')




#ģ��
load('mofang-qc.Rdata')
#���ݱ���
load('infiltrate-qc.Rdata')
load('newqc.Rdata')

#old
load('seurat-infiltrate-wanmei.Rdata')

save(a,am,am2,am3,qcmeta,qcdata,qcdata2,
     newmetadata,
     outdata,
     file = 'zs-infil-ori.Rdata')

save(am3,qcmeta,qcdata2,
     pbmc,pbmc1,pbmc2,pbmc3,
     pbmc.markers,
     file = 'zs-infil.Rdata')

save(am3,
     pbmc2,pbmc3,
     pbmc.markers,
     clu_ann,clu_markers,
     dui3,
     cancer,
     S1.markers,S2.markers,S4.markers,S6.markers,
     core.markers,Periphery.markers,
     file = 'zs-infil-anno2.Rdata')

save(pbmc3,
     clu_ann,clu_markers,pbmc.markers,
     cancer,cancer.markers,
     S1.markers,S2.markers,S4.markers,S6.markers,
     core.markers,Periphery.markers,
     C0.markers,C1.markers,C2.markers,C3.markers,C4.markers,C5.markers,C6.markers,C7.markers,
     #cds ��cancer
     cds,
     can_state1.markers,can_state2.markers,can_state3.markers,can_state4.markers,can_state5.markers,
     can_state1a.markers,can_state2a.markers,can_state3a.markers,can_state4a.markers,can_state5a.markers,
     #cds_all
     cds_all,
     state2.markers,state3.markers,state4.markers,
     #patient 2+4
     p24,p24_cds,
     p24_9_state9.markers,p24_9_state3.markers,
     file = 'cancer_ana.Rdata')



###################################################################################################################
#Ŀǰqc
#ȥ������23465 -> 19468
#am ����meta��Ϣ
#ȥ��ϸ����

##################################################################################################################

## a 23465  3589
## dat 19456  3589
#�������ϸ����С��3�Ļ���
dat=a[apply(a,1, function(x) sum(x>1) > 3),] 
dim(dat)

############################################################################################
#����meta
library(openxlsx)
newmetadata = read.xlsx('infiltrate-metadata.xlsx')
newmetadata2 <- t(newmetadata)

colnames(newmetadata2) <- colnames(a)
am <- rbind(dat,newmetadata2)


############################################################################################

#���������С��50��ϸ������
am2 <- t(am)
dat=am2[apply(am2,1, function(x) sum(x>1) > 65),] 
outdat=am2[apply(am2,1, function(x) sum(x>1) < 65),] 
outdata <- t(outdat)


dim(dat)
##3580 19468

am3 <- t(dat)
#am3Ϊȫ��������data��meta
qcmeta <- am3[19457:19468,]

qcdata <- am3[1:19456,]

rownamelinshi <- rownames(qcdata)
qcdata2=apply(qcdata,2,as.numeric)
rownames(qcdata2) <-rownamelinshi


############################################################################################
############################################################################################
#��׼����
pbmc.data <- qcdata2
dim(pbmc.data)
pbmc.data <- log2(1 + sweep(pbmc.data, 2, median(colSums(pbmc.data))/colSums(pbmc.data), '*'))
head(colnames(pbmc.data))
pbmc<- CreateSeuratObject(counts = pbmc.data, project = "10X pbmc", min.cells = 1, min.features = 0)
pbmc
head(pbmc@meta.data)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc@meta.data %>% head()

VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#plot3 <-FeatureScatter(pbmc,feature1 = 'nCount_RNA',feature2 = 'percent.mt')
#plot4 <-FeatureScatter(pbmc,feature1 = 'nCount_RNA',feature2 = 'nFeature_RNA')
#CombinePlots(plots = list(plot3,plot3))

#���mtRNAС��5%
pbmc_mt <- subset(pbmc,subset = percent.mt<5)


pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(pbmc), 10)
top10


plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2), ncol =1)


all.genes <- rownames(pbmc)
pbmc1 <- ScaleData(pbmc, features = all.genes)
## ���ǰ�����ݵ�����
#### raw counts, same as pbmc@assays$RNA@counts[1:6, 1:6]
pbmc1[["RNA"]]@counts[1:6, 1:6]
### library size normalized and log transformed data
pbmc1[["RNA"]]@data[1:6, 1:6]
### scaled data
pbmc1[["RNA"]]@scale.data[1:6, 1:6]
#��û��������һ�����Ҫ����Ϊ��Ҫ��PCAǰ�� ������DNAռ����һ��ȥ����������Ǽ�����PCA���
pbmc1 <- ScaleData(pbmc1, vars.to.regress = "percent.mt")
pbmc1 <- RunPCA(pbmc1, features = VariableFeatures(object = pbmc1), verbose = FALSE)
DimPlot(pbmc1, reduction = "pca")


#�鿴PC��Pֵ���ж�PC����
pbmc2 <- JackStraw(pbmc1, num.replicate = 100, dims = 50)
pbmc2 <- ScoreJackStraw(pbmc2, dims = 1:50)
JackStrawPlot(pbmc2, dims = 1:20)
ElbowPlot(pbmc2, ndims = 50)


pbmc3 <- FindNeighbors(pbmc2, dims = 1:20)

a
# Look at cluster IDs of the first 5 cells
head(Idents(pbmc3), 5)
#�����Խ�ά
pbmc3 <- RunUMAP(pbmc3, dims = 1:20)
pbmc3<- RunTSNE(pbmc3, dims = 1:20)
## after we run UMAP and TSNE, there are more entries in the reduction slot
str(pbmc3@reductions)
DimPlot(pbmc3, reduction = "umap", label = TRUE)
DimPlot(pbmc3, reduction = "tsne", label = TRUE)
DimPlot(pbmc3, reduction = "pca")

pbmc3@meta.data$sccelltype[which(pbmc3@meta.data$sccelltype=='MacrophageAndNeuron')] <-'Macrophage'
pbmc3@meta.data$sccelltype[which(pbmc3@meta.data$sccelltype=='cancer_cell')] <-'Cancer_cell'
pbmc3@meta.data$sccelltype[which(pbmc3@meta.data$sccelltype=='other_immune_cell')] <-'Other_immune_cell'
Idents(object = pbmc3) <- 'sccelltype'

#p = 
  DimPlot(pbmc3, reduction = "tsne", label = TRUE,
              label.size = 4,pt.size = 1)
ggsave('tsne_type.jpg',plot = p,dpi = 300,width = 15,height = 7.5,units = 'cm',
       path = 'C:/Users/Administrator/Desktop')


Idents(object = pbmc3) <- 'sccelltype'

table(pbmc3@meta.data$sccelltype)


Idents(object = pbmc3) <- 'seurat_clusters'
DimPlot(pbmc3, reduction = "umap", label = TRUE,pt.size = 0.7)
DimPlot(pbmc3, reduction = "tsne", label = TRUE,pt.size = 0.7)
DimPlot(pbmc3, reduction = "pca")





Idents(object = pbmc3) <- 'characteristics:.Patient.ID'
  p1 = DimPlot(pbmc3, reduction = "tsne", label = F,pt.size = 2,label.size = 8)
Idents(object = pbmc3) <- 'sccelltype'
p2 = DimPlot(pbmc3, reduction = "tsne", label = F,pt.size = 2,label.size = 8)
Idents(object = pbmc3) <- 'characteristics:.Tissue'
p3 = DimPlot(pbmc3, reduction = "tsne", label = F,pt.size = 2,label.size = 8)

jpeg(file='test.jpeg', height=1000, width=1000,bg = "white",res = 600)
p1
dev.off()


Idents(object = cancer) <- 'characteristics:.Tissue'
DimPlot(cancer, reduction = "tsne", label = F,pt.size = 2.5,label.size = 8)










# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc3, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
# ��log2FC
VlnPlot(pbmc3, features = c("GFAP"))
VlnPlot(pbmc3, features = c("LIF", 'LIFR', "IL6ST",'OSM','OSMR','TGFb'), slot = "counts", log = TRUE)
# MS4A1��CD79A ��Bϸ����־
# �ɶ������Ա�

FeaturePlot(pbmc3,features = c("GFAP"))
FeaturePlot(pbmc3, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))
# ����ֲ�ͼ + �����ֲ�ƴͼ

top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(pbmc3, features = top10$gene) + NoLegend()
# top10��ͼ

head(colnames(pbmc.data))

pbmc3 <- FindNeighbors(pbmc2, dims = 1:20)
for (res in c(0.05, 0.1, 0.2, 0.3,0.4, 0.5,0.6,0.8,1)) {
        pbmc3 <- FindClusters(pbmc3, graph.name = "RNA_snn", resolution = res, algorithm = 1)
}
library(clustree)
clustree(pbmc3@meta.data, prefix = "RNA_snn_res.") 







#########################################################################################
#�ڶ�����meta
qcmeta <- t(qcmeta)
pbmc3@meta.data <-cbind(pbmc3@meta.data,qcmeta)
view(pbmc3@meta.data)
#����plate��Ϣ
plate2=str_split(pbmc3@meta.data[,11],'000',simplify = T)[,2] 
plate=str_split(plate2,'[.]',simplify = T)[,1] 
table(plate2)
pbmc3@meta.data <-cbind(pbmc3@meta.data,plate)


#��ͼ������Ⱥ
library(ggplot2) 
genes_to_check = c('PTPRC', 'CD3D', 'CD3E','CD4', 'CD8A', ## Tcells
                   'CD68', 'CD163', 'CD14',  'CD86', 'CCL22','S100A4', ## DC(belong to monocyte)
                   'CD68',  'CD163','MRC1','MSR1' ,'S100A8', ## Macrophage (belong to monocyte)
                   'COL1A1','ACTA2','PECAM1', 'VWF' ,'PROX1', ## Fibroblasts,Endothelial
                   'EGFR','PROM1','DCN', ## epi or tumor
                   'MOG','AGXT2L1','GPR17','STMN2'
)
library(stringr)  
p <- DotPlot(pbmc3, features = unique(genes_to_check),
             assay='RNA' ,group.by = 'seurat_clusters' )  + coord_flip()

p
ggsave(plot=p, filename="check_marker_by_seurat_cluster.png")
################################################################################################
#markers��Ϣע��
annotations <- read.csv("annotation.csv")
test=pbmc.markers
test=data.frame(test)
test=  rownames_to_column(test)
dif=test[,4]-test[,5]
dif=data.frame(dif)
test1=cbind(test,dif)

testwithdisc <- test1 %>% 
        left_join(y = unique(annotations[, c("gene_name", "description")]),
                  by = c("rowname" = "gene_name"))
pbmc.markers<-testwithdisc

#####################################################################################################
#ϸ������
table(pbmc3@meta.data$seurat_clusters)
Cell_Num=data.frame(table(pbmc3@meta.data$seurat_clusters))
View(Cell_Num)
library(ggplot2)

mytheme <- theme(plot.title=element_text(face="bold", size="14", color="black"), #ָ��ͼ�ı���Ӧ��Ϊ��б����ɫ14��
                 axis.title=element_text(face="bold", size=10, color="black"),#��ı���Ϊ��б�����ɫ10
                 axis.text=element_text(face="bold", size=9, color="black"),#���ǩΪ���������ɫ9��
                 panel.background=element_rect(fill="white",color="darkblue"),#ͼƬ�����а�ɫ����������ɫ�ı߿�
                 panel.grid.major.y=element_line(color="grey",linetype=1),#��ˮƽ����Ӧ���ǻ�ɫ��ʵ��
                 panel.grid.minor.y=element_line(color="grey", linetype=2),#��ˮƽ����Ӧ���ǻ�ɫ������
                 panel.grid.minor.x=element_blank(), #��ֱ�������
                 legend.position="top",#ͼ��չʾ�ڶ���
                 axis.text.x=element_text(angle=50,hjust=0.5, vjust=0.5))

ggplot(data = Cell_Num, mapping = aes(x = Cell_Num$Var1, y = Cell_Num$Freq, fill = Cell_Num$Var1)) + geom_bar(stat = 'identity', position = 'identity')  + xlab('Cell Types')+geom_text(mapping = aes(label = Cell_Num$Freq))+mytheme


###################################################################################################################
##singleRע��
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
write.table(result_main_hpca, file = "HPCA_Main.txt", sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE) #���������������Ժ����

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

########################################################################################
#����
#����
#����
#����
#����
#����
########################################################################################
###################################################################################################################
#����12��
linshi <-pbmc3[,pbmc3@meta.data$seurat_clusters%in%c(12)]
view(linshi@meta.data)


ng=apply(pbmc.data,2,function(x)sum(x>1))
table(ng)
ng=data.frame(ng)

pbmctest[["ng.ident"]] <-   ng


###################################################################################################################
#����
##################################################################################################################
#��ԭԭ��QC am1 ȡ�������ĵ�ϸ��
am1<- t(am)
am2 = subset(am1,am1[,19463] == 'Tumor') #�����ڲ�ϸ����Ϊ2343,������3589
am1 <- t(am2)
rm(am2)
############################################################################################
#��qc
#amm ����meta��Ϣ,����meta������ϸ��+����
amm <- rbind(a,newmetadata2)

#amm1 ȡ�����ڲ�ϸ��
amm1<- t(amm)
amm2ls = subset(amm1,amm1[,23472] == 'Tumor') #�����ڲ�ϸ����Ϊ2343,������3589
amm1 <- t(amm2ls)

amm2ls <- t(amm1)

#ȥ��ϸ����2353 -> 2334
amm2=amm2ls[apply(amm2ls,1, function(x) sum(x>1) > 100),] 
dim(amm2)

#ȥ������23465 -> 19468

save(a,
     newmetadata,
     amm,amm1,
     file = 'newqc.Rdata')
############################################################################################################

###############С��ĳ���

library(SingleR)
pbmc3_for_SingleR <- GetAssayData(pbmc3, slot="data")
clusters=pbmc3@meta.data$seurat_clusters
mouseImmu <- ImmGenData()
pred.mouseImmu <- SingleR(test = pbmc3_for_SingleR, ref = mouseImmu, labels = mouseImmu$label.main,
                          method = "cluster", clusters = clusters, 
                          assay.type.test = "logcounts", assay.type.ref = "logcounts")

mouseRNA <- MouseRNAseqData()
pred.mouseRNA <- SingleR(test = pbmc3_for_SingleR, ref = mouseRNA, labels = mouseRNA$label.fine ,
                         method = "cluster", clusters = clusters, 
                         assay.type.test = "logcounts", assay.type.ref = "logcounts")

cellType=data.frame(ClusterID=levels(pbmc3@meta.data$seurat_clusters),
                    mouseImmu=pred.mouseImmu$labels,
                    mouseRNA=pred.mouseRNA$labels )
head(cellType)
pbmc3@meta.data$singleR=cellType[match(clusters,cellType$ClusterID),'mouseRNA']

pro='first_anno'
DimPlot(pbmc3,reduction = "umap",label=T, group.by = 'singleR')
ggplot2::ggsave(filename = paste0(pro,'_tSNE_singleR.pdf'))
DimPlot(pbmc3,reduction = "umap",label=T,split.by ='orig.ident',group.by = 'singleR')
ggplot2::ggsave(filename = paste0(pro,'_tSNE_singleR_by_orig.ident.pdf'))

save(pbmc3,file = 'last_pbmc3.Rdata') 



########################################################################################################################
#infercnv

install.packages("rjags")
if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
BiocManager::install("infercnv")

library(infercnv)

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=system.file("extdata", "oligodendroglioma_expression_downsampled.counts.matrix.gz", package = "infercnv"),
                                    annotations_file=system.file("extdata", "oligodendroglioma_annotations_downsampled.txt", package = "infercnv"),
                                    delim="\t",
                                    gene_order_file=system.file("extdata", "gencode_downsampled.EXAMPLE_ONLY_DONT_REUSE.txt", package = "infercnv"),
                                    ref_group_names=c("Microglia/Macrophage","Oligodendrocytes (non-malignant)")) 

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=tempfile(), 
                             cluster_by_groups=TRUE, 
                             denoise=TRUE,
                             HMM=TRUE)

########################################################################################################################
#ע��
Idents(object = pbmc3) <- 'seurat_clusters'
Idents(object = pbmc3) <- 'characteristics:.Patient.ID'


library(scCATCH)
clu_markers <- findmarkergenes(object = pbmc3,
                               species = 'Human',
                               cluster = 'All',
                               match_CellMatch = TRUE,
                               cancer = c("Glioma","Glioblastoma"),
                               tissue = c("Brain","Blood","Blood vessel"),
                               cell_min_pct = 0.25,
                               logfc = 0.25,
                               pvalue = 0.05)

clu_ann <- scCATCH(clu_markers$clu_markers,
                   species = "Human",
                   cancer = c("Glioma","Glioblastoma"),
                   tissue = c("Brain","Blood","Blood vessel"))

scnew.cluster.ids <- clu_ann$cell_type
names(scnew.cluster.ids) <- clu_ann$cluster
pbmc3 <- RenameIdents(pbmc3, scnew.cluster.ids)    #
##################Idents(object = pbmc3) <- "seurat_clusters"
DimPlot(pbmc3, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
##########################################################################################
Macrophage=c(0)
Microglial_Cell=c(2)
cancer_cell=c(1,5,6,7)
other_immune_cell =c(8)
Astrocyte=c(3,10,11)
Oligodendrocyte_Precursor_Cell=c(4)
Oligodendrocyte=c(9)
MacrophageAndNeuron=c(12)

current.cluster.ids <- c(MacrophageAndNeuron,
                         Astrocyte,
                         Microglial_Cell,
                         cancer_cell,
                         Oligodendrocyte_Precursor_Cell,
                         Oligodendrocyte,
                         other_immune_cell,
                         Macrophage)

new.cluster.ids <- c(
                     rep("MacrophageAndNeuron",length(MacrophageAndNeuron)),
                     rep("Astrocyte",length(Astrocyte)),
                     rep("Microglial_Cell",length(Microglial_Cell)),
                     rep("cancer_cell",length(cancer_cell)),
                     rep("Oligodendrocyte_Precursor_Cell",length(Oligodendrocyte_Precursor_Cell)),
                     rep("Oligodendrocyte",length(Oligodendrocyte)),
                     rep("other_immune_cell",length(other_immune_cell)),
                     rep("Macrophage",length(Macrophage))
                     )

pbmc3@meta.data$sccelltype <- plyr::mapvalues(x = as.integer(as.character(pbmc3@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)

table(pbmc3@meta.data$sccelltype)

####ע�͵ıȶ�
table(pbmc3@meta.data$`characteristics:.cell.type`)


meta = pbmc3@meta.data

table(meta[,'sccelltype'])
meta[,"sccelltype"][meta[,"sccelltype"] == 'cancer_cell'] = 'Neoplastic'
meta[,"sccelltype"][meta[,"sccelltype"] == 'Macrophage'] = 'Immune cell'
meta[,"sccelltype"][meta[,"sccelltype"] == 'Microglial_Cell'] = 'Immune cell'
meta[,"sccelltype"][meta[,"sccelltype"] == 'other_immune_cell'] = 'Immune cell'



compare = meta[,c(16,19)]
colnames(compare)[1] = "old"
colnames(compare)[2] = "new"

cancercells  = compare[,2][compare[,1]=="Neoplastic"]
cancercells=data.frame(cancercells)
1085/1119

immune =data.frame( compare[,2][compare[,1]=="Immune cell"])
table(immune[,1])
1385/1389
############################################################################################################
#���������Ժͻ���������
#ͳ��ϸ������
meta3 = meta[,c(13,14,19)]
meta3cancer =data.frame( meta3[,1][meta3[,3] == 'Neoplastic'])
table(meta3cancer[,1])
#79 1040
meta3immune =data.frame( meta3[,1][meta3[,3] == "Immune cell"])
table(meta3immune[,1])
#640  749 

meta3astrocyte =data.frame( meta3[,1][meta3[,3] == "Astrocyte"])
table(meta3astrocyte[,1])
#640  749 

##################
##��������ϸ������
cancer<-pbmc3[,pbmc3@meta.data$`sccelltype` %in%"cancer_cell"]

Idents(object = cancer) <- 'characteristics:.Tissue'
Idents(object = cancer) <- 'characteristics:.Patient.ID'

DimPlot(cancer, reduction = "tsne", label = F,pt.size = 2.5,label.size = 8)










pbmc3@meta.data = meta
immunecells<-pbmc3[,pbmc3@meta.data$`sccelltype` %in%"Immune cell"]



Idents(object = immunecells) <- 'characteristics:.Tissue'
Idents(object = immunecells) <- 'characteristics:.Patient.ID'

DimPlot(immunecells, reduction = "tsne", label = F,pt.size = 2.5,label.size = 8)











##########################################################################################################

Idents(object = pbmc3) <- 'sccelltype'
Idents(object = pbmc3) <- 'tissue'

DimPlot(pbmc3,reduction = "umap",label=T)

view(pbmc3@meta.data)

##pbmc3@meta.data = pbmc3@meta.data[,-19]



###################################################################################################################
#�����ض�ϸ������
cancer<-pbmc3[,pbmc3@meta.data$`sccelltype` %in%"cancer_cell"]
cancer<-pbmc3[,pbmc3@meta.data$`all_celltype` %in%"tumor_cells"]
Idents(object = cancer) <- 'position'
VlnPlot(cancer, features = c("LIF", 'LIFR', "IL6ST",'OSM','OSMR','TGFb'), slot = "counts", log = TRUE)


view(cancer@meta.data)

Idents(object = cancer) <- 'characteristics:.Tissue'

#���ܺ�����
core.markers <- FindMarkers(cancer, ident.1 = 'Tumor', ident.2 = 'Periphery',min.pct = 0.25)
Periphery.markers<- FindMarkers(cancer, ident.1 = 'Periphery', ident.2 = 'Tumor',min.pct = 0.25)

#���߼�
Idents(object = cancer) <- 'characteristics:.Patient.ID'
S1.markers <- FindMarkers(cancer, ident.1 = 'BT_S1', min.pct = 0.25)
S2.markers <- FindMarkers(cancer, ident.1 = 'BT_S2', min.pct = 0.25)
S4.markers <- FindMarkers(cancer, ident.1 = 'BT_S4', min.pct = 0.25)
S6.markers <- FindMarkers(cancer, ident.1 = 'BT_S6', min.pct = 0.25)

#���ӻ���ע��
annotations <- read.csv("annotation.csv")

test=C5.markers
test=data.frame(test)
test=  rownames_to_column(test)
dif=test[,4]-test[,5]
dif=data.frame(dif)
test1=cbind(test,dif)

testwithdisc <- test1 %>% 
        left_join(y = unique(annotations[, c("gene_name", "description")]),
                  by = c("rowname" = "gene_name"))
C5.markers<-testwithdisc



# Stash cell identity classes ��cluster��Ϣ
cancer[["cancer.ident"]] <-    Idents(object = cancer)
cancer <- StashIdent(object = cancer, save.name = "cancer.ident")

#����˫�ط��飬�ظ�tissue + patient.id == pat_tis
Idents(object = cancer) <- "characteristics:.Patient.ID"
cancer[["linshi"]]<-    Idents(object = cancer)
linshi<- cancer[["linshi"]]

Idents(object = cancer) <- "characteristics:.Tissue"
cancer[["linshi2"]]<-    Idents(object = cancer)
linshi2<- cancer[["linshi2"]]

linshi3 <- cbind(linshi,linshi2)
linshi3 <- data.frame(linshi3)
linshi4 <- tidyr::unite(linshi3, "pat_tis.ident", sep = "_")

cancer[["pat_tis.ident"]] <-   linshi4
Idents(object = cancer) <- "pat_tis.ident"
cancer <- StashIdent(object = cancer, save.name = "pat_tis.ident")

table(cancer@meta.data$pat_tis.ident)






















