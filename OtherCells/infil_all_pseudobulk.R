library(Seurat)
library(rmarkdown) 
library(tidyverse)
setwd('E:/R/glioma-infiltrate')

load('infil_all_pseudobulk.Rdata')


library(scater)
library(Seurat)
library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(apeglm)
library(DESeq2)
library(edgeR)
library(dplyr)
library(magrittr)
library(Matrix)
library(purrr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(SingleCellExperiment)
library(pheatmap)
library(png)
library(RColorBrewer)


load('zs-infil.Rdata')
load('cancer_ana.Rdata')

save(qcdata2,pbmc3,
     dds,res,res_tbl,sig_res,allmetadata,
     gathered_top20_sig,ei_meta,   res_table_thres,
     file='infil_all_pseudobulk.Rdata')
load('infil_all_pseudobulk.Rdata')




## dds
cluster_counts = qcdata2
#meta
colnames(pbmc3@meta.data)[13] = 'tissue'

allmetadata <- pbmc3@meta.data %>%
  dplyr::select(Sample.name,sccelltype, tissue) 
head(allmetadata)

table(allmetadata$tissue)
tissue = factor(c(rep('Periphery',1246),rep('Tumor',2334)),levels = c('Periphery','Tumor'))
tissue

dds <- DESeqDataSetFromMatrix(cluster_counts, 
                              allmetadata, 
                              design = ~tissue)



# Transform counts for data visualization
##小样本  rld <- rlog(dds, blind=TRUE)

rld <- vst(dds, blind=TRUE)


# Plot PCA 
DESeq2::plotPCA(rld, intgroup = "tissue")

#层次聚类
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)

# Plot heatmap
pheatmap(rld_cor, annotation = cluster_metadata[, c("tissue"), drop=F])



#运行deseq
dds <- DESeq(dds)





#离散估计图
plotDispEsts(dds)




levels(allmetadata$tissue)[2] = 'Tumor'
levels(allmetadata$tissue)[1] = 'Periphery'

contrast <- c("tissue", levels(allmetadata$tissue)[2], levels(allmetadata$tissue)[1])

# resultsNames(dds)
res <- results(dds, 
               contrast = contrast,
               alpha = 0.05)

res <- lfcShrink(dds, 
                 contrast =  contrast,
                 type="normal",
                 res=res)

res_tbl <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>%
  as_tibble()

# Check results output
res_tbl

# Write all results to file
write.csv(res_tbl,
          file = "all_all_genes.csv",
          quote = FALSE, 
          row.names = FALSE)



##checkpoint
save(qcdata2,pbmc3,rld_mat,rld_cor,rld,allmetadata,res,
     dds,sig_res,
     file='infil_all_pseudobulk.Rdata')
load('infil_all_pseudobulk.Rdata')






padj_cutoff <- 0.05

# Subset the significant results
sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
  dplyr::arrange(padj)

# Check significant genes output
sig_res

# Write significant results to file
write.csv(sig_res,
          #paste0("results", clusters[1], "_", levels(cluster_metadata$sample)[2], "_vs_", levels(cluster_metadata$sample)[1], "_sig_genes.csv"),
          file = "0.05_all_genes.csv",
          quote = FALSE, 
          row.names = FALSE)


###前 20 个最重要基因的标准化表达散点图
normalized_counts <- counts(dds, 
                            normalized = TRUE)

## Order results by padj values
top20_sig_genes <- sig_res %>%
  dplyr::arrange(padj) %>%
  dplyr::pull(gene) %>%
  head(n=20)


top20_sig_norm <- data.frame(normalized_counts) %>%
  rownames_to_column(var = "gene") %>%
  dplyr::filter(gene %in% top20_sig_genes)

gathered_top20_sig <- top20_sig_norm %>%
  gather(colnames(top20_sig_norm)[2:length(colnames(top20_sig_norm))], key = "samplename", value = "normalized_counts")






##deseq inner_join 作用？？

ei_meta = allmetadata[, c("sccelltype", "tissue" )]
ei_meta[,3]=rownames(ei_meta)
colnames(ei_meta)[3] = 'cell_names'
colnames(ei_meta)[2] = 'tissue'
colnames(ei_meta)[1] = 'sccelltype'





gathered_top20_sig <- inner_join(ei_meta[, c("cell_names", "tissue" )], gathered_top20_sig, by = c("cell_names" = "samplename"))

## plot using ggplot2
ggplot(gathered_top20_sig) +
  geom_point(aes(x = gene, 
                 y = normalized_counts, 
                 color = tissue), 
             position=position_jitter(w=0.1,h=0)) +
  scale_y_log10() +
  xlab("Genes") +
  ylab("log10 Normalized Counts") +
  ggtitle("Top 20 Significant DE Genes") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(plot.title = element_text(hjust = 0.5))


####热图

sig_norm <- data.frame(normalized_counts) %>%
  rownames_to_column(var = "gene") %>%
  dplyr::filter(gene %in% sig_res$gene)

# Set a color palette
heat_colors <- brewer.pal(6, "YlOrRd")

# Run pheatmap using the metadata data frame for the annotation
pheatmap(sig_norm[ , 2:length(colnames(sig_norm))], 
         color = heat_colors, 
         cluster_rows = T, 
         show_rownames = F,
         annotation = allmetadata[, c("sccelltype", "tissue")], 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20)        


#火山图
res_table_thres <- res_tbl %>% 
  dplyr::filter(abs(padj) < 0.05 & abs(log2FoldChange) >= 1)    #>1

threshold = abs(padj) < 0.05 & abs(log2FoldChange) >= 1

## Volcano plot
ggplot(res_table_thres) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = (abs(padj) < 0.05 & abs(log2FoldChange) >= 1))) +
  ggtitle("Volcano plot of stimulated B cells relative to control") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))             




#富集
library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)
library(GSEABase)
library(enrichplot)


##sig_res  应该是这个
sig_gene=filter(sig_res, padj<0.05)
rownames(sig_gene) = data.frame(sig_gene[,1])[,1]

#res_table_thres
sig_gene=filter(res_table_thres, padj<0.05)
rownames(sig_gene) = data.frame(sig_gene[,1])[,1]




geneList = sig_gene$log2FoldChange
names(geneList) = toupper(rownames(sig_gene))
geneList=sort(geneList,decreasing = T)


gmtfile = read.gmt('h.all.v7.4.symbols .gmt')  #（癌症）特征基因集合  #0
gmtfile = read.gmt('c6.all.v7.4.symbols.gmt')  #癌症特征基因  #1
gmtfile = read.gmt('c5.go.v7.4.symbols.gmt')  #GeneOntology
gmtfile = read.gmt('c8.all.v2022.1.Hs.symbols.gmt')   #0

gmtfile = read.gmt('c2.all.v2022.1.Hs.symbols.gmt')   #精选基因集
gmtfile = read.gmt('c3.all.v2022.1.Hs.symbols.gmt')    #调控目标基因集
gmtfile = read.gmt('c4.all.v2022.1.Hs.symbols.gmt')    #计算基因集
gmtfile = read.gmt('c7.all.v2022.1.Hs.symbols.gmt')    #免疫特征基因集    #0

gmtfile = read.gmt('GOBP_MICROGLIA_DIFFERENTIATION.v2022.1.Hs.gmt')  #0
gmtfile = read.gmt('KEGG_GLIOMA.v2022.1.Hs.gmt')   #1
gmtfile = read.gmt('COATES_MACROPHAGE_M1_VS_M2_UP.v2022.1.Hs.gmt' ) 
gmtfile = read.gmt('GOBP_CELLULAR_RESPONSE_TO_ACIDIC_PH.v2022.1.Hs.gmt' ) 
gmtfile = read.gmt('CUI_GLUCOSE_DEPRIVATION.v2023.1.Hs.gmt' ) 



{
egmt = GSEA(geneList=geneList ,TERM2GENE = gmtfile,minGSSize = 1,pvalueCutoff = 0.05,verbose = F )

gsea_result = egmt@result
rownames(gsea_result)

#gseaplot2(egmt,geneSetID = 'REACTOME_CELLULAR_RESPONSE_TO_STARVATION',pvalue_table=T)
#gseaplot2(egmt,geneSetID = 'HALLMARK_PEROXISOME',pvalue_table=T) 


down_kegg = gsea_result[gsea_result$pvalue<0.05 & gsea_result$enrichmentScore < 0.3,];down_kegg$group=-1
up_kegg = gsea_result[gsea_result$pvalue<0.05 & gsea_result$enrichmentScore > 0.3,];up_kegg$group=1

}




write.table(up_kegg,file='kegg.txt',sep = '\t' )
write.table(up_kegg,file='GO_up_kegg.txt',sep = '\t' )
write.table(down_kegg,file='GO_down_kegg.txt',sep = '\t' )

check = read.table(file='GO_up_kegg.txt',sep ='\t' )


#富集结果图
up_fig = up_kegg[,c(1,4,7)]


###    up_fig$p.adjust<-(-log10(up_fig$p.adjust))  #log10 p.adjust
up_fig<-up_fig[order(as.numeric(up_fig[,2]),decreasing=T),]
up_fig$enrichmentScore<-as.numeric(up_fig$enrichmentScore)

p = ggplot(up_fig,aes(enrichmentScore,ID))
p=p + geom_point(aes(color=-1*log10(p.adjust)))
p = p + theme_bw()
p = p+ geom_point(aes(size=-1*log10(p.adjust)))
p =p+labs(color=expression(-log[10](p.adjust)),size="size",  
          x="enrichmentScore",y="Pathway name",title="Pathway enrichment")
p






#收费内容 
#循环输出图
write.table(as.data.frame(down_kegg), file="pseudo down_kegg.xls", sep="\t", row.names=F)
write.table(as.data.frame(up_kegg), file="pseudo up_kegg.xls", sep="\t", row.names=F)

setwd('E:/R/glioma-infiltrate/kegg pseudo')

pro = 'test'

lapply(1:nrow(down_kegg),function(i){
  gseaplot2(egmt,down_kegg$ID[i],
            title=down_kegg$Description[i],pvalue_table = T)  
  ggsave(paste0(pro,'_down_kegg_',
                gsub('/','-',down_kegg$Description[i])
                ,'.pdf'))
}
)

lapply(1:nrow(up_kegg),function(i){
  gseaplot2(egmt,up_kegg$ID[i],
            title=up_kegg$Description[i],pvalue_table = T)  
  ggsave(paste0(pro,'_up_kegg_',
                gsub('/','-',up_kegg$Description[i])
                ,'.pdf'))
}
)




#提取目的基因分析表达情况
genelist = c('OSM','OSMR','LIF','LIFR','IL6ST')

geneexpression = subset(res_tbl,res_tbl[,"gene"] == genelist)



num = c(12169,12170,8583,8584,7664)

geneexpression = res_tbl[num,]


geneexpression[,8] = geneexpression[,3]/geneexpression[,5]
colnames(geneexpression)[8] = 'SD'



y <- rnorm(10,mean = -0.5638976,sd = 0.1494479)
y = data.frame(y)






