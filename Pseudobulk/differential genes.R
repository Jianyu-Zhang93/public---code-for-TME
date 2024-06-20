library(Seurat)
library(rmarkdown) 
library(tidyverse)
setwd('E:/R/glioma-infiltrate')

load('infil_diferential_genes.Rdata')


save(immune_all,all_cells,microglia,macrophage,astrocyte,cancer_Dgenes,mi_ma,
     
     file='infil_diferential_genes.Rdata')


immune_all = sig_res
all_cells = sig_res
microglia = sig_res
macrophage =sig_res
astrocyte = sig_res
cancer_Dgenes = sig_res

mi_ma = sig_res

sig_res = cancer_Dgenes


library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)
library(GSEABase)
library(enrichplot)


cancer_Dgenes
sig_res = mi_ma



##sig_res  应该是这个
sig_gene=filter(sig_res, padj<0.05)
rownames(sig_gene) = data.frame(sig_gene[,1])[,1]

#res_table_thres
#sig_gene=filter(res_table_thres, padj<0.05)
#rownames(sig_gene) = data.frame(sig_gene[,1])[,1]




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
gmtfile = read.gmt('COATES_MACROPHAGE_M1_VS_M2_UP.v2022.1.Hs.gmt')




egmt = GSEA(geneList=geneList ,TERM2GENE = gmtfile,minGSSize = 1,pvalueCutoff = 0.05,verbose = FALSE )

gsea_result = egmt@result
rownames(gsea_result)

#gseaplot2(egmt,geneSetID = 'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',pvalue_table=T)
#gseaplot2(egmt,geneSetID = 'HALLMARK_PEROXISOME',pvalue_table=T) 


down_kegg = gsea_result[gsea_result$pvalue<0.05 & gsea_result$enrichmentScore < 0.3,];down_kegg$group=-1
up_kegg = gsea_result[gsea_result$pvalue<0.05 & gsea_result$enrichmentScore > 0.3,];up_kegg$group=1

write.table(up_kegg,file='kegg.txt',sep = '\t' )
write.table(up_kegg,file='GO_up_kegg microglia.txt',sep = '\t' )
write.table(down_kegg,file='GO_down_kegg microglia.txt',sep = '\t' )

write.table(up_kegg,file='microglia_C7_up.txt',sep = '\t' )



check = read.table(file='microglia_C7_up.txt',sep ='\t' )


#富集结果图
up_fig = up_kegg[,c(1,4,7)]


###    up_fig$p.adjust<-(-log10(up_fig$p.adjust))  #log10 p.adjust
up_fig<-up_fig[order(as.numeric(up_fig[,2]),decreasing=T),]
up_fig$enrichmentScore<-as.numeric(up_fig$enrichmentScore)



p = ggplot(up_fig,aes(enrichmentScore,ID))
p=p + geom_point(aes(colour=-p.adjust,low = 'yellow',high='red'),size = 5)
p = p + theme_bw()
##p = p+ geom_point(aes(size=-1*log10(p.adjust)))

p = p +scale_fill_gradient(low='red', high='yellow')
p =p+labs(size="size",  color='adjust.p',
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


