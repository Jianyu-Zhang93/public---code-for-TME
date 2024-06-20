library(CellChat)
library(patchwork)
library(ggalluvial)
options(stringsAsFactors = FALSE)
setwd('E:/R/glioma-infiltrate')
load('infiltrate_cellchat_all.Rdata')

library(Seurat)
library(dplyr)
library(SeuratData)
library(patchwork) #最强大的拼图包
library(ggplot2)
library(CellChat)
library(ggalluvial)
library(svglite)

##合并
load('new-infiltrate_cellchat_peri222.Rdata')
cellchat_peri = cellchat
load('infiltrate_cellchat_core222.Rdata')
cellchat_core = cellchat

save(cellchat_core,cellchat_peri,cellchat,object.list,
     pbmc3,
     file='infiltrate_cellchat_all222.Rdata')
##

object.list <- list(core = cellchat_core, peri = cellchat_peri)

cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#> Merge the following slots: 'data.signaling','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.
cellchat

#交互总数和交互强度
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

##不同细胞群之间的相互作用数量或相互作用强度不同
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

#################################
gg1 <- netVisual_heatmap(cellchat_core,color.heatmap = c("white", "red"))
gg2 <- netVisual_heatmap(cellchat_peri,color.heatmap = c("white", "red"))
gg1 + gg2

#################################


gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4))
group.cellType <- factor(group.cellType, levels = c("FIB", "DC", "TC"))
object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#> Merge the following slots: 'data.signaling','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.



library(reticulate)
py_available()
#Win64OpenSSL_Light-1_1_1s
reticulate::py_install(packages = 'umap-learn')



cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
#> 2D visualization of signaling networks from datasets 1 2
netVisual_embeddingPairwiseZoomIn(cellchat, type = "functional", nCol = 2)





cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "structural")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "structural")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)
#> 2D visualization of signaling networks from datasets 1 2

netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)




rankSimilarity(cellchat, type = "functional")

pathway.show = c('CX3C','LIFR','OSM')
#"count": comparing the number of interactions;
gg1 <- rankNet(cellchat, measure = "count", mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, measure = "count", mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2



#"weight": comparing the total interaction weights (strength);
#gg1 <- rankNet(cellchat,  mode = "comparison", stacked = T, do.stat = TRUE)
#gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
#gg1 + gg2



#########
library(ComplexHeatmap)

i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 15)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 15)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))




#bubble图
#sources.use : targets.use = (1:7) : (1:7)(共7种细胞)  细胞A通讯至细胞B
#全部细胞
netVisual_bubble(cellchat, sources.use = c(1:7), targets.use = c(1:7),  comparison = c(1, 2), angle.x = 45)
#> Comparing communications on a merged object

##根据目的出图
pdf("./bubble.pdf",width = 16,height =30)
netVisual_bubble(cellchat, sources.use = c(1:7), targets.use = c(1:7),  comparison = c(1, 2), angle.x = 45, signaling = c("CCL","CXCL"))
dev.off()



pathways.show <- c('CCL') 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
#> 
#> Do heatmap based on a single object
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))



# Chord diagram
pathways.show <- c("OSM") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}




cellchat@meta$tissue = factor(cellchat@meta$tissue, levels = c("peri", "core")) # set factor level
plotGeneExpression(cellchat, signaling = "OSM", split.by = "datasets", colors.ggplot = T)



p1 = netVisual_heatmap(object.list[[2]], signaling = "OSM", color.heatmap = "Reds",title.name = paste("OSM", "signaling ",names(object.list)[i]))
p2 = netVisual_heatmap(object.list[[2]], signaling = "LIFR", color.heatmap = "Reds",title.name = paste("LIFR", "signaling ",names(object.list)[i]))
p1+p2





#####通路分类
genelist = read.table(file = 'genelist-immune.txt',sep = '\t')
genelist = read.table(file = 'genelist-onco.txt',sep = '\t')


genelist_immune = read.table(file = 'genelist-immune.txt',sep = '\t')
genelist_onco = read.table(file = 'genelist-onco.txt',sep = '\t')


##pathway.show = c('CX3C','LIFR','OSM')
pathway.show = genelist[,1]

pathway.show_immune = genelist_immune[,1]
pathway.show_onco = genelist_onco[,1]



#"count": comparing the number of interactions;
gg1 <- rankNet(cellchat, measure = "count", mode = "comparison", stacked = T, do.stat = TRUE,signaling = pathway.show)
gg2 <- rankNet(cellchat, measure = "count", mode = "comparison", stacked = F, do.stat = TRUE,signaling = pathway.show)
gg1 + gg2

#"weight": comparing the total interaction weights (strength);
gg1 <- rankNet(cellchat,  mode = "comparison", stacked = T, do.stat = TRUE,signaling = pathway.show)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE,signaling = pathway.show)
gg1 + gg2

##对比
gg1 <- rankNet(cellchat, measure = "count", mode = "comparison", stacked = T, do.stat = TRUE,signaling = pathway.show_onco, title='tumor-promoting pathways')
gg2 <- rankNet(cellchat, measure = "count", mode = "comparison", stacked = T, do.stat = TRUE,signaling = pathway.show_immune, title='immune-promoting pathways')
gg1 + gg2








