###
# 2021-10-14
###

setwd('E:/R/glioma-infiltrate')
load('infiltrate_cellchat_official.Rdata')

save(cellchat,groupSize,nPatterns,
     file='infiltrate_cellchat_official.Rdata')


#Install packages
#devtools::install_github("sqjin/CellChat")
#Library packages
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

#My data
load(file = '../../infiltrate_cellchat.Rdata') # This is a combined data from two biological conditions: normal and diseases
data.input = pbmc3@assays$RNA@data # normalized data matrix
data.input[1:3,1:3]
meta = pbmc3@meta.data # a dataframe with rownames containing cell mata data
meta[1:3,1:3]


# Prepare input data for CelChat analysis

# meta = data.frame(sccelltype = meta$sccelltype[cell.use], row.names = colnames(data.input)) # manually create a dataframe consisting of the cell sccelltype
unique(meta$sccelltype) # check the cell sccelltype
#>  [1] Inflam. FIB  FBN1+ FIB    APOE+ FIB    COL11A1+ FIB cDC2        
#>  [6] LC           Inflam. DC   cDC1         CD40LG+ TC   Inflam. TC  
#> [11] TC           NKT         
#> 12 Levels: APOE+ FIB FBN1+ FIB COL11A1+ FIB Inflam. FIB cDC1 cDC2 ... NKT

######################################################################################################
#p24
load('p24_cellchat_official.Rdata')
save(cellchat,groupSize,
     p24,pbmc3,
     file='p24_cellchat_official.Rdata')

p24<-pbmc3[,pbmc3@meta.data$`characteristics:.Patient.ID` %in% c("BT_S2","BT_S4")]
data.input = p24@assays$RNA@data # normalized data matrix
meta = p24@meta.data # a dataframe with rownames containing cell mata data
identity = data.frame(celltype =meta$sccelltype   , row.names = rownames(meta)) # 

unique(meta$sccelltype)

######################################################################################################


cellchat <- createCellChat(data.input)
cellchat <- addMeta(cellchat, meta = identity, meta.name = "celltype")


#> Create a CellChat object from a data matrix
#> Set cell identities for the new CellChat object
#> The cell groups used for CellChat analysis are  APOE+ FIB FBN1+ FIB COL11A1+ FIB Inflam. FIB cDC1 cDC2 LC Inflam. DC TC Inflam. TC CD40LG+ TC NKT

#cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "celltype") # set "sccelltype" as default cell identity

levels(cellchat@idents) # show factor levels of the cell sccelltype
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group



#save






CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)
#> Rows: 1,939
#> Columns: 11
#> $ interaction_name   <chr> "TGFB1_TGFBR1_TGFBR2", "TGFB2_TGFBR1_TGFBR2", "TGF鈥?
#> $ pathway_name       <chr> "TGFb", "TGFb", "TGFb", "TGFb", "TGFb", "TGFb", "T鈥?
#> $ ligand             <chr> "TGFB1", "TGFB2", "TGFB3", "TGFB1", "TGFB1", "TGFB鈥?
#> $ receptor           <chr> "TGFbR1_R2", "TGFbR1_R2", "TGFbR1_R2", "ACVR1B_TGF鈥?
#> $ agonist            <chr> "TGFb agonist", "TGFb agonist", "TGFb agonist", "T鈥?
#> $ antagonist         <chr> "TGFb antagonist", "TGFb antagonist", "TGFb antago鈥?
#> $ co_A_receptor      <chr> "", "", "", "", "", "", "", "", "", "", "", "", ""鈥?
#> $ co_I_receptor      <chr> "TGFb inhibition receptor", "TGFb inhibition recep鈥?
#> $ evidence           <chr> "KEGG: hsa04350", "KEGG: hsa04350", "KEGG: hsa0435鈥?
#> $ annotation         <chr> "Secreted Signaling", "Secreted Signaling", "Secre鈥?
#> $ interaction_name_2 <chr> "TGFB1 - (TGFBR1+TGFBR2)", "TGFB2 - (TGFBR1+TGFBR2鈥?

# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 4) # do parallel
#> Warning: [ONE-TIME WARNING] Forked processing ('multicore') is disabled
#> in future (>= 1.13.0) when running R from RStudio, because it is
#> considered unstable. Because of this, plan("multicore") will fall
#> back to plan("sequential"), and plan("multiprocess") will fall back to
#> plan("multisession") - not plan("multicore") as in the past. For more details,
#> how to control forked processing or not, and how to silence this warning in
#> future R sessions, see ?future::supportsMulticore
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)


#save
load('cellchat_TEST.Rdata')
save(cellchat,groupSize,
     pbmc3,
     file='cellchat_TEST.Rdata')


# project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.human)

cellchat <- computeCommunProb(cellchat,type='truncatedMean',trim=0.1)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))


linshi = data.frame(cellchat@netP$pathways)
write.table(file = 'test_infil_cellchat.txt',linshi,sep = '\t')


#save






par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}









#指定通路开始分析




#查看通路
cellchat@netP$pathways
### Select pathway
pathways.show <- c("LIFR")
###

# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver,layout = 'hierarchy')

# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")


# Chord diagram
#par(mfrow=c(1,1))
#netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.


# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
#> Do heatmap based on a single object

# Chord diagram
#group.cellType <- c(rep("Astrocyte", 1), rep("cancer_cell", 1),rep("Macrophage", 2),
#                  rep("Microglial_Cell", 1),rep("Oligodendrocyte", 2),
#                 rep("other_immune_cell", 1)) # grouping cell clusters into fibroblast, DC and TC cells
#names(group.cellType) <- levels(cellchat@idents)
#netVisual_chord_cell(cellchat, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))


#> Plot the aggregated cell-cell communication network at the signaling pathway level
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.



netAnalysis_contribution(cellchat, signaling = pathways.show)
pairLR.PTN <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.PTN[1,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver,layout = 'hierarchy')


# Circle plot
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")


# Chord diagram
#netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.






















##循环出图

# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
vertex.receiver = seq(1,4)

setwd('E:/R/glioma-infiltrate')

dir.create('test-infil-all-cellchat_fig')
setwd('E:/R/glioma-infiltrate/test-infil-all-cellchat_fig')

for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "circle")
  netVisual_heatmap(cellchat, signaling = pathways.show.all[i], color.heatmap = "Reds")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.png"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}


















# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), remove.isolate = FALSE)
#> Comparing communications on a single object


# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), signaling = c("CCL","CXCL"), remove.isolate = FALSE)
#> Comparing communications on a single object


# show all the significant interactions (L-R pairs) based on user's input (defined by `pairLR.use`)
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("CCL","CXCL","FGF"))
netVisual_bubble(cellchat, sources.use = c(3,4), targets.use = c(5:8), pairLR.use = pairLR.use, remove.isolate = TRUE)
#> Comparing communications on a single object


# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
# show all the interactions sending from Inflam.FIB
netVisual_chord_gene(cellchat, sources.use = 4, targets.use = c(5:11), lab.cex = 0.5,legend.pos.y = 30)



plotGeneExpression(cellchat, signaling = "EGF")
#> Registered S3 method overwritten by 'spatstat':
#>   method     from
#>   print.boxx cli
#> Scale for 'y' is already present. Adding another scale for 'y', which will
#> replace the existing scale.
#> Scale for 'y' is already present. Adding another scale for 'y', which will
#> replace the existing scale.
#> Scale for 'y' is already present. Adding another scale for 'y', which will
#> replace the existing scale.


plotGeneExpression(cellchat, signaling = "CXCL", enriched.only = FALSE)



# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)


# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("CXCL", "CCL"))
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2


















