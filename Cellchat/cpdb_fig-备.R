library(Seurat)
library(rmarkdown)
library(tidyverse)
setwd('E:/R/glioma-infiltrate')
load('cpdb_fig_sigmeans.Rdata')
save(pldf,peri_meansdf,core_meansdf,peri_pvalsdf,core_pvalsdf,all_meansdf,all_pvalsdf,
     file = 'cpdb_fig_sigmeans.Rdata')



peri = 'E:/R/glioma-infiltrate/peri cpdb/peri/'
core = 'E:/R/glioma-infiltrate/core cpdb/core/'


core_pvals <- read.delim(paste0(core,"pvalues.txt"), check.names = FALSE)
core_means <- read.delim(paste0(core,"cpdb_core.txt"), check.names = FALSE)

peri_pvals <- read.delim(paste0(peri,"pvalues.txt"), check.names = FALSE)
peri_means <- read.delim(paste0(peri,"cpdb_peri.txt"), check.names = FALSE)

#ºÏ²¢
core_pvals$pairname = paste0('core',"_",core_pvals$interacting_pair)
peri_pvals$pairname = paste0('peri',"_",peri_pvals$interacting_pair)

core_means$pairname = paste0('core',"_",core_means$interacting_pair)
peri_means$pairname = paste0('peri',"_",peri_means$interacting_pair)
##########################################################################################

invase_gene<- grep("OSM|OSMR|LIF|LIFR", core_means$interacting_pair,value = T)


##core
core_means %>% dplyr::filter(interacting_pair %in% invase_gene)%>%
  dplyr::select("interacting_pair",starts_with("cancer_cell"),ends_with("cancer_cell"))  %>%  
  reshape2::melt() -> core_meansdf
colnames(core_meansdf)<- c("interacting_pair","CC","means")

core_pvals %>% dplyr::filter(interacting_pair %in% invase_gene)%>%
  dplyr::select("interacting_pair",starts_with("cancer_cell"),ends_with("cancer_cell"))%>%  
  reshape2::melt()-> core_pvalsdf
colnames(core_pvalsdf)<- c("interacting_pair","CC","pvals")


core_meansdf$pairname = paste0('core',"_",core_meansdf$interacting_pair)
core_pvalsdf$pairname = paste0('core',"_",core_pvalsdf$interacting_pair)

##peri

peri_means %>% dplyr::filter(interacting_pair %in% invase_gene)%>%
  dplyr::select("interacting_pair",starts_with("cancer_cell"),ends_with("cancer_cell"))  %>%  
  reshape2::melt() -> peri_meansdf
colnames(peri_meansdf)<- c("interacting_pair","CC","means")

peri_pvals %>% dplyr::filter(interacting_pair %in% invase_gene)%>%
  dplyr::select("interacting_pair",starts_with("cancer_cell"),ends_with("cancer_cell"))%>%  
  reshape2::melt()-> peri_pvalsdf
colnames(peri_pvalsdf)<- c("interacting_pair","CC","pvals")


peri_meansdf$pairname = paste0('peri',"_",peri_meansdf$interacting_pair)
peri_pvalsdf$pairname = paste0('peri',"_",peri_pvalsdf$interacting_pair)

all_meansdf = rbind(peri_meansdf,core_meansdf)
all_pvalsdf = rbind(peri_pvalsdf,core_pvalsdf)



all_pvalsdf$joinlab<- paste0(all_pvalsdf$pairname,"_",all_pvalsdf$CC)
all_meansdf$joinlab<- paste0(all_meansdf$pairname,"_",all_meansdf$CC)
pldf <- merge(all_pvalsdf,all_meansdf,by = "joinlab")

##pldf$pair_fin = paste0(pldf$interacting_pair.x,"_",pldf$CC)



summary((filter(pldf,means >0))$means)

pldf%>% filter(means >0) %>% 
  ggplot(aes(CC.x,pairname.x) )+ 
  geom_point(aes(color=means,size=-log10(pvals+0.0000001)-(-log10(1+0.0000001)) )) +
  scale_size_continuous(range = c(0,10))+
  scale_color_gradient2(high="red",mid = "yellow",low ="darkblue",midpoint = 0.4  )+ theme_bw()+ 
  theme(axis.text.x = element_text(angle = -45,hjust = -0.1,vjust = 0.8)) 





pldf[is.na(pldf)==TRUE] <- 0



#######################################################################################################

{
-log10(pvals+0.0000001)-(-log10(1+0.0000001))


pldf%>% filter(means >0) %>% 
  ggplot(aes(CC.x,pairname.x) )+ 
  geom_point(aes(color=means,size=-log10(pvals+0.00001) )) +
  scale_size_continuous(range = c(0,10))+
  scale_color_gradient2(high="red",mid = "yellow",low ="darkblue",midpoint = 0.4  )+ theme_bw()+ 
  theme(axis.text.x = element_text(angle = -45,hjust = -0.1,vjust = 0.8)) 
  }

{
  chemokines <- grep("^CXC|CCL|CCR|CX3|XCL|XCR", mymeans$interacting_pair,value = T)
  chemokines <- grep("^CXC|CCL|CCR|CX3|XCL|XCR", mymeans$interacting_pair,value = T)
  th1 <- grep("IL2|IL12|IL18|IL27|IFNG|IL10|TNF$|TNF |LTA|LTB|STAT1|CCR5|CXCR3|IL12RB1|IFNGR1|TBX21|STAT4", 
              mymeans$interacting_pair,value = T)
  th2 <- grep("IL4|IL5|IL25|IL10|IL13|AREG|STAT6|GATA3|IL4R", 
              mymeans$interacting_pair,value = T)
  th17 <- grep("IL21|IL22|IL24|IL26|IL17A|IL17A|IL17F|IL17RA|IL10|RORC|RORA|STAT3|CCR4|CCR6|IL23RA|TGFB", 
               mymeans$interacting_pair,value = T)
  treg <- grep("IL35|IL10|FOXP3|IL2RA|TGFB", mymeans$interacting_pair,value = T)
  costimulatory <- grep("CD86|CD80|CD48|LILRB2|LILRB4|TNF|CD2|ICAM|SLAM|LT[AB]|NECTIN2|CD40|CD70|CD27|CD28|CD58|TSLP|PVR|CD44|CD55|CD[1-9]", 
                        mymeans$interacting_pair,value = T)
  coinhibitory <- grep("SIRP|CD47|ICOS|TIGIT|CTLA4|PDCD1|CD274|LAG3|HAVCR|VSIR", 
                       mymeans$interacting_pair,value = T)
  niche <- grep("CSF", mymeans$interacting_pair,value = T)
  
  mymeans %>% dplyr::filter(interacting_pair %in% invase_gene)%>%
    dplyr::select("interacting_pair",starts_with("cancer_cell"),ends_with("cancer_cell"))  %>%  
    reshape2::melt() -> meansdf
  
  colnames(meansdf)<- c("interacting_pair","CC","means")
  
  mypvals %>% dplyr::filter(interacting_pair %in% invase_gene)%>%
    dplyr::select("interacting_pair",starts_with("cancer_cell"),ends_with("cancer_cell"))%>%  
    reshape2::melt()-> pvalsdf
  
  colnames(pvalsdf)<- c("interacting_pair","CC","pvals")
  }









