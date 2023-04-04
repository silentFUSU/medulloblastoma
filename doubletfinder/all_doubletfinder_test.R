set.seed(1) 
myPaths <- .libPaths()
new <- c('/storage/zhangyanxiaoLab/suzhuojie/miniconda3/envs/R4.2/lib/R/library','/storage/zhangyanxiaoLab/qihongjian/R/x86_64-redhat-linux-gnu-library/4.2')
myPaths <- c(myPaths, new) 
.libPaths(myPaths)
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/medulloblastoma/")
library(Seurat)
library(vegan)
library(ggplot2)
library(stringr)
library(DoubletFinder)
sample_name <- c("s02_pri","s02_rec","s04_pri","s04_rec","s05_pri","s05_met",
                 "s06_pri","s06_met","s07_pri","s07_met","s08_pri","s08_met",
                 "s09_pri","s09_met","s10_pri","s10_met","s11_pri","s11_met",
                 "s12_pri","s12_rec","s13_pri","s13_met")
pc.num <- 1:30
s08 <- readRDS(paste0("data/merge_data/sample_merged_rds/",name,"_rna_tmp_base.rds"))
for(i in 1:22){
  name<-sample_name[i]
  se <- readRDS(paste0("data/merge_data/sample_merged_rds/",name,"_rna_tmp_base.rds"))
  se <- RunPCA(se, verbose = F)
  se <- RunUMAP(se, dims = pc.num)
  se <- FindNeighbors(se, dims = pc.num) %>% FindClusters(resolution = 0.5)
  sweep.res.list <- paramSweep_v3(se, PCs = pc.num, sct = F)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE) 
  bcmvn <- find.pK(sweep.stats)
  pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
  doubletrate = ncol(se) * 3.1 * 1e-6 
  homotypic.prop <- modelHomotypic(se$seurat_clusters)
  nExp <- round(doubletrate * ncol(se))
  nExp.adj <- round(nExp * (1-homotypic.prop))
  se <- doubletFinder_v3(se, PCs = pc.num, pN = 0.25, pK = pK_bcmvn, nExp = nExp.adj, sct = F)
  DF.name = colnames(se@meta.data)[grepl("DF.classification2", colnames(se@meta.data))]
}
#s08met 
DimPlot(se, label = TRUE, repel = TRUE, group.by = "DF.classifications_0.25_0.03_32")
se <- subset(se, cells=rownames(se@meta.data)[which(se@meta.data$DF.classifications_0.25_0.03_24=='Singlet')])
DimPlot(se, label = TRUE, repel = TRUE)
saveRDS(se,"data/merge_data/sample_merged_rds/s08_met_rna_doubletfinder.rds")

rm(list=ls())
load("data/merge_data/sample_merged_rds/all_the_samples_splited.RData")
s08.met<-readRDS("data/merge_data/sample_merged_rds/s08_met_rna_doubletfinder.rds")
pair<-merge(s02.pri,c(s02.rec,s04.pri,s04.rec,s05.pri,s05.met,s06.pri,s06.met,s07.pri,s07.met,s08.pri,
                      s08.met,s09.pri,s09.met,s10.pri,s10.met,s12.pri,s12.rec,s13.pri,s13.met),
            c("s02_pri","s02_rec","s04_pri","s04_rec","s05_pri","s05_met",
              "s06_pri","s06_met","s07_pri","s07_met","s08_pri","s08_met",
              "s09_pri","s09_met","s10_pri","s10_met",
              "s12_pri","s12_rec","s13_pri","s13_met"))
pair <- NormalizeData(pair)
pair <- FindVariableFeatures(pair, selection.method = "vst", nfeatures = 2000)
pair <- ScaleData(pair)
pair <- RunPCA(pair, verbose = F)
pair <- RunUMAP(pair, dims = 1:30)
pair <- FindNeighbors(pair, reduction = "pca", dims = 1:30)
pair <- FindClusters(pair, resolution = 2.5)
p1 <- DimPlot(pair, label = TRUE, repel = TRUE, raster=FALSE)
ggsave("result/merge_data/annotation/s08_met_doubletfind_resolution_2.5/pair_ori2.5.png",p1,limitsize = FALSE,width = 15,height = 10)
p2 <- DimPlot(pair, label = TRUE, repel = TRUE,raster=FALSE,group.by = "orig.ident")
ggsave("result/merge_data/annotation/s08_met_doubletfind_resolution_2.5/pair_orig.ident.png",p2,limitsize = FALSE,width = 15,height = 10)

pair <- FindClusters(pair, resolution = 5.0)
p3 <- DimPlot(pair, label = TRUE, repel = TRUE, raster=FALSE)
ggsave("result/merge_data/annotation/s08_met_doubletfind_resolution_5/pair_ori5.png",p3,limitsize = FALSE,width = 15,height = 10)

pair@meta.data$highlight<-as.character(pair@meta.data$RNA_snn_res.5)
#64 74 79 86 97 101  
pair@meta.data$highlight[which(pair@meta.data$highlight %in% c(0:63,65:73,75:78,80:85,87:96,98:100,102:108))]<-"NA"
DimPlot(pair, label = TRUE, repel = TRUE,group.by = "highlight")
pair@meta.data$highlight_ori<-as.character(pair@meta.data$orig.ident)
pair@meta.data$highlight_ori[which(pair@meta.data$highlight_ori %in% c("S02-primary","S02-recurrent","S04-primary", "S04-recurrent","S05-metastatic", "S05-primary", "S06-metastatic","S06-primary",
                                                                       "S07-metastatic","S07-primary", "S08-metastatic","S09-primary", "S10-metastatic","S10-primary","S12-recurrent", 
                                                                      "S13-metastatic", "S13-primary" ))]<-"NA"
highlight_ori<-DimPlot(pair, label = TRUE, repel = TRUE,raster=FALSE,group.by = "highlight_ori")
ggsave("result/merge_data/after_doubletfinder/s08_met_doubletfind_resolution_2.5/pair_highlight_ori.png",highlight_ori,limitsize = FALSE,width = 15,height = 10)
 
for(i in 1:20){
  name<-sample_name[i]
  se <- readRDS(paste0("data/merge_data/sample_merged_rds/",name,"_rna_tmp_base.rds"))
  se <- RunPCA(se, verbose = F)
  se <- RunUMAP(se, dims = pc.num)
  se <- FindNeighbors(se, dims = pc.num) %>% FindClusters(resolution = 0.5)
  ori <- DimPlot(se,label = TRUE, repel = TRUE, raster=FALSE)
  ggsave(paste0("result/merge_data/before_doubletfinder/ori/",name,"_ori.png"),ori,limitsize = FALSE,width = 15,height = 10)
  saveRDS(se,paste0("data/merge_data/before_doubletfinder/cluster_resolution0.5_rds/",name,"_rna_cluster.rds"))
}

