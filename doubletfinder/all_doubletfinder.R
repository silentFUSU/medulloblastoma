rm(list=ls())
set.seed(1)
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/miniconda3/envs/R4.2/lib/R/library","/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-redhat-linux-gnu-library/4.2"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/medulloblastoma/")
source("code/library.R")
library(Seurat)
library(vegan)
library(ggplot2)
library(stringr)
library(DoubletFinder)
sample_name <- c("s02_pri","s02_rec","s04_pri","s04_rec","s05_pri","s05_met",
                 "s06_pri","s06_met","s07_pri","s07_met","s08_pri","s08_met",
                 "s09_pri","s09_met","s10_pri","s10_met",
                 "s12_pri","s12_rec","s13_pri","s13_met")
pc.num <- 1:30
for(i in 3:20){
  name<-sample_name[i]
  se <- readRDS(paste0("data/merge_data/before_doubletfinder/cluster_resolution0.5_rds/",name,"_rna_cluster.rds"))
  #se <- RunPCA(se, verbose = F)
  #se <- RunUMAP(se, dims = pc.num)
  #se <- FindNeighbors(se, dims = pc.num) %>% FindClusters(resolution = 0.5)
  sweep.res.list <- paramSweep_v3(se, PCs = pc.num, sct = F)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE) 
  bcmvn <- find.pK(sweep.stats)
  pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
  if(ncol(se)>=1000 & ncol(se)<2000){
    rate <- 0.8
  }else if(ncol(se)>=2000 & ncol(se)<3000){
    rate <- 1.6
  }else if(ncol(se)>=3000 & ncol(se)<4000){
    rate <- 2.3  
  }else if (ncol(se) >= 4000 & ncol(se) < 5000) {
    rate <- 3.9
  }else if (ncol(se) >= 5000 & ncol(se) < 6000) {
    rate <- 4.6
  }else if (ncol(se) >= 6000 & ncol(se) < 7000) {
    rate <- 5.4
  }else if (ncol(se) >= 7000 & ncol(se) < 8000) {
    rate <- 6.1
  }else if (ncol(se) >= 9000 & ncol(se) < 10000) {
    rate <- 6.9
  }else {
    rate <- 7.6
  }
  doubletrate = ncol(se) * rate * 1e-6 
  homotypic.prop <- modelHomotypic(se$seurat_clusters)
  nExp <- round(doubletrate * ncol(se))
  nExp.adj <- round(nExp * (1-homotypic.prop))
  se <- doubletFinder_v3(se, PCs = pc.num, pN = 0.25, pK = pK_bcmvn, nExp = nExp.adj, sct = F)
  DF.name = colnames(se@meta.data)[grepl("DF.classification", colnames(se@meta.data))]
  singlet<-DimPlot(se, label = TRUE, repel = TRUE, group.by = DF.name)
  ggsave(paste0("result/merge_data/after_doubletfinder/",name,"_singlet.png"),singlet,limitsize = FALSE,width = 15,height = 10)
  se <- subset(se, cells=rownames(se@meta.data)[which(se@meta.data[,10]=='Singlet')])
  se <- NormalizeData(se)
  se <- FindVariableFeatures(se, selection.method = "vst", nfeatures = 2000)
  se <- ScaleData(se)
  se <- RunPCA(se, verbose = F)
  se <- RunUMAP(se, dims = pc.num)
  se <- FindNeighbors(se, dims = pc.num) %>% FindClusters(resolution = 0.5)
  saveRDS(se,paste0("data/merge_data/after_doubletfinder/",name,"_rna_doublet.rds"))
}

for(i in 1:20){
  name<-sample_name[i]
  se <- readRDS(paste0("data/merge_data/after_doubletfinder/",name,"_rna_doublet.rds"))
  ori <- DimPlot(se,label = TRUE, repel = TRUE, raster=FALSE)
  ggsave(paste0("result/merge_data/after_doubletfinder/ori/",name,"_ori.png"),ori,limitsize = FALSE,width = 15,height = 10)
}
rm(list=ls())
s02.pri<-readRDS(paste0("data/merge_data/after_doubletfinder/s02_pri_rna_doublet.rds"))
s02.rec<-readRDS(paste0("data/merge_data/after_doubletfinder/s02_rec_rna_doublet.rds"))
s04.pri<-readRDS(paste0("data/merge_data/after_doubletfinder/s04_pri_rna_doublet.rds"))
s04.rec<-readRDS(paste0("data/merge_data/after_doubletfinder/s04_rec_rna_doublet.rds"))
s05.pri<-readRDS(paste0("data/merge_data/after_doubletfinder/s05_pri_rna_doublet.rds"))
s05.met<-readRDS(paste0("data/merge_data/after_doubletfinder/s05_met_rna_doublet.rds"))
s06.pri<-readRDS(paste0("data/merge_data/after_doubletfinder/s06_pri_rna_doublet.rds"))
s06.met<-readRDS(paste0("data/merge_data/after_doubletfinder/s06_met_rna_doublet.rds"))
s07.pri<-readRDS(paste0("data/merge_data/after_doubletfinder/s07_pri_rna_doublet.rds"))
s07.met<-readRDS(paste0("data/merge_data/after_doubletfinder/s07_met_rna_doublet.rds"))
s08.pri<-readRDS(paste0("data/merge_data/after_doubletfinder/s08_pri_rna_doublet.rds"))
s08.met<-readRDS(paste0("data/merge_data/after_doubletfinder/s08_met_rna_doublet.rds"))
s09.pri<-readRDS(paste0("data/merge_data/after_doubletfinder/s09_pri_rna_doublet.rds"))
s09.met<-readRDS(paste0("data/merge_data/after_doubletfinder/s09_met_rna_doublet.rds"))
s10.pri<-readRDS(paste0("data/merge_data/after_doubletfinder/s10_pri_rna_doublet.rds"))
s10.met<-readRDS(paste0("data/merge_data/after_doubletfinder/s10_met_rna_doublet.rds"))
s12.pri<-readRDS(paste0("data/merge_data/after_doubletfinder/s12_pri_rna_doublet.rds"))
s12.rec<-readRDS(paste0("data/merge_data/after_doubletfinder/s12_rec_rna_doublet.rds"))
s13.pri<-readRDS(paste0("data/merge_data/after_doubletfinder/s13_pri_rna_doublet.rds"))
s13.met<-readRDS(paste0("data/merge_data/after_doubletfinder/s13_met_rna_doublet.rds"))
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
saveRDS(pair,"data/merge_data/after_doubletfinder/pair.rds")
p1 <- DimPlot(pair, label = TRUE, repel = TRUE, raster=FALSE)
ggsave("result/merge_data/after_doubletfinder/ori/pair_ori2.5_2.png",p1,limitsize = FALSE,width = 15,height = 10)
p2 <- DimPlot(pair, label = TRUE, repel = TRUE,raster=FALSE,group.by = "orig.ident")
ggsave("result/merge_data/after_doubletfinder/ori/pair_orig_2.ident.png",p2,limitsize = FALSE,width = 15,height = 10)
