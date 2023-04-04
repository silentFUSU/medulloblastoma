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
for(i in 1:20){
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
}
