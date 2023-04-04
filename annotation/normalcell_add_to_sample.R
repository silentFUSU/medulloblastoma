set.seed(1) 
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-redhat-linux-gnu-library/4.2","/storage/zhangyanxiaoLab/suzhuojie/miniconda3/envs/R4.2/lib/R/library"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/medulloblastoma/")
library(Seurat)
library(vegan)
library(ggplot2)
library(stringr)
library(dplyr)
pair <- readRDS("data/merge_data/before_doubletfinder/annotation_rds/pair_annotation.rds")
DimPlot(pair,label = TRUE, repel = TRUE,raster=FALSE,group.by = "annotation")
normalcell <- subset(pair,cells = rownames(pair@meta.data)[which(pair$annotation %in% c("Astrocyte","Endothelial","Fibroblast","Microglia","Oligodendrocyte","T cell"))])
DimPlot(normalcell,label = TRUE, repel = TRUE,raster=FALSE,group.by = "annotation")
normalcell <- FindVariableFeatures(normalcell, selection.method = "vst", nfeatures = 2000)
normalcell <- ScaleData(normalcell)
normalcell <- RunPCA(normalcell, features = VariableFeatures(object = normalcell))
normalcell <- RunUMAP(normalcell, dims = 1:30)
normalcell <- FindNeighbors(normalcell, dims = 1:30)
normalcell <- FindClusters(normalcell, resolution = 0.5)
DimPlot(pair,group.by = "orig.ident",raster=F,label = T)
samples <- c("s02_pri","s02_rec","s04_pri","s04_rec","s05_pri","s05_met","s06_pri","s06_met","s07_pri",
             "s07_met","s08_pri","s08_met","s09_pri","s09_met","s10_pri","s10_met","s12_pri","s12_rec",
             "s13_pri","s13_met")
for(i in 2:20){
  sample <- samples[i]
  se <- readRDS(paste0("data/merge_data/before_doubletfinder/annotation_rds/annotation_final_with_neurons/",sample,"_annotation.rds"))
 
  se <- subset(se,cells = rownames(se@meta.data)[which(!se$annotation %in% c("Astrocyte","Endothelial","Fibroblast","Microglia","Oligodendrocyte","T_cell"))])
  se <- merge(se,normalcell,add.cell.ids = c("cancer","normal_cell"))
  se <- NormalizeData(se, normalization.method = "LogNormalize", scale.factor = 10000)
  se <- FindVariableFeatures(se, selection.method = "vst", nfeatures = 2000)
  se <- ScaleData(se)
  se <- RunPCA(se, features = VariableFeatures(object = se))
  se <- RunUMAP(se, dims = 1:30)
  se <- FindNeighbors(se, dims = 1:30)
  se <- FindClusters(se, resolution = 0.5)
  DimPlot(se,label = TRUE, repel = TRUE,raster=FALSE,group.by = "annotation")
  ggsave(paste0("result/merge_data/before_doubletfinder/annotation/add_all_normal/",sample,".png"),width = 15,height = 10)
  saveRDS(se,paste0("data/merge_data/before_doubletfinder/annotation_rds/add_all_normal/",sample,".rds"))
}

# DimPlot(se,label = TRUE, repel = TRUE,raster=FALSE,group.by = "sub_cell_type")
# 
# pri <- readRDS("data/merge_data/before_doubletfinder/annotation_rds/add_all_normal/s12_pri.rds")
# rec <- readRDS("data/merge_data/before_doubletfinder/annotation_rds/add_all_normal/s12_rec.rds")
# 
# pri <- readRDS("data/merge_data/before_doubletfinder/annotation_rds/add_all_normal/s09_pri.rds")
# met <- readRDS("data/merge_data/before_doubletfinder/annotation_rds/add_all_normal/s09_met.rds")
# cluster <- pri$annotation
# sub_cell_type = case_when(
#   cluster %in% c(0)~"A",
#   cluster %in% c(1)~"B",
#   cluster %in% c(2)~"C",
#   cluster %in% c(3)~"D",
#   cluster %in% c(4)~"E",
#   cluster %in% c(5)~"F",
#   cluster %in% c(6)~"G",
#   cluster %in% c(7)~"H",
#   cluster %in% c(8)~"I",
#   cluster %in% c(9)~"J",
#   cluster %in% c(10)~"K",
#   cluster %in% c(11)~"L",
#   cluster %in% c(12)~"M",
#   cluster %in% c(13)~"N",
#   cluster %in% c(14)~"O",
#   cluster %in% c(15)~"P",
#   cluster %in% c(16)~"Q",
#   cluster %in% c(17)~"R",
#   cluster %in% c(18)~"S",
#   cluster %in% c("putative neurons")~"Z_putative neurons",
#   TRUE ~ as.character(cluster))
# pri$subcluster <- pri$annotation
# pri@meta.data$sub_cell_type= sub_cell_type
# pri$subcluster[which(pri$sub_cell_type %in% c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","Z_putative neurons"))] <- "s04_pri_sub" 
# DimPlot(s02_pri,group.by ="subcluster",label = T)
# 
# 
# 
# 
# cluster <- rec$annotation
# sub_cell_type = case_when(
#   cluster %in% c(0)~"A",
#   cluster %in% c(1)~"B",
#   cluster %in% c(2)~"C",
#   cluster %in% c(3)~"D",
#   cluster %in% c(4)~"E",
#   cluster %in% c(5)~"F",
#   cluster %in% c(6)~"G",
#   cluster %in% c(7)~"H",
#   cluster %in% c(8)~"I",
#   cluster %in% c(9)~"J",
#   cluster %in% c(10)~"K",
#   cluster %in% c(11)~"L",
#   cluster %in% c(12)~"M",
#   cluster %in% c(13)~"N",
#   cluster %in% c(14)~"O",
#   cluster %in% c(15)~"P",
#   cluster %in% c(16)~"Q",
#   cluster %in% c(17)~"R",
#   cluster %in% c(18)~"S",
#   cluster %in% c("putative neurons")~"Z_putative neurons",
#   TRUE ~ as.character(cluster))
# pri$subcluster <- pri$annotation
# rec@meta.data$sub_cell_type= sub_cell_type
# rec<-subset(rec,cells = rownames(rec@meta.data)[which(rec$sub_cell_type %in% c("A","B","C","D","E","F","G","I","J","K","L","M","N","O","P","Q","R","S","Z_putative neutons"))])
# rec$subcluster <- "s04_rec_sub1"
# rec$subcluster[which(rec$sub_cell_type %in% c("H","I","G"))] <- "s04_rec_sub2"
# se <- merge(pri,rec,add.cell.ids = c("pri","rec"))
# se <- NormalizeData(se, normalization.method = "LogNormalize", scale.factor = 10000)
# se <- FindVariableFeatures(se, selection.method = "vst", nfeatures = 2000)
# se <- ScaleData(se)
# se <- RunPCA(se, features = VariableFeatures(object = se))
# se <- RunUMAP(se, dims = 1:30)
# se <- FindNeighbors(se, dims = 1:30)
# se <- FindClusters(se, resolution = 0.5)
# DimPlot(se,group.by = "subcluster",label = T)
# ggsave("result/merge_data/before_doubletfinder/annotation/add_all_normal/s04_pri_rec_merge.png",width = 15,height = 10)
# saveRDS(se,"data/merge_data/before_doubletfinder/annotation_rds/add_all_normal/s04_pri_rec_merge.rds")
# 
# s09_pri<-readRDS("data/merge_data/before_doubletfinder/annotation_rds/add_all_normal/s09_pri.rds")
# s09_met<-readRDS("data/merge_data/before_doubletfinder/annotation_rds/add_all_normal/s09_met.rds")
# cluster <- s09_pri$annotation
# sub_cell_type = case_when(
#   cluster %in% c(0)~"A",
#   cluster %in% c(1)~"B",
#   cluster %in% c(2)~"C",
#   cluster %in% c(3)~"D",
#   cluster %in% c(4)~"E",
#   cluster %in% c(5)~"F",
#   cluster %in% c(6)~"G",
#   cluster %in% c(7)~"H",
#   cluster %in% c(8)~"I",
#   cluster %in% c(9)~"J",
#   cluster %in% c(10)~"K",
#   cluster %in% c(11)~"L",
#   cluster %in% c(12)~"M",
#   cluster %in% c(13)~"N",
#   cluster %in% c(14)~"O",
#   cluster %in% c(15)~"P",
#   cluster %in% c(16)~"Q",
#   cluster %in% c(17)~"R",
#   cluster %in% c(18)~"S",
#   cluster %in% c("putative neurons")~"Z_putative neurons",
#   TRUE ~ as.character(cluster))
# s09_pri$subcluster <- s09_pri$annotation
# s09_pri@meta.data$sub_cell_type= sub_cell_type
# s09_pri$subcluster[which(s09_pri$sub_cell_type =="A")]<-"s09_pri_sub1"
# s09_pri$subcluster[which(s09_pri$sub_cell_type =="F")]<-"s09_pri_sub2"
# s09_pri$subcluster[which(s09_pri$sub_cell_type %in% c("B","C","D","E","G","H","I","J","L"))]<-"s09_pri_sub3"
# cluster <- s09_met$annotation
# sub_cell_type = case_when(
#   cluster %in% c(0)~"A",
#   cluster %in% c(1)~"B",
#   cluster %in% c(2)~"C",
#   cluster %in% c(3)~"D",
#   cluster %in% c(4)~"E",
#   cluster %in% c(5)~"F",
#   cluster %in% c(6)~"G",
#   cluster %in% c(7)~"H",
#   cluster %in% c(8)~"I",
#   cluster %in% c(9)~"J",
#   cluster %in% c(10)~"K",
#   cluster %in% c(11)~"L",
#   cluster %in% c(12)~"M",
#   cluster %in% c(13)~"N",
#   cluster %in% c(14)~"O",
#   cluster %in% c(15)~"P",
#   cluster %in% c(16)~"Q",
#   cluster %in% c(17)~"R",
#   cluster %in% c(18)~"S",
#   cluster %in% c("putative neurons")~"Z_putative neurons",
#   TRUE ~ as.character(cluster))
# s09_met@meta.data$sub_cell_type= sub_cell_type
# s09_met<-subset(s09_met,cells = rownames(s09_met@meta.data)[which(s09_met$sub_cell_type %in% c("A","B","C","D","F","G","H","I","K","L","R"))])
# s09_met$subcluster<-"s09_met_sub1"
# s09_met$subcluster[which(s09_met$sub_cell_type %in% c("C","D"))]<-"s09_met_sub2"
# s09_met$subcluster[which(s09_met$sub_cell_type %in% c("F","G","H","I","K","L","R"))]<-"s09_met_sub3"
# se <- merge(s09_pri,s09_met,add.cell.ids = c("s09_pri","s09_met"))
# se <- NormalizeData(se, normalization.method = "LogNormalize", scale.factor = 10000)
# se <- FindVariableFeatures(se, selection.method = "vst", nfeatures = 2000)
# se <- ScaleData(se)
# se <- RunPCA(se, features = VariableFeatures(object = se))
# se <- RunUMAP(se, dims = 1:30)
# se <- FindNeighbors(se, dims = 1:30)
# se <- FindClusters(se, resolution = 0.5)
# DimPlot(se,group.by = "subcluster",label = T)
# ggsave("result/merge_data/before_doubletfinder/annotation/add_all_normal/s09_pri_met_merge.png",width = 15,height = 10)
# saveRDS(se,"data/merge_data/before_doubletfinder/annotation_rds/add_all_normal/s09_pri_met_merge.rds")
