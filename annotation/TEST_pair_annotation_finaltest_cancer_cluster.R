set.seed(1) 
myPaths <- .libPaths()
new <- c('/storage/zhangyanxiaoLab/suzhuojie/miniconda3/envs/R4.2/lib/R/library')
myPaths <- c(myPaths, new) 
.libPaths(myPaths)
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/medulloblastoma/")
library(Seurat)
library(vegan)
library(ggplot2)
library(stringr)
library(data.table)
library(harmony)
library(infercnv)
pair <- readRDS("data/merge_data/before_doubletfinder/annotation_rds/pair_annotation.rds")
clu75 <- rownames(pair@meta.data)[which(pair$RNA_snn_res.2.5=="75")]
clu57 <- rownames(pair@meta.data)[which(pair$RNA_snn_res.2.5=="57")]
clu70 <- rownames(pair@meta.data)[which(pair$RNA_snn_res.2.5=="70")]
clu76 <- rownames(pair@meta.data)[which(pair$RNA_snn_res.2.5=="76")]

se <- readRDS("data/merge_data/before_doubletfinder/annotation_rds/s02_pri_annotation.rds")
putative_cancer<-rownames(se@meta.data)[which(se$annotation=="putative cancer")]
se$annotation[which(rownames(se@meta.data) %in% putative_cancer)] <- se$seurat_clusters[which(rownames(se@meta.data) %in% putative_cancer)]
se$annotation <- factor(se$annotation,levels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14",
                                               "Astrocyte","Endothelial","Fibroblast","Microglia","Oligodendrocyte","T_cell"))
DimPlot(se,label = TRUE, repel = TRUE,raster=FALSE,group.by = "annotation")
ggsave("result/merge_data/before_doubletfinder/annotation/annotation_final_test_putative_neurons/s02_pri.png",
       width = 15,height = 10)
saveRDS(se,"data/merge_data/before_doubletfinder/annotation_rds/final_test_putative_neurons/s02_pri_anntation.rds")

se <- readRDS("data/merge_data/before_doubletfinder/annotation_rds/s13_met_annotation.rds")
se@meta.data<-se@meta.data[,-c(10)]
colnames(se@meta.data)[9]<-"annotation"
clu70 <- grep("S04_pri",clu70,value = T)
clu70<-as.data.frame(clu70)
colnames(clu70)[1]<-"clu"
clu70$clu<-str_split_fixed(clu70$clu, "_", 3)
clu70<-clu70$clu[,3]
se$annotation[which(rownames(se@meta.data) %in% clu70)]<-"putative neurons"
se$annotation[which(se$RNA_snn_res.0.5=="5")]<-"putative neurons"
se$annotation[which(se$annotation =="clu57")]<-"putative cancer"
putative_cancer<-rownames(se@meta.data)[which(se$annotation=="putative cancer")]

se$annotation[which(rownames(se@meta.data) %in% putative_cancer)] <- se$seurat_clusters[which(rownames(se@meta.data) %in% putative_cancer)]
se$annotation[which(se$annotation %in% c("clu76","clu57","old_clu38","clu70","clu75"))] <- "putative neurons"
se$annotation <- factor(se$annotation,levels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17",
                                               "Astrocyte","Endothelial","Fibroblast","Microglia","Oligodendrocyte","T_cell", "putative neurons"))
DimPlot(se,label = TRUE, repel = TRUE,raster=FALSE,group.by = "annotation")
ggsave("result/merge_data/before_doubletfinder/annotation/annotation_final_test_putative_neurons/s13_met.png",
       width = 15,height = 10)
saveRDS(se,"data/merge_data/before_doubletfinder/annotation_rds/final_test_putative_neurons/s13_met_anntation.rds")
se <- readRDS("data/merge_data/before_doubletfinder/annotation_rds/final_test_putative_neurons/s04_rec_anntation.rds")
