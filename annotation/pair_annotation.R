set.seed(1) 
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-redhat-linux-gnu-library/4.2","/storage/zhangyanxiaoLab/suzhuojie/miniconda3/envs/R4.2/lib/R/library"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/medulloblastoma/")
library(Seurat)
library(vegan)
library(ggplot2)
library(stringr)

pair <- readRDS("data/merge_data/before_doubletfinder/sample_merged_rds/all_paired_sample_merged_rna_tmp_base_clustering.rds")
pair <- FindClusters(pair, resolution = 2.5)
p1<-DimPlot(pair, label = TRUE, repel = TRUE,raster=FALSE)
ggsave("result/merge_data/annotation/pair_ori2.5.png",p1,width = 15,height = 10)
p2<-DimPlot(pair, label = TRUE, repel = TRUE,raster=FALSE, group.by = "orig.ident")
ggsave("result/merge_data/annotation/pair_orig.ident.png",p2,width = 15,height = 10)

FeaturePlot(pair, label = TRUE, repel = TRUE, feature = "CD96") #T cell 69
VlnPlot(pair, features = "CD96")

FeaturePlot(pair, label = TRUE, repel = TRUE, feature = "PECAM1") #Endothelial 54
VlnPlot(pair, features = "PECAM1")

FeaturePlot(pair, label = TRUE, repel = TRUE, feature = "DCN") #Fibroblast 47
VlnPlot(pair, features = "DCN")

FeaturePlot(pair, label = TRUE, repel = TRUE, feature = "CSF1R") #Microglia 44
VlnPlot(pair, features = "CSF1R")

FeaturePlot(pair, label = TRUE, repel = TRUE, feature = "AQP4") #Astrocyte 63
VlnPlot(pair, features = "AQP4")

FeaturePlot(pair, label = TRUE, repel = TRUE, feature = "SLC1A3")#Astrocyte
VlnPlot(pair, features = "SLC1A3")

FeaturePlot(pair, label = TRUE, repel = TRUE, feature = "MBP")#Oligodendrocyte 67
VlnPlot(pair, features = "MBP")

FeaturePlot(pair, label = TRUE, repel = TRUE, feature = "CLDN11")#Oligodendrocyte 67
VlnPlot(pair, features = "CLDN11")

FeaturePlot(pair, label = TRUE, repel = TRUE, feature = "PDGFRA")#OPC 67
VlnPlot(pair, features = "PDGFRA")



p3<-FeaturePlot(pair,label = TRUE, raster=FALSE,feature = c("CD96","PECAM1","DCN","CSF1R","AQP4","MBP","CLDN11","PDGFRA"),ncol=4)
ggsave("result/merge_data/annotation/pair_featureplot.png",p3,limitsize = FALSE,width = 60,height = 20)
p4<-VlnPlot(pair, raster=FALSE,feature = c("CD96","PECAM1","DCN","CSF1R","AQP4","MBP","CLDN11","PDGFRA"),ncol=4)
ggsave("result/merge_data/annotation/pair_Vlnplot.png",p4,limitsize = FALSE,width = 60,height = 20)

pair@meta.data$annotation<-as.character(pair@meta.data$RNA_snn_res.2.5)
pair@meta.data$annotation[which(pair@meta.data$annotation=='69')]<-"T cell"
pair@meta.data$annotation[which(pair@meta.data$annotation=='54')]<-"Endothelial"
pair@meta.data$annotation[which(pair@meta.data$annotation=='47')]<-"Fibroblast"
pair@meta.data$annotation[which(pair@meta.data$annotation=='44')]<-"Microglia"
pair@meta.data$annotation[which(pair@meta.data$annotation=='63')]<-"Astrocyte"
pair@meta.data$annotation[which(pair@meta.data$annotation=='67')]<-"Oligodendrocyte"
pair@meta.data$annotation[which(pair@meta.data$annotation %in% c(0:43,45:46,48:53,55:62,64:66,68,70:78))]<-"putative cancer"
p5 <- DimPlot(pair, label = TRUE, raster=FALSE,repel = TRUE,group.by = "annotation")
ggsave("result/merge_data/annotation/pair_annotation.png",p5,limitsize = FALSE,width = 15,height = 10)
saveRDS(pair,"data/merge_data/before_doubletfinder/annotation_rds/pair_annotation.rds")
#find doublet, this step should be done on each sample, if on merged data, could cause some issue
# pair_dou<-pair
# sweep.pair_dou <- paramSweep_v3(pair_dou, PCs = 1:30, sct = FALSE)
# sweep.stats_pair_dou <- summarizeSweep(sweep.pair_dou, GT = FALSE)
# bcmvn_pair_dou <- find.pK(sweep.stats_pair_dou)
# pK_bcmvn <- bcmvn_pair_dou$pK[which.max(bcmvn_pair_dou$BCmetric)] %>% as.character() %>% as.numeric()
# homotypic.prop <- modelHomotypic(pair_dou$seurat_clusters)
# nExp_poi <- round(0.075*nrow(pair_dou@meta.data))
# nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
# pair_dou <- doubletFinder_v3(pair_dou, PCs = 1:30, pN = 0.25, pK = pK_bcmvn, nExp = nExp.adj, sct = F)


# pair <- FindClusters(pair, resolution = 5)
# p6 <- DimPlot(pair, label = TRUE, repel = TRUE,raster=FALSE)
# ggsave("result/merge_data/annotation/resolution_5/pair_pair_ori5.png",p6,limitsize = FALSE,width = 15,height = 10)
# pair@meta.data$highlight<-as.character(pair@meta.data$RNA_snn_res.5)
# pair@meta.data$highlight[which(pair@meta.data$highlight %in% c(0:59,61:71,73,75:85,87:92,94,96:110))]<-"NA"
# DimPlot(pair, label = TRUE, repel = TRUE,group.by = "highlight")

# rm(list = ls())
# #after doublet finder
# pair <- readRDS("data/merge_data/after_doubletfinder/pair.rds")
# pair <- FindClusters(pair, resolution = 5)
# ori<-DimPlot(pair,label = TRUE, repel = TRUE,raster=FALSE)
# ggsave("result/merge_data/after_doubletfinder/ori/pair_pair_ori5.png",ori,limitsize = FALSE,width = 15,height = 10)
