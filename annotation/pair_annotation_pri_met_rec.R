set.seed(1) 
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-redhat-linux-gnu-library/4.2","/storage/zhangyanxiaoLab/suzhuojie/miniconda3/envs/R4.2/lib/R/library"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/medulloblastoma/")
library(Seurat)
library(vegan)
library(ggplot2)
library(stringr)
pair <- readRDS("data/merge_data/before_doubletfinder/annotation_rds/pair_annotation.rds")
pair <- FindClusters(pair, resolution = 2.5)
DimPlot(pair,label = TRUE, repel = TRUE,raster=FALSE,group.by = "RNA_snn_res.2.5")
pair$annotation[which(pair$RNA_snn_res.2.5==77)]<-"Astrocyte"
pair$annotation[which(pair$RNA_snn_res.2.5==74)]<-"Fibroblast"
DimPlot(pair,label = TRUE, repel = TRUE,raster=FALSE,group.by = "annotation")
ggsave(paste0("result/merge_data/before_doubletfinder/annotation/resolution_2.5/pair_annotation.png"),width = 15, height = 10) 
FeaturePlot(pair, label = TRUE, repel = TRUE, raster=FALSE,feature = "PAX6")#Neural Progenitor cell
VlnPlot(pair, features = "PAX6")

#将正常细胞的注释回溯到每个样本中
Astrocyte <- rownames(pair@meta.data)[which(pair@meta.data$annotation=="Astrocyte")]
T_cell <- rownames(pair@meta.data)[which(pair@meta.data$annotation=="T cell")]
Endothelial <- rownames(pair@meta.data)[which(pair@meta.data$annotation=="Endothelial")]
Microglia <- rownames(pair@meta.data)[which(pair@meta.data$annotation=="Microglia")]
Fibroblast <- rownames(pair@meta.data)[which(pair@meta.data$annotation=="Fibroblast")]
Oligodendrocyte <- rownames(pair@meta.data)[which(pair@meta.data$annotation=="Oligodendrocyte")]
neurons<-rownames(pair@meta.data)[which(pair@meta.data$annotation=="putative neurons")]
samplename<-c("s02_pri","s04_pri","s05_pri","s06_pri","s07_pri","s08_pri","s09_pri","s10_pri","s12_pri","s13_pri",
              "s05_met","s06_met","s07_met","s08_met","s09_met","s10_met","s13_met",
              "s02_rec","s04_rec","s12_rec")
samplename2<-c("S02_pri","S04_pri","S05_pri","S06_pri","S07_pri","S08_pri","S09_pri","S10_pri","S12_pri","S13_pri",
               "S05_met","S06_met","S07_met","S08_met","S09_met","S10_met","S13_met",
               "S02_rec","S04_rec","S12_rec")
annotation <- list(Astrocyte,T_cell,Endothelial,Microglia,Fibroblast,Oligodendrocyte,neurons)
celltype<-c("Astrocyte","T_cell","Endothelial","Microglia","Fibroblast","Oligodendrocyte","putative neurons")
for(i in 1:20){
  se <- readRDS(paste0("data/merge_data/before_doubletfinder/cluster_resolution0.5_rds/",samplename[i],"_rna_cluster.rds"))
  se@meta.data$annotation<-as.character(se$RNA_snn_res.0.5)
  for(j in 1:7){
    clu <- grep(samplename2[i],annotation[[j]],value = T)
    clu<-as.data.frame(clu)
    colnames(clu)[1]<-"clu"
    clu$clu<-str_split_fixed(clu$clu, "_", 3)
    clu<-clu$clu[,3]
    se@meta.data$annotation[which(rownames(se@meta.data)%in%clu)]<-celltype[j]
  }
  saveRDS(se,paste0("data/merge_data/before_doubletfinder/annotation_rds/annotation_final_with_neurons/",samplename[i],"_annotation.rds"))
}