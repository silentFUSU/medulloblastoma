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
pair <- readRDS("data/merge_data/sample_merged_rds/all_paired_sample_merged_rna_tmp_base_clustering.rds")
pair <- FindClusters(pair, resolution = 2.5)
#DimPlot(pair,raster=FALSE,label = T,group.by = "RNA_snn_res.0.8")
pairano<-readRDS("data/merge_data/before_doubletfinder/annotation_rds/pair_annotation_finaltest.rds")
#DimPlot(pairano,raster=FALSE,label = T,group.by = "annotaion_final")
pair$annotation <- pairano$annotaion_final
DimPlot(pair,raster=FALSE,label = T,group.by = "annotation")
# check the difference of definition of normal cell 
paircancer <- subset(pair,cells = rownames(pair@meta.data)[which(pair$annotation=="putative cancer")])
paircancer <- rownames(paircancer@meta.data)
paircancer2 <- subset(pair,cells = rownames(pair@meta.data)[which(!(pair$RNA_snn_res.0.8 %in% c("28","30","36","33","34","40","42")))])
paircancer2 <- rownames(paircancer2@meta.data)
sect <- intersect(paircancer,paircancer2)
paircancer2dif <- paircancer2[which(!paircancer2%in%sect)]
pair$annotation[which(rownames(pair@meta.data) %in% paircancer2dif)] <-"putative cancer"
pair <-subset(pair,cells = rownames(pair@meta.data)[which(!rownames(pair@meta.data)%in%paircancer2dif)])
DimPlot(pair,raster=FALSE,label = T,cells = paircancer2dif,group.by = "annotation")
DimPlot(pair,raster=FALSE,label = T,group.by = "annotation")

#delete the normal cells which is closer to cancer cells
Tcell_cancer <- rownames(pair@meta.data)[which(pair$annotation=="T cell" & pair@reductions[["umap"]]@cell.embeddings[,2] > 10)]
Tcell_cancer2 <- rownames(pair@meta.data)[which(pair$annotation=="T cell" & pair@reductions[["umap"]]@cell.embeddings[,2] > 0)]
Endothelial_cancer <-  rownames(pair@meta.data)[which(pair$annotation=="Endothelial" & pair@reductions[["umap"]]@cell.embeddings[,1] > 0)]
neurons_cancer <- rownames(pair@meta.data)[which(pair$annotation=="putative neurons" & pair@reductions[["umap"]]@cell.embeddings[,1] > -5)]
pair <- subset(pair, cells = rownames(pair@meta.data)[which(!rownames(pair@meta.data) %in% Endothelial_cancer)])
pair <- subset(pair, cells = rownames(pair@meta.data)[which(!rownames(pair@meta.data) %in% Tcell_cancer)])
pair <- subset(pair, cells = rownames(pair@meta.data)[which(!rownames(pair@meta.data) %in% Tcell_cancer2)])
pair <- subset(pair, cells = rownames(pair@meta.data)[which(!rownames(pair@meta.data) %in% neurons_cancer)])
DimPlot(pair,raster=FALSE,label = T,group.by = "annotation")
ggsave("result/merge_data/before_doubletfinder/annotation/clean_annotation/pair.png",width = 15,height = 10)
saveRDS(pair,"data/merge_data/before_doubletfinder/annotation_rds/pair_annotation_final.rds")

#each sample annotation
samplename<-c("s02","s04","s05","s06","s07","s08","s09","s10","s12","s13")
samplename2<-c("S02","S04","S05","S06","S07","S08","S09","S10","S12","S13")
Astrocyte <- rownames(pair@meta.data)[which(pair@meta.data$annotation=="Astrocyte")]
T_cell <- rownames(pair@meta.data)[which(pair@meta.data$annotation=="T cell")]
Endothelial <- rownames(pair@meta.data)[which(pair@meta.data$annotation=="Endothelial")]
Microglia <- rownames(pair@meta.data)[which(pair@meta.data$annotation=="Microglia")]
Fibroblast <- rownames(pair@meta.data)[which(pair@meta.data$annotation=="Fibroblast")]
Oligodendrocyte <- rownames(pair@meta.data)[which(pair@meta.data$annotation=="Oligodendrocyte")]
neurons <- rownames(pair@meta.data)[which(pair@meta.data$annotation=="putative neurons")]
annotation <- list(Astrocyte,T_cell,Endothelial,Microglia,Fibroblast,Oligodendrocyte,neurons)
celltype<-c("Astrocyte","T_cell","Endothelial","Microglia","Fibroblast","Oligodendrocyte","putative neurons")
pair_barcode <- rownames(pair@meta.data)
i=1
  sample<-samplename[i]
  sample2<-samplename2[i]
  se <- readRDS(paste0("data/merge_data/before_doubletfinder/cluster_resolution0.5_rds/",sample,"_pri_rna_cluster.rds"))
  se_barcode <- grep(paste0(sample2,"_pri"),pair_barcode,value = T)
  se_barcode<-as.data.frame(se_barcode)
  colnames(se_barcode)[1]<-"barcode"
  se_barcode$barcode<-str_split_fixed(se_barcode$barcode, "_", 3)
  se_barcode<-se_barcode$barcode[,3]
  se <- subset(se,cells = se_barcode)
  se$annotation <- se$RNA_snn_res.0.5
  for(j in 1:7){
    clu <- grep(paste0(sample2,"_pri"),annotation[[j]],value = T)
    clu<-as.data.frame(clu)
    colnames(clu)[1]<-"clu"
    clu$clu<-str_split_fixed(clu$clu, "_", 3)
    clu<-clu$clu[,3]
    se$annotation <- as.character(se$annotation)
    se@meta.data$annotation[which(rownames(se@meta.data)%in%clu)]<-celltype[j]
    se$annotation <- as.factor(se$annotation)
  }
 se$annotation <- factor(se$annotation,levels=c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17",
                                                 "Astrocyte","Endothelial","Fibroblast","Microglia","Oligodendrocyte","T_cell", "putative neurons"))
 DimPlot(se,label=T,repel = T,group.by = "annotation")
 i=i+1
#s02_pri
  table(se$annotation)
  se <- subset(se,cells = rownames(se@meta.data)[which(se$annotation !="13")])
  DimPlot(se,label=T,repel = T,group.by = "annotation")
  ggsave(paste0("result/merge_data/before_doubletfinder/annotation/clean_annotation/",sample,"_pri.png"),width = 15,height = 10) 
  saveRDS(se,paste0("data/merge_data/before_doubletfinder/annotation_rds/clean_annotation/",sample,"_pri.rds"))
  
#s04_pri
  table(se$annotation)
  se <- subset(se,cells = rownames(se@meta.data)[which(se$annotation !="11")])
  DimPlot(se,label=T,repel = T,group.by = "annotation")
  ggsave(paste0("result/merge_data/before_doubletfinder/annotation/clean_annotation/",sample,"_pri.png"),width = 15,height = 10) 
  saveRDS(se,paste0("data/merge_data/before_doubletfinder/annotation_rds/clean_annotation/",sample,"_pri.rds"))

#s05_pri
  table(se$annotation)
  DimPlot(se,label=T,repel = T,group.by = "annotation",cols = c("putative neurons"="red"))
  se <- subset(se,cells = rownames(se@meta.data)[which(se$annotation !="putative neurons")])
  DimPlot(se,label=T,repel = T,group.by = "annotation")
  ggsave(paste0("result/merge_data/before_doubletfinder/annotation/clean_annotation/",sample,"_pri.png"),width = 15,height = 10) 
  saveRDS(se,paste0("data/merge_data/before_doubletfinder/annotation_rds/clean_annotation/",sample,"_pri.rds"))
  
#s06_pri
  table(se$annotation)
  se <- subset(se,cells = rownames(se@meta.data)[which(se$annotation !="7")])
  se <- subset(se,cells = rownames(se@meta.data)[which(se$annotation !="T_cell")])
  DimPlot(se,label=T,repel = T,group.by = "annotation",cols = c("0"="red"))
  DimPlot(se,label=T,repel = T,group.by = "annotation")
  
  ggsave(paste0("result/merge_data/before_doubletfinder/annotation/clean_annotation/",sample,"_pri.png"),width = 15,height = 10) 
  saveRDS(se,paste0("data/merge_data/before_doubletfinder/annotation_rds/clean_annotation/",sample,"_pri.rds"))
  
#s07_pri
  table(se$annotation)
  se <- subset(se,cells = rownames(se@meta.data)[which(se$annotation !="Oligodendrocyte")])
  se <- subset(se,cells = rownames(se@meta.data)[which(se$annotation !="T_cell")])
  se <- subset(se,cells = rownames(se@meta.data)[which(se$annotation !="9")])
  se <- subset(se,cells = rownames(se@meta.data)[which(se$annotation !="putative neurons")])
  DimPlot(se,label=T,repel = T,group.by = "annotation",cols = c("0"="red"))
  DimPlot(se,label=T,repel = T,group.by = "annotation")
  
  ggsave(paste0("result/merge_data/before_doubletfinder/annotation/clean_annotation/",sample,"_pri.png"),width = 15,height = 10) 
  saveRDS(se,paste0("data/merge_data/before_doubletfinder/annotation_rds/clean_annotation/",sample,"_pri.rds"))
  
#s08_pri
  table(se$annotation)
  se <- subset(se,cells = rownames(se@meta.data)[which(se$annotation !="9")])
  se <- subset(se,cells = rownames(se@meta.data)[which(se$annotation !="8")])
  se <- subset(se,cells = rownames(se@meta.data)[which(se$annotation !="putative neurons")])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="5" & se@reductions[["umap"]]@cell.embeddings[,1] < -10))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="3" & se@reductions[["umap"]]@cell.embeddings[,2] < -5))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="4" & se@reductions[["umap"]]@cell.embeddings[,2] < -5))])
  DimPlot(se,label=T,repel = T,group.by = "annotation",cols = c("4"="red"))
  DimPlot(se,label=T,repel = T,group.by = "annotation")
  ggsave(paste0("result/merge_data/before_doubletfinder/annotation/clean_annotation/",sample,"_pri.png"),width = 15,height = 10) 
  saveRDS(se,paste0("data/merge_data/before_doubletfinder/annotation_rds/clean_annotation/",sample,"_pri.rds"))

#s09_pri
  table(se$annotation)
  se <- subset(se,cells = rownames(se@meta.data)[which(se$annotation !="10")])
  se <- subset(se,cells = rownames(se@meta.data)[which(se$annotation !="12")])  
  se <- subset(se,cells = rownames(se@meta.data)[which(se$annotation !="putative neurons")])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="8" & se@reductions[["umap"]]@cell.embeddings[,1] < -10))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="5" & se@reductions[["umap"]]@cell.embeddings[,1] < -10))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="5" & se@reductions[["umap"]]@cell.embeddings[,1] < -2))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="5" & se@reductions[["umap"]]@cell.embeddings[,1] < -1.5 & se@reductions[["umap"]]@cell.embeddings[,2] < -6))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="7" & se@reductions[["umap"]]@cell.embeddings[,1] < -3 & se@reductions[["umap"]]@cell.embeddings[,2] < -5))])
  DimPlot(se,label=T,repel = T,group.by = "annotation")
  DimPlot(se,label=T,repel = T,group.by = "annotation",cols = c("5"="red"))
  ggsave(paste0("result/merge_data/before_doubletfinder/annotation/clean_annotation/",sample,"_pri.png"),width = 15,height = 10) 
  saveRDS(se,paste0("data/merge_data/before_doubletfinder/annotation_rds/clean_annotation/",sample,"_pri.rds"))

#s10_pri
  table(se$annotation)
  se <- subset(se,cells = rownames(se@meta.data)[which(se$annotation !="10")])
  se <- subset(se,cells = rownames(se@meta.data)[which(se$annotation !="12")])  
  se <- subset(se,cells = rownames(se@meta.data)[which(se$annotation !="putative neurons")])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="5" & se@reductions[["umap"]]@cell.embeddings[,2] < -10))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="8" & se@reductions[["umap"]]@cell.embeddings[,1] < -8))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="9" & se@reductions[["umap"]]@cell.embeddings[,2] < -11))])
  DimPlot(se,label=T,repel = T,group.by = "annotation")
  DimPlot(se,label=T,repel = T,group.by = "annotation",cols = c("9"="red"))
  ggsave(paste0("result/merge_data/before_doubletfinder/annotation/clean_annotation/",sample,"_pri.png"),width = 15,height = 10) 
  saveRDS(se,paste0("data/merge_data/before_doubletfinder/annotation_rds/clean_annotation/",sample,"_pri.rds"))
  
#s12_pri
  table(se$annotation)
  se <- subset(se,cells = rownames(se@meta.data)[which(se$annotation !="9")])
  se <- subset(se,cells = rownames(se@meta.data)[which(se$annotation !="13")])  
  se <- subset(se,cells = rownames(se@meta.data)[which(se$annotation !="14")])
  se <- subset(se,cells = rownames(se@meta.data)[which(se$annotation !="Astrocyte")])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="1" & se@reductions[["umap"]]@cell.embeddings[,2] < -9))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="0" & se@reductions[["umap"]]@cell.embeddings[,1] > 7))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="6" & se@reductions[["umap"]]@cell.embeddings[,1] > 7))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="3" & se@reductions[["umap"]]@cell.embeddings[,1] > 6))])
  DimPlot(se,label=T,repel = T,group.by = "annotation")
  DimPlot(se,label=T,repel = T,group.by = "annotation",cols = c("3"="red"))
  ggsave(paste0("result/merge_data/before_doubletfinder/annotation/clean_annotation/",sample,"_pri.png"),width = 15,height = 10) 
  saveRDS(se,paste0("data/merge_data/before_doubletfinder/annotation_rds/clean_annotation/",sample,"_pri.rds"))

#s13_pri
  table(se$annotation)
  se <- subset(se,cells = rownames(se@meta.data)[which(se$annotation !="10")])
  se <- subset(se,cells = rownames(se@meta.data)[which(se$annotation !="11")])  
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="6" & se@reductions[["umap"]]@cell.embeddings[,2] < -10))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="6" & se@reductions[["umap"]]@cell.embeddings[,1] < -10))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="7" & se@reductions[["umap"]]@cell.embeddings[,1] < -10))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="9" & se@reductions[["umap"]]@cell.embeddings[,1] < -10))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="5" & se@reductions[["umap"]]@cell.embeddings[,2] < -10))])
  DimPlot(se,label=T,repel = T,group.by = "annotation")
  DimPlot(se,label=T,repel = T,group.by = "annotation",cols = c("5"="red"))
  ggsave(paste0("result/merge_data/before_doubletfinder/annotation/clean_annotation/",sample,"_pri.png"),width = 15,height = 10) 
  saveRDS(se,paste0("data/merge_data/before_doubletfinder/annotation_rds/clean_annotation/",sample,"_pri.rds"))

samplename<-c("s05","s06","s07","s08","s09","s10","s13")
samplename2<-c("S05","S06","S07","S08","S09","S10","S13")  
i=1
  sample<-samplename[i]
  sample2<-samplename2[i]
  se <- readRDS(paste0("data/merge_data/before_doubletfinder/cluster_resolution0.5_rds/",sample,"_met_rna_cluster.rds"))
  se_barcode <- grep(paste0(sample2,"_met"),pair_barcode,value = T)
  se_barcode<-as.data.frame(se_barcode)
  colnames(se_barcode)[1]<-"barcode"
  se_barcode$barcode<-str_split_fixed(se_barcode$barcode, "_", 3)
  se_barcode<-se_barcode$barcode[,3]
  se <- subset(se,cells = se_barcode)
  se$annotation <- se$RNA_snn_res.0.5
  for(j in 1:7){
    clu <- grep(paste0(sample2,"_met"),annotation[[j]],value = T)
    clu<-as.data.frame(clu)
    colnames(clu)[1]<-"clu"
    clu$clu<-str_split_fixed(clu$clu, "_", 3)
    clu<-clu$clu[,3]
    se$annotation <- as.character(se$annotation)
    se@meta.data$annotation[which(rownames(se@meta.data)%in%clu)]<-celltype[j]
    se$annotation <- as.factor(se$annotation)
  }
  se$annotation <- factor(se$annotation,levels=c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17",
                                                 "Astrocyte","Endothelial","Fibroblast","Microglia","Oligodendrocyte","T_cell", "putative neurons"))
  DimPlot(se,label=T,repel = T,group.by = "annotation")
i=i+1

#s05_met
  table(se$annotation)
  se <- subset(se,cells = rownames(se@meta.data)[which(se$annotation !="14")])
  se <- subset(se,cells = rownames(se@meta.data)[which(se$annotation !="Astrocyte")])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="2" & se@reductions[["umap"]]@cell.embeddings[,1] < -10))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="4" & se@reductions[["umap"]]@cell.embeddings[,1] < -10))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="9" & se@reductions[["umap"]]@cell.embeddings[,1] < -9))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="10" & se@reductions[["umap"]]@cell.embeddings[,1] < -9))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="11" & se@reductions[["umap"]]@cell.embeddings[,1] < -13))])
  DimPlot(se,label=T,repel = T,group.by = "annotation")
  DimPlot(se,label=T,repel = T,group.by = "annotation",cols = c("15"="red"))
  ggsave(paste0("result/merge_data/before_doubletfinder/annotation/clean_annotation/",sample,"_met.png"),width = 15,height = 10) 
  saveRDS(se,paste0("data/merge_data/before_doubletfinder/annotation_rds/clean_annotation/",sample,"_met.rds"))
  
#s06_met
  table(se$annotation)
  se <- subset(se,cells = rownames(se@meta.data)[which(se$annotation !="5")])
  se <- subset(se,cells = rownames(se@meta.data)[which(se$annotation !="8")])
  se <- subset(se,cells = rownames(se@meta.data)[which(se$annotation !="T_cell")])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="0" & se@reductions[["umap"]]@cell.embeddings[,1] < -10))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="1" & se@reductions[["umap"]]@cell.embeddings[,2] > 5))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="3" & se@reductions[["umap"]]@cell.embeddings[,1] < -10))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="3" & se@reductions[["umap"]]@cell.embeddings[,2] > 5))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="7" & se@reductions[["umap"]]@cell.embeddings[,1] < -9))])
  DimPlot(se,label=T,repel = T,group.by = "annotation")
  DimPlot(se,label=T,repel = T,group.by = "annotation",cols = c("Oligodendrocyte"="red"))
  ggsave(paste0("result/merge_data/before_doubletfinder/annotation/clean_annotation/",sample,"_met.png"),width = 15,height = 10) 
  saveRDS(se,paste0("data/merge_data/before_doubletfinder/annotation_rds/clean_annotation/",sample,"_met.rds"))

#s07_met
  table(se$annotation)
  se <- subset(se,cells = rownames(se@meta.data)[which(se$annotation !="8")])
  se <- subset(se,cells = rownames(se@meta.data)[which(se$annotation !="Oligodendrocyte")])
  se <- subset(se,cells = rownames(se@meta.data)[which(se$annotation !="T_cell")])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="0" & se@reductions[["umap"]]@cell.embeddings[,2] > 9))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="2" & se@reductions[["umap"]]@cell.embeddings[,2] > 9))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="3" & se@reductions[["umap"]]@cell.embeddings[,2] > 9))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="1" & se@reductions[["umap"]]@cell.embeddings[,2] < -10))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="4" & se@reductions[["umap"]]@cell.embeddings[,2] < -10))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="5" & se@reductions[["umap"]]@cell.embeddings[,2] < -10))])
  DimPlot(se,label=T,repel = T,group.by = "annotation")
  DimPlot(se,label=T,repel = T,group.by = "annotation",cols = c("1"="red"))
  ggsave(paste0("result/merge_data/before_doubletfinder/annotation/clean_annotation/",sample,"_met.png"),width = 15,height = 10) 
  saveRDS(se,paste0("data/merge_data/before_doubletfinder/annotation_rds/clean_annotation/",sample,"_met.rds"))

#s08_met
  table(se$annotation)
  se <- subset(se,cells = rownames(se@meta.data)[which(se$annotation !="6")])
  se <- subset(se,cells = rownames(se@meta.data)[which(se$annotation !="7")])
  se <- subset(se,cells = rownames(se@meta.data)[which(se$annotation !="9")])
  se <- subset(se,cells = rownames(se@meta.data)[which(se$annotation !="10")])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="0" & se@reductions[["umap"]]@cell.embeddings[,2] > 10))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="0" & se@reductions[["umap"]]@cell.embeddings[,1] > 5))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="2" & se@reductions[["umap"]]@cell.embeddings[,2] > 10))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="2" & se@reductions[["umap"]]@cell.embeddings[,1] > 3))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="3" & se@reductions[["umap"]]@cell.embeddings[,2] > 10))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="3" & se@reductions[["umap"]]@cell.embeddings[,1] > 10))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="4" & se@reductions[["umap"]]@cell.embeddings[,1] > 3.5))])
  DimPlot(se,label=T,repel = T,group.by = "annotation")
  DimPlot(se,label=T,repel = T,group.by = "annotation",cols = c("6"="red"))
  ggsave(paste0("result/merge_data/before_doubletfinder/annotation/clean_annotation/",sample,"_met.png"),width = 15,height = 10) 
  saveRDS(se,paste0("data/merge_data/before_doubletfinder/annotation_rds/clean_annotation/",sample,"_met.rds"))
  
#s09_met
  table(se$annotation)
  se <- subset(se,cells = rownames(se@meta.data)[which(se$annotation !="7")])
  se <- subset(se,cells = rownames(se@meta.data)[which(se$annotation !="10")])
  se <- subset(se,cells = rownames(se@meta.data)[which(se$annotation !="Microglia")])

  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="0" & se@reductions[["umap"]]@cell.embeddings[,2] > 10))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="0" & se@reductions[["umap"]]@cell.embeddings[,1] > 5))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="2" & se@reductions[["umap"]]@cell.embeddings[,2] > 10))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="2" & se@reductions[["umap"]]@cell.embeddings[,1] > 3))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="3" & se@reductions[["umap"]]@cell.embeddings[,2] > 10))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="3" & se@reductions[["umap"]]@cell.embeddings[,1] > 10))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="4" & se@reductions[["umap"]]@cell.embeddings[,1] > 3.5))])
  DimPlot(se,label=T,repel = T,group.by = "annotation")
  DimPlot(se,label=T,repel = T,group.by = "annotation",cols = c("10"="red"))
  ggsave(paste0("result/merge_data/before_doubletfinder/annotation/clean_annotation/",sample,"_met.png"),width = 15,height = 10) 
  saveRDS(se,paste0("data/merge_data/before_doubletfinder/annotation_rds/clean_annotation/",sample,"_met.rds"))

#s10_met
  table(se$annotation)
  se <- subset(se,cells = rownames(se@meta.data)[which(se$annotation !="7")])
  se <- subset(se,cells = rownames(se@meta.data)[which(se$annotation !="Astrocyte")])

  
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="0" & se@reductions[["umap"]]@cell.embeddings[,2] > 6))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="0" & se@reductions[["umap"]]@cell.embeddings[,1] < -10))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="2" & se@reductions[["umap"]]@cell.embeddings[,2] > 6))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="3" & se@reductions[["umap"]]@cell.embeddings[,1] < -10))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="5" & se@reductions[["umap"]]@cell.embeddings[,2] > 6))])
  DimPlot(se,label=T,repel = T,group.by = "annotation") 
  DimPlot(se,label=T,repel = T,group.by = "annotation",cols = c("6"="red"))
  ggsave(paste0("result/merge_data/before_doubletfinder/annotation/clean_annotation/",sample,"_met.png"),width = 15,height = 10) 
  saveRDS(se,paste0("data/merge_data/before_doubletfinder/annotation_rds/clean_annotation/",sample,"_met.rds"))

#s13_met
  table(se$annotation)
  se <- subset(se,cells = rownames(se@meta.data)[which(se$annotation !="8")])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="0" & se@reductions[["umap"]]@cell.embeddings[,1] < -10))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="1" & se@reductions[["umap"]]@cell.embeddings[,1] < -10))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="2" & se@reductions[["umap"]]@cell.embeddings[,2] > 7))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="6" & se@reductions[["umap"]]@cell.embeddings[,2] > 7))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="5" & se@reductions[["umap"]]@cell.embeddings[,2] > 6))])
  DimPlot(se,label=T,repel = T,group.by = "annotation") 
  DimPlot(se,label=T,repel = T,group.by = "annotation",cols = c("10"="red"))
  ggsave(paste0("result/merge_data/before_doubletfinder/annotation/clean_annotation/",sample,"_met.png"),width = 15,height = 10) 
  saveRDS(se,paste0("data/merge_data/before_doubletfinder/annotation_rds/clean_annotation/",sample,"_met.rds"))

samplename<-c("s02","s04","s12")
samplename2<-c("S02","S04","S12")  
i=1
  sample<-samplename[i]
  sample2<-samplename2[i]
  se <- readRDS(paste0("data/merge_data/before_doubletfinder/cluster_resolution0.5_rds/",sample,"_rec_rna_cluster.rds"))
  se_barcode <- grep(paste0(sample2,"_rec"),pair_barcode,value = T)
  se_barcode<-as.data.frame(se_barcode)
  colnames(se_barcode)[1]<-"barcode"
  se_barcode$barcode<-str_split_fixed(se_barcode$barcode, "_", 3)
  se_barcode<-se_barcode$barcode[,3]
  se <- subset(se,cells = se_barcode)
  se$annotation <- se$RNA_snn_res.0.5
  for(j in 1:7){
    clu <- grep(paste0(sample2,"_rec"),annotation[[j]],value = T)
    clu<-as.data.frame(clu)
    colnames(clu)[1]<-"clu"
    clu$clu<-str_split_fixed(clu$clu, "_", 3)
    clu<-clu$clu[,3]
    se$annotation <- as.character(se$annotation)
    se@meta.data$annotation[which(rownames(se@meta.data)%in%clu)]<-celltype[j]
    se$annotation <- as.factor(se$annotation)
  }
  se$annotation <- factor(se$annotation,levels=c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17",
                                                 "Astrocyte","Endothelial","Fibroblast","Microglia","Oligodendrocyte","T_cell", "putative neurons"))
  DimPlot(se,label=T,repel = T,group.by = "annotation")
i=i+1

#s02_rec
  table(se$annotation)
  se <- subset(se,cells = rownames(se@meta.data)[which(se$annotation !="10")])
  se <- subset(se,cells = rownames(se@meta.data)[which(se$annotation !="11")])
  se <- subset(se,cells = rownames(se@meta.data)[which(se$annotation !="Endothelial")])
  se <- subset(se,cells = rownames(se@meta.data)[which(se$annotation !="Fibroblast")])
  se <- subset(se,cells = rownames(se@meta.data)[which(se$annotation !="T_cell")])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="2" & se@reductions[["umap"]]@cell.embeddings[,1] < -10))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="3" & se@reductions[["umap"]]@cell.embeddings[,1] < -10))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="4" & se@reductions[["umap"]]@cell.embeddings[,1] < -10))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="5" & se@reductions[["umap"]]@cell.embeddings[,2] > 10))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="6" & se@reductions[["umap"]]@cell.embeddings[,1] < -10))])
  DimPlot(se,label=T,repel = T,group.by = "annotation") 
  DimPlot(se,label=T,repel = T,group.by = "annotation",cols = c("6"="red"))
  ggsave(paste0("result/merge_data/before_doubletfinder/annotation/clean_annotation/",sample,"_rec.png"),width = 15,height = 10) 
  saveRDS(se,paste0("data/merge_data/before_doubletfinder/annotation_rds/clean_annotation/",sample,"_rec.rds"))
  
#s04_rec
  table(se$annotation)
  se <- subset(se,cells = rownames(se@meta.data)[which(se$annotation !="Fibroblast")])
  se <- subset(se,cells = rownames(se@meta.data)[which(se$annotation !="11")])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="6" & se@reductions[["umap"]]@cell.embeddings[,1] > 10))])
  DimPlot(se,label=T,repel = T,group.by = "annotation") 
  DimPlot(se,label=T,repel = T,group.by = "annotation",cols = c("6"="red"))
  ggsave(paste0("result/merge_data/before_doubletfinder/annotation/clean_annotation/",sample,"_rec.png"),width = 15,height = 10) 
  saveRDS(se,paste0("data/merge_data/before_doubletfinder/annotation_rds/clean_annotation/",sample,"_rec.rds"))
  
#s12_rec
  table(se$annotation)
  se <- subset(se,cells = rownames(se@meta.data)[which(se$annotation !="10")])
  se <- subset(se,cells = rownames(se@meta.data)[which(se$annotation !="11")])
  se <- subset(se,cells = rownames(se@meta.data)[which(se$annotation !="putative neurons")])
  se <- subset(se,cells = rownames(se@meta.data)[which(se$annotation !="Astrocyte")])
  se <- subset(se,cells = rownames(se@meta.data)[which(se$annotation !="Endothelial")])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="0" & se@reductions[["umap"]]@cell.embeddings[,1] < 1))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="1" & se@reductions[["umap"]]@cell.embeddings[,1] < -10))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="2" & se@reductions[["umap"]]@cell.embeddings[,1] < 1))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="3" & se@reductions[["umap"]]@cell.embeddings[,1] < -15))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="4" & se@reductions[["umap"]]@cell.embeddings[,1] < 1 & se@reductions[["umap"]]@cell.embeddings[,2] > 3))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="5" & se@reductions[["umap"]]@cell.embeddings[,1] < -9))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="5" & se@reductions[["umap"]]@cell.embeddings[,1] > 0))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="6" & se@reductions[["umap"]]@cell.embeddings[,1] < -9))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="8" & se@reductions[["umap"]]@cell.embeddings[,1] < -19))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="13" & se@reductions[["umap"]]@cell.embeddings[,1] < 0))])
  se <- subset(se,cells = rownames(se@meta.data)[which(!(se$annotation=="T_cell" & se@reductions[["umap"]]@cell.embeddings[,1] < -5))])
  DimPlot(se,label=T,repel = T,group.by = "annotation") 
  DimPlot(se,label=T,repel = T,group.by = "annotation",cols = c("T_cell"="red"))
  ggsave(paste0("result/merge_data/before_doubletfinder/annotation/clean_annotation/",sample,"_rec.png"),width = 15,height = 10) 
  saveRDS(se,paste0("data/merge_data/before_doubletfinder/annotation_rds/clean_annotation/",sample,"_rec.rds"))
  