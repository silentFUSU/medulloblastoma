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
pair <- readRDS("data/merge_data/before_doubletfinder/annotation_rds/pair_annotation.rds")
pair_old <- readRDS("data/old_data/after_doubletfinder/rna_anno_rds/7pair_annotation.rds")
clu75 <- as.data.frame(table(pair$orig.ident[which(pair$RNA_snn_res.2.5=="75")]))
clu75 <- rownames(pair@meta.data)[which(pair$RNA_snn_res.2.5=="75")]
s06_met <- readRDS("data/merge_data/before_doubletfinder/annotation_rds/s06_met_annotation.rds")
s06_met_old <- readRDS("data/old_data/after_doubletfinder/rna_anno_rds/s06_met_annotation.rds")
DimPlot(s06_met_old,label = TRUE, repel = TRUE,raster=FALSE,group.by = "test_clus")

clu75 <- grep("S06_met",clu75,value = T)
clu75<-as.data.frame(clu75)
colnames(clu75)[1]<-"clu"
clu75$clu<-str_split_fixed(clu75$clu, "_", 3)
clu75<-clu75$clu[,3]
s06_met$new_clu75<-s06_met$annotation
s06_met$new_clu75[which(rownames(s06_met@meta.data) %in% clu75)]<-"clu75"
DimPlot(s06_met,label = TRUE, repel = TRUE,raster=FALSE,group.by = "new_clu75")
ggsave("result/merge_data/before_doubletfinder/annotation/s06_met_new_clu75.png",width = 15,height = 10)

s06_met_old$new_clu75 <- s06_met_old$annotation
s06_met_old$new_clu75[which(rownames(s06_met_old@meta.data) %in% clu75)]<-"clu75"
DimPlot(s06_met_old,label = TRUE, repel = TRUE,raster=FALSE,group.by = "new_clu75")
ggsave("result/old_data/after_doubletfinder/annotation/s06_met_new_clu75.png",width = 15,height = 10)

old_clu38 <- rownames(pair_old@meta.data)[which(pair_old$RNA_snn_res.1.5==38)]
old_clu38 <- grep("S06_met",old_clu38,value = T)
old_clu38<-as.data.frame(old_clu38)
colnames(old_clu38)[1]<-"clu"
old_clu38$clu<-str_split_fixed(old_clu38$clu, "_", 3)
old_clu38<-old_clu38$clu[,3]
s06_met$new_clu75_old_clu38 <- s06_met$new_clu75
s06_met$new_clu75_old_clu38[which(rownames(s06_met@meta.data) %in% old_clu38)]<-"old_clu38"
DimPlot(s06_met,label = TRUE, repel = TRUE,raster=FALSE,group.by = "new_clu75_old_clu38")
ggsave("result/merge_data/before_doubletfinder/annotation/s06_met_new_clu75_old_clu38.png",width = 15,height = 10)

saveRDS(s06_met,"data/merge_data/before_doubletfinder/annotation_rds/s06_met_annotation.rds")
saveRDS(s06_met_old,"data/old_data/after_doubletfinder/rna_anno_rds/s06_met_annotation.rds")
pair$clu75<-"NA"
pair$clu75[which(pair$RNA_snn_res.2.5=="75")]<-"clu75"
DimPlot(pair,label = TRUE, repel = TRUE,raster=FALSE,group.by = "clu75",cols = c("NA"="#eeeeee","clu75"="#fc5185"))
ggsave("result/merge_data/before_doubletfinder/annotation/pair_clu75.png",width = 15,height = 10)

pair$clu57<-"NA"
pair$clu57[which(pair$RNA_snn_res.2.5=="57")]<-"clu57"
DimPlot(pair,label = TRUE, repel = TRUE,raster=FALSE,group.by = "clu57",cols = c("NA"="#eeeeee","clu57"="#fc5185"))
ggsave("result/merge_data/before_doubletfinder/annotation/pair_clu57.png",width = 15,height = 10)

clu57 <- as.data.frame(table(pair$orig.ident[which(pair$RNA_snn_res.2.5=="57")]))
sample <- c("S02_rec","S06_met","S05_met","S04_rec","S08_pri","S07_pri","S04_pri","S05_pri","S08_met","S09_pri","S12_rec")
sample2 <- c("s02_rec","s06_met","s05_met","s04_rec","s08_pri","s07_pri","s04_pri","s05_pri","s08_met","s09_pri","s12_rec")

for(i in 9:11){
  name <- sample[i]
  name2 <- sample2[i]
  se<-readRDS(paste0("data/merge_data/before_doubletfinder/annotation_rds/",name2,"_annotation.rds"))
  se$new_clu57 <- se$annotation
  clu57 <- rownames(pair@meta.data)[which(pair$RNA_snn_res.2.5=="57")]
  clu57 <- grep(name,clu57,value = T)
  clu57<-as.data.frame(clu57)
  colnames(clu57)[1]<-"clu"
  clu57$clu<-str_split_fixed(clu57$clu, "_", 3)
  clu57<-clu57$clu[,3]
  se$new_clu57[which(rownames(se@meta.data) %in% clu57)]<-"clu57"
  ano<-DimPlot(se, label = TRUE, repel = TRUE,group.by = "new_clu57")
  ggsave(paste0("result/merge_data/before_doubletfinder/annotation/",name2,"_new_clu57.png"),ano,width = 15, height = 10) 
  saveRDS(se,paste0("data/merge_data/before_doubletfinder/annotation_rds/",name2,"_annotation.rds"))
}

se<-readRDS(paste0("data/merge_data/before_doubletfinder/annotation_rds/s02_rec_annotation.rds"))
se$new_clu57_old_clu38 <- se$new_clu57
old_clu38 <- rownames(pair_old@meta.data)[which(pair_old$RNA_snn_res.1.5==38)]
old_clu38 <- grep("S02_rec",old_clu38,value = T)
old_clu38<-as.data.frame(old_clu38)
colnames(old_clu38)[1]<-"clu"
old_clu38$clu<-str_split_fixed(old_clu38$clu, "_", 3)
old_clu38<-old_clu38$clu[,3]
se$new_clu57_old_clu38[which(rownames(se@meta.data) %in% old_clu38)]<-"old_clu38"
ano<-DimPlot(se, label = TRUE, repel = TRUE,group.by = "new_clu57_old_clu38")
ggsave(paste0("result/merge_data/before_doubletfinder/annotation/s02_rec_new_clu57_old_clu38.png"),ano,width = 15, height = 10) 
saveRDS(se,paste0("data/merge_data/before_doubletfinder/annotation_rds/s02_rec_annotation.rds"))

se<-readRDS(paste0("data/merge_data/before_doubletfinder/annotation_rds/s05_met_annotation.rds"))
se$new_clu57_old_clu38 <- se$new_clu57
old_clu38 <- rownames(pair_old@meta.data)[which(pair_old$RNA_snn_res.1.5==38)]
old_clu38 <- grep("S05_met",old_clu38,value = T)
old_clu38<-as.data.frame(old_clu38)
colnames(old_clu38)[1]<-"clu"
old_clu38$clu<-str_split_fixed(old_clu38$clu, "_", 3)
old_clu38<-old_clu38$clu[,3]
se$new_clu57_old_clu38[which(rownames(se@meta.data) %in% old_clu38)]<-"old_clu38"
ano<-DimPlot(se, label = TRUE, repel = TRUE,group.by = "new_clu57_old_clu38")
ggsave(paste0("result/merge_data/before_doubletfinder/annotation/s05_met_new_clu57_old_clu38.png"),ano,width = 15, height = 10) 
saveRDS(se,paste0("data/merge_data/before_doubletfinder/annotation_rds/s05_met_annotation.rds"))

pair$clu70<-"NA"
pair$clu70[which(pair$RNA_snn_res.2.5=="70")]<-"clu70"
DimPlot(pair,label = TRUE, repel = TRUE,raster=FALSE,group.by = "clu70",cols = c("NA"="#eeeeee","clu70"="#fc5185"))
ggsave("result/merge_data/before_doubletfinder/annotation/pair_clu70.png",width = 15,height = 10)
clu70 <- as.data.frame(table(pair$orig.ident[which(pair$RNA_snn_res.2.5=="70")]))
se<-readRDS(paste0("data/merge_data/before_doubletfinder/annotation_rds/s08_met_annotation.rds"))
se$new_clu57_new_clu70 <- se$new_clu57
clu70 <- rownames(pair@meta.data)[which(pair$RNA_snn_res.2.5=="70")]
clu70 <- grep("S08_met",clu70,value = T)
clu70<-as.data.frame(clu70)
colnames(clu70)[1]<-"clu"
clu70$clu<-str_split_fixed(clu70$clu, "_", 3)
clu70<-clu70$clu[,3]
se$new_clu57_new_clu70[which(rownames(se@meta.data) %in% clu70)]<-"clu70"
ano<-DimPlot(se, label = TRUE, repel = TRUE,group.by = "new_clu57_new_clu70")
ggsave(paste0("result/merge_data/before_doubletfinder/annotation/s08_met_new_clu57_new_clu70.png"),ano,width = 15, height = 10) 
saveRDS(se,paste0("data/merge_data/before_doubletfinder/annotation_rds/s08_met_annotation.rds"))
old_clu31 <- rownames(pair_old@meta.data)[which(pair_old$RNA_snn_res.1.5==31)]
old_clu31 <- grep("S08_met",old_clu31,value = T)
old_clu31<-as.data.frame(old_clu31)
colnames(old_clu31)[1]<-"clu"
old_clu31$clu<-str_split_fixed(old_clu31$clu, "_", 3)
old_clu31<-old_clu31$clu[,3]
se$new_clu57_new_clu70_old_clu31<-se$new_clu57_new_clu70
se$new_clu57_new_clu70_old_clu31[which(rownames(se@meta.data) %in% old_clu31)]<-"old_clu31"
ano<-DimPlot(se, label = TRUE, repel = TRUE,group.by = "new_clu57_new_clu70_old_clu31")
ggsave(paste0("result/merge_data/before_doubletfinder/annotation/s08_met_new_clu57_new_clu70_old_clu31.png"),ano,width = 15, height = 10) 
saveRDS(se,paste0("data/merge_data/before_doubletfinder/annotation_rds/s08_met_annotation.rds"))

pair$clu76<-"NA"
pair$clu76[which(pair$RNA_snn_res.2.5=="76")]<-"clu76"
DimPlot(pair,label = TRUE, repel = TRUE,raster=FALSE,group.by = "clu76",cols = c("NA"="#eeeeee","clu76"="#fc5185"))
ggsave("result/merge_data/before_doubletfinder/annotation/pair_clu76.png",width = 15,height = 10)
clu76 <- as.data.frame(table(pair$orig.ident[which(pair$RNA_snn_res.2.5=="76")]))
se <- readRDS(paste0("data/merge_data/before_doubletfinder/annotation_rds/s02_rec_annotation.rds"))
se$new_clu57_old_clu38_new_clu76 <- se$new_clu57_old_clu38
clu76 <- rownames(pair@meta.data)[which(pair$RNA_snn_res.2.5=="76")]
clu76 <- grep("S02_rec",clu76,value = T)
clu76<-as.data.frame(clu76)
colnames(clu76)[1]<-"clu"
clu76$clu<-str_split_fixed(clu76$clu, "_", 3)
clu76<-clu76$clu[,3]
se$new_clu57_old_clu38_new_clu76[which(rownames(se@meta.data) %in% clu76)]<-"clu76"
ano <- DimPlot(se,label = TRUE, repel = TRUE,raster=FALSE,group.by = "new_clu57_old_clu38_new_clu76")
ggsave(paste0("result/merge_data/before_doubletfinder/annotation/s02_rec_new_clu57_old_clu38_new_clu76.png"),ano,width = 15, height = 10) 


s02_rec_clu5 <- rownames(se@meta.data)[which(se$RNA_snn_res.0.5=="5")]
s02_rec_clu5 <- paste0("S02_rec_",s02_rec_clu5)
pair$s02_rec_clu5<-"NA"
pair$s02_rec_clu5[which(rownames(pair@meta.data) %in% s02_rec_clu5)] <-"s02_rec_clu5"
DimPlot(pair,label = TRUE, repel = TRUE, raster=FALSE,group.by = "s02_rec_clu5", cols = c("NA"="#eeeeee","s02_rec_clu5"="#fc5185"))
ggsave(paste0("result/merge_data/before_doubletfinder/annotation/pair_s02_rec_clu5.png"),width = 15, height = 10)
se$clu5<-"NA"
se$clu5[which(se$RNA_snn_res.0.5=="5")]<-"clu5"
DimPlot(se,label = TRUE, repel = TRUE, raster=FALSE,group.by = "clu5")
saveRDS(se,"data/merge_data/before_doubletfinder/annotation_rds/s02_rec_annotation.rds")

se$new_clu57_old_clu38_new_clu76_new_clu70<-se$new_clu57_old_clu38_new_clu76
clu70 <- rownames(pair@meta.data)[which(pair$RNA_snn_res.2.5=="70")]
clu70 <- grep("S02_rec",clu70,value = T)
clu70<-as.data.frame(clu70)
colnames(clu70)[1]<-"clu"
clu70$clu<-str_split_fixed(clu70$clu, "_", 3)
clu70<-clu70$clu[,3]
se$new_clu57_old_clu38_new_clu76_new_clu70[which(rownames(se@meta.data) %in% clu70)]<-"clu70"
DimPlot(se,label = TRUE, repel = TRUE,raster=FALSE,group.by = "new_clu57_old_clu38_new_clu76_new_clu70")
ggsave(paste0("result/merge_data/before_doubletfinder/annotation/s02_rec_new_clu57_old_clu38_new_clu76_new_clu70.png"),width = 15, height = 10)
saveRDS(se,"data/merge_data/before_doubletfinder/annotation_rds/s02_rec_annotation.rds")

se<-readRDS(paste0("data/merge_data/before_doubletfinder/annotation_rds/s04_pri_annotation.rds"))
clu70 <- rownames(pair@meta.data)[which(pair$RNA_snn_res.2.5=="70")]
clu70 <- grep("S04_pri",clu70,value = T)
clu70<-as.data.frame(clu70)
colnames(clu70)[1]<-"clu"
clu70$clu<-str_split_fixed(clu70$clu, "_", 3)
clu70<-clu70$clu[,3]
se$new_clu57_new_clu70 <- se$new_clu57
se$new_clu57_new_clu70[which(rownames(se@meta.data) %in% clu70)]<-"clu70"
DimPlot(se,label = TRUE, repel = TRUE,raster=FALSE,group.by = "new_clu57_new_clu70")
ggsave(paste0("result/merge_data/before_doubletfinder/annotation/s04_pri_new_clu57_new_clu70.png"),width = 15, height = 10)


DimPlot(pair,label = TRUE, repel = TRUE,raster=FALSE,group.by = "annotation2")
pair$annotaion_final<-pair$annotation2
pair$annotaion_final[which(pair$RNA_snn_res.2.5 %in% c("57","70","75","76"))]<-"putative neurons"
DimPlot(pair,label = TRUE, repel = TRUE,raster=FALSE,group.by = "annotaion_final")
ggsave(paste0("result/merge_data/before_doubletfinder/annotation/pair_annotation_final.png"),width = 15, height = 10)

se <- readRDS("data/merge_data/before_doubletfinder/annotation_rds/s02_pri_annotation.rds")
DimPlot(se,label = TRUE, repel = TRUE,raster=FALSE,group.by = "RNA_snn_res.0.5")
se$clu6_clu12<-"NA"
se$clu6_clu12[which(se$RNA_snn_res.0.5=="6")]<-"clu6"
se$clu6_clu12[which(se$RNA_snn_res.0.5=="12")]<-"clu12"
DimPlot(se,label = TRUE, repel = TRUE,raster=FALSE,group.by = "clu6_clu12")
ggsave(paste0("result/merge_data/before_doubletfinder/annotation/s02_pri_clu6_clu12.png"),width = 15, height = 10)
s02_pri_clu6<-rownames(se@meta.data)[which(se$RNA_snn_res.0.5=="6")]
s02_pri_clu12<-rownames(se@meta.data)[which(se$RNA_snn_res.0.5=="12")]
s02_pri_clu6<-paste0("S02_pri_",s02_pri_clu6)
s02_pri_clu12<-paste0("S02_pri_",s02_pri_clu12)
pair$s02_pri_clu6_clu12<-"NA"
pair$s02_pri_clu6_clu12[which(rownames(pair@meta.data)%in%s02_pri_clu6)]<-"s02_pri_clu6"
pair$s02_pri_clu6_clu12[which(rownames(pair@meta.data)%in%s02_pri_clu12)]<-"s02_pri_clu12"
DimPlot(pair,label = TRUE, repel = TRUE,raster=FALSE,group.by = "s02_pri_clu6_clu12",cols = c("NA"="#eeeeee","s02_pri_clu12"="#fc5185","s02_pri_clu6"="#3fc1c9"))
ggsave(paste0("result/merge_data/before_doubletfinder/annotation/pair_s02_pri_clu6_clu12.png"),width = 15, height = 10)

se <- readRDS("data/merge_data/before_doubletfinder/annotation_rds/s04_rec_annotation.rds")
DimPlot(se,label = TRUE, repel = TRUE,raster=FALSE,group.by = "RNA_snn_res.0.5")
se$clu12<-"NA"
se$clu12[which(se$RNA_snn_res.0.5=="12")]<-"clu12"
DimPlot(se,label = TRUE, repel = TRUE,raster=FALSE,group.by = "clu12")
ggsave(paste0("result/merge_data/before_doubletfinder/annotation/s04_rec_clu12.png"),width = 15, height = 10)
s04_rec_clu12<-rownames(se@meta.data)[which(se$RNA_snn_res.0.5=="12")]
s04_rec_clu12<-paste0("S04_rec_",s04_rec_clu12)

pair$s04_rec_clu12<-"NA"
pair$s04_rec_clu12[which(rownames(pair@meta.data)%in%s04_rec_clu12)]<-"s04_rec_clu12"
DimPlot(pair,label = TRUE, repel = TRUE,raster=FALSE,group.by = "s04_rec_clu12",cols = c("NA"="#eeeeee","s04_rec_clu12"="#fc5185"))
ggsave("result/merge_data/before_doubletfinder/annotation/pair_s04_rec_clu12.png",width = 15, height = 10)

se <- readRDS("data/merge_data/before_doubletfinder/annotation_rds/s07_met_annotation.rds")
s07_met_clu9 <- rownames(se@meta.data)[which(se$RNA_snn_res.0.5=="9")]
s07_met_clu9 <- paste0("S07_met_",s07_met_clu9)
pair$s07_met_clu9 <- "NA"
pair$s07_met_clu9[which(rownames(pair@meta.data)%in%s07_met_clu9)]<-"s07_met_clu9"
DimPlot(pair,label = TRUE, repel = TRUE,raster=FALSE,group.by = "s07_met_clu9",cols = c("NA"="#eeeeee","s07_met_clu9"="#fc5185"))
ggsave("result/merge_data/before_doubletfinder/annotation/pair_s07_met_clu9.png",width = 15, height = 10)

saveRDS(pair,"data/merge_data/before_doubletfinder/annotation_rds/pair_annotation_finaltest.rds")
