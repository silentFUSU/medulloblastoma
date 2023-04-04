rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/miniconda3/envs/R4.2/lib/R/library","/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-redhat-linux-gnu-library/4.2"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/medulloblastoma/")
source("code/library.R")
library(infercnv)
library(ComplexHeatmap)
library(AnnoProbe)
library(dplyr)
samples <- c("s02_pri","s02_rec","s04_pri","s04_rec","s05_pri","s05_met","s06_pri","s06_met","s07_pri",
             "s07_met","s08_pri","s08_met","s09_pri","s09_met","s10_pri","s10_met","s12_pri","s12_rec",
             "s13_pri","s13_met")

for(i in c(1,3:20)){
  sample<-samples[i]
  se <- readRDS(paste0("data/merge_data/before_doubletfinder/annotation_rds/add_all_normal/",sample,".rds"))
  dat <- GetAssayData(se,assay = "RNA",slot = "counts")
  dat <- as.data.frame(dat)
  geneInfor=annoGene(rownames(dat),"SYMBOL",'human')
  geneInfor=geneInfor[c(1,4:6)] 
  geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
  dat=dat[match(geneInfor[,1], rownames(dat)),] 
  rownames(geneInfor) <- geneInfor$SYMBOL
  geneInfor <- geneInfor[,-1]
  # ano_table <- as.data.frame(table(se$annotation))
  # ano_table <- ano_table[-which(ano_table[,2]<=1),]
  # ano_table <- as.character(ano_table$Var1)
  # se<-subset(se,cells = rownames(se@meta.data)[which(se$annotation %in% ano_table)])
  cluster <- se$annotation
  sub_cell_type = case_when(
    cluster %in% c(0)~"A",
    cluster %in% c(1)~"B",
    cluster %in% c(2)~"C",
    cluster %in% c(3)~"D",
    cluster %in% c(4)~"E",
    cluster %in% c(5)~"F",
    cluster %in% c(6)~"G",
    cluster %in% c(7)~"H",
    cluster %in% c(8)~"I",
    cluster %in% c(9)~"J",
    cluster %in% c(10)~"K",
    cluster %in% c(11)~"L",
    cluster %in% c(12)~"M",
    cluster %in% c(13)~"N",
    cluster %in% c(14)~"O",
    cluster %in% c(15)~"P",
    cluster %in% c(16)~"Q",
    cluster %in% c(17)~"R",
    cluster %in% c(18)~"S",
    cluster %in% c("putative neurons")~"Z_putative neurons",
    TRUE ~ as.character(cluster))
  se@meta.data$sub_cell_type= sub_cell_type
  meta <- subset(se@meta.data, select = c("sub_cell_type"))
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix=dat,
                                      annotations_file=meta,
                                      delim="\t",
                                      gene_order_file=geneInfor,
                                      ref_group_name=c("Astrocyte","Endothelial","Fibroblast","Microglia","Oligodendrocyte","T cell"))
  infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                               out_dir=paste0("result/merge_data/before_doubletfinder/infercnv/",sample,"add_all_normal/"), 
                               cluster_by_groups=T,  # 选择TRUE是按样本分组 改为FALSE会进行按另一个参数k_obs_groups给出的分组数（默认为1）进行分组
                               denoise=T,     #是否去噪
                               num_threads=16,
                               plot_chr_scale = T,
                               k_nn=30,  # 1.9.1 default param
                               leiden_resolution = 1,  # 1.9.1 default param
                               leiden_method = "simple",  # 1.9.1 default param
                               leiden_function = "modularity",
                               HMM=F) 
  se$sub_cell_type <- factor(se$sub_cell_type,levels=c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S",
                                                       "Astrocyte","Endothelial","Fibroblast","Microglia","Oligodendrocyte","T cell", "Z_putative neurons"))
  ori<-DimPlot(se,label = TRUE, repel = TRUE,raster=FALSE,group.by = "sub_cell_type")
  ggsave(paste0("result/merge_data/before_doubletfinder/annotation/add_all_normal/",sample,"_ABC.png"),ori,width = 15,height = 10)
  saveRDS(infercnv_obj,paste0("data/merge_data/before_doubletfinder/infercnv_obj/",sample,"_add_all_normal_obj.rds"))
}

samples <- c("s02_pri","s02_rec","s04_pri","s04_rec","s05_pri","s05_met","s06_pri","s06_met","s07_pri",
             "s07_met","s08_pri","s08_met","s09_pri","s09_met","s10_pri","s10_met","s12_pri","s12_rec",
             "s13_pri","s13_met")

for(i in c(1:20)){
  sample<-samples[i]
  se <- readRDS(paste0("data/merge_data/before_doubletfinder/annotation_rds/add_all_normal/",sample,".rds"))
  cluster <- se$annotation
  sub_cell_type = case_when(
    cluster %in% c(0)~"A",
    cluster %in% c(1)~"B",
    cluster %in% c(2)~"C",
    cluster %in% c(3)~"D",
    cluster %in% c(4)~"E",
    cluster %in% c(5)~"F",
    cluster %in% c(6)~"G",
    cluster %in% c(7)~"H",
    cluster %in% c(8)~"I",
    cluster %in% c(9)~"J",
    cluster %in% c(10)~"K",
    cluster %in% c(11)~"L",
    cluster %in% c(12)~"M",
    cluster %in% c(13)~"N",
    cluster %in% c(14)~"O",
    cluster %in% c(15)~"P",
    cluster %in% c(16)~"Q",
    cluster %in% c(17)~"R",
    cluster %in% c(18)~"S",
    cluster %in% c("putative neurons")~"Z_putative neurons",
    TRUE ~ as.character(cluster))
  se@meta.data$sub_cell_type= sub_cell_type
  meta <- subset(se@meta.data, select = c("sub_cell_type"))
  se$sub_cell_type <- factor(se$sub_cell_type,levels=c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S",
                                                       "Astrocyte","Endothelial","Fibroblast","Microglia","Oligodendrocyte","T cell", "Z_putative neurons"))
  ori<-DimPlot(se,label = TRUE, repel = TRUE,raster=FALSE,group.by = "sub_cell_type")
  ggsave(paste0("result/merge_data/before_doubletfinder/annotation/add_all_normal/",sample,"_ABC.png"),ori,width = 15,height = 10)
}
