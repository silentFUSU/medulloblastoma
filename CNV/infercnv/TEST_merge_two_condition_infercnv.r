rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-redhat-linux-gnu-library/4.2","/storage/zhangyanxiaoLab/suzhuojie/miniconda3/envs/R4.2/lib/R/library"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/medulloblastoma/")
source("code/library.R")
library(infercnv)
library(ComplexHeatmap)
library(AnnoProbe)
library(dplyr)
se <- readRDS(paste0("data/merge_data/before_doubletfinder/annotation_rds/add_all_normal/s04_pri_rec_merge.rds"))
dat <- GetAssayData(se,assay = "RNA",slot = "counts")
dat <- as.data.frame(dat)
geneInfor=annoGene(rownames(dat),"SYMBOL",'human')
geneInfor=geneInfor[c(1,4:6)] 
geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
dat=dat[match(geneInfor[,1], rownames(dat)),] 
rownames(geneInfor) <- geneInfor$SYMBOL
geneInfor <- geneInfor[,-1]
meta <- subset(se@meta.data, select = c("subcluster"))
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=dat,
                                    annotations_file=meta,
                                    delim="\t",
                                    gene_order_file=geneInfor,
                                    ref_group_name=c("Astrocyte","Endothelial","Fibroblast","Microglia","Oligodendrocyte","T cell"))
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=paste0("result/merge_data/before_doubletfinder/infercnv/s04_pri_rec_merge/"), 
                             cluster_by_groups=F,  # 选择TRUE是按样本分组 改为FALSE会进行按另一个参数k_obs_groups给出的分组数（默认为1）进行分组
                             denoise=T,     #是否去噪
                             num_threads=16,
                             plot_chr_scale = T,
                             k_nn=30,  # 1.9.1 default param
                             leiden_resolution = 1,  # 1.9.1 default param
                             leiden_method = "simple",  # 1.9.1 default param
                             leiden_function = "modularity",
                             HMM=F) 