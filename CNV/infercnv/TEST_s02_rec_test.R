rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/miniconda3/envs/R4.2/lib/R/library","/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-redhat-linux-gnu-library/4.2"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/medulloblastoma/")
source("code/library.R")
library(infercnv)
se <- readRDS("data/merge_data/before_doubletfinder/annotation_rds/final_test_putative_neurons/s02_rec_anntation.rds")
expFile = system.file("extdata","oligodendroglioma_expression_downsampled.counts.matrix.gz", package = "infercnv")
geneFile=system.file("extdata","oligodendroglioma_annotations_downsampled.txt", package ="infercnv")
groupFiles = system.file("extdata","gencode_downsampled.EXAMPLE_ONLY_DONT_REUSE.txt",package ="infercnv")
library(AnnoProbe)
dat <- GetAssayData(se,assay = "RNA",slot = "counts")
dat <- as.data.frame(dat)
geneInfor=annoGene(rownames(dat),"SYMBOL",'human')
geneInfor=geneInfor[c(1,4:6)] 
geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
dat=dat[match(geneInfor[,1], rownames(dat)),] 
rownames(geneInfor) <- geneInfor$SYMBOL
geneInfor <- geneInfor[,-1]
meta <- subset(se@meta.data, select = c("annotation"))
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=dat,
                                    annotations_file=meta,
                                    delim="\t",
                                    gene_order_file=geneInfor,
                                    ref_group_name=c("Astrocyte","Endothelial","Fibroblast","Oligodendrocyte","Microglia","T_cell"))
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=paste0("result/merge_data/before_doubletfinder/infercnv/s02_rec3/"), 
                             cluster_by_groups=T,  # 选择TRUE是按样本分组 改为FALSE会进行按另一个参数k_obs_groups给出的分组数（默认为1）进行分组
                             denoise=T,     #是否去噪
                             num_threads=16,
                             HMM=F) 
infercnv::plot_cnv(infercnv_obj, #上两步得到的infercnv对象
                   plot_chr_scale = T, #画染色体全长，默认只画出（分析用到的）基因
                   output_filename = "s02_rec_better_plot3",output_format = "png", #保存为pdf文件
                   ) #改颜色
