rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-redhat-linux-gnu-library/4.2","/storage/zhangyanxiaoLab/suzhuojie/miniconda3/envs/R/lib/R/library"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/medulloblastoma/")
#source("code/library.R")
# library(AnnoProbe)
library(vegan)
library(tibble)
#library(infercnv)
library(stringr)
library(ggplot2)
library(vcfR)
library(VennDiagram)
# library(Rsamtools)
library(Seurat)
library(ggplot2)
library(pheatmap)
library(ComplexHeatmap)
library(tidyr)
sample1<-"S07_pri"
sample2<-"S07_met"
sample <- "s07"
sample_se1<-"s07_pri"
sample_se2<-"s07_met"
pri<-read.vcfR(paste0('/storage/zhangyanxiaoLab/suzhuojie/software/varCA/out_othersample_gatk_varscan/classify/',sample1,'_snp/final.vcf.gz'), verbose = FALSE )
other<-read.vcfR(paste0('/storage/zhangyanxiaoLab/suzhuojie/software/varCA/out_othersample_gatk_varscan/classify/',sample2,'_snp/final.vcf.gz'), verbose = FALSE )

filter <- read.table(file=paste0('/storage/zhangyanxiaoLab/suzhuojie/software/varCA/out_othersample_gatk_varscan/merged_snp/',sample1,'/filter.tsv.gz'),header = T)
filter2<-read.table(file=paste0('/storage/zhangyanxiaoLab/suzhuojie/software/varCA/out_othersample_gatk_varscan/merged_snp/',sample2,'/filter.tsv.gz'),header = T)

#pri_bed <- as.data.frame(read.table(paste0('/storage/zhangyanxiaoLab/suzhuojie/software/varCA/out_othersample_gatk_varscan/peaks/',sample1,'/peaks.bed'),header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
#other_bed <- as.data.frame(read.table(paste0('/storage/zhangyanxiaoLab/suzhuojie/software/varCA/out_othersample_gatk_varscan/peaks/',sample2,'/peaks.bed'),header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))

pri_snp<-as.data.frame(pri@fix)
pri_snp$label<-paste0(pri_snp$CHROM,"-",pri_snp$POS)
#pri_bed$label<-paste0(pri_bed$V1,"-",pri_bed$V3)
filter$label<-paste0(filter$CHROM,"-",filter$POS)
#pri_snp<-pri_snp[which(pri_snp$label %in% pri_bed$label),]
pri_filter<-filter[which(filter$label %in% pri_snp$label),]
pri_filter<-pri_filter[-which(pri_filter$varscan.snp.RD==0),]
pri_filter<-pri_filter[-which(pri_filter$varscan.snp.AD==0),]
pri_filter<-pri_filter[-which(pri_filter$varscan.snp.RD<10),]
pri_snp<-pri_snp[which(pri_snp$label%in%pri_filter$label),]


other_snp<-as.data.frame(other@fix)
other_snp$label<-paste0(other_snp$CHROM,"-",other_snp$POS)
#other_bed$label<-paste0(other_bed$V1,"-",other_bed$V3)
#other_snp<-other_snp[which(other_snp$POS %in% other_bed$V3),]
filter2$label<-paste0(filter2$CHROM,"-",filter2$POS)
other_filter<-filter2[which(filter2$label %in% other_snp$label),]
other_filter<-other_filter[-which(other_filter$varscan.snp.RD==0),]
other_filter<-other_filter[-which(other_filter$varscan.snp.AD==0),]
other_filter<-other_filter[-which(other_filter$varscan.snp.RD<10),]
other_snp<-other_snp[which(other_snp$label%in%other_filter$label),]

other_snp$symbol<-paste0(other_snp$label,"-",other_snp$REF,"-",other_snp$ALT)
pri_snp$symbol<-paste0(pri_snp$label,"-",pri_snp$REF,"-",pri_snp$ALT)
overlap <- intersect(pri_snp$symbol,other_snp$symbol)
pri_snp2<-pri_snp[which(pri_snp$symbol %in% overlap),]
other_snp2<-other_snp[which(other_snp$symbol %in% overlap),]



# overlap<-read.table("result/merge_data/before_doubletfinder/snp/s04_pri_snp_position.txt",header = F,row.names = NULL)

pri_se<-readRDS(paste0("data/merge_data/before_doubletfinder/annotation_rds/clean_annotation/",sample_se1,".rds"))
other_se<-readRDS(paste0("data/merge_data/before_doubletfinder/annotation_rds/clean_annotation/",sample_se2,".rds"))

pri_se_name<-rownames(pri_se@meta.data)[which(!pri_se$annotation %in% c("Astrocyte","Endothelial","Fibroblast","Microglia","Oligodendrocyte","T cell"))]
pri_se_name<-paste0("pri_",pri_se_name)
other_se_name<-rownames(other_se@meta.data)[which(!other_se$annotation %in% c("Astrocyte","Endothelial","Fibroblast","Microglia","Oligodendrocyte","T cell"))]
other_se_name<-paste0("met_",other_se_name)
overlap<-as.data.frame(overlap)
overlap<-separate(overlap,overlap,c("V1","V2","V3","V4"),"-")
overlap<-overlap[-which(!overlap$V4%in%c("A","T","C","G") ),]
write.table(overlap,paste0("result/merge_data/before_doubletfinder/snp/",sample,"_snp_position.txt"),quote = F,row.names = F,col.names = F)
overlap<-paste0(overlap$V1,"-",overlap$V2,"-",overlap$V3,"-",overlap$V4)
pri_se<-subset(pri_se,cells=rownames(pri_se@meta.data)[which(!pri_se$annotation %in% c("Astrocyte","Endothelial","Fibroblast","Microglia","Oligodendrocyte","T cell"))])
pri_se <- NormalizeData(pri_se, normalization.method = "LogNormalize", scale.factor = 10000)
pri_se <- FindVariableFeatures(pri_se, selection.method = "vst", nfeatures = 2000)
pri_se <- ScaleData(pri_se)
pri_se <- RunPCA(pri_se, features = VariableFeatures(object = pri_se))
pri_se <- RunUMAP(pri_se, dims = 1:30)
pri_se <- FindNeighbors(pri_se, dims = 1:30)
pri_se <- FindClusters(pri_se, resolution = 0.5)

other_se<-subset(other_se,cells=rownames(other_se@meta.data)[which(!other_se$annotation %in% c("Astrocyte","Endothelial","Fibroblast","Microglia","Oligodendrocyte","T cell"))])
other_se <- NormalizeData(other_se, normalization.method = "LogNormalize", scale.factor = 10000)
other_se <- FindVariableFeatures(other_se, selection.method = "vst", nfeatures = 2000)
other_se <- ScaleData(other_se)
other_se <- RunPCA(other_se, features = VariableFeatures(object = other_se))
other_se <- RunUMAP(other_se, dims = 1:30)
other_se <- FindNeighbors(other_se, dims = 1:30)
other_se <- FindClusters(other_se, resolution = 0.5)
DimPlot(pri_se,label = T)
DimPlot(other_se,label = T)

pri_se$seurat_clusters_label<-paste0("pri_",pri_se$seurat_clusters)
other_se$seurat_clusters_label<-paste0("met_",other_se$seurat_clusters)
merge<-merge(pri_se,other_se, add.cell.ids=c("pri","met"),project="merge")
merge <- NormalizeData(merge, normalization.method = "LogNormalize", scale.factor = 10000)
merge <- FindVariableFeatures(merge, selection.method = "vst", nfeatures = 2000)
merge <- ScaleData(merge)
merge <- RunPCA(merge, features = VariableFeatures(object = merge))
merge <- RunUMAP(merge, dims = 1:30)
merge <- FindNeighbors(merge, dims = 1:30)
merge <- FindClusters(merge, resolution = 0.5)

save.image(paste0("data/merge_data/before_doubletfinder/snp_Rdata/snp_",sample,".rdata"))
