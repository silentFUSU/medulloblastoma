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
old_clu38 <- rownames(pair_old@meta.data)[which(pair_old$RNA_snn_res.1.5==38)]
old_clu31 <- rownames(pair_old@meta.data)[which(pair_old$RNA_snn_res.1.5==31)]
pair$old_clu <- "NA"

pair$old_clu[which(rownames(pair@meta.data) %in% old_clu31)] <- "clu31"
DimPlot(pair,label = TRUE, repel = TRUE,raster=FALSE,group.by = "old_clu",
        cols = c("NA"="#eeeeee","clu31"="#ff5722","clu38"="#fecea8"))
ggsave(paste0("result/merge_data/before_doubletfinder/annotation/pair_oldclu31_38.png"),width = 15, height = 10) 

pair$new_clu_neurons <- "NA"
pair$new_clu_neurons[which(pair$RNA_snn_res.2.5 %in% c("57","75","70","76"))] <- "putative neurons"
DimPlot(pair,label = TRUE, repel = TRUE,raster=FALSE,group.by = "new_clu_neurons",
        cols = c("NA"="#eeeeee","putative neurons"="#ff5722"))
ggsave(paste0("result/merge_data/before_doubletfinder/annotation/pair_new_clu_neurons.png"),width = 15, height = 10) 

putative_neurons_clu <- as.data.frame(table(pair$orig.ident[which(pair$new_clu_neurons=="putative neurons")]))
myLabel = as.vector(putative_neurons_clu[,1]) 
myLabel = paste(myLabel, "(", round(putative_neurons_clu[,2] / sum(putative_neurons_clu[,2]) * 100, 2), "%)        ", sep = "") 
p = ggplot(putative_neurons_clu, aes(x = "", y = putative_neurons_clu[,2], fill = putative_neurons_clu[,1])) + 
  geom_bar(stat = "identity", width = 1) +    
  coord_polar(theta = "y") + 
  labs(x = "", y = "", title = "") + 
  theme(axis.ticks = element_blank()) + 
  theme(legend.title = element_blank(), legend.position = "top") + 
  scale_fill_discrete(breaks = putative_neurons_clu[,1], labels = myLabel)
p
ggsave(paste0("result/merge_data/before_doubletfinder/annotation/putative_neurons_pie.png"),p,width = 15, height = 10)
pair$new_4clu <- pair$RNA_snn_res.2.5
pair$new_4clu[!(pair$RNA_snn_res.2.5 %in% c("57","75","70","76"))] <- "NA"
DimPlot(pair,label = TRUE, repel = TRUE,raster=FALSE,group.by = "new_4clu")
ggsave(paste0("result/merge_data/before_doubletfinder/annotation/pair_new_4clu.png"),width = 15, height = 10)

library(enrichplot)
library(clusterProfiler)
#library(GOplot)
library(DOSE)
library(ggnewscale)
#library(topGO)
library(circlize)
library(ComplexHeatmap)
library(org.Hs.eg.db)
pair$new_clu_neurons <- as.character(pair$RNA_snn_res.2.5)
pair$new_clu_neurons[which(pair$RNA_snn_res.2.5 %in% c("57","75","70","76"))] <- "putative neurons"
DimPlot(pair,label = TRUE, repel = TRUE,raster=FALSE,group.by = "new_clu_neurons")

neuron_markers <- FindMarkers(pair,ident.1 = 'putative neurons', group.by = "new_clu_neurons")
Go_database <- org.Hs.eg.db
gene <- bitr(rownames(neuron_markers),fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = Go_database)#6.33% of input gene IDs are fail to map

GO<-enrichGO( gene$ENTREZID,#GO富集分析
              OrgDb = Go_database,
              keyType = "ENTREZID",#设定读取的gene ID类型
              ont = "ALL",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
              pvalueCutoff = 0.05,#设定p值阈值
              qvalueCutoff = 0.05,#设定q值阈值
              readable = T)
barplot(GO, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")
ggsave(paste0("result/merge_data/before_doubletfinder/annotation/neurons_GO.png"),width = 15, height = 25)
