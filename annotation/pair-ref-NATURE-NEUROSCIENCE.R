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
ref <- readRDS("data/nature2021_cerebellar_development.rds")
DefaultAssay(ref) <- "RNA"
ref2 <- DietSeurat(ref, assays = "RNA")
ref2 <- UpdateSeuratObject(ref2)
pair <- readRDS("data/merge_data/before_doubletfinder/annotation_rds/pair_annotation_finaltest.rds")

pair_ref <- merge(pair,ref,c("exp","ref"))
DefaultAssay(ref) <- "RNA"
ref <- DietSeurat(ref, assays = "RNA")
ref = UpdateSeuratObject(object = ref)
DimPlot(ref,label = TRUE, repel = TRUE,raster=FALSE,group.by = "figure_clusters")

# mat <- fread("data/cbl-dev/exprMatrix.tsv.gz")
# meta <- read.table("data/cbl-dev/meta.tsv", header=T, sep="\t", as.is=T, row.names=1)
# genes = mat[,1][[1]]
# genes = gsub(".+[|]", "", genes)
# mat = data.frame(mat[,-1], row.names=genes)
# pbmc <- CreateSeuratObject(counts = mat, project = "ref", meta.data=meta)

#so <- NormalizeData(so)
#so <- ScaleData(so)
#so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = 2000)
#so <- RunPCA(so,npcs = 30,features = VariableFeatures(object = so),eval = FALSE)
# pbmc <- PercentageFeatureSet(pbmc, pattern = "^MT-", col.name = "percent.mt")
# pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt",method = 'glmGamPoi')

pbmc<- ref2
pbmc<- NormalizeData(pbmc, verbose = FALSE)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc, verbose = FALSE)

pbmc <-RunPCA(pbmc)
pbmc <- RunUMAP(pbmc, reduction = "pca", dims = 1:30)
saveRDS(pbmc,"data/nature2021_cerebellar_development2.rds")
DimPlot(pbmc,label = TRUE, repel = TRUE,raster=FALSE,group.by = "figure_clusters")
ggsave("result/merge_data/before_doubletfinder/annotation/ref_figure_clusters.png",width = 15,height = 10)
pair_normal<-subset(pair,cells = rownames(pair@meta.data)[which(!(pair$annotaion_final %in% "putative cancer"))])
DimPlot(pair_normal,label = TRUE, repel = TRUE,raster=FALSE,group.by = "annotaion_final")
pair_ref <- merge(pair_normal,pbmc,c("exp","ref"))
pair_ref<- NormalizeData(pair_ref, verbose = T)
pair_ref <- FindVariableFeatures(pair_ref, selection.method = "vst", nfeatures = 2000)
pair_ref <- ScaleData(pair_ref, verbose = T)
pair_ref <-RunPCA(pair_ref)
pair_ref <- RunUMAP(pair_ref, reduction = "pca", dims = 1:30)
annotation <- c(pair_ref$annotaion_final[1:5196],pair_ref$fig_cell_type[5197:74370])
pair_ref$annotation_final<-annotation
DimPlot(pair_ref,label = TRUE, repel = TRUE,raster=FALSE,group.by = "annotation_final")
ggsave("result/merge_data/before_doubletfinder/annotation/pair_ref_before_harmony.png",width = 15,height = 10)

pair_ref2 <- merge(pair_normal,pbmc,c("exp","ref"))
pair_ref2 <- NormalizeData(pair_ref2, verbose = T)
pair_ref2 <- FindVariableFeatures(pair_ref2, selection.method = "vst", nfeatures = 2000)
pair_ref2 <- ScaleData(pair_ref2, verbose = T)
pair_ref2$annotation_final <- annotation
pair_ref2 <- RunHarmony(pair_ref2, group.by.vars="annotation_final")

pair_ref2 <- RunUMAP(pair_ref2, reduction = "harmony",dims = 1:30)
DimPlot(pair_ref2,label = TRUE, repel = TRUE,raster=FALSE,group.by = "annotation_final")
ggsave("result/merge_data/before_doubletfinder/annotation/pair_ref_after_harmony_annotation_final.png",width = 15,height = 10)
pair_ref2 <- RunHarmony(pair_ref2, group.by.vars="orig.ident")
pair_ref2 <- RunUMAP(pair_ref2, reduction = "harmony",dims = 1:30)
DimPlot(pair_ref2, label = TRUE, repel = TRUE,raster=FALSE,group.by = "annotation_final",cols = c("H-GN"="#eeeeee","H-GCP"="#eeeeee","H-PIP"="#eeeeee","H-MLI"="#eeeeee","H-iCN"="#eeeeee","H-PC"="#eeeeee","H-Choroid"="#eeeeee",
                                                                                                  "H-iCN"="#eeeeee","H-BS Choroid/Ependymal"="#eeeeee","H-RL"="#eeeeee","H-Pericytes"="#eeeeee","H-Glia"="#eeeeee","H-BG"="#eeeeee",
                                                                                                  "H-Ast/Ependymal"="#eeeeee","H-Ast"="#eeeeee","H-Meninges"="#eeeeee","T cell"="#eeeeee","Fibroblast"="#eeeeee","Astrocyte"="#eeeeee",
                                                                                                  "H-Brainstem"="#eeeeee","H-OPC"="#eeeeee",
                                                                                                  "H-Microglia"="#14ffec","Microglia"="#e84545",
                                                                                                  "H-Endothelial"="#3ec1d3","Endothelial"="#ff5722",
                                                                                                  "H-eCN/UBC"="#7effdb","putative neurons"="#e84a5f",
                                                                                                  "H-Committed OPC"="#005691","Oligodendrocyte"="#d72323"))
ggsave("result/merge_data/before_doubletfinder/annotation/pair_ref_after_harmony_oriident.png",width = 15,height = 10)
