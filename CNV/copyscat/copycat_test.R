set.seed(1) 
myPaths <- .libPaths()
new <- c('/storage/zhangyanxiaoLab/suzhuojie/miniconda3/envs/R4.2/lib/R/library')
myPaths <- c(myPaths, new) 
.libPaths(myPaths)
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/medulloblastoma/")
library(Seurat)
library(CopyscAT)
library(Signac)
library(pheatmap)
library(vegan)
files=c('S02-primary','S02-recurrent')
files2=c('s02_pri','s02_rec')
initialiseEnvironment(genomeFile="data/hg38_chrom_sizes.tsv",
                      cytobandFile="data/hg38_1e+06_cytoband_densities_granges.tsv",
                      cpgFile="data/hg38_1e+06_cpg_densities.tsv",
                      binSize=1e6,
                      minFrags=1e4,
                      cellSuffix=c("-1","-2"),
                      lowerTrim=0.5,
                      upperTrim=0.8)
chrom_order <- c("chr1p","chr1q","chr2p","chr2q","chr3p","chr3q","chr4p","chr4q","chr5p","chr5q","chr6p","chr6q",
                 "chr7p","chr7q","chr8p","chr8q","chr9p","chr9q","chr10p","chr10q","chr11p","chr11q","chr12p","chr12q","chr13p",
                 "chr13q","chr14p","chr14q","chr15p","chr15q","chr16p","chr16q","chr17p","chr17q","chr18p","chr18q","chr19p",
                 "chr19q","chr20p","chr20q","chr21q","chr22q","chrXp","chrXq")
for(i in 1:1){
  setOutputFile(paste0("result/merge_data/before_doubletfinder/copyscat/",files[i],"_atac/"),paste0(files[i],"_atac_"))
  dir.create(paste0("result/merge_data/before_doubletfinder/copyscat/",files[i],"_atac/"))
  #part1 INITIAL DATA 
  scData<-readInputTable(paste0("data/merge_data/copyscat/",files[i],"_fragments_process.tsv"), sep ="\t")
  pbmc.atac<-readRDS(paste0("data/merge_data/before_doubletfinder/cluster_resolution0.5_rds/",files2[i],"_rna_cluster.rds"))
  pbmc.atac.barcode<-rownames(pbmc.atac@meta.data)
  scData<-subset(scData,rownames(scData)%in%pbmc.atac.barcode)#seurat_obj subset
  scData_k_norm <- normalizeMatrixN(scData,
                                    logNorm = FALSE,
                                    maxZero=2000,
                                    imputeZeros = FALSE,
                                    blacklistProp = 0.8,
                                    blacklistCutoff=0,
                                    dividingFactor=1,
                                    upperFilterQuantile = 1)#1 QC
  summaryFunction<-cutAverage
  scData_collapse<-collapseChrom3N(scData_k_norm,
                                   summaryFunction=summaryFunction,
                                   binExpand = 1,
                                   minimumChromValue = 100,
                                   logTrans = FALSE,
                                   tssEnrich = 1,
                                   logBase=2,
                                   minCPG=300,
                                   powVal=0.73) 
  #scData_collapse_check<-filterCells(scData_collapse,minimumSegments = 40,minDensity = 0.1)#2 QC
  
  graphCNVDistribution(scData_collapse,outputSuffix = "violinsn")
  median_iqr <- computeCenters(scData_collapse,summaryFunction=summaryFunction)
  
  #part2 ASSESMENT OF CHROMOSOME-LEVEL CNVs
  #option1 using all cells to genertate 'normal'
  candidate_cnvs<-identifyCNVClusters(scData_collapse,
                                      median_iqr,
                                      useDummyCells = TRUE,
                                      propDummy=0.25,
                                      minMix=0.01,
                                      deltaMean = 0.03,
                                      deltaBIC2 = 0.25,
                                      bicMinimum = 0.1, 
                                      subsetSize=600,
                                      fakeCellSD = 0.08, 
                                      uncertaintyCutoff = 0.55,
                                      summaryFunction=summaryFunction,
                                      maxClust = 4,
                                      mergeCutoff = 3,
                                      IQRCutoff= 0.2,
                                      medianQuantileCutoff = 0.4)
  candidate_cnvs_clean<-clusterCNV(initialResultList = candidate_cnvs,medianIQR = candidate_cnvs[[3]],minDiff=1.5)
  #warning:In xtfrm.data.frame(x) : cannot xtfrm data frames
  final_cnv_list<-annotateCNV4(candidate_cnvs_clean, 
                               saveOutput=TRUE,
                               outputSuffix = "clean_cnv",
                               sdCNV = 0.5,
                               filterResults=F,
                               filterRange=-5)
  #"chr1p"  "chr1q"  "chr2p"  "chr2q"  "chr3p"  "chr3q"  "chr4p"  "chr4q"  "chr5p"  "chr5q"  "chr6p"  "chr6q"  
  #"chr7p"  "chr7q" "chr8p"  "chr8q"  "chr9p"  "chr9q"  "chr10p" "chr10q" "chr11p" "chr11q" "chr12p" "chr12q" 
  #"chr13q" "chr14q" "chr15q" "chr16p""chr16q" "chr17p" "chr17q" "chr18p" "chr18q" "chr19p" "chr19q" "chr20p" 
  #"chr20q" "chr21q" "chr22q" "chrXp"  "chrXq"  counts: 41
  
  #按照结果做个heatmap，但是好像发现chrom的顺序有点不对，重新做一个表
  final_cnv_list_matrix<-final_cnv_list[[3]]
  sample_barcode <- final_cnv_list_matrix[,1]  
  rownames(final_cnv_list_matrix) <- sample_barcode  
  final_cnv_list_matrix <- final_cnv_list_matrix[,-1]
  #final_cnv_list_matrix2 <- as.data.frame(final_cnv_list_matrix[,which(colnames(final_cnv_list_matrix)==chrom_order[1])])
  #rownames(final_cnv_list_matrix2) <- rownames(final_cnv_list_matrix)
  #colnames(final_cnv_list_matrix2)[1]<-chrom_order[1]
  #tmpm <- list()
  #final_cnv_list_matrix2<-reorder(final_cnv_list_matrix,chrom_order)
  #chrom_oder_subset<-subset(chrom_order,)
  tmpm_base <- as.data.frame(matrix(nrow=length(sample_barcode),ncol=1))
  rownames(tmpm_base) <- sample_barcode
  
  for(chr in 1:length(chrom_order)){
    tmpm <- as.data.frame(final_cnv_list_matrix[,which(colnames(final_cnv_list_matrix)==chrom_order[chr])])
    if(!isEmpty(tmpm)){
      colnames(tmpm) <- chrom_order[chr]
      tmpm_base <- cbind(tmpm_base,tmpm,stringsAsFactors=T)
    }
  }
  final_cnv_list_matrix<-tmpm_base[,-1]
  bk <- c(seq(-0.5,1.9,by=0.01),seq(2,3,by=0.01))
  heatmap<-pheatmap(final_cnv_list_matrix,cluster_cols = F,show_rownames = F,
                    color = c(colorRampPalette(colors = c("#3f72af", "white"))(length(bk)/3*2),
                              colorRampPalette(colors = c("white","#d72323"))(length(bk)/3))
                   )
  
  ggsave(paste0("result/merge_data/before_doubletfinder/copyscat/",files[i],"_atac/",files[i],"_heatmap.png"),heatmap,width = 10, height = 20)
  
  #option2 AUTOMATICALLY IDENTIFY NON-NEOPLASTIC CELLS AND USE THESE FOR CONTROL
  #work best if tumor cellularity is <90%
  #run these after filtercells
  # nmf_results<-identifyNonNeoplastic(scData_collapse,estimatedCellularity = 0.90,methodHclust = "ward.D")
  # #warning:Recycling array of length 1 in array-vector arithmetic is deprecated Use c() or as.vector() instead.
  # #estimatedCellularity:expected cellularity of tumour
  # write.table(x=rownames_to_column(data.frame(nmf_results$cellAssigns),var="Barcode"),
  #             file=str_c(scCNVCaller$locPrefix,scCNVCaller$outPrefix,"nmf_clusters.csv"),
  #             quote=FALSE,row.names = FALSE,sep=",")
  # print(paste("Normal cluster is: ",nmf_results$clusterNormal))
  
  #part2B SmoothingCNV calls with clusters
  ##smoothedCNVList <- smoothClusters(scDataSampClusters, inputCNVList = final_cnv_list[[3]],percentPositive = 0.4,removeEmpty = F)
  #percentPositive: minimum percentage of a cluster to use to call positive for CNV, default is 0.5
  #`funs()` was deprecated in dplyr 0.8.0.
  
  #part3 : identify double minutes / amplification
  # library(compiler)
  # dmRead <- cmpfun(identifyDoubleMinutes)
  # dm_candidates<-dmRead(scData_k_norm,minCells=100,qualityCutoff2 = 100,minThreshold = 3) 
  # #minthread is a time-saving option that doesn't call changepoints on any cell with a maximum Z score less than the number you set.
  # 
  # cluster <- read.csv(paste0('result/after_doubletfinder/copyscat/',files[i],'_atac/',files[i],'_atac_nmf_clusters.csv'))
  # cluster<-cluster[which(cluster$nmf_results.cellAssigns==nmf_results$clusterNormal),]
  # cluster<-cluster$Barcode
  # mallignant<-data.frame(pbmc.atac.barcode,"mallignant")
  # colnames(mallignant)[2]<-"mallignant"
  # for(cluster_nums in 1:length(cluster)){
  #   mallignant$mallignant[which(mallignant$pbmc.atac.barcode == cluster[cluster_nums])]<-"normal"
  # }
  # rownames(mallignant)<-mallignant$pbmc.atac.barcode
  # pbmc.atac<-AddMetaData(pbmc.atac,metadata = mallignant)
  # p1 <- DimPlot(pbmc.atac, group.by = "mallignant", label = T)
  # p1
  # ggsave(paste0("result/after_doubletfinder/copyscat/",files[i],"_atac/",files[i],"_mallignant.png"),p1,width = 10, height = 10)
  # p2 <- DimPlot(pbmc.atac, group.by = "predicted.id", label = T)
  # p2
  # ggsave(paste0("result/after_doubletfinder/copyscat/",files[i],"_atac/",files[i],"_ori.png"),p2,width = 10, height = 10)
}
chr <- read.csv("result/merge_data/before_doubletfinder/copyscat/S02-primary_atac/S02-primary_atac_clean_cnv_cnv_scores.csv",header = T, sep = ",")
chr17q <- data.frame(chr$rowname,chr$chr17q)
chr17q_1<-as.data.frame(chr17q[,-1])
rownames(chr17q_1)<-chr17q$chr.rowname
pbmc.atac$chr17q <- chr17q_1
DimPlot(pbmc.atac,group.by="chr17q",cols = c('white','white',"blue","red"))
DimPlot(pbmc.atac,group.by="chr17q",cols = c("blue","red",'white','white'))
