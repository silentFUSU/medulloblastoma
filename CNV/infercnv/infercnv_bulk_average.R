rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-redhat-linux-gnu-library/4.2","/storage/zhangyanxiaoLab/suzhuojie/miniconda3/envs/R/lib/R/library"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/medulloblastoma/")
#source("code/library.R")
library(AnnoProbe)
library(vegan)
library(tibble)
library(infercnv)
library(stringr)
library(ggplot2)

file=c("s02_pri","s02_rec","s04_pri","s04_rec","s05_pri","s05_met","s06_pri","s06_met","s07_pri","s07_met",
       "s08_pri","s08_met","s09_pri","s09_met","s10_pri","s10_met","s12_pri","s12_rec","s13_pri","s13_met")


for(i in 1:20){
  filename<-file[i]
  da<-read.csv(paste0("result/merge_data/before_doubletfinder/infercnv/",filename,"_expr.csv"),header=TRUE,row.names = 1)
  geneInfor=annoGene(rownames(da),"SYMBOL",'human')
  geneInfor=geneInfor[c(1,4:6)] 
  geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
  row_mean = apply(da,1,mean)
  da<-add_column(da,chr=geneInfor$chr)
  da<-add_column(da,mean=row_mean)
  da<-add_column(da,rowname=rownames(da))
  rowname<-rownames(da)
  da$rowname<-factor(da$rowname,levels = rowname)
  da$chr<-factor(da$chr,levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12",
                                 "chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22"))
  p1<-ggplot(data = da,aes(x=rowname,y=mean,colour=chr))+
    geom_point(size=1)+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+
    theme(axis.title.y=element_text(size=20))+ylim(0.8,1.2)+
    geom_hline(aes(yintercept=1))
  ggsave(paste0("result/merge_data/before_doubletfinder/infercnv/",filename,"_chr_annotation.png"),p1,width = 30, height = 5)
}

# file=c("s01_pri_rna","s02_pri_rna","s02_rec_rna","s03_pri_rna","s04_pri_rna",
#        "s04_rec_rna","s05_met_rna","s05_pri_rna","s06_met_rna","s06_pri_rna",
#        "s07_met_rna","s07_pri_rna","s08_met_rna","s08_pri_rna","s09_pri_rna",
#        "s10_pri_rna","s12_pri_rna","s12_rec_rna")
#file2=c("s04_rec_rna","s06_met_rna","s06_pri_rna","s07_met_rna","s09_pri_rna","s08_met_rna","s08_pri_rna","s12_pri_rna")
# for(i in 2:20){
#   filename<-file[i]
#   infercnv_obj<- readRDS(paste0("result/merge_data/before_doubletfinder/infercnv/",filename,"add_all_normal/preliminary.infercnv_obj"))
#   da<-infercnv_obj@expr.data
#   da<-da[,which(str_detect(colnames(da),"cancer"))]
#   write.csv(da,paste0("result/merge_data/before_doubletfinder/infercnv/",filename,"_expr.csv"))
# }
# s02_pri<-read.csv(paste0("result/merge_data/before_doubletfinder/infercnv/s02_pri_expr.csv"),header=TRUE,row.names = 1)
# s02_rec<-read.csv(paste0("result/merge_data/before_doubletfinder/infercnv/s02_rec_expr.csv"),header=TRUE,row.names = 1)
# s04_pri<-read.csv(paste0("result/merge_data/before_doubletfinder/infercnv/s04_pri_expr.csv"),header=TRUE,row.names = 1)
# s04_rec<-read.csv(paste0("result/merge_data/before_doubletfinder/infercnv/s04_rec_expr.csv"),header=TRUE,row.names = 1)
# s05_pri<-read.csv(paste0("result/merge_data/before_doubletfinder/infercnv/s05_pri_expr.csv"),header=TRUE,row.names = 1)
# s05_met<-read.csv(paste0("result/merge_data/before_doubletfinder/infercnv/s05_met_expr.csv"),header=TRUE,row.names = 1)
# s06_pri<-read.csv(paste0("result/merge_data/before_doubletfinder/infercnv/s06_pri_expr.csv"),header=TRUE,row.names = 1)
# s06_met<-read.csv(paste0("result/merge_data/before_doubletfinder/infercnv/s06_met_expr.csv"),header=TRUE,row.names = 1)
# s07_pri<-read.csv(paste0("result/merge_data/before_doubletfinder/infercnv/s07_pri_expr.csv"),header=TRUE,row.names = 1)
# s07_met<-read.csv(paste0("result/merge_data/before_doubletfinder/infercnv/s07_met_expr.csv"),header=TRUE,row.names = 1)
# s08_pri<-read.csv(paste0("result/merge_data/before_doubletfinder/infercnv/s08_pri_expr.csv"),header=TRUE,row.names = 1)
# s08_met<-read.csv(paste0("result/merge_data/before_doubletfinder/infercnv/s08_met_expr.csv"),header=TRUE,row.names = 1)
# s09_pri<-read.csv(paste0("result/merge_data/before_doubletfinder/infercnv/s09_pri_expr.csv"),header=TRUE,row.names = 1)
# s09_met<-read.csv(paste0("result/merge_data/before_doubletfinder/infercnv/s09_met_expr.csv"),header=TRUE,row.names = 1)
# s10_pri<-read.csv(paste0("result/merge_data/before_doubletfinder/infercnv/s10_pri_expr.csv"),header=TRUE,row.names = 1)
# s10_met<-read.csv(paste0("result/merge_data/before_doubletfinder/infercnv/s10_met_expr.csv"),header=TRUE,row.names = 1)
# s12_pri<-read.csv(paste0("result/merge_data/before_doubletfinder/infercnv/s12_pri_expr.csv"),header=TRUE,row.names = 1)
# s12_rec<-read.csv(paste0("result/merge_data/before_doubletfinder/infercnv/s12_rec_expr.csv"),header=TRUE,row.names = 1)
# s13_pri<-read.csv(paste0("result/merge_data/before_doubletfinder/infercnv/s13_pri_expr.csv"),header=TRUE,row.names = 1)
# s13_met<-read.csv(paste0("result/merge_data/before_doubletfinder/infercnv/s13_met_expr.csv"),header=TRUE,row.names = 1)
# row_mean = apply(s02_pri,1,mean)
# s02_pri<-add_column(s02_pri,rowname=rownames(s02_pri))
# s02_pri<-add_column(s02_pri,s02_pri=row_mean)
# s02_pri_mean<-as.data.frame(s02_pri[,c((ncol(s02_pri)-1):ncol(s02_pri))])
# 
# row_mean = apply(s02_rec,1,mean)
# s02_rec<-add_column(s02_rec,rowname=rownames(s02_rec))
# s02_rec<-add_column(s02_rec,s02_rec=row_mean)
# s02_rec_mean<-as.data.frame(s02_rec[,c((ncol(s02_rec)-1):ncol(s02_rec))])
# 
# row_mean = apply(s04_pri,1,mean)
# s04_pri<-add_column(s04_pri,rowname=rownames(s04_pri))
# s04_pri<-add_column(s04_pri,s04_pri=row_mean)
# s04_pri_mean<-as.data.frame(s04_pri[,c((ncol(s04_pri)-1):ncol(s04_pri))])
# 
# row_mean = apply(s04_rec,1,mean)
# s04_rec<-add_column(s04_rec,rowname=rownames(s04_rec))
# s04_rec<-add_column(s04_rec,s04_rec=row_mean)
# s04_rec_mean<-as.data.frame(s04_rec[,c((ncol(s04_rec)-1):ncol(s04_rec))])
# 
# row_mean = apply(s05_pri,1,mean)
# s05_pri<-add_column(s05_pri,rowname=rownames(s05_pri))
# s05_pri<-add_column(s05_pri,s05_pri=row_mean)
# s05_pri_mean<-as.data.frame(s05_pri[,c((ncol(s05_pri)-1):ncol(s05_pri))])
# 
# row_mean = apply(s05_met,1,mean)
# s05_met<-add_column(s05_met,rowname=rownames(s05_met))
# s05_met<-add_column(s05_met,s05_met=row_mean)
# s05_met_mean<-as.data.frame(s05_met[,c((ncol(s05_met)-1):ncol(s05_met))])
# 
# s02_s05<-merge(s02_pri_mean,c(s02_rec_mean,s04_pri_mean,s04_rec_mean,s05_pri_mean,s05_met_mean),by='rowname',all=F)
# da<-merge(s02_pri_mean,s02_rec_mean,by='rowname',all=F)
# da<-merge(da,s04_pri_mean,by='rowname',all=F)
# da<-merge(da,s04_rec_mean,by='rowname',all=F)
