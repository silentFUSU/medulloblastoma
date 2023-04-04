rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-redhat-linux-gnu-library/4.2","/storage/zhangyanxiaoLab/suzhuojie/miniconda3/envs/R/lib/R/library"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/medulloblastoma/")
#source("code/library.R")
library(AnnoProbe)
library(vegan)
library(tibble)
#library(infercnv)
library(stringr)
library(ggplot2)
library(corrplot)
file=c("s02_pri","s02_rec","s04_pri","s04_rec","s05_pri","s05_met","s06_pri","s06_met","s07_pri","s07_met",
       "s08_pri","s08_met","s09_pri","s09_met","s10_pri","s10_met","s12_pri","s12_rec","s13_pri","s13_met")
for(i in 1:20){
  filename<-file[i]
  da<-read.csv(paste0("result/merge_data/before_doubletfinder/infercnv/",filename,"_expr.csv"),header=TRUE,row.names = 1)
  row_mean = apply(da,1,mean)
  da<-add_column(da,mean=row_mean)
  da<-add_column(da,rowname=rownames(da))
  da_mean<-da[,(ncol(da)-1):ncol(da)]
  colnames(da)[1]<-filename
  write.csv(da_mean,paste0("result/merge_data/before_doubletfinder/infercnv/",filename,"_expr_mean.csv"))
}

s02_pri<-read.csv(paste0("result/merge_data/before_doubletfinder/infercnv/s02_pri_expr_mean.csv"),header=TRUE,row.names = 1)
s02_rec<-read.csv(paste0("result/merge_data/before_doubletfinder/infercnv/s02_rec_expr_mean.csv"),header=TRUE,row.names = 1)
s04_pri<-read.csv(paste0("result/merge_data/before_doubletfinder/infercnv/s04_pri_expr_mean.csv"),header=TRUE,row.names = 1)
s04_rec<-read.csv(paste0("result/merge_data/before_doubletfinder/infercnv/s04_rec_expr_mean.csv"),header=TRUE,row.names = 1)
s05_pri<-read.csv(paste0("result/merge_data/before_doubletfinder/infercnv/s05_pri_expr_mean.csv"),header=TRUE,row.names = 1)
s05_met<-read.csv(paste0("result/merge_data/before_doubletfinder/infercnv/s05_met_expr_mean.csv"),header=TRUE,row.names = 1)
s06_pri<-read.csv(paste0("result/merge_data/before_doubletfinder/infercnv/s06_pri_expr_mean.csv"),header=TRUE,row.names = 1)
s06_met<-read.csv(paste0("result/merge_data/before_doubletfinder/infercnv/s06_met_expr_mean.csv"),header=TRUE,row.names = 1)
s07_pri<-read.csv(paste0("result/merge_data/before_doubletfinder/infercnv/s07_pri_expr_mean.csv"),header=TRUE,row.names = 1)
s07_met<-read.csv(paste0("result/merge_data/before_doubletfinder/infercnv/s07_met_expr_mean.csv"),header=TRUE,row.names = 1)
s08_pri<-read.csv(paste0("result/merge_data/before_doubletfinder/infercnv/s08_pri_expr_mean.csv"),header=TRUE,row.names = 1)
s08_met<-read.csv(paste0("result/merge_data/before_doubletfinder/infercnv/s08_met_expr_mean.csv"),header=TRUE,row.names = 1)
s09_pri<-read.csv(paste0("result/merge_data/before_doubletfinder/infercnv/s09_pri_expr_mean.csv"),header=TRUE,row.names = 1)
s09_met<-read.csv(paste0("result/merge_data/before_doubletfinder/infercnv/s09_met_expr_mean.csv"),header=TRUE,row.names = 1)
s10_pri<-read.csv(paste0("result/merge_data/before_doubletfinder/infercnv/s10_pri_expr_mean.csv"),header=TRUE,row.names = 1)
s10_met<-read.csv(paste0("result/merge_data/before_doubletfinder/infercnv/s10_met_expr_mean.csv"),header=TRUE,row.names = 1)
s12_pri<-read.csv(paste0("result/merge_data/before_doubletfinder/infercnv/s12_pri_expr_mean.csv"),header=TRUE,row.names = 1)
s12_rec<-read.csv(paste0("result/merge_data/before_doubletfinder/infercnv/s12_rec_expr_mean.csv"),header=TRUE,row.names = 1)
s13_pri<-read.csv(paste0("result/merge_data/before_doubletfinder/infercnv/s13_pri_expr_mean.csv"),header=TRUE,row.names = 1)
s13_met<-read.csv(paste0("result/merge_data/before_doubletfinder/infercnv/s13_met_expr_mean.csv"),header=TRUE,row.names = 1)

merge<-Reduce(function(x, y) merge(x, y, by="rowname",all=F), list(s05_pri,s05_met,s06_pri,s06_met,s07_pri,s07_met,
                                                  s08_pri,s08_met,s09_pri,s09_met,s10_pri,s10_met,s13_pri,s13_met,s02_pri,s02_rec,s04_pri,s04_rec,s12_pri,s12_rec))
colnames(merge) <- c("rowname","s05_pri","s05_met","s06_pri","s06_met","s07_pri","s07_met",
                     "s08_pri","s08_met","s09_pri","s09_met","s10_pri","s10_met","s13_pri","s13_met","s02_pri","s02_rec","s04_pri","s04_rec","s12_pri","s12_rec")
rownames(merge) <- merge$rowname
merge<-merge[,-1]
cor_data <- cor(merge,method="spearman")
corrplot(cor_data, method="color",order = 'hclust')
corrplot(cor_data, method="color")
ggsave("result/merge_data/before_doubletfinder/infercnv/bulk_infercnv_correlation.png",width = 10,height = 10)
