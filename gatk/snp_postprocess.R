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
library(ggrepel)
library(corrplot)
library(pheatmap)
sample1<-"S07_pri"
sample2<-"S07_met"
sample <- "s07"
sample_se1<-"s07_pri"
sample_se2<-"s07_met"
load(paste0("data/merge_data/before_doubletfinder/snp_Rdata/snp_",sample,".rdata"))
pri_da<-as.data.frame(pri_se$seurat_clusters)
rownames(pri_da)<-paste0("pri_",rownames(pri_da))
pri_da$rownames<-rownames(pri_da)
pri_da$`pri_se$seurat_clusters`<-paste0("pri_",pri_da$`pri_se$seurat_clusters`)
pri_da<-as.data.frame(pri_da[which(rownames(pri_da)%in%pri_se_name),])
colnames(pri_da)[1]<-"seurat_cluster"

other_da<-as.data.frame(other_se$seurat_clusters)
rownames(other_da)<-paste0("met_",rownames(other_da))
other_da$rownames<-rownames(other_da)
other_da$`other_se$seurat_clusters`<-paste0("met_",other_da$`other_se$seurat_clusters`)
other_da<-as.data.frame(other_da[which(rownames(other_da)%in%other_se_name),])
colnames(other_da)[1]<-"seurat_cluster"
da<-rbind(pri_da,other_da)
df.empty <- data.frame(matrix(ncol = length(overlap), nrow = (length(pri_se_name)+length(other_se_name))))
rownames(df.empty)<-rownames(da)
df.empty <- cbind(seurat_cluster=da[,1],df.empty)
colnames(df.empty)[2:ncol(df.empty)]<-overlap
df.empty[1:5,1:5]
for(i in 1:length(overlap)){
  for(j in 1:2){
    if(file.info(paste0("result/merge_data/before_doubletfinder/snp/snp_",sample_se1,"/",overlap[i],"-pri1.txt"))$size>0 &
       file.info(paste0("result/merge_data/before_doubletfinder/snp/snp_",sample_se1,"/",overlap[i],"-pri2.txt"))$size>0){
      barcode<- read.table(paste0("result/merge_data/before_doubletfinder/snp/snp_",sample_se1,"/",overlap[i],"-pri",j,".txt"))
      barcode$V1<-paste0("pri_",barcode$V1)
      df.empty[which(rownames(df.empty)%in%barcode$V1),i+1]<-(3-j)
    }
  }
}

for(i in 1:length(overlap)){
  for(j in 1:2){
    if(file.info(paste0("result/merge_data/before_doubletfinder/snp/snp_",sample_se2,"/",overlap[i],"-met1.txt"))$size>0 &
       file.info(paste0("result/merge_data/before_doubletfinder/snp/snp_",sample_se2,"/",overlap[i],"-met2.txt"))$size>0){
      barcode<- read.table(paste0("result/merge_data/before_doubletfinder/snp/snp_",sample_se2,"/",overlap[i],"-met",j,".txt"))
      barcode$V1<-paste0("met_",barcode$V1)
      df.empty[which(rownames(df.empty)%in%barcode$V1),i+1]<-(3-j)
    }
  }
}
#Ref=2 Alt=1 other =0

df.empty[is.na(df.empty)]<-0

saveRDS(df.empty,paste0("result/merge_data/before_doubletfinder/snp/",sample,"_snp.rds"))
df.empty<-readRDS(paste0("result/merge_data/before_doubletfinder/snp/",sample,"_snp.rds"))
df.empty[1:5,1:5]
# df.empty2<-df.empty[which(rowSums(df.empty[,2:ncol(df.empty)]) > 0),]
#df.empty$group<-"pri"
#df.empty$group[(length(pri_se_name)+1):(length(pri_se_name)+length(rec_se_name))]<-"rec"

# pca1 <- prcomp(df.empty[,2:ncol(df.empty)],center = F,scale = F)
# saveRDS(pca1,"result/merge_data/before_doubletfinder/snp/s04_19037_pca.rds")
# pca1<-readRDS("result/merge_data/before_doubletfinder/snp/s04_19037_pca.rds")
# df1 <- pca1$x
# df1 <- as.data.frame(df1)
# summ1 <- summary(pca1)
# xlab1 <- paste0("PC1(",round(summ1$importance[2,1]*100,2),"%)")
# ylab1 <- paste0("PC2(",round(summ1$importance[2,2]*100,2),"%)")
# df1$group<-df.empty$seurat_cluster
# ggplot(data = df1,aes(x = PC1,y = PC2,color = group))+
#  # 添加置信椭圆
#   geom_point(size = 1)+
#   labs(x = xlab1,y = ylab1,color = "Condition",title = "PCA Scores Plot")+
#   guides(fill = "none")+
#   theme_bw()+
#   #scale_colour_manual(values = c("purple","orange","pink"))+
#   theme(plot.title = element_text(hjust = 0.5,size = 15),
#         axis.text = element_text(size = 11),axis.title = element_text(size = 13),
#         legend.text = element_text(size = 11),legend.title = element_text(size = 13),
#         plot.margin = unit(c(0.4,0.4,0.4,0.4),'cm'))
# ggsave(filename = "result/merge_data/before_doubletfinder/snp/s04_19037_PCA.png",width = 13,height = 10)
# 
# 
# annotation_col<-as.data.frame(df.empty[,19038])
# colnames(annotation_col)[1]<-"condition"
# rownames(annotation_col)<-rownames(df.empty)
# df.empty<-as.data.frame(t(df.empty))
# df.empty<-as.matrix(df.empty)
# df.empty2<-as.data.frame(df.empty[1:200,])
# df.empty2[1:200,]<-as.numeric(df.empty2[1:200,])
# df.empty2<-as.matrix(df.empty2)
# df.empty<-df.empty[-19038,]
# df.empty[1:19037,]<-as.numeric(df.empty[1:19037,])
# pdf("result/merge_data/before_doubletfinder/snp/s04_19037.pdf", width = 30,height = 50)
# Heatmap(
#   df.empty,cluster_columns = T,cluster_rows = T,col = c("blue","red","black","white")
# )
# dev.off()
rownames(merge@meta.data)
rownames<-as.data.frame(rownames(df.empty))
DimPlot(merge,group.by = "orig.ident",cells.highlight = rownames)

rownames <- as.data.frame(rownames(merge@meta.data))
colnames(rownames)[1]<-"barcode"
write.table(rownames,paste0("result/merge_data/before_doubletfinder/gatk/barcode_cancer/",sample,".txt"),row.names = F,col.names = F,quote = FALSE)


# df <- data.frame(matrix(ncol =length(overlap), nrow = nrow(table(df.empty$seurat_cluster))))
# rownames<-as.data.frame(table(df.empty$seurat_cluster))
# rownames(df)<-rownames$Var1
# colnames(df)<-overlap
# for(i in 1:length(overlap)){
#   for(j in 1:length(rownames(df))){
#     da<-data.frame(c0=0,c1=0,c2=0)
#     tmp<-as.data.frame(table(df.empty[,i+1][which(df.empty$seurat_cluster==rownames(df)[j])]))
#     tmp$Var1<-as.numeric(as.character(tmp$Var1))
#     for(k in 1:length(rownames(tmp))){
#       da[,tmp[k,1]+1]<-tmp[k,2]
#     }
#     da=da$c1/(da$c0+da$c1+da$c2)
#     df[j,i]=da
#   }
# }#Ref=2 Alt=1 other =0
# saveRDS(df,"result/merge_data/before_doubletfinder/snp/s04_19037_seurat_cluster.rds")#19037删去了11个snp（这11个snp没有alt）故最后只剩19026
# df<-readRDS("result/merge_data/before_doubletfinder/snp/s04_19037_seurat_cluster.rds")
#由原先的alt/（alt+ref+other）改为alt/（alt+ref）进行聚类
df <- data.frame(matrix(ncol =length(overlap), nrow = nrow(table(df.empty$seurat_cluster))))
rownames<-as.data.frame(table(df.empty$seurat_cluster))
rownames(df)<-rownames$Var1
colnames(df)<-overlap
for(i in 1:length(overlap)){
  for(j in 1:length(rownames(df))){
    da<-data.frame(c0=0,c1=0,c2=0)
    tmp<-as.data.frame(table(df.empty[,i+1][which(df.empty$seurat_cluster==rownames(df)[j])]))
    tmp$Var1<-as.numeric(as.character(tmp$Var1))
    for(k in 1:length(rownames(tmp))){
      da[,tmp[k,1]+1]<-tmp[k,2]
    }
    da=da$c1/(da$c1+da$c2)#alt/(ref+alt)
    df[j,i]=da
  }
}
saveRDS(df,paste0("result/merge_data/before_doubletfinder/snp/",sample,"_seurat_cluster_alt+ref.rds"))
df<-readRDS(paste0("result/merge_data/before_doubletfinder/snp/",sample,"_seurat_cluster_alt+ref.rds"))
df[is.na(df)]<-0
pca <- prcomp(df,center = F,scale = F)
df1 <- pca$x
df1 <- as.data.frame(df1)
summ1 <- summary(pca)
xlab1 <- paste0("PC1(",round(summ1$importance[2,1]*100,2),"%)")
ylab1 <- paste0("PC2(",round(summ1$importance[2,2]*100,2),"%)")
df1$group<-rownames(df)
ggplot(data = df1,aes(x = PC1,y = PC2,color = group))+
  # 添加置信椭圆
  geom_point(size = 1)+
labs(x = xlab1,y = ylab1,color = "Condition",title = "PCA Scores Plot")+
  guides(fill = "none")+
  theme_bw()+
  geom_text_repel(aes(PC1, PC2, label = group))+
  #scale_colour_manual(values = c("purple","orange","pink"))+
  theme(plot.title = element_text(hjust = 0.5,size = 15),
        axis.text = element_text(size = 11),axis.title = element_text(size = 13),
        legend.text = element_text(size = 11),legend.title = element_text(size = 13),
        plot.margin = unit(c(0.4,0.4,0.4,0.4),'cm'))
ggsave(filename = paste0("result/merge_data/before_doubletfinder/snp/",sample,"_PCA_seurat_clusteralt+ref.png"),width = 13,height = 10)

df_t<-as.data.frame(t(df))
cor_data <- cor(df_t,method="spearman")
corrplot(cor_data, method="color",order = 'hclust',addrect=4,rect.col="red")
#"pri_4","pri_5","pri_3","pri_1","pri_0","pri_2","rec_1","rec_4","rec_0","rec_2","rec_8","rec_3","rec_5","rec_7","rec_6","rec_9"
DimPlot(merge,group.by = "seurat_clusters_label",label = T,repel = T,
        cells.highlight = rownames(merge@meta.data)[which(merge$seurat_clusters_label%in%c("pri_7"))],
        sizes.highlight =0.5
)
DimPlot(merge,group.by = "seurat_clusters_label",label = T,repel = T,
        sizes.highlight =0.5
)
pheatmap(df_t, cluster_rows = FALSE,show_rownames = F, clustering_distance_cols = "manhattan",color = c("white","#71c9ce"))

percent<- data.frame(matrix(nrow=1,ncol=ncol(df.empty)))
colnames(percent)<-colnames(df)
for(i in 1:ncol(percent)-1){
  da<-data.frame(c0=0,c1=0,c2=0)
  tmp<-as.data.frame(table(df.empty[,i+1]))
  tmp$Var1<-as.numeric(as.character(tmp$Var1))
  for(k in 1:length(rownames(tmp))){
    da[,tmp[k,1]+1]<-tmp[k,2]
  }
  da=da$c1/(da$c1+da$c2)#alt/(ref+alt)
  percent[1,i]=da
}
percent<-as.data.frame(t(percent))
ggplot(percent, aes(x = V1)) +geom_histogram(fill = "lightblue", colour = "black")

# pca <- prcomp(df,center = F,scale = F)
# df1 <- pca$x
# df1 <- as.data.frame(df1)
# summ1 <- summary(pca)
# xlab1 <- paste0("PC1(",round(summ1$importance[2,1]*100,2),"%)")
# ylab1 <- paste0("PC2(",round(summ1$importance[2,2]*100,2),"%)")
# df1$group<-rownames(df2)
# ggplot(data = df1,aes(x = PC1,y = PC2,color = group))+
#   # 添加置信椭圆
#   geom_point(size = 1)+
#   labs(x = xlab1,y = ylab1,color = "Condition",title = "PCA Scores Plot")+
#   guides(fill = "none")+
#   theme_bw()+
#   geom_text_repel(aes(PC1, PC2, label = group))+
#   #scale_colour_manual(values = c("purple","orange","pink"))+
#   theme(plot.title = element_text(hjust = 0.5,size = 15),
#         axis.text = element_text(size = 11),axis.title = element_text(size = 13),
#         legend.text = element_text(size = 11),legend.title = element_text(size = 13),
#         plot.margin = unit(c(0.4,0.4,0.4,0.4),'cm'))
#   ggsave(filename = "result/merge_data/before_doubletfinder/snp/s04_19037_PCA_seurat_cluster.png",width = 13,height = 10)
# 
#   

# #Ref=2 Alt=1 other =0
# percent<-as.data.frame(t(percent))
# ggplot(percent, aes(x = V1)) +geom_histogram(fill = "lightblue", colour = "black")
# 
# df_t<-as.data.frame(t(df))
# cor_data <- cor(df_t,method="spearman")
# corrplot(cor_data, method="color",order = 'hclust')
# pheatmap(df_t, cluster_rows = FALSE,show_rownames = F, clustering_distance_cols = "manhattan")

