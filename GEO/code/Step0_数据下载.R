### ---
### Title: "Untitled"
### Create: "Yuansh"
### Date: "5/01/2020"
### Email: yuansh3354@163.com
### output: html_document
### ---


### step0 准备

##### 1. 设置镜像
if(T){
  options()$repos 
  options()$BioC_mirror
  options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")
  options("repos" = c(CRAN="http://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
  #setwd('/Volumes/Lexar/ZG/')
  #BiocManager::install('randomForestSRC')
  #install.packages('包')
}

##### 2. 导入包
if(T){
  library(limma)
  library(GEOquery)
  library(org.Hs.eg.db)
  library(RColorBrewer)
  library(AnnotationDbi)
  library(affy)
  library(gcrma)
  library(stringr)
  library(hgu133plus2.db )
  library(org.Hs.eg.db)
  library(GenomicFeatures)
  library(rtracklayer)
  library(biomaRt)
  library(glmnet)
  library(survival)
  library(Hmisc)
  library(clusterProfiler)
}

# 训练
file = 'GSE19188'
setwd('/Volumes/Lexar/ZG/GSE19188分析结果')
if(T){
  gset <- getGEO(file, destdir=".",
                 AnnotGPL = F,     ## 注释文件
                 getGPL = F)
  a=gset[[1]]
  tab <- select(hgu133plus2.db, keys = keys(hgu133plus2.db), columns = c("ENTREZID"))
  e <- exprs(a)
  geneExpr <- t(sapply(split(tab[,1], tab[,2]), function(ids){
    colMeans(e[ids,,drop=FALSE])
  }))
  pd=pData(a) 
  # 1. 将样本信息导出
  write.csv(pd,file = paste(file,'_clinic.csv',sep = ''))
  write.csv(geneExpr,file = paste(file,'.csv',sep = ''))
  rm(a,gset)
}
paste(file,'_clinic.csv',sep = '')
pd = read.csv(paste(file,'_clinic.csv',sep = ''),header = T, row.names = 1)
df_expr = read.csv(paste(file,'.csv',sep = ''),header = T, row.names = 1)

# 验证集1 
file = 'GSE50081' # GSE50081 GSE31210
if(T){
  gset <- getGEO(file, destdir=".",
                 AnnotGPL = F,     ## 注释文件
                 getGPL = F)
  
  a=gset[[1]]
  tab <- select(hgu133plus2.db, keys = keys(hgu133plus2.db), columns = c("ENTREZID"))
  e <- exprs(a)
  e = na.omit(e)
  geneExpr <- t(sapply(split(tab[,1], tab[,2]), function(ids){
    colMeans(e[ids,,drop=FALSE])
  }))
  #gse31210需要归一化处理
  # geneExpr=geneExpr[apply(geneExpr,1,sd)>0,]
  # geneExpr=log2(geneExpr+1)
  # geneExpr=normalizeBetweenArrays(geneExpr)
  # pd=pData(a) 
  # 1. 将样本信息导出
  write.csv(pd,file = paste(file,'_clinic.csv',sep = ''))
  write.csv(geneExpr,file = paste(file,'.csv',sep = ''))
  rm(a,gset)
}
paste(file,'_clinic.csv',sep = '')
pd = read.csv(paste(file,'_clinic.csv',sep = ''),header = T, row.names = 1)
df_expr = read.csv(paste(file,'.csv',sep = ''),header = T, row.names = 1)

# 验证集2
file = 'GSE31210' # GSE50081 GSE31210
if(T){
  gset <- getGEO(file, destdir=".",
                 AnnotGPL = F,     ## 注释文件
                 getGPL = F)
  
  a=gset[[1]]
  tab <- select(hgu133plus2.db, keys = keys(hgu133plus2.db), columns = c("ENTREZID"))
  e <- exprs(a)
  e = na.omit(e)
  geneExpr <- t(sapply(split(tab[,1], tab[,2]), function(ids){
    colMeans(e[ids,,drop=FALSE])
  }))
  gse31210需要归一化处理
  geneExpr=geneExpr[apply(geneExpr,1,sd)>0,]
  geneExpr=log2(geneExpr+1)
  geneExpr=normalizeBetweenArrays(geneExpr)
  pd=pData(a) 
  #1. 将样本信息导出
  write.csv(pd,file = paste(file,'_clinic.csv',sep = ''))
  write.csv(geneExpr,file = paste(file,'.csv',sep = ''))
  rm(a,gset)
}
paste(file,'_clinic.csv',sep = '')
pd = read.csv(paste(file,'_clinic.csv',sep = ''),header = T, row.names = 1)
df_expr = read.csv(paste(file,'.csv',sep = ''),header = T, row.names = 1)



