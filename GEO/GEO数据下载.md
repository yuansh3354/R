##### 1. 设置镜像

```R
if(T){
  options()$repos 
  options()$BioC_mirror
  options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")
  options("repos" = c(CRAN="http://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
  #setwd('/Volumes/Lexar/ZG/')
  #BiocManager::install('randomForestSRC')
  #install.packages('包')
  #library(devtools)
  #install_github("jmzeng1314/AnnoProbe")
}
```

##### 2. 导入包

```R
if(T){
  library(AnnoProbe)
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
```

```R
file = 'GSE14520' # 要下载的数据
library(hthgu133a)
gset <- getGEO(file, destdir=".",
               AnnotGPL = F,     ## 注释文件
               getGPL = F)
               
a=gset[[1]]
gpl='GPL3921'
tab=idmap(gpl,type = 'pipe') 
e <- exprs(a)
e = as.data.frame(e)

# 这时候把e打出来看一下是否需要进行归一化等等的处理
geneExpr <- t(sapply(split(tab[,1], tab[,2]), function(ids){
  colMeans(e[ids,,drop=FALSE])
}))

geneExpr = na.omit(geneExpr)
pd=pData(a) 

### 根据需求决定是否需要归一化
#geneExpr=geneExpr[apply(geneExpr,1,sd)>0,]
#geneExpr=log2(geneExpr+1)
#geneExpr=normalizeBetweenArrays(geneExpr)

# 1. 将样本信息导出
write.csv(pd,file = paste(file,'_clinic.csv',sep = ''))
write.csv(geneExpr,file = paste(file,'_exp.csv',sep = ''))
paste(file,'_clinic.csv',sep = '')
pd = read.csv(paste(file,'_clinic.csv',sep = ''),header = T, row.names = 1)
df_expr = read.csv(paste(file,'.csv',sep = ''),header = T, row.names = 1)
```
