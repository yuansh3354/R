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
}
```

##### 2. 导入包

```R
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
```

### 下载不需要归一化的数据
 不同的平台导入的R包不一样
|平台|注释包|
|:---:|:---:|
| GPL201  |  hgfocus|
| GPL96   | hgu133a|
| GPL571  |  hgu133a2|
| GPL97   | hgu133b|
| GPL570  |  hgu133plus2|
| GPL13667|    hgu219|
| GPL8300 |   hgu95av2|
| GPL91   | hgu95av2|
| GPL92   | hgu95b|
| GPL93   | hgu95c|
| GPL94   | hgu95d|
| GPL95   | hgu95e|
| GPL887  |  hgug4110b| 
| GPL886  |  hgug4111a | 
| GPL1708 |   hgug4112a |
| GPL13497|    HsAgilentDesign026652 |
| GPL6244 |   hugene10sttranscriptcluster |
| GPL11532|    hugene11sttranscriptcluster |
| GPL6097 |   illuminaHumanv1 |
| GPL6102 |   illuminaHumanv2|
| GPL6947 |   illuminaHumanv3 |
| GPL10558|    illuminaHumanv4|
| GPL6885 |   illuminaMousev2 |
| GPL81   | mgu74av2 |
| GPL82   | mgu74bv2 |
| GPL83   | mgu74cv2 |
| GPL339  |  moe430a|
| GPL6246 |   mogene10sttranscriptcluster |
| GPL340  |  mouse4302|
| GPL1261 |   mouse430a2|
| GPL8321 |   mouse430a2|

```R
file = 'GSEXXXXX' # 要下载的数据
setwd('/Volumes/Lexar/ZG/GSE19188分析结果')
if(T){
  gset <- getGEO(file, destdir=".",
                 AnnotGPL = F,     ## 注释文件
                 getGPL = F)
  a=gset[[1]]
  tab <- select(XXXXX.db, keys = keys(XXXXX.db), columns = c("ENTREZID")) # 这里注意一下不同的平台信息要用不同的包
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
  write.csv(geneExpr,file = paste(file,'.csv',sep = ''))
}
paste(file,'_clinic.csv',sep = '')
pd = read.csv(paste(file,'_clinic.csv',sep = ''),header = T, row.names = 1)
df_expr = read.csv(paste(file,'.csv',sep = ''),header = T, row.names = 1)
```
