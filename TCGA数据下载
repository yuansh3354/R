
```R
rm(list = ls()) 
getwd()
setwd('/Volumes/Lexar/ff_20200421')
tof = TRUE # 运行if代码
tof = FALSE # 注释if代码
```

### 设置下载R包的环境变量

```R
tof = T
if(tof){
  #清空当前工作空间变量  
  options()$repos  
  #查看当前工作空间默认的下载包路径
  options()$BioC_mirror 
  #查看使用BioCManager下载包的默认路径
  options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") 
  # 指定使用BioCManager下载的路径
  options("repos" = c(CRAN="http://mirrors.tuna.tsinghua.edu.cn/CRAN/")) 
  options("repos" = c(CRAN="http://mirrors.cloud.tencent.com/CRAN/")
  # 指定使用install.packages下载包的路径
  options()$repos 
  options()$BioC_mirror
  # https://bioconductor.org/packages/release/bioc/html/GEOquery.html
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager") 
  #判断是否存在BiocManger包，没有的话下载该包
} 
```

### 下载TCGA数据分析必要的包
### 如果已经有了就不要下载了
```R
tof = F
if(tof){
  if(!require("KEGG.db")) BiocManager::install("KEGG.db",ask = F,update = F)
  if(!require("GSEABase")) BiocManager::install("GSEABase",ask = F,update = F)
  if(!require("GSVA")) BiocManager::install("GSVA",ask = F,update = F)
  if(!require("clusterProfiler")) BiocManager::install("clusterProfiler",ask = F,update = F)
  if(!require("GEOquery")) BiocManager::install("GEOquery",ask = F,update = F)
  if(!require("limma")) BiocManager::install("limma",ask = F,update = F)
  if(!require("impute")) BiocManager::install("impute",ask = F,update = F)
  if(!require("genefu")) BiocManager::install("genefu",ask = F,update = F)
  if(!require("org.Hs.eg.db")) BiocManager::install("org.Hs.eg.db",ask = F,update = F)
  if(!require("hgu133plus2.db")) BiocManager::install("hgu133plus2.db",ask = F,update = F)
  if(!require("ConsensusClusterPlus")) BiocManager::install("ConsensusClusterPlus",ask = F,update = F)
  
  BiocManager::install(c('airway','DESeq2','edgeR','limma'),
                       ask = F,update = F)
  
  if(! require("maftools")) BiocManager::install("maftools",ask = F,update = F)
  if(! require("genefilter")) BiocManager::install("genefilter",ask = F,update = F)
  
  if(! require("CLL")) BiocManager::install("CLL",ask = F,update = F)
  if(! require("org.Hs.eg.db")) BiocManager::install('org.Hs.eg.db',ask = F,update = F)
  
  if(! require("maftools")) BiocManager::install("maftools",ask = F,update = F)
  if(! require("RTCGA")) BiocManager::install("RTCGA",ask = F,update = F)
  if(! require("RTCGA.clinical")) BiocManager::install("RTCGA.clinical",ask = F,update = F)
  # https://bioconductor.org/packages/3.6/data/experiment/src/contrib/RTCGA.clinical_20151101.8.0.tar.gzn
  if(! require("RTCGA.miRNASeq")) BiocManager::install("RTCGA.miRNASeq",ask = F,update = F)
  
} 
```

### 自定义批量下载R包
### 输入自己需要下载的包

```R
packages_names = c('TCGAbiolinks')
BiocManager::install(packages_names,
                     ask = F,update = F)
```

### 导入R包

```R
tof = T
if(tof){
  library(TCGAbiolinks)
}
```

### 查看数据

```R
tof = T
if(tof){
  # 1.project； 要下载的疾病种类
  # 2.data.category；数据类别
  # 3.data.type； 数据类型
  # 4.workflow.type； 数据处理方式
  # 5.legacy = FALSE；false的话是用参考基因组hg19 true是hg38
  # 6.access；一般用open
  # 7.platform；
  # 8.file.type；
  # 9.barcode；只下载单个样本的时候用
  # 10.experimental.strategy；
  # 11.sample.type
  
  project_list = TCGAbiolinks:::getGDCprojects()$project_id # 首先查看一下都有哪些数据
  project = project_list[grep('LIHC',project_list)] # 用grep获取目标数据
  project
  TCGAbiolinks:::getProjectSummary(project) # 查看数据类型
}
```
### 下载表达谱数据
```R
tof = F
if(tof){
  query <- GDCquery(project = project, 
                    legacy = FALSE, 
                    experimental.strategy = "RNA-Seq", 
                    data.category = "Transcriptome Profiling", 
                    data.type = "Gene Expression Quantification", 
                    workflow.type = "HTSeq - Counts")
  
  #再使用命令GDCdownload(）下载
  GDCdownload(query)
}
```

### 下载临床数据

```R
tof = T
if(tof){
  clinical <- GDCquery_clinic(project =project, type = "clinical")
  head(clinical)
  write.csv(clinical,file= paste(project,'clinical.csv',sep = ""))
}
```


