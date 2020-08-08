# 因为我自己用的是 mac 电脑，所以如果有安装包的话首选直接使用终端+梯子下载各种数据和 R 包
# 除非代理不能用，不然不使用镜像下载
# 这里有个问题注意一下，mac 上的终端下载数据还是 R 包速度都是相当的块，比 rstudio 快很多且稳定

rm(list = ls())   
options()$repos 
options()$BioC_mirror
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options()$repos 
options()$BioC_mirror

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

##  先用这个，这个是 R 包最通用的下载方法
## 一般常规来说，如果能挂梯子就不要用这种方法，太麻烦了
## repos 是指定镜像
install.packages("geomnet", repos='https://mran.microsoft.com/snapshot/2019-02-01/')

## 生物学分析包用这个下载
BiocManager::install("KEGG.db",ask = F,update = F)

## 如果是存放在 github 上的包就用这个
## 这个就是强烈建议用终端下载
library(devtools)
devtools::install_github("OHDSI/Achilles")
