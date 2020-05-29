rm(list = ls())   
options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
options("repos" = c(CRAN="http://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

# https://bioconductor.org/packages/release/bioc/html/GEOquery.html
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
options()$repos
##  先用这个
install.packages('gdata')
# 不行用这个
BiocManager::install("KEGG.db",ask = F,update = F)
#再不行最后用这个
library(devtools)
devtools::install_github("OHDSI/Achilles")
