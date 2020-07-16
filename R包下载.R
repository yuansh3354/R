# 因为我自己用的是 mac 电脑，所以如果有安装包的话首选直接使用终端+梯子下载各种数据和 R 包
# 除非代理不能用，不然不使用镜像下载

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

##  先用这个，这个是 R 的官方库，可能需要代理才比较好进去
install.packages("geomnet", repos='https://mran.microsoft.com/snapshot/2019-02-01/')

# 不行用这个
BiocManager::install("KEGG.db",ask = F,update = F)

# 如果是存放在 github 上的包就用这个
library(devtools)
devtools::install_github("OHDSI/Achilles")
