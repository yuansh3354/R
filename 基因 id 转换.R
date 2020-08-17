library("clusterProfiler")
gene.df <- bitr(gene, fromType = "ENTREZID", #fromType是指你的数据ID类型是属于哪一类的
                toType = c("ENSEMBL", "SYMBOL"), #toType是指你要转换成哪种ID类型，可以写多种，也可以只写一种
                OrgDb = org.Hs.eg.db)#Orgdb是指对应的注释包是哪个
head(gene.df)
