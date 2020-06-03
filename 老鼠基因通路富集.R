options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
options("repos" = c(CRAN="http://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
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
  library(org.Mm.eg.db)
  library(GenomicFeatures)
  library(rtracklayer)
  library(biomaRt)
  library(glmnet)
  library(survival)
  library(Hmisc)
  library(clusterProfiler)
  library(KEGG.db)
}
setwd("/Users/yuansh/Downloads")
top = read.csv('OB_EDGER.csv',header = T)
gene = top[,1]
gene.df <- bitr(gene, fromType = "SYMBOL", #fromType是指你的数据ID类型是属于哪一类的
                toType = c("ENTREZID", "SYMBOL"), #toType是指你要转换成哪种ID类型，可以写多种，也可以只写一种
                OrgDb = org.Mm.eg.db)#Orgdb是指对应的注释包是哪个


if(T){
  logFC_t=1
  top$g= 0
  top[which(top$logFC<= -1),]$g = 'DOWN'
  top[which(top$logFC>= 1),]$g = 'UP'
  
  top = merge(top,gene.df)  
  gene_up= top[top$g == 'UP','ENTREZID'] 
  gene_down=top[top$g == 'DOWN','ENTREZID'] 
  gene_diff=c(gene_up,gene_down)
}
if(T){
  if(T){
    kk.up <- enrichKEGG(gene         = gene_up,
                        organism     = 'mmu',
                        universe     = gene_diff,
                        pvalueCutoff = 0.9,
                        qvalueCutoff =0.9,
                        use_internal_data = T)
    head(kk.up)[,1:6]
    barplot(kk.up )
    ggsave('kk.up.barplot.png')
    
    kk.down <- enrichKEGG(gene         =  gene_down,
                          organism     = 'mmu',
                          universe     = gene_diff,
                          pvalueCutoff = 0.9,
                          qvalueCutoff =0.9,
                          use_internal_data = T)
    head(kk.down)[,1:6]
    barplot(kk.down )
    ggsave('kk.down.barplot.png')
    
    kk.diff <- enrichKEGG(gene         = gene_diff,
                          organism     = 'mmu',
                          pvalueCutoff = 0.05,
                          use_internal_data = T)
    head(kk.diff)[,1:6]
    barplot(kk.diff )
    ggsave('kk.diff.barplot.png')
    
    kegg_diff_dt <- as.data.frame(kk.diff)
    kegg_down_dt <- as.data.frame(kk.down)
    kegg_up_dt <- as.data.frame(kk.up)
    down_kegg<-kegg_down_dt[kegg_down_dt$pvalue<0.05,];down_kegg$group=-1
    up_kegg<-kegg_up_dt[kegg_up_dt$pvalue<0.05,];up_kegg$group=1
    g_kegg=kegg_plot(up_kegg,down_kegg)
    print(g_kegg)
    
    ggsave(g_kegg,filename = 'kegg_up_down.png')
  }#kegg
}
write.csv(kegg_down_dt,'OB_EDGER_down.csv')
write.csv(kegg_up_dt,'OB_EDGER_up.csv')
write.csv(kegg_diff_dt,'OB_EDGER_diff.csv')
