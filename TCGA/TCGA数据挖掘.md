### 导入R包
```R
rm(list=ls())
library("FactoMineR")
library("factoextra")
library(GSEABase)
library(GSVA)
library(clusterProfiler)
library(genefu)
library(ggplot2)
library(ggpubr)
library(hgu133plus2.db)
library(limma)
library(org.Hs.eg.db)
library(pheatmap)
library(stringr)
library(devtools) 
library(TCGAbiolinks)
library(survival)
library(survminer)
```
### 导入数据
```
expset=as.data.frame(read.table('COADREAD_exp.txt',header = T, row.names = 1))
```
如果直接导入R，则样本名中的“_”会被“.”所替换
所以需要列名处理一下
```
t=str_split(colnames(expset),'\\.',simplify = T) 
colnames(expset)=paste(t[,1],t[,2],t[,3],t[,4],sep = '-')
head(expset)
clinic=as.data.frame(read.csv('clinic.csv',header = T))

survival_table=as.data.frame(read.csv('COADREAD_survival.csv',header = T))
cat('表达谱样本数: ',dim(expset)[2],sep='\n')
cat('临床信息样本数: ',dim(clinic)[1],sep='\n')
cat('生存信息样本数: ',dim(survival_table)[1],sep='\n')


a = colnames(expset)
b = clinic$ID
c = intersect(b,a)
ix = c
length(c)
rownames(clinic)=clinic$ID
rownames(survival_table)=survival_table$sample

```

### check一下
```
dat=expset
exprSet=dat
pheatmap::pheatmap(cor(exprSet)) 
# 组内的样本的相似性应该是要高于组间的！
colD=data.frame(group_list=group_list)
rownames(colD)=colnames(exprSet)
pheatmap::pheatmap(cor(exprSet),
                   annotation_col = colD,
                   show_rownames = F,
                   filename = 'cor_all.png')
dim(exprSet)
exprSet=exprSet[apply(exprSet,1, function(x) sum(x>1) > 5),]
dim(exprSet)
​
exprSet=log(edgeR::cpm(exprSet)+1)
dim(exprSet)
exprSet=exprSet[names(sort(apply(exprSet, 1,mad),decreasing = T)[1:500]),]
dim(exprSet)
M=cor(log2(exprSet+1)) 
pheatmap::pheatmap(M,annotation_col = colD)
pheatmap::pheatmap(M,
                   show_rownames = F,
                   annotation_col = colD,
                   filename = 'cor_top500.png')
                   
                   
group_list = clinic$pathologic_stage
table(group_list)
# 每次都要检测数据
dat[1:4,1:4]
## 下面是画PCA的必须操作，需要看说明书。
dat=t(dat)#画PCA图时要求是行名时样本名，列名时探针名，因此此时需要转换
dat=as.data.frame(dat)#将matrix转换为data.frame
dat=cbind(dat,group_list) #cbind横向追加，即将分组信息追加到最后一列
library("FactoMineR")#画主成分分析图需要加载这两个包
library("factoextra") 
# The variable group_list (index = 54676) is removed
# before PCA analysis
dat.pca <- PCA(dat[,-ncol(dat)], graph = FALSE)#现在dat最后一列是group_list，需要重新赋值给一个dat.pca,这个矩阵是不含有分组信息的

fviz_pca_ind(dat.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = dat$group_list, # color by groups
             # palette = c("#00AFBB", "#E7B800"),
             addEllipses = F, # Concentration ellipses
             legend.title = "Groups"
)
ggsave('all_samples_PCA.png')
```

### limma 分析
```
bp=function(g){         #定义一个函数g，函数为{}里的内容
  library(ggpubr)
  df=data.frame(gene=g,stage=group_list)
  p <- ggboxplot(df, x = "stage", y = "gene",
                 color = "stage", palette = "jco",
                 add = "jitter")
  #  Add p-value
  p + stat_compare_means()
}
#limma
group_list=clinic$sample_type
design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
head(design)
dat=expset
exprSet=dat
rownames(design)=colnames(exprSet)
colnames(design)=c('normal','tumor')
contrast.matrix<-makeContrasts("tumor-normal",
                               levels = design)
contrast.matrix ##这个矩阵声明，我们要把 Tumor 组跟 Normal 进行差异分析比较
​
deg = function(exprSet,design,contrast.matrix){
  ##step1
  fit <- lmFit(exprSet,design)
  ##step2
  fit2 <- contrasts.fit(fit, contrast.matrix) 
  ##这一步很重要，大家可以自行看看效果
  
  fit2 <- eBayes(fit2)  ## default no trend !!!
  ##eBayes() with trend=TRUE
  ##step3
  tempOutput = topTable(fit2, coef=1, n=Inf)
  nrDEG = na.omit(tempOutput) 
  #write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
  head(nrDEG)
  return(nrDEG)
}
​
deg = deg(exprSet,design,contrast.matrix)
​
head(deg)
## for volcano 
if(T){
  nrDEG=deg
  head(nrDEG)
  attach(nrDEG)
  plot(logFC,-log10(P.Value))
  library(ggpubr)
  df=nrDEG
  df$v= -log10(P.Value) #df新增加一列'v',值为-log10(P.Value)
  ggscatter(df, x = "logFC", y = "v",size=0.5)
  
  df$g=ifelse(df$P.Value>0.01,'stable', #if 判断：如果这一基因的P.Value>0.01，则为stable基因
              ifelse( df$logFC >2,'up', #接上句else 否则：接下来开始判断那些P.Value<0.01的基因，再if 判断：如果logFC >1.5,则为up（上调）基因
                      ifelse( df$logFC < -2,'down','stable') )#接上句else 否则：接下来开始判断那些logFC <1.5 的基因，再if 判断：如果logFC <1.5，则为down（下调）基因，否则为stable基因
  )
  table(df$g)
  df$name=rownames(df)
  head(df)
  ggscatter(df, x = "logFC", y = "v",size=0.5,color = 'g')
  ggscatter(df, x = "logFC", y = "v", color = "g",size = 0.5,
            label = "name", repel = T,
            #label.select = rownames(df)[df$g != 'stable'] ,
            label.select = c('TTC9', 'AQP3', 'CXCL11','PTGS2'), #挑选一些基因在图中显示出来
            palette = c("#00AFBB", "#E7B800", "#FC4E07") )
  ggsave('volcano.png')
  
  ggscatter(df, x = "AveExpr", y = "logFC",size = 0.2)
  df$p_c = ifelse(df$P.Value<0.001,'p<0.001',
                  ifelse(df$P.Value<0.01,'0.001<p<0.01','p>0.01'))
  table(df$p_c )
  ggscatter(df,x = "AveExpr", y = "logFC", color = "p_c",size=0.2, 
            palette = c("green", "red", "black") )
  ggsave('MA.png')
  
  
}
​
## for heatmap 
if(T){ 
  load(file = 'step1-output.Rdata')
  # 每次都要检测数据
  dat[1:4,1:4]
  table(group_list)
  x=deg$logFC #deg取logFC这列并将其重新赋值给x
  names(x)=rownames(deg) #deg取probe_id这列，并将其作为名字给x
  cg=c(names(head(sort(x),100)),#对x进行从小到大排列，取前100及后100，并取其对应的探针名，作为向量赋值给cg
       names(tail(sort(x),100)))
  library(pheatmap)
  pheatmap(dat[cg,],show_colnames =F,show_rownames = F) #对dat按照cg取行，所得到的矩阵来画热图
  n=t(scale(t(dat[cg,])))
  #通过“scale”对log-ratio数值进行归一化，现在的dat是行名为探针，列名为样本名，由于scale这个函数应用在不同组数据间存在差异时，需要行名为样本，因此需要用t(dat[cg,])来转换，最后再转换回来
  
  n[n>2]=2
  n[n< -2]= -2
  n[1:4,1:4]
  pheatmap(n,show_colnames =F,show_rownames = F)
  ac=data.frame(g=group_list)
  rownames(ac)=colnames(n) #将ac的行名也就分组信息（是‘no TNBC’还是‘TNBC’）给到n的列名，即热图中位于上方的分组信息
  pheatmap(n,show_colnames =F,
           show_rownames = F,
          cluster_cols = F, 
           annotation_col=ac,filename = 'heatmap_top200_DEG.png') #列名注释信息为ac即分组信息
  
  
}

write.csv(deg,file = 'deg.csv')
```

### 生存分析
```
survival_table=as.data.frame(read.csv('COADREAD_survival.csv',header = T))
rownames(survival_table)=survival_table$sample
#信息矩阵制作
meta = na.omit(survival_table[tumor,1:4])
meta = meta[,-2]
colnames(meta) = c('ID','event',"time")

meta$time_m = meta$time/30
dat = expset[,as.character(meta$ID)]
a = dat['APOC1',]
c = apply(a,1,mean)
b = ifelse(t(a)[,1]>c,1,0)
meta$up_APOC1 = b
sfit <- survfit(Surv(time, event)~up_APOC1, data=meta)
sfit
ggsurvplot(sfit, conf.int=F, pval=TRUE)
## more complicate figures.
ggsurvplot(sfit,palette = c("#E7B800", "#2E9FDF"),
           risk.table =TRUE,pval =TRUE,
           conf.int =TRUE,xlab ="Time in months", 
           ggtheme =theme_light(), 
           ncensor.plot = TRUE)
```
