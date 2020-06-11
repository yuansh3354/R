### ---
### Title: "Untitled"
### Create: "Yuansh"
### Date: "5/01/2020"
### Email: yuansh3354@163.com
### output: html_document
### ---

rm(list = ls())

### step4 识别风险样本
# 训练集风险识别GSE31210
#setwd("/Volumes/Lexar/ZG/GSE19188分析结果/GSE31210")
result = read.csv('GSE19188-gene-pair.csv')
if(T){
  df = read.csv('GSE31210.csv',row.names = 1)
  pd = read.csv('GSE31210_clinic.csv',row.names = 1)
  # 过滤正常样本
  sit = -which(pd[,10]=='tissue: normal lung')
  pd = pd[sit,]
  df = df[,sit]
  rownames(pd) == colnames(df)
} # GSE31210

if(T){
  df = read.csv('GSE50081.csv',row.names = 1)
  pd = read.csv('GSE50081_clinic.csv',row.names = 1)
  rownames(pd) == colnames(df)
} # GSE50081


if(T){
  pair = result[,1:2]
  pair[,1] = as.character(pair[,1])
  pair[,2] = as.character(pair[,2])
  
  n = dim(pair)[1]
  label = NULL
  for(i in 1:n){
    t = df[pair[i,1],] - df[pair[i,2],]
    label = rbind(label,t)
  }
  label = ifelse(label>0,1,0)
  label = apply(label, 2, sum)
  # 因为是正常和癌症找差异,所以,在全为癌症的样本中区分高低风险,按70%的投票规则
  label1 = ifelse(label>9,'low' ,'hig')
  table(label1)
} # 识别风险

if(T){
  group_list = as.character(label1)
  design <- model.matrix(~0+factor(group_list))
  colnames(design)=levels(factor(group_list))
  rownames(design)=colnames(df)
  print(paste(colnames(design)[1],colnames(design)[2],sep = '-'))
  contrast.matrix<-makeContrasts(paste(colnames(design)[1],colnames(design)[2],sep = '-'),
                                 levels = design)
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
  deg = deg(df,design,contrast.matrix)
  deg = deg[order(deg$adj.P.Val),]
  top = deg
  top = top[which(top$adj.P.Val<0.05 & abs(top$logFC)>1),]
  write.csv(top,paste('tarin','-deg.csv',sep = ''))
  } # 定义label
if(T){
  top = deg
  top = top[which(top$adj.P.Val<0.05 & abs(top$logFC)>1),]
  top = top[order(top$adj.P.Val),]
  int_gene = rownames(top[which(top$adj.P.Val<0.05 & abs(top$logFC)>1),])
  id = rownames(int_gene)
  write.csv(top,paste('tarin','-deg.csv',sep = ''))
} # 差异基因
# 取p值前100的进行聚类
if(T){
  library(pheatmap)
  int = int_gene[1:100]
  n = t(scale(t(df[int,])))
  n[n>2] = 2
  n[n< -2] = -2
  ac = data.frame(g=group_list)
  ac$names = colnames(n) 
  ac = ac[order(ac[,1]),]
  rownames(ac) = ac$names
  a = as.data.frame(ac[,1])
  colnames(a) = 'Type'
  rownames(a) = rownames(ac)
  pheatmap(n,show_colnames =F,
           show_rownames = F,
           annotation_col=a,filename = 'train基因聚类.png')
}
### 5.通路富集
##### 其中go是功能富集,kegg是通路富集
if(T){
  logFC_t=1
  top$g= 0
  top[which(top$logFC<= -1),]$g = 'DOWN'
  top[which(top$logFC>= 1),]$g = 'UP'
  
  top$ENTREZID=rownames(top)
  
  gene_up= top[top$g == 'UP','ENTREZID'] 
  gene_down=top[top$g == 'DOWN','ENTREZID'] 
  gene_diff=c(gene_up,gene_down)
  source('kegg_and_go_up_and_down.R')
  run_go(gene_up,gene_down,pro='NORMAL-TUMOR')
  run_kegg(gene_up,gene_down,pro='NORMAL-TUMOR')
}

