### ---
### Title: "Untitled"
### Create: "Yuansh"
### Date: "5/01/2020"
### Email: yuansh3354@163.com
### output: html_document
### ---

### step1 差异基因

##### 1.提取grouplist 并构建分组矩阵
if(T){
  group_list = as.character(pd[,44])
  design <- model.matrix(~0+factor(group_list))
  colnames(design)=levels(factor(group_list))
  rownames(design)=colnames(df_expr)
  print(paste(colnames(design)[1],colnames(design)[2],sep = '-'))
  contrast.matrix<-makeContrasts(paste(colnames(design)[1],colnames(design)[2],sep = '-'),
                                 levels = design)
}

##### 2. 寻找差异基因
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
deg = deg(df_expr,design,contrast.matrix)
deg = deg[order(deg$adj.P.Val),]
# 输出文件夹中的GSE13911-deg.csv就是差异基因
top = deg
top = top[which(top$adj.P.Val<0.05 & abs(top$logFC)>1),]
write.csv(top,paste(file,'-deg.csv',sep = ''))

##### 3. 提取满足条件的差异基因
if(T){
  top = deg
  top = top[which(top$adj.P.Val<0.05 & abs(top$logFC)>1),]
  top = top[order(top$adj.P.Val),]
  top = top[1:100,]
  int_gene = rownames(top[which(top$adj.P.Val<0.05 & abs(top$logFC)>1),])
  id = rownames(int_gene)
  write.csv(top,paste(file,'-top-deg.csv',sep = ''))
}


##### 4.聚类
if(T){
  library(pheatmap)
  n = t(scale(t(df_expr[int_gene,])))
  # n[n>2] = 2
  # n[n< -2] = -2
  ac = data.frame(g=group_list)
  ac$names = colnames(n) 
  ac = ac[order(ac[,1]),]
  rownames(ac) = ac$names
  a = as.data.frame(ac[,1])
  colnames(a) = 'Type'
  rownames(a) = rownames(ac)
  pheatmap(n,show_colnames =F,
           show_rownames = F,
           annotation_col=a,filename = '符合条件的基因聚类.png')
}
