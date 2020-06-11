#全局变量设置#
if(T){
	Sys.setenv(R_MAX_NUM_DLLS=999) ##Sys.setenv修改环境设置，R的namespace是有上限的，如果导入包时超过这个上次就会报错,R_MAX_NUM_DLLS可以修改这个上限
	options(stringsAsFactors = F) ##options:允许用户对工作空间进行全局设置，stringsAsFactors防止R自动把字符串string的列辨认成factor
	options()$repos  ## 查看使用install.packages安装时的默认镜像
	options()$BioC_mirror ##查看使用bioconductor的默认镜像
	options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") ##指定镜像，这个是中国科技大学镜像
	options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")) ##指定install.packages安装镜像，这个是清华镜像
	options()$repos 
	options()$BioC_mirror
	if (!requireNamespace("BiocManager", quietly = TRUE)) 
	  install.packages("BiocManager") ##判断是否存在BiocManager包，不存在的话安装
}

########################################################################
#安装包#
#BiocManager::install("pathview",'DESeq2','reshap','qrqc','org.Mm.eg.db','DOSE','clusterProfiler','topGO','pathview')
#BiocManager::install("topGO")
########################################################################
#编辑环境#
rm(list=ls())
if(T){
	library(reshape2)
	library(edgeR)
	library(DESeq2)
	library(reshape)
	library(ggplot2)
	library(qrqc)
	library(stringr)
	library(optparse)
	library(org.Mm.eg.db)
	library(DOSE)
	library(clusterProfiler)
	library(pathview)
	library(topGO)
}

########################################################################
#分析前准备
	setwd("/Volumes/Lexar/project/R")
	getwd()
	dir()
	#数据预处理#
	a=as.data.frame(read.table("all.id.txt",stringsAsFactors = F, header = T, row.names = 1))
	rownames(a)=str_split(rownames(a),'\\.',simplify = T) [,1]
	dat=a[apply(a[,6:dim(a)[2]],1,sum) > 0,6:dim(a)[2]] 
	plate=str_split(colnames(dat),'\\.',simplify = T) [,6:7]#取列名，以'_'号分割，提取第三列。
	colnames(a)[6:dim(a)[2]]=paste(plate[,1],plate[,2],sep = '_')
	colnames(dat)=paste(plate[,1],plate[,2],sep = '_')
	#原始表达谱的清洗和分离
	#将列名中的一些奇怪字符串去掉
	#删去没有counts的基因，根据具体需求设置阈值，这里设置只要有一个reads匹配到就行>0

	dat=log2(edgeR::cpm(dat)+1) 
	hc=hclust(dist(t(dat))) 
	png("hc.png")
	plot(hc)
	dev.off()
	##归一化的一种选择，这里是CPM(count-per-million，每百万碱基中每个转录本的count值)
	###CPM只对read count相对总reads数做了数量的均一化，去除文库大小差异。
	# 可以看到dist函数计算样本直接距离和cor函数计算样本直接相关性，是完全不同的概念。虽然我都没有调它们两个函数的默认的参数。
	# 总结：
	# - dist函数计算行与行（样本）之间的距离
	# - cor函数计算列与列（样本）之间的相关性
	# - scale函数默认对每一列（样本）内部归一化
	#层次聚类
	##样本间层次聚类
	# 原始表达矩阵转置后，细胞在行，所以计算的是细胞与细胞之间的距离。
	## statquest 有详细讲解背后的统计学原理。
	#t:矩阵转置，行转列，列转行
	#分类时常常需要估算不同样本之间的相似性(Similarity Measurement)
	# 这时通常采用的方法就是计算样本间”距离”(Distance)。
	#dist函数是R语言计算距离的主要函数。dist函数可以计算行与行两两间的距离。
	# 所以之前的矩阵里面行是基因，转置后行是样本，因为我们要计算样本与样本之间的距离。
	# dist()函数计算变量间距离
	#hclust函数用来层次聚类

	clus = cutree(hc, 3) #对hclust()函数的聚类结果进行剪枝，即选择输出指定类别数的系谱聚类结果。
	group_list= as.factor(clus) ##转换为因子属性
	table(group_list) ##统计频数

	#提取批次信息
	plate=str_split(colnames(dat),'_',simplify = T)[,1] #取列名，以'_'号分割，提取第三列。
	table(plate)
	n_g = apply(a[6:dim(a)[2]],2,function(x) sum(x>1))
	df=data.frame(g=group_list,plate=plate,n_g=n_g) 
	#统计每个样本有表达的有多少行（基因）
	# 这里我们定义， reads数量大于1的那些基因为有表达，一般来说单细胞转录组过半数的基因是不会表达的。
	#新建数据框(细胞的属性信息)
	##(样本为行名，列分别为：样本分类信息，样本分组，样本表达的基因数【注意：不是表达量的和，而是种类数或者说个数】)

	df$all='all' #添加列，列名为"all"，没事意思，就是后面有需要
	metadata=df
	save(a,dat,df,file = 'input.Rdata') #保存a,dat,df这变量到上级目录的input.Rdata
	# 因为另外一个项目也需要使用这个数据集，所以保存到了上级目录。
	#PCA#
	rm(list = ls())  ## 魔幻操作，一键清空~
	options(stringsAsFactors = F)
	load(file = 'input.Rdata')
	plate=df$plate
	## 下面是画PCA的必须操作，需要看不同做PCA的包的说明书。
	dat_back=dat
	dat=dat_back
	dat=t(dat)
	dat=as.data.frame(dat)
	dat=cbind(dat,plate ) #cbind根据列进行合并，即叠加所有列 #矩阵添加批次信息
	table(dat$plate)
	library("FactoMineR")
	library("factoextra") 
	# The variable plate (index = ) is removed
	# before PCA analysis
	dat.pca <- PCA(dat[,-ncol(dat)], graph = FALSE)
	fviz_pca_ind(dat.pca,#repel =T,
	             geom.ind = "point", # show points only (nbut not "text")
	             col.ind = dat$plate, # color by groups
	             #palette = c("#00AFBB", "#E7B800"),
	             addEllipses = TRUE, # Concentration ellipses
	             legend.title = "Groups"
	) 
	ggsave('all_cells_PCA_by_plate.png') 




#差异分析#

##控制组标签==1，疾病组==2
setwd("/Volumes/Lexar/project/R")
rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)
load(file = 'input.Rdata')
##（多变量）
dir.create("7小时处理差异")
setwd('./7小时处理差异')

run_go = function(res,contrast,ordb){
  tmp=select(ordb, keys=res$ENSEMBL, columns="GO", keytype="ENSEMBL")
  ensembl_go=unlist(tapply(tmp[,2],as.factor(tmp[,1]),function(x) paste(x,collapse ="|"),simplify =F))
  res$go=ensembl_go[res$ENSEMBL]
  diff_gene_deseq2 <-subset(res,FDR < 0.05)
  write.csv(diff_gene_deseq2,paste0("diff_gene_",contrast,".csv"),quote = F,row.names = F)
  write.csv(res,paste0(contrast,"_clean",".csv"),quote = F,row.names = F)
}
run_edgr = function(exprSet,group_list){
  suppressMessages(library(edgeR)) 
  dgelist <- DGEList(counts = exprSet, group = group_list)
  keep <- rowSums(cpm(dgelist) > 1 ) >= 2
  dgelist <- dgelist[keep,keep.lib.sizes = FALSE]
  dgelist_norm <- calcNormFactors(dgelist, method = 'TMM')
  cls = ifelse(group_list=='C1','red','blue')
  png("样品间的相似性.png")
  plotMDS(dgelist_norm, col = cls,dim = c(1, 2))
  dev.off()
  
  design <- model.matrix(~group_list)    #构建分组矩阵
  #design 查看分组信息（attr(,"assign")）分组0相对于分组1的比较
  dge <- estimateDisp(dgelist_norm, design, robust = TRUE) #估算离散值
  png("离散度.png")
  plotBCV(dge)
  dev.off()
  
  fit <- glmFit(dge, design, robust = TRUE)     #拟合模型
  lrt <- glmLRT(fit)   #统计检验
  topTags(lrt)
  write.csv(topTags(lrt, n = nrow(dgelist$counts)), 'glmLRT.csv', quote = FALSE)        #输出主要结果
  dge_de <- decideTestsDGE(lrt, adjust.method = 'fdr', p.value = 0.05)  #查看默认方法获得的差异基因
  summary(dge_de)
  png("glmLRT.png")
  plotMD(lrt, status = dge_de, values = c(1, -1), col =c('blue', 'red'))     #作图观测
  abline(h = c(-1, 1), col = 'gray', lty = 2)
  dev.off()
  
  fit <- glmQLFit(dge, design, robust = TRUE)        #拟合模型
  lqt <- glmQLFTest(fit)    #统计检验
  topTags(lqt)
  write.csv(topTags(lqt, n = nrow(dgelist$counts)), 'glmQLFTest.csv', quote = FALSE)        #输出主要结果
  dge_de <- decideTestsDGE(lqt, adjust.method = 'fdr', p.value = 0.05)  #查看默认方法获得的差异基因
  summary(dge_de)
  png("glmLQT.png")
  plotMD(lqt, status = dge_de, values = c(1, -1), col = c('blue', 'red'))     #作图观测
  abline(h = c(-1, 1), col = 'gray', lty = 2)
  dev.off()
  
  dge_et <- exactTest(dge) #检验
  topTags(dge_et)
  write.csv(topTags(dge_et, n = nrow(dgelist$counts)), 'exactTest.csv', quote = FALSE)        #输出主要结果
  dge_de <- decideTestsDGE(dge_et, adjust.method = 'fdr', p.value = 0.05)   #查看默认方法获得的差异基因
  summary(dge_de)
  detags <- rownames(dge)[as.logical(dge_de)]
  png("dge_de.png")
  plotSmear(dge_et, de.tags = detags, cex = 0.5)      #作图观测
  abline(h = c(-1, 1), col = 'gray', lty = 2)
  dev.off()
}
reu_DESeq= function(exprSet,colData,condition){
  suppressMessages(library(DESeq2)) 
  dds <- DESeqDataSetFromMatrix(countData = exprSet, colData = colData,
                                design= ~ condition)
  dds <- dds[ rowSums(counts(dds)) > 1, ]
  dds <- DESeq(dds)
  #~在R里面用于构建公式对象，~左边为因变量，右边为自变量。
  res <- results(dds)
  mcols(res, use.names = TRUE)
  summary(res)
  res$FDR=res$padj
  # 获取padj（p值经过多重校验校正后的值）小于0.05，表达倍数取以2为对数后大于1或者小于-1的差异表达基因。
  table(res$FDR<0.05)     #计算p值小于0.05的基因的个数
  res_deseq <- res[order(res$FDR),]   #根据res的padj值进行排序，并赋值给res_deseq
  diff_gene_deseq2 <- subset(res_deseq, FDR<0.05)    #将res_deseq过滤：要求是p值小于0.05且log2FoldChange在1~-1范围内的
  diff_gene_deseq2 <- row.names(diff_gene_deseq2)
  res_diff_data <- merge(as.data.frame(res),as.data.frame(counts(dds,normalize=TRUE)),by="row.names",sort=FALSE)
  write.csv(res_diff_data,file = "DESeq_diff.csv",row.names = F)
}
t1=c(9,10,11) ##删除非目标样本列
t2=c(4,5,6)  #删除非目标样本行
ordb = org.Mm.eg.db
ors ='mmu'
#差异分析
if(T){
  a =a[,-t1]
  dat = dat[,-t2]
  df = df[-t2,]
  exprSet=a[,6:dim(a)[2]]
  condition = factor(ifelse(df$g==1,'contral','KD'))
  colData = as.data.frame(condition,row.names = rownames(df))
  group_list=ifelse(df$g==1,'C1','C2')
  contrast=paste('C1','C2',sep="-")
  print("run_edgr")
  run_edgr(exprSet,group_list)
  print("run_DESeq")
  reu_DESeq(exprSet,colData,condition) 


  ext=as.data.frame(read.table("exactTest.csv",stringsAsFactors = F, header = T, row.names = 1,sep = ","))
  head(ext)
  deg=as.data.frame(read.table("DESeq_diff.csv",stringsAsFactors = F, header = T, row.names = 1,sep = ","))
  head(deg)
  lrt=as.data.frame(read.table("glmLRT.csv",stringsAsFactors = F, header = T, row.names = 1,sep = ","))
  lqt=as.data.frame(read.table("glmQLFTest.csv" ,stringsAsFactors = F, header = T, row.names = 1,sep = ","))
  ext$Row.names=rownames(ext)
  deg$Row.names=rownames(deg)
  lrt$Row.names=rownames(lrt)
  lqt$Row.names=rownames(lqt)
  ##还原基因名称 
  print("还原基因名称")
  k=keys(ordb,keytype = "ENSEMBL")
  list=select(ordb,keys=k,columns = c("ENTREZID","SYMBOL"), keytype="ENSEMBL")
  ID_list=na.omit(list[match(rownames(a),list[,"ENSEMBL"]),])
  ##去空
  res_ext=na.omit(merge(ID_list,ext,by.x="ENSEMBL",by.y= "Row.names",all=TRUE))
  res_ext$up_down=ifelse(res_ext$logFC>0,1,-1)
  res_lrt=na.omit(merge(ID_list,lrt,by.x="ENSEMBL",by.y= "Row.names",all=TRUE))
  res_lrt$up_down=ifelse(res_lrt$logFC>0,1,-1)
  res_lqt=na.omit(merge(ID_list,lqt,by.x="ENSEMBL",by.y= "Row.names",all=TRUE))
  res_lqt$up_down=ifelse(res_lqt$logFC>0,1,-1)
  res_deg=na.omit(merge(ID_list,deg,by.x="ENSEMBL",by.y= "Row.names",all=TRUE))
  res_deg$up_down=ifelse(res_deg$log2FoldChange>0,1,-1)
  
  ##GO注释
  print("注释")
  res = res_deg
  contrast = "res_deg"
  run_go(res,contrast,ordb)
  
  res = res_lrt
  contrast = "res_lrt"
  run_go(res,contrast,ordb)
  
  res = res_lqt
  contrast = "res_lqt"
  run_go(res,contrast,ordb)
  
  res = res_ext
  contrast = "res_ext"
  run_go(res,contrast,ordb)  
}
#plot
  print("富集")
  ext=as.data.frame(read.table("diff_gene_res_ext.csv",stringsAsFactors = F, header = T, row.names = 1,sep = ","))
  deg=as.data.frame(read.table("diff_gene_res_deg.csv",stringsAsFactors = F, header = T, row.names = 1,sep = ","))
  lrt=as.data.frame(read.table("diff_gene_res_lrt.csv",stringsAsFactors = F, header = T, row.names = 1,sep = ","))
  lqt=as.data.frame(read.table("diff_gene_res_lqt.csv" ,stringsAsFactors = F, header = T, row.names = 1,sep = ","))
  print(1)
  difgen = ext
  contrast = 'ext'
  diff_up = difgen[difgen$up_down == 1,]$ENTREZID
  diff_down= difgen[difgen$up_down == -1,]$ENTREZID
  ego_up <- enrichGO(gene=diff_up,ordb,ont="CC",pvalueCutoff=0.05)
  ego_down <- enrichGO(gene=diff_down,ordb,ont="CC",pvalueCutoff=0.05)
  ekk_up<- enrichKEGG(gene=diff_up,organism=ors,pvalueCutoff=0.05)
  ekk_down <- enrichKEGG(gene=diff_down,organism=ors,pvalueCutoff=0.05)
  write.csv(ekk_up@result,paste0("ekk_up_",contrast,".csv"),quote = F,row.names = F)
  write.csv(ekk_down@result,paste0("ekk_down_",contrast,".csv"),quote = F,row.names = F)
  write.csv(ego_up@result,paste0("ego_up_",contrast,".csv"),quote = F,row.names = F)
  write.csv(ego_down@result,paste0("ego_down_",contrast,".csv"),quote = F,row.names = F)  
  png(paste0("ego_up_",contrast,'.png'))
  barplot(ego_up,title="EnrichmentGO_up")
  dev.off()
  png(paste0("ego_down_",contrast,'.png'))
  barplot(ego_down,title="EnrichmentGO_down")
  dev.off()
  png(paste0("ekk_up_",contrast,'.png'))
  barplot(ekk_up,title="EnrichmentKEGG_up")
  dev.off()
  png(paste0("ekk_down_",contrast,'.png'))
  barplot(ekk_down,title="EnrichmentKEGG_down")
  dev.off()
  print(2)
  difgen = deg
  contrast = 'deg'
  diff_up = difgen[difgen$up_down == 1,]$ENTREZID
  diff_down= difgen[difgen$up_down == -1,]$ENTREZID
  ego_up <- enrichGO(gene=diff_up,ordb,ont="CC",pvalueCutoff=0.05)
  ego_down <- enrichGO(gene=diff_down,ordb,ont="CC",pvalueCutoff=0.05)
  ekk_up<- enrichKEGG(gene=diff_up,organism=ors,pvalueCutoff=0.05)
  ekk_down <- enrichKEGG(gene=diff_down,organism=ors,pvalueCutoff=0.05)
  write.csv(ekk_up@result,paste0("ekk_up_",contrast,".csv"),quote = F,row.names = F)
  write.csv(ekk_down@result,paste0("ekk_down_",contrast,".csv"),quote = F,row.names = F)
  write.csv(ego_up@result,paste0("ego_up_",contrast,".csv"),quote = F,row.names = F)
  write.csv(ego_down@result,paste0("ego_down_",contrast,".csv"),quote = F,row.names = F)  
  png(paste0("ego_up_",contrast,'.png'))
  barplot(ego_up,title="EnrichmentGO_up")
  dev.off()
  png(paste0("ego_down_",contrast,'.png'))
  barplot(ego_down,title="EnrichmentGO_down")
  dev.off()
  png(paste0("ekk_up_",contrast,'.png'))
  barplot(ekk_up,title="EnrichmentKEGG_up")
  dev.off()
  png(paste0("ekk_down_",contrast,'.png'))
  barplot(ekk_down,title="EnrichmentKEGG_down")
  dev.off()
  print(3)
  difgen = lrt
  contrast = 'lrt'
  diff_up = difgen[difgen$up_down == 1,]$ENTREZID
  diff_down= difgen[difgen$up_down == -1,]$ENTREZID
  ego_up <- enrichGO(gene=diff_up,ordb,ont="CC",pvalueCutoff=0.05)
  ego_down <- enrichGO(gene=diff_down,ordb,ont="CC",pvalueCutoff=0.05)
  ekk_up<- enrichKEGG(gene=diff_up,organism=ors,pvalueCutoff=0.05)
  ekk_down <- enrichKEGG(gene=diff_down,organism=ors,pvalueCutoff=0.05)
  write.csv(ekk_up@result,paste0("ekk_up_",contrast,".csv"),quote = F,row.names = F)
  write.csv(ekk_down@result,paste0("ekk_down_",contrast,".csv"),quote = F,row.names = F)
  write.csv(ego_up@result,paste0("ego_up_",contrast,".csv"),quote = F,row.names = F)
  write.csv(ego_down@result,paste0("ego_down_",contrast,".csv"),quote = F,row.names = F)  
  png(paste0("ego_up_",contrast,'.png'))
  barplot(ego_up,title="EnrichmentGO_up")
  dev.off()
  png(paste0("ego_down_",contrast,'.png'))
  barplot(ego_down,title="EnrichmentGO_down")
  dev.off()
  png(paste0("ekk_up_",contrast,'.png'))
  barplot(ekk_up,title="EnrichmentKEGG_up")
  dev.off()
  png(paste0("ekk_down_",contrast,'.png'))
  barplot(ekk_down,title="EnrichmentKEGG_down")
  dev.off()
  print(4)
  difgen = lqt
  contrast = 'lqt'
  diff_up = difgen[difgen$up_down == 1,]$ENTREZID
  diff_down= difgen[difgen$up_down == -1,]$ENTREZID
  ego_up <- enrichGO(gene=diff_up,ordb,ont="CC",pvalueCutoff=0.05)
  ego_down <- enrichGO(gene=diff_down,ordb,ont="CC",pvalueCutoff=0.05)
  ekk_up<- enrichKEGG(gene=diff_up,organism=ors,pvalueCutoff=0.05)
  ekk_down <- enrichKEGG(gene=diff_down,organism=ors,pvalueCutoff=0.05)
  write.csv(ekk_up@result,paste0("ekk_up_",contrast,".csv"),quote = F,row.names = F)
  write.csv(ekk_down@result,paste0("ekk_down_",contrast,".csv"),quote = F,row.names = F)
  write.csv(ego_up@result,paste0("ego_up_",contrast,".csv"),quote = F,row.names = F)
  write.csv(ego_down@result,paste0("ego_down_",contrast,".csv"),quote = F,row.names = F)  
  png(paste0("ego_up_",contrast,'.png'))
  barplot(ego_up,title="EnrichmentGO_up")
  dev.off()
  png(paste0("ego_down_",contrast,'.png'))
  barplot(ego_down,title="EnrichmentGO_down")
  dev.off()
  png(paste0("ekk_up_",contrast,'.png'))
  barplot(ekk_up,title="EnrichmentKEGG_up")
  dev.off()
  png(paste0("ekk_down_",contrast,'.png'))
  barplot(ekk_down,title="EnrichmentKEGG_down")
  dev.off()





#可视化

plotGOgraph(ego_up)
plotGOgraph(ego_down)
path<- pathview(gene.data = res$ENTREZID,pathway.id = "mmu04110", species = ors,imit = list(gene=max(abs(geneList)), cpd=1))







#热图
rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)
load(file = 'DEG_diff_gene.Rdata')
load(file = 'input.Rdata')
library(pheatmap)
genes = diff_gene_deseq2$ENSEMBL
ac=data.frame(group=df$g)
rownames(ac)=colnames(dat)
mat = dat[genes,]
n=t(scale(t(mat)))
pheatmap(n,show_rownames = T,show_colnames = F, 
         cluster_rows = F,cluster_cols = F,
         annotation_col = ac)



# 斑马鱼：Bioconductor - org.Dr.eg.db - /packages/release/data/annotation/html/org.Dr.eg.db.html

# Details biocViews AnnotationData , Danio_rerio , OrgDb Version 3

# 拟南芥：Bioconductor - org.At.tair.db - /packages/release/data/annotation/html/org.At.tair.db.html

# Details biocViews AnnotationData , Arabidopsis_thaliana , OrgDb Version 3

# 小鼠：Bioconductor - org.Mm.eg.db - /packages/release/data/annotation/html/org.Mm.eg.db.html

# Details biocViews AnnotationData , Mus_musculus , OrgDb , mouseLLMappings Version 3

# 人类：Bioconductor - org.Hs.eg.db - /packages/release/data/annotation/html/org.Hs.eg.db.html



#可视化
barplot(ego_up, showCategory=20,title="EnrichmentGO_up")
barplot(ego_down, showCategory=20,title="EnrichmentGO_down")
barplot(ekk_up, showCategory=20,title="EnrichmentKEGG_up")
barplot(ekk_down, showCategory=20,title="EnrichmentKEGG_down")
plotGOgraph(ego_up)
plotGOgraph(ego_down)
path<- pathview(gene.data = res$ENTREZID,pathway.id = "mmu04110", species = ors,imit = list(gene=max(abs(geneList)), cpd=1))








