#全局变量设置#
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

	########################################################################
#安装包#
	BiocManager::install("multicore")
	BiocManager::install("pathview")
	BiocManager::install("topGO")
	########################################################################
#编辑环境#
	rm(list=ls())
	library(reshape2)
	library(edgeR)
	library(DESeq2)
	library(reshape)
	library(ggplot2)
	library(multicore)
	library(qrqc)
	library(stringr)
	library(stringr)
	library(optparse)
	library(org.Mm.eg.db)
	library(DOSE)
	library(clusterProfiler)
	library(pathview)
	library(topGO)
	########################################################################
#工作目录设置#
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
	plot(hc)
	hc=hclust(dist(t(dat))) 
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
	#全局设置
		rm(list = ls())  ## 魔幻操作，一键清空~
		options(stringsAsFactors = F)
		load(file = 'input.Rdata')
		setwd("/Volumes/Lexar/project/R")
		## 载入第0步准备好的表达矩阵，及细胞的一些属性（hclust分群，plate批次，检测到的细胞数量）
		# 注意 变量a是原始的counts矩阵，变量 dat是logCPM后的表达量矩阵。
		getwd()
		dir()
		load('input.Rdata')
		#多重比较全局变量设置（少用）
			# exprSet=a[,6:dim(a)[2]]
			# condition <- as.factor(df[,2])
			# colData = as.data.frame(condition,row.names = rownames(df))
			# group_list=df$g
			# contrast=paste(paste('group',unique(group_list),sep="_"),collapse = "-")
			# group_list=factor(group_list,labels=paste('group',unique(group_list),sep="_"))
		#两类样本比较全局变量设置
			exprSet=a[,6:dim(a)[2]]
			condition = ifelse(df$g==1,'me','other')
			colData = as.data.frame(condition,row.names = rownames(df))
			group_list=df$g
			contrast=paste(paste('group',unique(group_list),sep="_"),collapse = "-")
			group_list=factor(group_list,labels=paste('group',unique(group_list),sep="_"))
	#limma

		suppressPackageStartupMessages(library(limma))
		design=model.matrix(~0+factor(group_list))
		colnames(design)=levels(factor(group_list))
		fit=lmFit(exprSet,design)
		cont.matrix=makeContrasts(contrasts=contrast ,levels = design)
		fit2=contrasts.fit(fit,cont.matrix)
		fit2=eBayes(fit2)
		results=topTable(fit2,adjust='BH',n=Inf)
		nr_results=na.omit(results)
		write.table(results,paste0("limma_DEG.",contrast,".txt"),quote = F,sep = '\t',row.names = T)
		write.table(nr_results,paste0("nr_limma_DEG.",contrast,".txt"),quote = F,sep = '\t',row.names = T)
	#DESeq
		#结果解读
			#——————————————————————————————————————————————————————————————————————
			#｜baseMean			｜mean of normalized counts for all samples。     
			#｜log2FoldChange 	｜log2 fold change (MLE): condition control vs KD 
			#｜lfcSE	            ｜standard error: condition control vs KD		  
			#｜stat	            ｜Wald statistic: condition control vs KD         
			#｜pvalue			    ｜Wald test p-value: condition control vs KD
			#｜padj	            ｜BH adjusted p-values
	        # ——————————————————————————————————————————————————————————————————————
		suppressMessages(library(DESeq2)) 
	                               
		dds <- DESeqDataSetFromMatrix(countData = exprSet, colData = colData,
	                                design= ~ condition)
	  	#~在R里面用于构建公式对象，~左边为因变量，右边为自变量。
	  	dds <- DESeq(dds)

  		png("qc_dispersions.png", 10000, 10000, pointsize=20)
  		plotDispEsts(dds, main="Dispersion plot")
 		dev.off()
 		##标准化基因count曲线
 		##黑色的点为估计出来的标准化后的基因count均值
 		##红色是拟合曲线
 		##蓝色的点代表标准后真实的基因count均值
	  	resultsNames(dds) # lists the coefficients
	  	res <- results(dds)
	  	summary(res)
	  	# 获取padj（p值经过多重校验校正后的值）小于0.05，表达倍数取以2为对数后大于1或者小于-1的差异表达基因。
	  	table(res$padj<0.05) #取P值小于0.05的结果
	  	res <- res[order(res$padj),]
		resdata <-  merge(as.data.frame(res),as.data.frame(counts(dds,normalize=TRUE)),by="row.names",sort=FALSE)
		# 得到csv格式的差异表达分析结果
		write.csv(resdata,file= "DESeq_res_raw.cvs",row.names = F)
		##还原基因名称 
		  k=keys(org.Mm.eg.db,keytype = "ENSEMBL")
		  list=select(org.Mm.eg.db,keys=k,columns = c("ENTREZID","SYMBOL"), keytype="ENSEMBL")
		  ID_list=na.omit(list[match(rownames(a),list[,"ENSEMBL"]),])
		##去空
		  res=na.omit(merge(ID_list,resdata,by.x="ENSEMBL",by.y= "Row.names",all=TRUE))
		  res$up_down=ifelse(res$log2FoldChange>0,1,-1)
		##GO注释
			tmp=select(org.Mm.eg.db, keys=res$ENSEMBL, columns="GO", keytype="ENSEMBL")
			ensembl_go=unlist(tapply(tmp[,2],as.factor(tmp[,1]),function(x) paste(x,collapse ="|"),simplify =F))
			res$go=ensembl_go[res$ENSEMBL]
		##输出
			diff_gene_deseq2 <-subset(res,padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
			diff_gene_deseq2=diff_gene_deseq2[order(diff_gene_deseq2$padj),]#
			write.csv(diff_gene_deseq2,file="diff_gene_deseq2.csv",row.names =F)
			write.csv(res,file= "DESeq_res_clean.cvs",row.names = F)
			save(diff_gene_deseq2,file = 'DEG_diff_gene.Rdata')
		##富集
			diff_up = diff_gene_deseq2[diff_gene_deseq2$up_down == 1,]$ENTREZID
			diff_down= diff_gene_deseq2[diff_gene_deseq2$up_down == -1,]$ENTREZID
			ego_up <- enrichGO(gene=diff_up,org.Mm.eg.db,ont="CC",pvalueCutoff=0.05,readable=TRUE)
			ego_down <- enrichGO(gene=diff_down,org.Mm.eg.db,ont="CC",pvalueCutoff=0.05,readable=TRUE)
			ekk_up<- enrichKEGG(gene=diff_up,organism="mmu",pvalueCutoff=0.05)
			ekk_down <- enrichKEGG(gene=diff_down,organism="mmu",pvalueCutoff=0.05)
			write.csv(summary(ekk_up),"KEGG_up-enrich.csv",row.names =F)
			write.csv(summary(ekk_down),"KEGG_down-enrich.csv",row.names =F)
			write.csv(summary(ego_up),"GO_up-enrich.csv",row.names =F)
			write.csv(summary(ego_down),"GO_down-enrich.csv",row.names =F)

			#可视化
			barplot(ego_up, showCategory=20,title="EnrichmentGO_up")
			barplot(ego_down, showCategory=20,title="EnrichmentGO_down")
			barplot(ekk_up, showCategory=20,title="EnrichmentKEGG_up")
			barplot(ekk_down, showCategory=20,title="EnrichmentKEGG_down")
			plotGOgraph(ego_up)
			plotGOgraph(ego_down)
			path<- pathview(gene.data = res$ENTREZID,pathway.id = "mmu04110", species = "mmu",imit = list(gene=max(abs(geneList)), cpd=1))
	#edegR
		suppressPackageStartupMessages(library(edgeR)) 
		dgelist = DGEList(exprSet,  group = condition)
		keep <- rowSums(cpm(dgelist) > 1 ) >= 2 #过滤
		dgelist <- dgelist[keep, ,keep.lib.sizes = FALSE]
		dgelist_norm <- calcNormFactors(dgelist, method = 'TMM')
		dev.new()
		plotMDS(dgelist_norm, col = rep(c('red', 'blue'), each = 5), dim = c(1, 2))
		design <- model.matrix( ~ condition)    #构建分组矩阵
		dge <- estimateDisp(dgelist_norm, design, robust = TRUE) #估算离散值
		plotBCV(dge) #作图查看
		dev.off()
		#负二项式广义对数线性模型
			fit <- glmFit(dge, design, robust = TRUE)     #拟合模型
			lrt <- glmLRT(fit)   #统计检验
			
			topTags(lrt)
			write.csv(topTags(lrt, n = nrow(dgelist$counts)), 'glmLRT.csv', quote = FALSE)        #输出主要结果
			 
			dge_de <- decideTestsDGE(lrt, adjust.method = 'fdr', p.value = 0.05)  #查看默认方法获得的差异基因
			summary(dge_de)
			#纵轴为log2 Fold Change值；横轴为log2 CPM值，反映了基因表达量信息；蓝色的点表示上调基因，红色的点表示下调基因，黑色的点为无差异基因。如果你想根据导出的数据，选择使用ggplot2自定义绘制火山图
			dev.new()
			plotMD(lrt, status = dge_de, values = c(1, -1), col = c('blue', 'red'))     #作图观测
			abline(h = c(-1, 1), col = 'gray', lty = 2)
			dev.off()
			res = as.data.frame(topTags(lrt, n = nrow(dgelist$counts)))
			res$up_down=ifelse(res$logFC>0,1,-1)
			diff_gene_deseq2 <-subset(res,FDR < 0.05)
			k=keys(org.Mm.eg.db,keytype = "ENSEMBL")
		 	list=select(org.Mm.eg.db,keys=k,columns = c("ENTREZID","SYMBOL"), keytype="ENSEMBL")
		  	ID_list=na.omit(list[match(rownames(diff_gene_deseq2),list[,"ENSEMBL"]),])
		  	diff_gene_deseq2$Row.names=rownames(diff_gene_deseq2)
		  	diff_gene_deseq2=na.omit(merge(ID_list,diff_gene_deseq2,by.x="ENSEMBL",by.y= "Row.names",all=TRUE))
		  	tmp=select(org.Mm.eg.db, keys=diff_gene_deseq2$ENSEMBL, columns="GO", keytype="ENSEMBL")
			ensembl_go=unlist(tapply(tmp[,2],as.factor(tmp[,1]),function(x) paste(x,collapse ="|"),simplify =F))
			diff_gene_deseq2$go=ensembl_go[diff_gene_deseq2$ENSEMBL]


			diff_up = diff_gene_deseq2[diff_gene_deseq2$up_down == 1,]$ENTREZID
			diff_down= diff_gene_deseq2[diff_gene_deseq2$up_down == -1,]$ENTREZID
			ego_up <- enrichGO(gene=diff_up,org.Mm.eg.db,ont="CC",pvalueCutoff=0.05,readable=TRUE)
			ego_down <- enrichGO(gene=diff_down,org.Mm.eg.db,ont="CC",pvalueCutoff=0.05,readable=TRUE)
			ekk_up<- enrichKEGG(gene=diff_up,organism="mmu",pvalueCutoff=0.05)
			ekk_down <- enrichKEGG(gene=diff_down,organism="mmu",pvalueCutoff=0.05)
			write.csv(summary(ekk_up),"irt_KEGG_up-enrich.csv",row.names =F)
			write.csv(summary(ekk_down),"irt_KEGG_down-enrich.csv",row.names =F)
			write.csv(summary(ego_up),"irt_GO_up-enrich.csv",row.names =F)
			write.csv(summary(ego_down),"irt_GO_down-enrich.csv",row.names =F)
		#类似然负二项式广义对数线性模型
			#quasi-likelihood negative binomial generalized log-linear model 拟合
			fit <- glmQLFit(dge, design, robust = TRUE)        #拟合模型
			lqt <- glmQLFTest(fit)    #统计检验
			 
			topTags(lqt)
			write.csv(topTags(lqt, n = nrow(dgelist$counts)), 'glmQLFTest.csv', quote = FALSE)        #输出主要结果
			 
			dge_de <- decideTestsDGE(lqt, adjust.method = 'fdr', p.value = 0.05)  #查看默认方法获得的差异基因
			summary(dge_de)
			dev.new()
			plotMD(lqt, status = dge_de, values = c(1, -1), col = c('blue', 'red'))     #作图观测
			abline(h = c(-1, 1), col = 'gray', lty = 2)
			dev.off()
			res = as.data.frame(topTags(lqt, n = nrow(dgelist$counts)))
			res$up_down=ifelse(res$logFC>0,1,-1)
			diff_gene_deseq2 <-subset(res,FDR < 0.05)
			k=keys(org.Mm.eg.db,keytype = "ENSEMBL")
		 	list=select(org.Mm.eg.db,keys=k,columns = c("ENTREZID","SYMBOL"), keytype="ENSEMBL")
		  	ID_list=na.omit(list[match(rownames(diff_gene_deseq2),list[,"ENSEMBL"]),])
		  	diff_gene_deseq2$Row.names=rownames(diff_gene_deseq2)
		  	diff_gene_deseq2=na.omit(merge(ID_list,diff_gene_deseq2,by.x="ENSEMBL",by.y= "Row.names",all=TRUE))
		  	tmp=select(org.Mm.eg.db, keys=diff_gene_deseq2$ENSEMBL, columns="GO", keytype="ENSEMBL")
			ensembl_go=unlist(tapply(tmp[,2],as.factor(tmp[,1]),function(x) paste(x,collapse ="|"),simplify =F))
			diff_gene_deseq2$go=ensembl_go[diff_gene_deseq2$ENSEMBL]


			diff_up = diff_gene_deseq2[diff_gene_deseq2$up_down == 1,]$ENTREZID
			diff_down= diff_gene_deseq2[diff_gene_deseq2$up_down == -1,]$ENTREZID
			ego_up <- enrichGO(gene=diff_up,org.Mm.eg.db,ont="CC",pvalueCutoff=0.05,readable=TRUE)
			ego_down <- enrichGO(gene=diff_down,org.Mm.eg.db,ont="CC",pvalueCutoff=0.05,readable=TRUE)
			ekk_up<- enrichKEGG(gene=diff_up,organism="mmu",pvalueCutoff=0.05)
			ekk_down <- enrichKEGG(gene=diff_down,organism="mmu",pvalueCutoff=0.05)
			write.csv(summary(ekk_up),"irt_KEGG_up-enrich.csv",row.names =F)
			write.csv(summary(ekk_down),"irt_KEGG_down-enrich.csv",row.names =F)
			write.csv(summary(ego_up),"irt_GO_up-enrich.csv",row.names =F)
			write.csv(summary(ego_down),"irt_GO_down-enrich.csv",row.names =F)
		#配对检验
		dge_et <- exactTest(dge) #检验
		 
		topTags(dge_et)
		write.csv(topTags(dge_et, n = nrow(dgelist$counts)), 'exactTest.csv', quote = FALSE)        #输出主要结果
		 
		dge_de <- decideTestsDGE(dge_et, adjust.method = 'fdr', p.value = 0.05)   #查看默认方法获得的差异基因
		summary(dge_de)
		 
		detags <- rownames(dge)[as.logical(dge_de)]
		plotSmear(dge_et, de.tags = detags, cex = 0.5)      #作图观测
		abline(h = c(-1, 1), col = 'gray', lty = 2)
		#voom 线性建模（limma）
			limma_voom <- voom(dgelist_norm, design, plot = TRUE)
			 
			fit <- lmFit(limma_voom, design)  #拟合
			fit <- eBayes(fit)
			 
			topTable(fit, coef = 2)
			write.csv(topTable(fit, coef = 2, number = nrow(dgelist$counts)), 'limma_voom.csv', quote = FALSE)      #输出主要结果
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
			path<- pathview(gene.data = res$ENTREZID,pathway.id = "mmu04110", species = "mmu",imit = list(gene=max(abs(geneList)), cpd=1))








