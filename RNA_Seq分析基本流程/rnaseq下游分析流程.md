```R
load(file = "expr_df.Rdata")
```

![img](https://pic4.zhimg.com/80/v2-58050097fedd1a21aee954520b23c6ab_hd.jpg)

第一列是ensemble的ID，后面是3V3的样本

---

### 构建metadata文件

> 使用DEseq分析时，需要制作一个分组指示信息

```R
metadata <- data.frame(sample_id = colnames(expr_df)[-1])
sample <- rep(c("con","treat"),each=3)
metadata$sample <- relevel(factor(sample),"con")
```

![img](https://pic2.zhimg.com/80/v2-aee3740450087c05bcfa44243f6b40f1_hd.jpg)

### **构建dds对象**

这一步由`DESeqDataSetFromMatrix`这个函数来完成，他需要输入我们的表达矩阵，制作好的metadata，还要制定分组的列，在这里是sample，最后一个tidy的意思是，我们第一列是基因ID，需要自动处理。

```R
library(DESeq2)
## 第一列有名称，所以tidy=TRUE
dds <-DESeqDataSetFromMatrix(countData=expr_df, 
                             colData=metadata, 
                             design=~sample,
                             tidy=TRUE)
```

### **数据过滤**

别看有这么多行，不是每一个就要都要表达的，如果一行基因在所有样本中的counts数小于等于1，我们就把他删掉，实际上，不做这一步，对差异分析的结果没有影响，可能会对GSEA的结果有影响。

```R
dds <- dds[rowSums(counts(dds))>1,]
```

### **样本聚类**

有时候，我们会把实验样本的的标签搞错，比如，明明是加药组，但是有一组没忘记加药，如果对他聚类，他会被划归为正常组，这个是需要我们知道的。

用vst来标化数据，实际上还有rlog方法，或者就是log2的方法，官网推荐< 30个样本用rlog，大于30个样本用vst，速度更快，这里我们不要计较那么多了，就用vst，因为真实的TCGA数据，样本往往大于30个。

```R
vsd <- vst(dds, blind = FALSE)
```

再用dist这个函数计算样本间的距离

```R
sampleDists <- dist(t(assay(vsd)))
```

用hclust来进行层次聚类

```R
hc <- hclust(sampleDists, method = "ward.D2")
```

然后画图

```R
plot(hc, hang = -1)
```

![img](https://pic3.zhimg.com/80/v2-7ecfda4346034efef3429fe6ea869316_hd.jpg)

---

```r
library(factoextra)
res <- hcut(sampleDists, k = 2, stand = TRUE)
# Visualize
fviz_dend(res, 
          # 加边框
          rect = TRUE,
          # 边框颜色
          rect_border="cluster",
          # 边框线条类型
          rect_lty=2,
          # 边框线条粗细
          lwd=1.2,
          # 边框填充
          rect_fill = T,
          # 字体大小
          cex = 1,
          # 字体颜色
          color_labels_by_k=T,
          # 平行放置
          horiz=T)
```

![img](https://pic1.zhimg.com/80/v2-0de0f153f6975acfc74bdf12eae02e0c_hd.jpg)

---

### **Deseq2 计算**

主程序是Deseq这个函数，里面顺序执行了一系列函数，每一步都可以单独运行。这一步，只有6个样本基本上就是10s以内，如果是1000个样本，小电脑跑不过去，跑过去也需要5个小时以上，很耗时间。

```text
dds <- DESeq(dds)
```

![img](https://pic2.zhimg.com/80/v2-979aed99304a67e3fbc04d66149357ad_hd.jpg)

得到dds之后，我们可以通过counts这个函数得到能作图的标注化后的counts数据，他矫正了样本间测序的深度，使得样本间可以直接比较。

```text
normalized_counts <- as.data.frame(counts(dds, normalized=TRUE))
```

![img](https://pic4.zhimg.com/80/v2-f48e48caa9020a3f2eb3766baa5ffd1b_hd.jpg)

Deseq2内置了一个画图函数，可以方便地制定基因作图

```text
plotCounts(dds, gene = "ENSG00000172183", intgroup=c("sample"))
```

![img](https://pic4.zhimg.com/80/v2-514e15ea9957151ef45e84b47ff621d3_hd.jpg)

这个功能本质上没有什么用，但是，可以提前确定实验的质量。比如，你的两组是敲减某个基因以及对照组，通过制定那个基因作图，就可以看出实验有没有成功，如果这个基因没有任何改变，也可以不用往下做了。回去重新做实验送样本吧。

我们说过，这种图自己看就可以，给老板看就算了，你得美化一下，那么他里面有个内置的参数returnData可以返回作图数据。
一旦返回数据，我们就可以用ggplot2自己简单画一下。

```text
plotdata <- plotCounts(dds, gene = "ENSG00000172183", intgroup=c("sample"),returnData = T)
library(ggplot2)
ggplot(plotdata,aes(x=sample,y=count,col=sample))+
  geom_point()+
  theme_bw()
```

![img](https://pic1.zhimg.com/80/v2-eafd8eeff1bf544a1150b10a7cc03ea4_hd.jpg)



现在看起来平淡无奇，如果样本多，可以画出点图配boxplot，如果是配对样本，那么还可以画出配对的图。

有一点要强调一下，vst，以及log2的标准化，跟标准化的counts不是一个概念。前者是为了以后的聚类，热土，PCA分析，比如，我们计算样本间的距离就是用的vst 标化的数据，而标准化的counts是为了差异作图，你看纵坐标就会发现，他的数值一般很大。

### **LogFC的矫正**

这一步，对于依赖logFC变化值的分析，很重要，比如GSEA分析。我们画一个MAplot图，看图说话

```text
contrast <- c("sample", "treat", "con")
dd1 <- results(dds, contrast=contrast, alpha = 0.05)
plotMA(dd1, ylim=c(-2,2))
```

![img](https://pic1.zhimg.com/80/v2-854b8f098d6bc11e3629327d544b77f0_hd.jpg)

MA-plot上的点是每一个基因。横坐标是标准化后的counts平均值，纵坐标是log后的变化值。红色的点是p.adjust的值小于0.05的基因，由counts函数中的alpha 参数指定。
我们发现在左侧，有很多counts很小的基因，发生了很大的变化，但是没有明显意义。类似于从(1,3,9)变成了（20,12,3）他们的counts很小，波动性很大，对logFC产生了很大的影响。
GSEA分析中，排序就是按照logFC来进行的。按照这个结果往下做，GSEA那里，富集不到任何条目。

那就需要矫正。用的函数是lfcShrink，有很多参数，我们只演示一种

```text
dd2 <- lfcShrink(dds, contrast=contrast, res=dd1)
plotMA(dd2, ylim=c(-2,2))
```

![img](https://pic3.zhimg.com/80/v2-3cbb44e1a5d41e3ea7aae2080e17a28e_hd.jpg)



这样，原先那些波动性很大的基因，就被矫正了。而此时有LogFC以，红色的点为主。

### **差异分析的结果**

这个结果实际上已经通过counts函数获得了，我们不在担心，处理组和对照组完全相反这种情况，因为他内置了参数设定比较组。比如，我们有5个处理组，我们不需要做5次Deseq，我们在results中指定即可。

用summary这个函数，可以看到差异分析的结果,高表达和低表达的比例。低丰度基因所占的比例。

```text
summary(dd2, alpha = 0.05)
```

![img](https://pic4.zhimg.com/80/v2-004699d33047825715c1535f9a324857_hd.jpg)

再把差异分析的结果转化成data.frame的格式

```text
library(dplyr)
library(tibble)
res <- dd2 %>% 
  data.frame() %>% 
  rownames_to_column("gene_id")
```

![img](https://pic3.zhimg.com/80/v2-52c15d9fc7b58a35e78277af00d0d6e6_hd.jpg)

### **基因ID转换**

以前我们从gtf文件转换的, 但是我们需要gtf文件来提取mRNA以及lncRNA，就顺手做了ID转换，而且，mRNA和lncRNA是分别做的Deseq2，这从原理上来讲，是有问题的，Deseq2矫正了测序的深度，而这个深度应该是所有基因算在一起的深度，不应该分开来算。
[TGCA数据的标准化以及差异分析](https://link.zhihu.com/?target=https%3A//mp.weixin.qq.com/s%3F__biz%3DMzIyMzA2MTcwMg%3D%3D%26mid%3D2650733328%26idx%3D1%26sn%3D0756a9a60a1d13478d43b11aab78cb72%26chksm%3Df029a8b9c75e21af023393f2a82204fb313fc991cbc96c31a735b513b9018569bba412dd7ea5%26token%3D1524745954%26lang%3Dzh_CN%23rd)
用两个包来转换，得到ENTREZID用于后续分析，得到SYMBOL便于识别。

```text
library(AnnotationDbi)
library(org.Hs.eg.db)
res$symbol <- mapIds(org.Hs.eg.db,
                     keys=res$gene_id,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
res$entrez <- mapIds(org.Hs.eg.db,
                     keys=res$gene_id,
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
```

![img](https://pic3.zhimg.com/80/v2-818b03d1013ad6e7a37495cae76e3c86_hd.jpg)

有了这个文件，里面有logFC，p值，还有基因名称，我们可以完成GO，KEGG，热图，火山图所有操作。这一部分内容参考
[最有诚意的GEO数据库教程](https://link.zhihu.com/?target=https%3A//mp.weixin.qq.com/s%3F__biz%3DMzIyMzA2MTcwMg%3D%3D%26mid%3D2650733210%26idx%3D1%26sn%3D76f5609282d92ba24644729fede24897%26chksm%3Df029a933c75e2025bc8c54a1ad046c3d5f837849d8f68969aa66bbbe58e194e3df0ed981f5d9%26token%3D1524745954%26lang%3Dzh_CN%23rd)
这个帖子目前已经有接近4000次访问，这在我们这样一个小号是不容易的，靠的是真诚。

### **制作geneList**

我们在那个帖子里面并没有讲GSEA分析，今天来展示一下。原理略过。我们这里还是用Y叔的神包clusterprofier，神包虽好，记得引用。

使用这个包做GSEA，要制作一个genelist，这个是一个向量，他的内容是排序后的logFC值，他的名称是ENTREZID，而这两个我们都是不缺的，在上一步得到的差异结果中。

```text
library(dplyr)
gene_df <- res %>% 
  dplyr::select(gene_id,log2FoldChange,symbol,entrez) %>% 
  ## 去掉NA
  filter(entrez!="NA") %>% 
  ## 去掉重复
  distinct(entrez,.keep_all = T)
```

![img](https://pic3.zhimg.com/80/v2-2202db8e1b6812f5ef8cbe9331e64246_hd.jpg)

制作genelist的三部曲

```text
## 1.获取基因logFC
geneList <- gene_df$log2FoldChange
## 2.命名
names(geneList) = gene_df$entrez
## 3.排序很重要
geneList = sort(geneList, decreasing = TRUE)
```

看一下这个genelist，增加感性的理解

```text
head(geneList)
```

![img](https://pic3.zhimg.com/80/v2-bf0a4c9802f8743a767c678c142b415a_hd.png)

### **GSEA分析**

完成KEGG的GSEA分析

```text
library(clusterProfiler)
## 没有富集到任何数据
gseaKEGG <- gseKEGG(geneList     = geneList,
                    organism     = 'hsa',
                    nPerm        = 1000,
                    minGSSize    = 20,
                    pvalueCutoff = 0.05,
                    verbose      = FALSE)
```

作图展示富集分布图

```text
library(ggplot2)
dotplot(gseaKEGG,showCategory=12,split=".sign")+facet_grid(~.sign)
```

![img](https://pic2.zhimg.com/80/v2-c6144ad5a6463f524844fccb0a18f919_hd.jpg)



这时候，我们看到有一些通路是被激活的，有一些通路是被抑制的。比如Cell cycle是被抑制的，我们可以选取单个通路来作图。



把富集的结果转换成data.frame,找到Cell cycle的通路ID是hsa04110

```text
gseaKEGG_results <- gseaKEGG@result
```

![img](https://pic2.zhimg.com/80/v2-f3cf0708188dbdf333124fc9bca4c941_hd.jpg)

使用gseaplot2把他画出来

```text
library(enrichplot)
pathway.id = "hsa04110"
gseaplot2(gseaKEGG, 
          color = "red",
          geneSetID = pathway.id,
          pvalue_table = T)
```

![img](https://pic2.zhimg.com/80/v2-4ab3f004456aef9767020ce208d72df5_hd.jpg)

也可以画出一个激活的

```text
pathway.id = "hsa04060"
gseaplot2(gseaKEGG, 
          color = "red",
          geneSetID = pathway.id,
          pvalue_table = T)
```

![img](https://pic4.zhimg.com/80/v2-37eb63076683e3db89e1d8748af42047_hd.jpg)

### **pathview 展示**

我们现在知道cell cycle是被抑制的，如果还想看一下这个通路里面的基因是如何变化的，应该怎么办呢，pathview 可以帮到我们。

```text
library(pathview)
pathway.id = "hsa04110"
pv.out <- pathview(gene.data  = geneList,
                   pathway.id = pathway.id,
                   species    = "hsa")
```

![img](https://pic1.zhimg.com/80/v2-9050509771304d6be14b00596f28f82c_hd.jpg)

改变一下参数，可以得到另外一种构图

```text
library(pathview)
pathway.id = "hsa04110"
pv.out <- pathview(gene.data  = geneList,
                   pathway.id = pathway.id,
                   species    = "hsa",
                   kegg.native = F)
```

![img](https://pic4.zhimg.com/80/v2-3891ba8013977d8a6eced457c376f6ff_hd.jpg)

一眼看过去，都是绿的，说明这个通路确实是被抑制了，还可以在图上缕一缕，哪些是核心分子，一般说来，越往上游越核心。

### **总结**

写到这里，GEO的分析，TCGA的基本分析，RNA-seq的基本分析都写完了。这几个帖子可以把大部分的培训班给搞定。里面的内容，只要会一点R语言，就可以重复。说到底，R语言的培训，最应该培训的是基本技能。

但是你有没有发现，虽然图做的这么好看，我们总觉得还欠缺了些什么。这也是目前生物信息分析和实验结合的痛点。

> 我如果开题，你这一通分析，还是没有让我确定要研究的核心基因。

而这个，是一篇帖子，或者普通培训班不能给予的。

> 如果未来的生信培训要出点什么彩，这个一定是个方向。

生信技术是通用的，优点就是可被重复，可以被写成教程，但是挖掘能做实验可发文章的核心基因，目前并没有系统教程，类似于玄学。这需要我们长年累月地积累，才能有点感觉，而我觉得这是科研人员做生信分析的核心技能。在这方面，我也是个初学者。

> 否则，我们就是个能做实验会做分析的技术人员，谈不上一个独立的科研人。