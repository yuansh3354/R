```bash
总的来说，Cell Ranger主要的流程有：拆分数据 mkfastq、细胞定量 count、定量组合 aggr、调参reanalyze，还有一些小工具比如mkref、mkgtf、upload、sitecheck、mat2csv、vdj、mkvdjref、testrun
```

### make fastq
***几乎用不到，并且跑一次需要三五天***
```bash
# 第一种
$ cellranger mkfastq --id=bcl \
                    --run=/path/to/bcl \
                    --samplesheet=samplesheet-1.2.0.csv
# 第二种
$ cellranger mkfastq --id=bcl \
                    --run=/path/to/bcl \
                    --csv=simple-1.2.0.csv
# 其中id指定输出目录的名称，run指的是下机的原始BCL文件目录
# 重要的就是测序lane、样本名称、index等信息
```

# count
```bash
ref=refdata-cellranger-GRCh38-1.2.0 # 参考基因组
fastqFilePath=~/fq/
cat SRR_List.txt |\
while read id; do \
(cellranger count --id=$id \ #前缀名字
                  --transcriptome=$ref \
                  --fastqs=$fastqFilePath \
                  --sample=$id \
                  --localcores=10\
                  --localmem=30);
done
# id指定输出文件存放目录名
# transcriptome指定与CellRanger兼容的参考基因组
# fastqs指定mkfastq或者自定义的测序文件
# sample要和fastq文件的前缀中的sample保持一致，作为软件识别的标志
# expect-cells指定复现的细胞数量，这个要和实验设计结合起来
# nosecondary 只获得表达矩阵，不进行后续的降维、聚类和可视化分析(因为后期会自行用R包去做)
```