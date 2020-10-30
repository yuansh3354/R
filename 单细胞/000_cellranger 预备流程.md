### 环境配置（这里的环境配置和 RNA 流程的环境配置一样）
```bash
conda install -y sra-tools
conda install -y trimmomatic
conda install -y cutadapt multiqc 
conda install -y trim-galore
conda install -y star hisat2 bowtie2
conda install -y subread tophat htseq bedtools deeptools
conda install -y salmon
source deactivate #注销当前的rna环境
```

### 下载 sra 数据
```bash
# 使用 prefetch 下载，这个下载及其的慢
cat SRR_Acc_List-2586-4.txt |while read i
do prefetch $i -O `pwd` && echo "** ${i}.sra done **"
done
# 这个是使用 ascp 下载
cat SRR_Acc_List.txt|while read id
do
x=$(echo $id | cut -b1-6)
y=$(echo $id | cut -b10-10)
echo $id
ascp -QT -l 300m -P33001  -i \
${wkd}/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/$x/00$y/$id/ ./
done
```

### sra 转换为 fastq
```bash
#--gzip将生成的结果fastq文件进行压缩
cd fq
 for i in $wkd/dataset/*sra
do
        echo $i
        time fastq-dump --gzip --split-files ./$i
done

# 修改 fastq 文件名称
cat SRR_Acc_List-9245-3.txt | while read i ;
do
    mv ${i}_1*.gz ${i}_S1_L001_I1_001.fastq.gz;\
    mv ${i}_2*.gz ${i}_S1_L001_R1_001.fastq.gz;\
    mv ${i}_3*.gz ${i}_S1_L001_R2_001.fastq.gz\
done
```

### 质控
```bash
# 以P2586-4为例
mkdir -p $wkd/qc
cd $wkd/qc
find $wkd/raw/P2586-4 -name '*R1*.gz'>P2586-4-id-1.txt
find $wkd/raw/P2586-4 -name '*R2*.gz'>P2586-4-id-2.txt
cat P2586-4-id-1.txt P2586-4-id-2.txt >P2586-4-id-all.txt

cat P2586-4-id-all.txt| xargs fastqc -t 20 -o ./
```

### cellranger 软件下载

> 要去下载新的版本

```bash
# 2.0版本下载(732M)
curl -o cellranger-2.0.2.tar.gz "http://cf.10xgenomics.com/releases/cell-exp/cellranger-2.0.2.tar.gz?Expires=1557256518&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cDovL2NmLjEweGdlbm9taWNzLmNvbS9yZWxlYXNlcy9jZWxsLWV4cC9jZWxscmFuZ2VyLTIuMC4yLnRhci5neiIsIkNvbmRpdGlvbiI6eyJEYXRlTGVzc1RoYW4iOnsiQVdTOkVwb2NoVGltZSI6MTU1NzI1NjUxOH19fV19&Signature=HoJUuPo4iTFdQgzFU1GH7uKf3uGitQxTjB6WOA9qGPlejf7tNcBPjO65WuSUZ~w8WWdeAvky-oV7XGfheY-bUr2b7QHr7jQEqc84cyU~PLvT~fYjkgC7cG7nlpbJOT~b7U~YH9amvR~SCLlyynp7scPDIA~9~keCYrIPgevTf2QyktybuSyjNTwugefOic~~XFkc9lrS~WQ9MNA1CLl4ExlQKsxWS77PEB6mwrMZXX65obDnZW9fIs3dIny6H5YoadbkgmsT52jmLien6PsG1g2jpAO90pPuHoru8LL64Q9gmB3I0nJAqi3EmrO3GKnUpHUhGb6doKmjSN6XccpmsA__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"
# 2.1版本下载
curl -o cellranger-2.1.1.tar.gz "http://cf.10xgenomics.com/releases/cell-exp/cellranger-2.1.1.tar.gz?Expires=1557260110&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cDovL2NmLjEweGdlbm9taWNzLmNvbS9yZWxlYXNlcy9jZWxsLWV4cC9jZWxscmFuZ2VyLTIuMS4xLnRhci5neiIsIkNvbmRpdGlvbiI6eyJEYXRlTGVzc1RoYW4iOnsiQVdTOkVwb2NoVGltZSI6MTU1NzI2MDExMH19fV19&Signature=RNQd-gTASTQhtnUSBfQWrnqo6Pyy2wDXtV5tlxkG97727GvoRhMqFXbEsz4gJl2BMckdVvW3S1tZRwRo5pmxPzmhq-8RKxf99pGqlzo84HYqhbIRkxXlIbLbj-u3PUJqo8cesWpbSVSKkS2TCNS-9GMFNieQswqMS2-DN4BqoBOJnWr7T4wlOMd9hypXWwOsW2P2fqaM-WP2ooMyo-oIxm3y9gDghXdDEP5lvHU7GCQcFGGexkdIrD6S5p8JPJ1DB5XieGrtEuP1YVp6tLMGXFoRWXS8dQLI1egWDYlOuRaiQgLIb3o5ZxBg5NpzLPP5kDHMAVzJFdBpf~~rkyNYTA__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"

# 解压 写入环境变量
tar zxvf cellranger-2.0.2.tar.gz
export PATH=/home/biosoft/cellranger-2.2.0:$PATH
```

### 软件检测
```bash

cellranger testrun --id=tiny
# 我使用了12个CPU，大约需要20分钟检查完
# 如果成功完整地安装的话，最后会给出这样一个报告：
cellranger testrun (2.0.2)
Copyright (c) 2017 10x Genomics, Inc. All rights reserved.
-------------------------------------------------------------------------------

Running Cell Ranger in test mode...

Martian Runtime - 2.0.2-2.2.2

Running preflight checks (please wait)...
[runtime] (ready)           ID.tiny.SC_RNA_COUNTER_CS.SC_RNA_COUNTER.SETUP_CHUNKS
[runtime] (split_complete) ID.tiny.SC_RNA_COUNTER_CS.SC_RNA_COUNTER.SETUP_CHUNKS
...

Pipestance completed successfully!
```

### 下载注释信息
```bash
curl -O http://cf.10xgenomics.com/supp/cell-exp/refdata-cellranger-GRCh38-1.2.0.tar.gz
# 然后解压
tar -xzvf refdata-cellranger-GRCh38-1.2.0.tar.gz
```