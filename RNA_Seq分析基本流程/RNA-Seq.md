# RNA-Seq

前期准备

> 设置代理

```shell
alias proxy='export all_proxy=socks5://127.0.0.1:1086'
alias unproxy='unset all_proxy'
proxy
curl cip.cc
```



> 环境创建

```shell
# 创建名为rna的软件安装环境
# 查看当前conda环境
# 激活conda的rna环境
# 注销当前的rna环境
conda info --envs 
source activate base 
#source deactivate 
```

> 软件安装

```shell
# conda search 搜索包
# conda install 装包

conda search sratools
conda install -y sra-tools
conda install -y trimmomatic
conda install -y cutadapt multiqc 
conda install -y trim-galore
conda install -y star hisat2 bowtie2
conda install -y subread tophat htseq bedtools deeptools
conda install -y salmon
```

> 下载数据

```shell
ascp -QT -l 300m -P33001  -i ${wkd}/asperaweb_id_dsa.openssh \
 era-fasp@fasp.sra.ebi.ac.uk:sra/sra-instant/reads/ByRun/sra/SRR/SRR358/SRR3589958/SRR3589958.sra  ./
wkd=/Volumes/Lexar/project #设置工作目录

mkdir dataset && cd dataset

cat SRR_Acc_List.txt | while read id; do (prefetch  ${id} );done

```

> Sra 数据转换
>
> 默认是双端测序
>
> ***若单端要查函数具体用法***

```shell
for i in $wkd/dataset/*sra
do
        echo $i
        fastq-dump --split-3 --skip-technical --clip --gzip $i  
done
```

> 检查测序质量

```shell
# -t 是线程，理论上是越大越好，之前设置24直接被kill所以小一点
wkd=$PWD #设置工作目录
mkdir qc_out
fastqc $wkd/dataset/*.gz -t 6 -o $wkd/qc_out/

multiqc $wkd/dataset/ -o $wkd/qc_out/

```

> 过滤和去街头
>
> 用时2个小时

```shell
#PH7D-10_S10_L001_R1_001.fastq.gz 
#PH7D-10_S10_L001_R2_001.fastq.gz

mkdir $wkd/clean && cd $wkd/clean
ls $wkd/seq/*_R1_001.fastq.gz >1
ls $wkd/seq/*_R2_001.fastq.gz >2
paste 1 2  > config
chmod 777 qc.sh 
```

---

> Qc 脚本 sh

```
# Global variable
#########################################################
wkd='/Volumes/Lexar/project'
source activate base
bin_trim_galore=trim_galore
dir=${wkd}/clean

#########################################################
# Main program 
#########################################################
start=$(date +%s) #脚本启始时间
echo '脚本开始执行'
#########################################################



cat $1 |while read id
do
        arr=(${id})
        fq1=${arr[0]}
        fq2=${arr[1]} 
 $bin_trim_galore -q 25 --phred33 --length 36 -e 0.1 --stringency 3 --paired -o $dir  $fq1 $fq2  
done 

source deactivate 



#########################################################
end=$(date +%s) #脚本终止时间
take=$(( end - start ))
h=$(( take / 3600 ))
m=$((( take - 3600 * h ) / 60 ))
s=$((( take - 3600 * h - 60 * m )))
echo 脚本运行完毕用时 ${h} 小时 ${m} 分钟 ${s} 秒
#########################################################


 
```

---

> 比对代码

```shell
cd $wkd/clean 

ls *gz|cut -d"_" -f 1 |sort -u |while read id;do
ls -lh ${id}_1_val_1.fq.gz   ${id}_2_val_2.fq.gz 
hisat2 -p 10 -x /public/reference/index/hisat/hg38/genome -1 ${id}_1_val_1.fq.gz -2 ${id}_2_val_2.fq.gz  -S ${id}.hisat.sam;

done 
```

> 转换

```shell
ls *.sam|while read id ;do (samtools sort -O bam -@ 5  -o $(basename ${id} ".sam").bam   ${id});done
rm *.sam 
```

> 简历索引和统计

```shell
ls *.bam |xargs -i samtools index {}
ls *.bam |while read id ;do ( samtools flagstat -@ 1 $id >  $(basename ${id} ".bam").flagstat  );done
```

> 整合

```shell
mkdir $wkd/align 
cd $wkd/align 
# 如果一个个样本单独计数，输出多个文件使用代码是：
for fn in {508..523}
do
featureCounts -T 5 -p -t exon -g gene_id  -a /public/reference/gtf/gencode/gencode.v25.annotation.gtf.gz -o $fn.counts.txt SRR1039$fn.bam
done
# 如果是批量样本的bam进行计数，使用代码是：
mkdir $wkd/align 
cd $wkd/align 
source activate rna
gtf="/public/reference/gtf/gencode/gencode.v25.annotation.gtf.gz"   
featureCounts -T 5 -p -t exon -g gene_id  -a $gtf -o  all.id.txt  *.bam  1>counts.id.log 2>&1 &
# 这样得到的  all.id.txt  文件就是表达矩阵啦，但是，这个 featureCounts有非常多的参数可以调整。
source deactivate 
```