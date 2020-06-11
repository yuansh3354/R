# Global variable
#########################################################
wkd='/Volumes/Lexar/project'
source activate rna
bin_trim_galore=trim_galore
R1=R1
R2=R2

cd ${wkd}
#########################################################
# Main program 
#########################################################
start=$(date +%s) #脚本启始时间
echo "脚本开始执行"
echo "开始时间："`date` 
#########################################################
#Sept1
ls $wkd/seq/*_${R1}.fastq.gz >1
ls $wkd/seq/*_${R2}.fastq.gz >2
paste 1 2  > config


cat $1 |while read id
do
        arr=(${id})
        fq1=${arr[0]}
        fq2=${arr[1]} 
	$bin_trim_galore -q 25 --phred33 --length 36 -e 0.1 --stringency 3 --paired -o ${wkd}/clean $fq1 $fq2 
done 

source deactivate 



#########################################################
end=$(date +%s) #脚本终止时间
take=$(( end - start ))
h=$(( take / 3600 ))
m=$((( take - 3600 * h ) / 60 ))
s=$((( take - 3600 * h - 60 * m )))
echo 脚本运行完毕，用时 ${h} 小时 ${m} 分钟 ${s} 秒
echo "结束时间："`date` 
say "The program has finished running"
#########################################################


