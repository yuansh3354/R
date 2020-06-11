# Global variable
#########################################################
wkd='/Volumes/Lexar/project'
source activate rna
#########################################################

start=$(date +%s) #脚本启始时间
echo "脚本开始执行"
echo "开始时间："`date` 
#########################################################
# Main program
cd $wkd/clean 
ls *gz|cut -d"_" -f 1 |sort -u |while read id;do
ls -lh ${id}_R1_val_1.fq.gz   ${id}_R2_val_2.fq.gz 
hisat2 -p 10 -x /Volumes/Elements/reference/index/hisat/grcm38/genome -1 ${id}_R1_val_1.fq.gz -2 ${id}_R2_val_2.fq.gz  -S ${id}.hisat.sam;
done 


conda deactivate 
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


