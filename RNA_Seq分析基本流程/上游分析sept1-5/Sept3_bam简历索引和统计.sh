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

ls ${wkd}/clean/*.sam|while read id ;do (samtools sort -O bam -@ 5  -o $(basename ${id} ".sam").bam   ${id});done
mv ${wkd}/clean/*.bam ${wkd}/bam/
rm ${wkd}/clean/*.sam
ls ${wkd}/bam/*.bam |xargs -i samtools index {}
ls ${wkd}/bam/*.bam |while read id ;do ( samtools flagstat -@ 1 $id >  $(basename ${id} ".bam").flagstat  );done




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


