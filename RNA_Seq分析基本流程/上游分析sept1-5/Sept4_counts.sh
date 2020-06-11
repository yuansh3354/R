# Global variable
#########################################################
wkd='/Volumes/Lexar/project'
cd $wkd/align 
gtf="/Volumes/Elements/reference/gtf/gencode/mouse/gencode.vM23.chr_patch_hapl_scaff.annotation.gtf.gz"   

source activate rna
#########################################################

start=$(date +%s) #脚本启始时间
echo "脚本开始执行"
echo "开始时间："`date` 
#########################################################
# Main program

featureCounts -T 6 -p -t exon -g gene_id \
  -a $gtf -o  all.id.txt  $wkd/bam/*.bam \
   1>counts.id.log 2>&1 





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


