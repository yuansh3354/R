# Global variable
#########################################################
wkd='/Volumes/Lexar/HBV_Rseq'
source activate base
bin_trim_galore=trim_galore
dir=${wkd}/clean

#########################################################
# Main program 
#########################################################
start=$(date +%s) #脚本启始时间
echo "脚本开始执行"
echo "开始时间："`date` 
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
echo 脚本运行完毕，用时 ${h} 小时 ${m} 分钟 ${s} 秒
echo "结束时间："`date` 
say "The program has finished running"
#########################################################


