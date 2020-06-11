# Global variable
#########################################################
wkd='/Volumes/Lexar/小鼠RNA-seq/bam'
biosoft="/Users/yuansh/opt/leafcutter"
source activate rna
cd $wkd
#########################################################

start=$(date +%s) #脚本启始时间
echo "脚本开始执行"
echo "开始时间："`date` 
#########################################################
# Main program
ls *.bam > bam_path.txt
cat bam_path.txt |while read id
do
    echo Converting $id to $id.junc
    sh ${biosoft}/scripts/bam2junc.sh $id $id.junc
    echo $id.junc >> test_juncfiles.txt
done
mv *junc ../junc 
mv test_juncfiles.txt ../junc
cd ../junc
python ${biosoft}/clustering/leafcutter_cluster.py -j test_juncfiles.txt -m 50 -o npc113 -l 500000


ls *.bam.junc > group_info.txt



${biosoft}/scripts/leafcutter_ds.R --num_threads 4 \
 --exon_file=${biosoft}/leafcutter/data/gencode19_exons.txt.gz \
npc113_perind_numers.counts.gz group_info.txt

${biosoft}/scripts/ds_plots.R -e  ${biosoft}/leafcutter/data/gencode19_exons.txt.gz npc113_perind_numers.counts.gz group_info.txt leafcutter_ds_cluster_significance.txt -f 0.05

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

