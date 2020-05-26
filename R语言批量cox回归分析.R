##### 1. 提取感兴趣的基因ID (IDs)
##### 2. 将IDs变成基因名(gene)
##### 3. 构建survival数据框,第一个为生存时间,第二列为状态(survival)
##### 4. 提取基因表达谱(sample)
##### 5. 转置基因表达谱,行为样本列为基因
##### 6. 合并survival和基因表达谱(cox_dat)
##### 7. 对基因和列名都要进行gsub转换,把-转为'-'


##### 主函数

top = deg
top = top[which(top$adj.P.Val<0.05 & abs(top$logFC)>1),]
top = top[order(top$adj.P.Val),]
IDs = rownames(top)
gene = paste0('ID_',IDs)
sample = df_expr[IDs,]
survival = pd[,42:43] 
survival = survival[-which(survival[,1] == 'Not available'),]
survival$status.ch1 = ifelse(survival$status.ch1 == 'alive',0,1)
sample = sample[,colnames(sample) %in% rownames(survival)]
gene = gsub(gene,pattern = '-', replacement = '_')
rownames(sample) = gene
colnames(survival) = c('time', 'status')
cox_dat = as.data.frame(cbind(survival,t(sample))) 
cox_dat[,1] = as.numeric(cox_dat[,1])
cox_dat[,2] = as.numeric(cox_dat[,2])
##### a是基因gene 注意这里如果报错 要看一下a里面有没有乱码比如@符号 这时候去掉这个基因就行


##### 函数包
library("survival")
library("survminer")
library(clusterProfiler)
cox_analy = function(gene,survival_info){
  uni_cox = function(single_gene){
    formula = as.formula(paste0('Surv(time,status)~',single_gene))
    surv_uni_cox = summary(coxph(formula, data = cox_dat)) 
    ph_hypothesis_p = cox.zph(coxph(formula,data = cox_dat))$table[1:3]
    if(surv_uni_cox$coefficients[,5]<0.05 & ph_hypothesis_p > 0.05){
      single_cox_report = data.frame(
        'ID'=single_gene,
        'beta' = surv_uni_cox$coefficients[,1],
        'Hazard_Ratio' = exp(surv_uni_cox$coefficients[,1]),
        'z_pvalue'=surv_uni_cox$coefficients[,5],
        'Wald_pvalue'= as.numeric(surv_uni_cox$waldtest[3]),
        'Likelihood_pvalue'=as.numeric(surv_uni_cox$logtest[3]))
      single_cox_report
    }
  }
  uni_cox_list = lapply(gene,uni_cox)
  do.call(rbind,uni_cox_list)
}
a = gene
uni_cox_df = cox_analy(a,cox_dat)
