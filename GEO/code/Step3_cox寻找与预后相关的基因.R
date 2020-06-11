### ---
### Title: "Untitled"
### Create: "Yuansh"
### Date: "5/01/2020"
### Email: yuansh3354@163.com
### output: html_document
### ---

### step2 cox分析
##### 1.构建数据:其中sample是所需表达谱,survival是所需生存信息
if(T){
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
  
} 
##### 2.构建cox模型 识别和预后相关的基因
if(T){
  library("survival")
  library("survminer")
  library(clusterProfiler)
  library(stringr)
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
  cox_IDs = str_split(uni_cox_df[,1],'_',simplify = T) 
  uni_cox_df[,1] = cox_IDs[,2]
  cox_IDs = cox_IDs[,2]
}
uni_cox_df$z_pvalue_adjust = p.adjust(uni_cox_df$z_pvalue ,method = "BH")
write.csv(uni_cox_df,paste(file,'-cox-gene.csv',sep = ''),row.names = F)
