library(stringr)
t=str_split(colnames(expset),'\\.',simplify = T) 
colnames(expset)=paste(t[,1],t[,2],t[,3],t[,4],sep = '-')
