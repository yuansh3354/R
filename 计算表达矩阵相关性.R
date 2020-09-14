library(magrittr)

##############################
##########   注意   ##########
##############################

# 1.a,b矩阵行为基因ID，列为样本ID
# 2. a,b矩阵列的个数必需一样，即样本个数一样
# 3. 行可以不一样

a = read.csv('data1.csv',head = T, row.names =1)
b = read.csv('data2.csv',head = T, row.names =1)

cor(t(a),t(b))
