# cell_type ： df中要排序的横坐标
# AD ：按Add值进行排序
# subgroup： 分组
df$cell_type = with(df, factor(cell_type, levels=cell_type[order(ave(AD, subgroup, FUN=min),AD)]))
