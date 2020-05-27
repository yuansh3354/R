##### 导入R包

```R

library(ggplot2)
library(patchwork)

```

##### 构建基本图形

```R

p1 <- ggplot(mtcars) + 
  geom_point(aes(mpg, disp)) + 
  ggtitle('Plot 1')
p2 <- ggplot(mtcars) + 
  geom_boxplot(aes(gear, disp, group = gear)) + 
  ggtitle('Plot 2')

p3 <- ggplot(mtcars) + 
  geom_point(aes(hp, wt, colour = mpg)) + 
  ggtitle('Plot 3')

p4 <- ggplot(mtcars) + 
  geom_bar(aes(gear)) + 
  facet_wrap(~cyl) + 
  ggtitle('Plot 4')

```

##### 组合

```R

# 普通组合
1 + p2 + p3 + p4 + plot_layout(nrow = 3, byrow = FALSE)

# 堆积组合
p1 / p2

# 并排组合
p1 | p2

# 打包组合
(p1 | (p2 / p3)) + 
  plot_annotation(title = 'The surprising story about mtcars')

# 图片编号
p1 + p2 + p3 + 
  plot_annotation(tag_levels = 'A')

```
