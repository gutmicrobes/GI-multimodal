mrna1
labels
rownames(labels)<-colnames(mrna1)
labels <-as.data.frame(labels)
# 使用mutate和recode函数替换取值
labels <- labels %>%
  mutate(value = recode(V1, `0` = "I", `1` = "II", `2` = "III"))
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
#Kruskal-Wallis
data <- t(下一步immune)
# 确保样本顺序一致
data <- data[rownames(labels),]
data <-as.data.frame(data)
data$Subtype <- labels[, 2]
# 确保所有列名是有效的R变量名
colnames(data) <- make.names(colnames(data))

# 执行Kruskal-Wallis测试
p_values <- sapply(colnames(data)[-ncol(data)], function(gene) {
  kruskal.test(as.formula(paste(gene, "~ Subtype")), data = data)$p.value
})

# 将结果转化为数据框
p_values_df <- data.frame(gene = names(p_values), p_value = p_values)

# 按p值排序
p_values_df <- p_values_df %>%
  arrange(p_value)
# 显示前20个差异显著的基因
head(p_values_df, 20)

# 选择前5个差异显著的基因
top_genes <- head(p_values_df$gene, 6)

# 将数据转换为长格式
data_long <- data %>%
  select(all_of(c(top_genes, "Subtype"))) %>%
  pivot_longer(-Subtype, names_to = "gene", values_to = "expression")
data_long$Subtype <- as.factor(data_long$Subtype)
#colnames(data_long)[colnames(data_long) == "label"] <- "Subtype"

# 绘制箱线图
ggplot(data_long, aes(x = Subtype, y = expression, fill = Subtype)) +
  geom_boxplot() +
  facet_wrap(~ gene, scales = "free_y") +
  theme_minimal() +
  labs(title = "Expression of Top Differential mRNAs", y = "Expression Level", x = "Sample Subtype")




# 提取前20个差异显著的基因
top_20_genes <- head(p_values_df$gene, 20)
# 提取前20个差异显著基因的表达数据
expression_data <- data[, top_20_genes]

# 创建样本标签颜色向量
annotation <- data.frame(Subtype = factor(data$Subtype))
rownames(annotation) <- rownames(data)

# 按标签排序数据
sorted_indices <- order(annotation$Subtype)
sorted_expression_data <- expression_data[sorted_indices, ]
sorted_annotation <- annotation[sorted_indices, , drop = FALSE]

# 确保行名和列名一致
rownames(sorted_expression_data) <- rownames(sorted_annotation)
#heatmap_colors <- colorRampPalette(c("navy", "white", "firebrick3"))(50)
heatmap_colors <- colorRampPalette(c("blue","blue","blue","white", "red","red","red"))(50)
#heatmap_colors <- colorRampPalette(c("blue","blue","white","red","red"))(500)

# 绘制热图
pheatmap(t(sorted_expression_data), 
         annotation_col = sorted_annotation, 
         show_rownames = TRUE,  # 显示行名
         show_colnames = FALSE,
         cluster_rows = TRUE,  # 聚类行
         cluster_cols = FALSE,  # 不聚类列
         treeheight_row = 0,  # 隐藏行的树状图
         color = heatmap_colors,  # 设置颜色
         scale = "row",
         main = "Heatmap of Top 20 Differential mRNAs")




#免疫细胞
top_genes <- head(p_values_df$gene, 20)

# 将数据转换为长格式
data_long <- data %>%
  select(all_of(c(top_genes, "Subtype"))) %>%
  pivot_longer(-Subtype, names_to = "gene", values_to = "expression")
data_long$Subtype <- as.factor(data_long$Subtype)
# 计算均值
df_mean <- data_long %>%
  group_by(gene, Subtype) %>%
  summarise(Mean_Expression = mean(expression))
# 转换为宽格式
df_wide <- df_mean %>%
  pivot_wider(names_from = Subtype, values_from = Mean_Expression)
# 转换为矩阵
matrix_data <- as.matrix(df_wide[,-1])
rownames(matrix_data) <- df_wide$gene
library(RColorBrewer)
# 绘制热图
pheatmap(matrix_data, 
         color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100), 
         scale = "row", 
         treeheight_row = 0,
         treeheight_col = 0,
         angle_col = 0,
         main = "Heatmap of Differential mRNAs")  # 旋转 x 轴标签


# 创建饼图数据框
df_pie <- df_mean %>%
  group_by(Subtype) %>%
  mutate(Percentage = Mean_Expression / sum(Mean_Expression) * 100)
library(RColorBrewer)
colors <- brewer.pal(n = 11, name = "Paired")  # Set3 是一个科研常用的配色方案，你可以根据需要选择其他方案

# 绘制饼图
ggplot(df_pie, aes(x = "", y = Percentage, fill = gene)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  facet_wrap(~ Subtype) +
  scale_fill_manual(values = colors) +
  theme_void() +
  labs(title = "Proportions of Immune Cell Subtypes in Different Categories", fill = "Immune Cell types")

# 绘制横向条形图
ggplot(df_pie, aes(x = Subtype, y = Percentage, fill = gene)) +
  geom_bar(stat = "identity", position = "fill",width = 0.5) +
  coord_flip() +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  labs(title = "Proportions of Immune Cells in Different Subtypes", x = "Subtype", y = "Percentage", fill = "Immune Cell types")
