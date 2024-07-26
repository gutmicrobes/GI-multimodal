#寻找差异基因

setwd("E:\\BaiduNetdiskDownload\\王瑜多_毕业论文代码数据\\dyang")
crc_mirna <-read.table("E:\\BaiduNetdiskDownload\\王瑜多_毕业论文代码数据\\王瑜多_毕业论文代码数据\\第三章\\data\\colon\\miRAN\\combined_miRNA_Expr.txt",header=TRUE)
stad_mirna <-read.table("E:\\BaiduNetdiskDownload\\王瑜多_毕业论文代码数据\\王瑜多_毕业论文代码数据\\第三章\\data\\stomach\\miRNA\\combined_miRNA_Expr.txt",header=TRUE)
stad_mrna <-read.csv('stad_gene_414_58387.csv',row.names = 1)
crc_mrna <-read.csv('crc_gene_525_58387.csv',row.names = 1)
#合并每种组学数据并区分癌症与正常人
mrna_all <- merge(stad_mrna, crc_mrna, by = "row.names")
rownames(mrna_all)<-mrna_all$Row.names
mrna_all <-mrna_all[,-1]
mirna_all <- merge(stad_mirna, crc_mirna, by = "row.names")
rownames(mirna_all)<-mirna_all$Row.names
mirna_all <-mirna_all[,-1]

gp <- substring(colnames(mrna_all), 14, 15)
table(gp)
mrna_tumor <- mrna_all[, which(as.numeric(gp) < 10)]
mrna_normal <- mrna_all[, which(as.numeric(gp) >= 10)]
gp <- substring(colnames(mirna_all), 14, 15)
table(gp)
mirna_tumor <- mirna_all[, which(as.numeric(gp) < 10)]
mirna_normal <- mirna_all[, which(as.numeric(gp) >= 10)]

#对每个样本创建患病否的标签
gp <- substring(colnames(mrna_all), 14, 15)
# 创建标签向量，如果 as.numeric(gp) < 10 则为1，否则为0
label <- ifelse(as.numeric(gp) < 10, 1, 0)
# 将标签向量的名字设置为mrna_all的行名
label <-as.data.frame(label)
rownames(label) <- colnames(mrna_all)
#免疫细胞计算
library(dplyr)
library(tidyverse)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(CIBERSORT)
#读取包自带的LM22文件（免疫细胞特征基因文件）
sig_matrix <- system.file("extdata", "LM22.txt", package = "CIBERSORT")
data(LM22)
LM22[1:4,1:4]
mixed_expr <-mrna_all
results <- cibersort(sig_matrix = LM22, mixture_file = mixed_expr,perm = 1000,QN = F)
res_cibersort <- results[,1:22]
#剔除表达很低的免疫细胞
ciber.res <- res_cibersort[,colSums(res_cibersort) > 0] %>% as.data.frame()
write.csv(ciber.res,'immune_allmrna.csv')


ciber.res <-immune_allmrna
#占比图
# 设定保存文件路径
output_file_tiff <- "E:/BaiduNetdiskDownload/王瑜多_毕业论文代码数据/dyang/pic/allmrna_immune.tif"
# 开始保存图形到高分辨率 TIFF 文件，设置图形大小和分辨率
tiff(output_file_tiff, width = 2800, height = 1800, res = 300)
# 重置绘图参数
#par(mfrow = c(1, 1))
# 创建彩虹色板（带70%透明度）
mycol <- ggplot2::alpha(rainbow(ncol(ciber.res)), 0.7)
# 设置绘图参数，调整边距以适应图形
par(bty = "n", mgp = c(2, 0.5, 0), mar = c(5, 4, 4, 10), tcl = -0.3, las = 1, xpd = F)
# 绘制柱状图
barplot(as.matrix(t(ciber.res)),
        border = NA, # 柱子无边框
        names.arg = rep("", nrow(ciber.res)), # 无横坐标样本名
        yaxt = "n", # 先不绘制y轴
        ylab = "  Relative percentage", # 修改y轴名称
        col = mycol,
        mgp = c(2.5, 1, 0)) # 采用彩虹色板
# 补齐y轴添加百分号
axis(side = 2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1),
     labels = c("0%", "20%", "40%", "60%", "80%", "100%"))
# 添加图例，使用 x 和 y 参数手动设置位置
legend("topright", 
       inset = c(-0.30, 0), 
       legend = colnames(ciber.res), 
       xpd = TRUE,
       fill = mycol,
       cex = 0.8, 
       border = NA, 
       y.intersp = 1,
       x.intersp = 0.5,
       bty = "n")

# 关闭设备，完成图形保存
dev.off()

# 打印消息
cat("图形已保存到文件:", output_file_tiff, "\n")

#相关图
pdf("cor.pdf",height=13,width=13)
par(oma=c(0.5,1,1,1.2))
M=cor(ciber.res)
library(corrplot)
corrplot(M,
         order = "hclust",  # 根据层次聚类排序
         method = "color",  # 使用颜色方法
         addCoef.col = "black",  # 相关系数文本的颜色，注意参数名的修改
         diag = TRUE,  # 在对角线上绘制对角线元素
         col = colorRampPalette(c("blue", "white", "red"))(200),
         tl.cex = 0.8# 设置颜色渐变
)
dev.off()

#小提琴图
immune_cells_data <- cbind(ciber.res, label)
immune_cells_data <- immune_cells_data %>% rownames_to_column("sample")
immune_cells_data_long <- melt(immune_cells_data, id.vars = c("sample", "label"), variable.name = "CellType", value.name = "Percentage")
immune_cells_data_long$label <- as.factor(immune_cells_data_long$label)
immune_cells_data_long$CellType <- as.factor(immune_cells_data_long$CellType)

pl=ggplot(immune_cells_data_long, aes(x = CellType, y = Percentage, fill = label)) +
  geom_boxplot() +
  scale_y_continuous(labels = scales::percent) +
  labs(y = "Relative percentage", x = "Immune Cell Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "right")
ggsave(filename = "E:/BaiduNetdiskDownload/王瑜多_毕业论文代码数据/dyang/pic/allmrna_immune.pdf", plot = pl, device = "pdf", width = 12, height = 6)


#选择有差异的列
3,4,6,8,9,11,13,14,15,16,19,20
new_ciber <-ciber.res[, c(3, 4, 6, 8, 9, 11, 13, 14, 15, 16, 19, 20)]
write.csv(new_ciber,'shaixuan_immune.csv')
#以上完成了第一步 即把mrna，mirna，免疫细胞计算出了
#下一步进行差异基因的筛选



#差异mrna和mirna
library(limma)
library(edgeR)
# limma 包建议使用对数转换后的数据，其实我们从 Xena 下载的直接就是log2(count+1)后的数据，
# 但是这里为了展示 limma 包内部的标准化方法的使用，我们这里还是使用原始计数矩阵作为输入数据。
library(dplyr)
merged_mrna_sample <- cbind(mrna_tumor, mrna_normal)
group <- c(rep('tumor', ncol(mrna_tumor)), rep('normal', ncol(mrna_normal)))
group <- factor(group, levels = c("normal", "tumor"))


merged_mirna_sample <- cbind(mirna_tumor, mirna_normal)
group <- c(rep('tumor', ncol(mirna_tumor)), rep('normal', ncol(mirna_normal)))
group <- factor(group, levels = c("normal", "tumor"))
#这里变换即可
exprSet <- merged_mirna_sample
# 创建设计矩阵，指定组别信息
design <- model.matrix(~0 + factor(group))
colnames(design) = levels(factor(group))
rownames(design) = colnames(exprSet)
# 创建 DGEList 对象
dge <- DGEList(counts = exprSet, group = group)
# 这里我们使用上面提到的 filterByExpr() 进行自动过滤，去除低表达基因
keep <- filterByExpr(dge)
table(keep)
dge <- dge[keep, , keep.lib.sizes = FALSE]
# 归一化，得到的归一化系数被用作文库大小的缩放系数
dge <- calcNormFactors(dge)
# 使用 voom 方法进行标准化,将计数数据转换为一个正态分布的尺度
dev.new()
v <- voom(dge, design, plot = TRUE, normalize = "quantile")
# 使用线性模型进行拟合
fit <- lmFit(v$E, design)

# 和上面两个包一样，需要说明是谁比谁
con <- paste(rev(levels(group)), collapse = "-")
con
# 创建对比矩阵
cont.matrix <- makeContrasts(contrasts = c(con), levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

# 获取差异表达基因结果
tempOutput <- topTable(fit2, coef = con, n = Inf)
DEG_limma_voom <- na.omit(tempOutput)
head(DEG_limma_voom)
# 将差异表达基因结果保存到 Rdata 文件中
#save(DEG_limma_voom, file = './last_DEG_limma_mirna.Rdata')

# 画图图 —— 火山图和热图

# 加change列,标记上下调基因，可根据需求设定阈值
logFC = 1
P.Value = 0.01
k1 <- (DEG_limma_voom$P.Value < P.Value) & (DEG_limma_voom$logFC < -logFC)
k2 <- (DEG_limma_voom$P.Value < P.Value) & (DEG_limma_voom$logFC > logFC)
DEG_limma_voom <- mutate(DEG_limma_voom, change = ifelse(k1, "down", ifelse(k2, "up", "stable")))
table(DEG_limma_voom$change)
# 假设 mat 是你的矩阵
chayimrna <- head(rownames(DEG_limma_voom), 1)
chayimrna_biaoda <-merged_mirna_sample[rownames(merged_mirna_sample) %in% chayimrna,]
rt <-as.data.frame(chayimrna_biaoda)
rt <-t(rt)
rt <-as.data.frame(rt)
getout <-c("CMTM5","FAM180B","KRT24","LGI1","WNT2")#"hsa-miR-21-5p",
getout <-c("hsa-miR-182-5p")
rtout <-rt[ ,!colnames(rt)%in%getout]
tenmrna <- as.data.frame(chayimrna)

rt <-t(rt)
# 定义一个函数来检查行名的第14和15个字符
label_row <- function(row_name) {
  # 提取第14和15个字符
  num <- as.numeric(substr(row_name, 14, 15))
  if (!is.na(num) && num < 10) {
    return("tumor")
  } else {
    return("normal")
  }
}
rt$label
# 应用函数到每个行名并创建新列
rt$label <- sapply(rownames(rt), label_row)
library(ggpubr)

ggviolin(rt, "label", "`hsa-miR-21-5p`", fill = "label",
        palette = c("#00AFBB", "#FC4E07"),add = "boxplot", add.params = list(fill = "white"))

911 53
tumor1 <-rt[1:911,]
normal1 <-rt[912:964,]
# 调整图形参数，特别是 mgp 参数
par(mgp = c(3.5, 0.8, 0))  # 调整 y 轴标题的位置

# 绘制小提琴图
vioplot(tumor1[,1], normal1[,1], names = c("Tumor", "Normal"), 
        col = c("red", "blue"), 
        ylab = "`hsa-miR-21-5p`", xlab = "Group")
# 添加均值线
abline(h = mean(tumor1[,1]), col = "darkred", lwd = 2, lty = 2)  # Tumor组的均值线
abline(h = mean(normal1[,1]), col = "darkblue", lwd = 2, lty = 2)  # Normal组的均值线

mrna down stable     up 
3360  23933   2484 

mirna down stable     up 
91    386     72 

p <- ggplot(data = DEG_limma_voom, 
            aes(x = logFC, 
                y = -log10(P.Value))) +
  geom_point(alpha = 0.4, size = 3.5, 
             aes(color = change)) +
  ylab("-log10(Pvalue)")+
  scale_color_manual(values = c("blue4", "grey", "red3"))+
  geom_vline(xintercept = c(-logFC, logFC), lty = 4, col = "black", lwd = 0.8) +
  geom_hline(yintercept = -log10(P.Value), lty = 4, col = "black", lwd = 0.8) +
  theme_bw()
p
# 差异基因热图
deg_opt <- DEG_limma_voom %>% filter(DEG_limma_voom$change != "stable")

exp_brca_heatmap <- exprSet %>% filter(rownames(exprSet) %in% rownames(deg_opt))
write.csv(exp_brca_heatmap,'mirna_chayi_0607.csv',row.names = TRUE)
data_tumor = exp_brca_heatmap[,which(substr(colnames(exp_brca_heatmap), 14, 15)<=10)]
write.csv(data_tumor,'mirna_chayi_tumor0607.csv',row.names = TRUE)

mrna_next <-data_tumor
mirna_next <-data_tumor

annotation_col <- data.frame(group = group)
rownames(annotation_col) <- colnames(exp_brca_heatmap) 
library('pheatmap')
p1 <- pheatmap(exp_brca_heatmap, show_colnames = F, show_rownames = F,
               scale = "row",
               cluster_cols = F,
               annotation_col = annotation_col,
               #cluster_rows = FALSE,  # 添加这行以确保不显示行聚类树
               breaks = seq(-3, 3, length.out = 100)) 
p1

#ggsave(filename = "./heatmap_mrna0607.pdf", plot = p1, device = "pdf", width = 5, height = 6)


#数据合并
#new_ciber
#mrna_tumor
#mirna_tumor
#pig
image <-read.csv("E:\\BaiduNetdiskDownload\\王瑜多_毕业论文代码数据\\王瑜多_毕业论文代码数据\\第四章\\pig_dimna.csv")
rownames(image) <-image$X
image <-image[,-1]

mrna_tumor <-read.csv("E:\\BaiduNetdiskDownload\\王瑜多_毕业论文代码数据\\dyang\\mrna_chayi_tumor0607.csv")
rownames(mrna_tumor) <-mrna_tumor$X
mrna_tumor <-mrna_tumor[,-1]

mirna_tumor <-read.csv("E:\\BaiduNetdiskDownload\\王瑜多_毕业论文代码数据\\dyang\\mirna_chayi_tumor0607.csv")
rownames(mirna_tumor) <-mirna_tumor$X
mirna_tumor <-mirna_tumor[,-1]



ciber.res <-read.csv("E:\\BaiduNetdiskDownload\\王瑜多_毕业论文代码数据\\dyang\\shaixuan_immune.csv")
rownames(ciber.res) <-ciber.res$X
ciber.res <-ciber.res[,-1]


mrna_tumor <-mrna_next
mirna_tumor <-mirna_next
ciber.res <-new_ciber

common_samples <- Reduce(intersect, list(names(image), rownames(ciber.res), colnames(mrna_tumor), names(mirna_tumor)))
mrna1 <- mrna_tumor[, common_samples, drop = FALSE]
mirna1 <- mirna_tumor[, common_samples, drop = FALSE]
immune <-ciber.res[ common_samples,, drop = FALSE]
immune1 <-t(immune)
image1 <-image[, common_samples, drop = FALSE]
write.csv(mrna1,"下一步mrna.csv")
write.csv(mirna1,"下一步mirna.csv")
write.csv(immune1,"下一步immune.csv")
write.csv(image1,"下一步image.csv")



library(dynamicTreeCut)
library(stats)
library(fastcluster)
#biocLite("GO.db")
library(WGCNA)
library(lattice)
library(latticeExtra)
#软阈值
image1 <- image1[2:(nrow(image1)-1), ]  # 去掉第一行和最后一行
# 假设datExpr是基因表达矩阵，行是基因，列是样本
image1 <-下一步image
#log2Data <- log2(mirna1 + 1)
rt=t(image1)
standardized_data <- scale(rt)

#rt=log2(rt + 1)
#rt <- rt[2:(nrow(rt)-1), ]  # 去掉第一行和最后一行
#rt <-as.matrix(rt)
#rt <- rt[rownames(rt) != "TCGA.AA.3664.01A", ]
rt<-image1
rt <-rt[-c(180),]

rt <-下一步immune

rt <-mirna1
rt <- rt[!(rownames(rt) %in% c("hsa-miR-21-5p", "hsa-miR-192-5p")), ]

rt <-mrna1
sampleTree = hclust(dist(rt), method = "average")
par(mfrow = c(1,1)); 
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
powers = c(c(1:10), seq(from = 12, to=60, by=2))
sft = pickSoftThreshold(rt, powerVector = powers, networkType = 'signed',verbose = 5)
sft$fitIndices
sizeGrWindow(9, 5)
par(mfrow = c(1,2)); 
cex1 = 0.8;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("mRNA Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("mRNA Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
softPower <- sft$powerEstimate

#image  3  
#mirna  14 
#immune 9
#mrna   8

softPower <- 3
sampleAdjacency <- adjacency(rt, power = softPower)
TOM = TOMsimilarity(sampleAdjacency);
dissTOM = 1-TOM
transformed_data <- log(dissTOM+1)  # 加1以避免对0取对数

#distanceMatrix <- exp(-TOM)


# 可视化邻接矩阵的热图（如果数据量不是太大）
library(pheatmap)
pheatmap(transformed_data)
pheatmap(dissTOM)

rownames(dissTOM)<-rownames(rt)
colnames(dissTOM)<-rownames(rt)


write.csv(dissTOM,'dis_immune0627.csv')
dis_immune <-dissTOM

write.csv(dissTOM,'dis_mirna0627.csv')
dis_mirna <-dissTOM

write.csv(dissTOM,'dis_mrna0627.csv')
dis_mrna <-dissTOM

write.csv(dissTOM,'dis_image0627.csv')
dis_image <-dissTOM
 
setwd("E:\\BaiduNetdiskDownload\\王瑜多_毕业论文代码数据\\dyang\\pre-data")

dis_immune <-read.csv("dis_immune.csv")
rownames(dis_immune) <-dis_immune[,1]
dis_immune <-dis_immune[,-1]

dis_mirna <-read.csv("dis_mirna.csv")
rownames(dis_mirna) <-dis_mirna[,1]
dis_mirna <-dis_mirna[,-1]

dis_mrna <-read.csv("dis_mrna.csv")
rownames(dis_mrna) <-dis_mrna[,1]
dis_mrna <-dis_mrna[,-1]


dis_image <-read.csv("dis_image1.csv")
rownames(dis_image) <-dis_image[,1]
dis_image <-dis_image[,-1]


install.packages("cluster")
install.packages("factoextra")
library(cluster)
library(factoextra)

dist_matrix <-dis_image


mrna <-下一步mrna
datExpr <-mrna
library(WGCNA)
# 选择一个合适的软阈值
powers <- c(1:30)
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
powers = c(c(1:10), seq(from = 12, to=40, by=2))
sft = pickSoftThreshold(rt, powerVector = powers, networkType = 'signed',verbose = 5)
sft$fitIndices
sizeGrWindow(9, 5)
par(mfrow = c(1,2)); 
cex1 = 0.8;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("mRNA Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("mRNA Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# 通常选择scaleFreeToplogyFitIndex高，且平均连通度不为零的最小阈值
softPower <- 16
adjacency <- adjacency(t(datExpr), power = softPower)
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1-TOM

# 模块检测
geneTree = hclust(as.dist(dissTOM), method = "average")
minModuleSize <- 30
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
table(dynamicMods)
