# 加载包

library(data.table)  # 用于高效处理大数据集的库
library(dplyr)       # 数据操作和转换的库
library(tidyverse)   # 数据处理和可视化的综合库
setwd("E:\\BaiduNetdiskDownload\\王瑜多_毕业论文代码数据\\dyang")

#导入数据
stad_mrna <-read.table("E:\\BaiduNetdiskDownload\\王瑜多_毕业论文代码数据\\王瑜多_毕业论文代码数据\\第三章\\data\\stomach\\mRNA\\combined_RNA_Expr.txt",header=TRUE)
stad_mirna <-read.table("E:\\BaiduNetdiskDownload\\王瑜多_毕业论文代码数据\\王瑜多_毕业论文代码数据\\第三章\\data\\stomach\\miRNA\\combined_miRNA_Expr.txt",header=TRUE)
#基因转换信息
dlbc.pro <- fread("E:\\BaiduNetdiskDownload\\王瑜多_毕业论文代码数据\\dyang\\gencode.v22.annotation.gene.probeMap", header = T, sep = '\t', data.table = F)
head(dlbc.pro)
# 提取前两列用于进行转换
dlbc.pro <- dlbc.pro[ , c(1, 2)]
head(dlbc.pro)
dlbc.pro$id <- gsub("\\..*", "", dlbc.pro$id)
# 将行名作为新的一列添加到数据框中，列名为 "id"
stad_mrna$id <- rownames(stad_mrna)
stad_mrna.pro <- merge(dlbc.pro, stad_mrna, by.y  = "id", by.x = "id" )
head(stad_mrna.pro)[1:5, 1:5]
stad_mrna.pro <- distinct(stad_mrna.pro,gene, .keep_all = T)
rownames(stad_mrna.pro) <- stad_mrna.pro$gene
stad_mrna.pro <- stad_mrna.pro[ , -c(1,2)]
write.csv(stad_mrna.pro,'stad_gene_414_58387.csv',row.names = TRUE)

crc_mirna <-read.table("E:\\BaiduNetdiskDownload\\王瑜多_毕业论文代码数据\\王瑜多_毕业论文代码数据\\第三章\\data\\colon\\miRAN\\combined_miRNA_Expr.txt",header=TRUE)
crc_mrna <-read.table("E:\\BaiduNetdiskDownload\\王瑜多_毕业论文代码数据\\王瑜多_毕业论文代码数据\\第三章\\data\\colon\\mRNA\\combined_RNA_Expr.txt",header=TRUE)
# 将行名作为新的一列添加到数据框中，列名为 "id"
crc_mrna$id <- rownames(crc_mrna)
crc_mrna.pro <- merge(dlbc.pro, crc_mrna, by.y  = "id", by.x = "id" )
head(crc_mrna.pro)[1:5, 1:5]
crc_mrna.pro <- distinct(crc_mrna.pro,gene, .keep_all = T)
rownames(crc_mrna.pro) <- crc_mrna.pro$gene
crc_mrna.pro <- crc_mrna.pro[ , -c(1,2)]
write.csv(crc_mrna.pro,'crc_gene_525_58387.csv',row.names = TRUE)


#获取分组信息
#gp_stad <- substring(colnames(stad_mrna.pro), 14, 15)
#table(gp_stad)
merged_mrna_sample <- merge(stad_mrna.pro, crc_mrna.pro, by = "row.names")
rownames(merged_mrna_sample)<-merged_mrna_sample$Row.names
merged_mrna_sample <-merged_mrna_sample[,-1]
gp <- substring(colnames(merged_mrna_sample), 14, 15)
table(gp)
mrna_tumor <- merged_mrna_sample[, as.numeric(gp) < 10]
mrna_normal <- merged_mrna_sample[, as.numeric(gp) >= 10]
merged_mrna_sample <- cbind(mrna_tumor, mrna_normal)
group <- c(rep('tumor', ncol(mrna_tumor)), rep('normal', ncol(mrna_normal)))
group <- factor(group, levels = c("normal", "tumor"))
table(group)

merged_mirna_sample <- merge(stad_mirna, crc_mirna, by = "row.names")
rownames(merged_mirna_sample)<-merged_mirna_sample$Row.names
merged_mirna_sample <-merged_mirna_sample[,-1]
gp <- substring(colnames(merged_mrna_sample), 14, 15)
table(gp)
mirna_tumor <- merged_mirna_sample[, as.numeric(gp) < 10]
mirna_normal <- merged_mrna_sample[, as.numeric(gp) >= 10]
merged_mirna_sample <- cbind(mirna_tumor, mirna_normal)
group <- c(rep('tumor', ncol(mirna_tumor)), rep('normal', ncol(mirna_normal)))
group <- factor(group, levels = c("normal", "tumor"))
table(group)

