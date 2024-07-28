# 加载必须的包并做参数设置
library(MASS)
library(class)
library(cluster)
library(impute)
library(Hmisc)
library(WGCNA)
options(stringsAsFactors = F)
library(DESeq2)
# rownames(labelall) <-colnames(dds)
# write.csv(labelall,"alllabel.csv")
datExpr <-mrna1
# dds <- log(datExpr+1)
# write.csv(dds,"logdata.csv")
# datExpr <-dds
# 筛选方法：mad>1且top5000
WGCNA_matrix = t(datExpr[order(apply(datExpr,1,mad),decreasing = T)[1:5000],])
datExpr_filted <- WGCNA_matrix
#datExpr_filted <- t(mrna1)

gsg = goodSamplesGenes(datExpr_filted, verbose = 3)
gsg$allOK
sampleTree = hclust(dist(datExpr_filted), method = "average")
par(mfrow = c(1,1));
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
save(datExpr_filted, file = "datExpr_filted_hclust.RData")
datExpr_filted <- datExpr_filted[!(rownames(datExpr_filted) %in% c("TCGA.BR.A4IV.01A", "TCGA.RD.A8MW.01A","TCGA.VQ.A91V.01A")), ]
#datExpr_filted <- datExpr_filted[!(rownames(datExpr_filted) %in% c("TCGA.CD.8531.01A")), ]


# Constructing a weighted gene network entails the choice of the soft thresholding power to which coexpression similarity is raised to calculate adjacency.
# Set up a bunch of power gradients(设定一些列power梯度)
powers = c(c(1:10), seq(from = 12, to=30, by=2))

sft = pickSoftThreshold(datExpr_filted, powerVector = powers, verbose = 5) #this step will take some time
# The "sft" object contains the network characteristics that calculated for each power value(在sft这个对象中保存了每个power值计算出来的网络的特???)
str(sft) 
# make figures to show the soft thresholding power
par(mfrow = c(1,2))
cex1 = 0.9
# x axis is the Soft threshold (power)，y axis is evaluation parameters for scale-free networks(纵轴是无标度网络的评估参???)
# The higher R^2 ,the more scale-free the network is.(数值越高，网络越符合无标度特征 (non-scale))
# "FitIndices" stores the characteristics of each network corresponding to its power. (fitIndices保存了每个power对应的网络的特征)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red")
# R^2 value(h) usually around 0.9. But if there's some big changes between your samples, R^2 value will be lower. 
# In my case, I use 0.7 to set the criteria.
abline(h=0.9,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# System will recommend the best power value to you. 
# But you can not totally trust it, depend on your data.
sft$powerEstimate  #I use power =3 建议3


#上一步估计的最佳 power 值
powers <-3

#获得 TOM 矩阵
adjacency <- adjacency(datExpr_filted, power = powers)
tom_sim <- TOMsimilarity(adjacency)
rownames(tom_sim) <- rownames(adjacency)
colnames(tom_sim) <- colnames(adjacency)
tom_sim[1:6,1:6]

#输出 TOM 矩阵
#write.table(tom_sim, 'TOMsimilarity.txt', sep = '\t', col.names = NA, quote = FALSE)

#TOM 相异度 = 1 – TOM 相似度
tom_dis  <- 1 - tom_sim

#层次聚类树，使用中值的非权重成对组法的平均聚合聚类
geneTree <- hclust(as.dist(tom_dis), method = 'average')
par(mfrow = c(1,1))
plot(geneTree, xlab = '', sub = '', main = 'Gene clustering on TOM-based dissimilarity',
     labels = FALSE, hang = 0.04)
k <- softConnectivity(datE=datExpr_filted,power=3)
head(k)
sizeGrWindow(10, 5)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k,main="Check Scale free topology\n")
# 连通性直方图：如果网络符合无标度拓扑特性，连通性直方图应显示大多数基因的连通性较低，而少数基因的连通性较高。
# 无标度拓扑特性检查：如果网络符合无标度拓扑特性，双对数图中的点应大致沿一条直线分布，表示连通性的频率分布遵循幂律分布。

#使用动态剪切树挖掘模块
minModuleSize <- 30  #模块基因数目
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = tom_dis,
                             deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize)

table(dynamicMods)
#模块颜色指代
dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)

plotDendroAndColors(geneTree, dynamicColors, 'Dynamic Tree Cut',
                    dendroLabels = FALSE, addGuide = TRUE, hang = 0.03, guideHang = 0.05,
                    main = 'Gene dendrogram and module colors')


# #基因表达聚类树和共表达拓扑热图
# plot_sim <- -(1-tom_sim)
# #plot_sim <- log(tom_sim)
# diag(plot_sim) <- NA
# TOMplot(plot_sim, geneTree, dynamicColors,
#         main = 'Network heatmap plot, selected genes')
# 
# 
# 
# library(pheatmap)
# # 自定义颜色映射函数
# color_palette <- colorRampPalette(c("white", "blue", "red"))(100)
# 
# # 使用 pheatmap 绘制热图
# pheatmap(plot_sim,
#          cluster_rows = geneTree,
#          cluster_cols = geneTree,
#          color = color_palette,
#          main = "Network heatmap plot, selected genes",
#          annotation_colors = list(module = dynamicColors),
#          show_rownames = FALSE,
#          show_colnames = FALSE)

#计算基因表达矩阵中模块的特征基因（第一主成分）
MEList <- moduleEigengenes(datExpr_filted, colors = dynamicColors)
#模块中每个样本的综合表达
MEs <- MEList$eigengenes
head(MEs)[1:6]

#输出模块特征基因矩阵
#write.table(MEs, 'moduleEigengenes.txt', sep = '\t', col.names = NA, quote = FALSE)

#通过模块特征基因计算模块间相关性，表征模块间相似度
ME_cor <- cor(MEs)
ME_cor[1:6,1:6]
# 自定义颜色映射函数
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
library(pheatmap)
# 绘制相关性热图
pheatmap(ME_cor, 
         color = color_palette, 
         main = "Module Eigengene Correlation Heatmap",cluster_rows = FALSE,
         cluster_cols = FALSE)
#绘制聚类树观察
METree <- hclust(as.dist(1-ME_cor), method = 'average')
par(mfrow=c(1,1))
plot(METree, main = 'Clustering of module eigengenes', xlab = '', sub = '')

#探索性分析，观察模块间的相似性
#height 值可代表模块间的相异度，并确定一个合适的阈值作为剪切高度
#以便为低相异度（高相似度）的模块合并提供依据
abline(h = 0.2, col = 'blue')
abline(h = 0.35, col = 'red')

#模块特征基因聚类树热图
plotEigengeneNetworks(MEs,"Eigengene adjacency heatmap", cex.lab = 0.8, xLabelsAngle= 90,plotDendrograms = FALSE,
                      marDendro = c(0, 4, 1, 2), marHeatmap = c(3, 4, 1, 2))


#相似模块合并，以 0.25 作为合并阈值（剪切高度），在此高度下的模块将合并
#近似理解为相关程度高于 0.75 的模块将合并到一起
merge_module <- mergeCloseModules(datExpr_filted, dynamicColors, cutHeight = 0.25, verbose = 3)
mergedColors <- merge_module$colors
table(mergedColors)

#基因表达和模块聚类树
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c('Dynamic Tree Cut', 'Merged dynamic'),
                    dendroLabels = FALSE, addGuide = TRUE, hang = 0.03, guideHang = 0.05)



datTraits <-labelall
rownames(datTraits) <-colnames(mrna1)
datTraits <- datTraits[!(rownames(datTraits) %in% c("TCGA.BR.A4IV.01A", "TCGA.RD.A8MW.01A","TCGA.VQ.A91V.01A")), ]

#患者的临床表型数据
trait <- datTraits
colnames(trait) <-c("Subtype I","Subtype II","Subtype III")

#使用上一步新组合的共表达模块的结果
module <- merge_module$newMEs

#患者基因共表达模块和临床表型的相关性分析
moduleTraitCor <- cor(module, trait, use = 'p')
moduleTraitCor[1:3,1:3]  #相关矩阵

#相关系数的 p 值矩阵
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(module))

#输出相关系数矩阵或 p 值矩阵
#write.table(moduleTraitCor, 'moduleTraitCor.txt', sep = '\t', col.names = NA, quote = FALSE)
#write.table(moduleTraitPvalue, 'moduleTraitPvalue.txt', sep = '\t', col.names = NA, quote = FALSE)
# 移除 grey 模块
moduleTraitCor <- moduleTraitCor[rownames(moduleTraitCor) != "MEgrey", ]
moduleTraitPvalue <- moduleTraitPvalue[rownames(moduleTraitPvalue) != "MEgrey", ]

#相关图绘制
textMatrix <- paste(signif(moduleTraitCor, 2), '\n(', signif(moduleTraitPvalue, 1), ')', sep = '')
dim(textMatrix) <- dim(moduleTraitCor)
par(mfrow = c(1,1))
#某个模块特征向量与某个表型呈负相关，表示该表型值增加时，该模块特征向量的值倾向于减少。
labeledHeatmap(Matrix = moduleTraitCor, main = paste('Module-trait relationships'),
               xLabels = names(trait), yLabels = rownames(moduleTraitPvalue), ySymbols = rownames(moduleTraitPvalue),
               colorLabels = FALSE,  xLabelsAngle = 0,  xLabelsAdj =0.5,
               colors = blueWhiteRed(45), cex.text = 0.7, zlim = c(-1,1),
               textMatrix = textMatrix, setStdMargins = FALSE)



# library(ggplot2)
# # 将数据转换为长格式
# library(reshape2)
# moduleTraitCor_long <- melt(moduleTraitCor)
# moduleTraitPvalue_long <- melt(moduleTraitPvalue)
# 
# # 合并相关性和 p 值数据
# moduleTrait_long <- data.frame(
#   Module = rep(rownames(moduleTraitCor), times = ncol(moduleTraitCor)),
#   Trait = rep(colnames(moduleTraitCor), each = nrow(moduleTraitCor)),
#   Correlation = moduleTraitCor_long$value,
#   Pvalue = moduleTraitPvalue_long$value
# )
# 
# # 绘制热图
# ggplot(moduleTrait_long, aes(x = Trait, y = Module, fill = Correlation)) +
#   geom_tile() +
#   scale_fill_gradient2(low = "#5573BC", mid = "white", high = "#ED1E24", midpoint = 0) +
#   geom_text(aes(label = sprintf("%.2f\n(%.1e)", Correlation, Pvalue)), size = 3) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +  # 设置 x 轴标签水平
#   labs(title = "Module-trait relationships", x = "", y = "")


#基因与模块的对应关系列表
gene_module <- data.frame(gene_name = colnames(datExpr_filted), module = mergedColors, stringsAsFactors = FALSE)
head(gene_module)

#“black”模块内的基因名称0
gene_module_select <- subset(gene_module, module == 'turquoise')$gene_name
#1
gene_module_select <- subset(gene_module, module == 'brown')$gene_name

gene_module_select <- subset(gene_module, module == 'black')$gene_name


#“black”模块内基因在各样本中的表达值矩阵（基因表达值矩阵的一个子集）
gene_select <- t(datExpr_filted[,gene_module_select])

#“black”模块内基因的共表达相似度（在 TOM 矩阵中提取子集）
tom_select <- tom_sim[gene_module_select,gene_module_select]

# kwithin0 <-rowSums(tom_select)-1
# kwithin0 <-as.data.frame(kwithin0)
# rownames(kwithin0) <-rownames(tom_select)
# 
# kwithin1 <-rowSums(tom_select)
# kwithin1 <-as.data.frame(kwithin1)
# rownames(kwithin1) <-rownames(tom_select)
# 
# 
# kwithin2 <-rowSums(tom_select)
# kwithin2 <-as.data.frame(kwithin2)
# rownames(kwithin2) <-rownames(tom_select)
# 
# colorh1<-mergedColors
# kbe=abs(cor(datExpr_filted,use="p"))^softPower 
# degrees=intramodularConnectivity(kbe, mergedColors) 
# kwithin0
# 定义一个简单的加法函数
calcukwi <- function(color) {
  moduleGenes <- which(mergedColors == color)
  moduleExpr <- datExpr_filted[, moduleGenes]
  
  # 计算模块内的加权邻接矩阵
  kbe_module <- abs(cor(moduleExpr, use = "p"))^powers
  
  # 计算模块内的连通性
  connectivity <- rowSums(kbe_module) - 1
  result <- as.data.frame(connectivity)
  return(result)
}
kwithin0 <-calcukwi("turquoise")
kwithin1 <-calcukwi("brown")
kwithin2 <-calcukwi("black")


HubGenes <- chooseTopHubInEachModule(datExpr_filted,mergedColors)
# black         blue        brown         cyan        green  greenyellow      magenta 
# "COL5A2"      "SPAG5" "AC010976.2"     "KCTD19"      "TACR1"      "P2RY8"      "SNTB1" 
# midnightblue         pink       purple          red       salmon          tan    turquoise 
# "PCLO"      "RANP1"       "APOB"       "IL16"     "FCGR1A"       "POMC"       "MPDZ" 
# yellow 
# "IGKC" 

# 为每个模块单独计算内部模块连通性
# uniqueColors <- unique(mergedColors)
# moduleConnectivity <- list()
# 
# for (color in uniqueColors) {
#   moduleGenes <- which(mergedColors == color)
#   moduleExpr <- datExpr_filted[, moduleGenes]
#   
#   # 计算模块内的加权邻接矩阵
#   kbe_module <- abs(cor(moduleExpr, use = "p"))^softPower
#   
#   # 计算模块内的连通性
#   connectivity <- rowSums(kbe_module) - 1
#   moduleConnectivity[[color]] <- connectivity
# }


#输出
#write.table(gene_select, 'gene_select.txt', sep = '\t', col.names = NA, quote = FALSE)
#write.table(tom_select, 'tom_select.txt', sep = '\t', col.names = NA, quote = FALSE)


#选择 black 模块内的基因
gene_black <- datExpr_filted[ ,gene_module_select]

#基因的模块成员度（module membership）计算
#即各基因表达值与相应模块特征基因的相关性，其衡量了基因在全局网络中的位置
geneModuleMembership <- signedKME(gene_black, module['MEblack'], outputColumnName = 'MM')
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nrow(module)))

#各基因表达值与临床表型的相关性分析
# 检查基因表达数据中的变异度

geneTraitSignificance <- as.data.frame(cor(gene_black, trait['Subtype III'], use = 'p'))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nrow(trait)))
GeneSignificance <- abs(geneTraitSignificance)
# 计算 black 模块中的模块成员系数（MM）

plot(geneModuleMembership[,1],GeneSignificance[,1] , 
     xlab = "Module Membership in black module", 
     ylab = "Gene significance for trait",
     main = "Module membership vs. gene significance\n",
     pch = 21, col = "black", bg = "white")
#选择显著（p<0.01）、高 black 模块成员度（MM>=0.8），与 TNM 表型高度相关（r>=0.8）的基因
geneModuleMembership[geneModuleMembership<0.5 | MMPvalue>0.01] <- 0
geneTraitSignificance[abs(geneTraitSignificance)<0.1 | GSPvalue>0.01] <- 0

select <- cbind(geneModuleMembership, geneTraitSignificance)


select0 <- subset(select, geneModuleMembership>=0.8 & geneTraitSignificance>=0.45)

# [1] "FOXN3"         "NR3C1"         "PLEKHO1"       "ATP8B2"        "TTC28"        
# [6] "GYPC"          "C20orf194"     "CLIP4"         "MAGI2-AS3"     "RBMS3"        
# [11] "CNRIP1"        "PDE1A"         "SLC9A9"        "NR2F2-AS1"     "RP11-875O11.1"
# [16] "RP11-730A19.9"
kwithin00 <- kwithin0[rownames(select0), , drop = FALSE]  # 取df2的第一列
select0$kwithin <- kwithin00[,1] 
write.csv(select0,"sub1_hub.csv")
plotNetworkHeatmap(datExpr_filted,
                   plotGenes = rownames(select0),
                   networkType = "unsigned",
                   useTOM = TRUE,
                   power=powers,
                   main="unsigned correlations")


select1 <- subset(select, geneModuleMembership>=0.8 & abs(geneTraitSignificance)>=0.1)
rownames(select1)
# [1] "MALAT1"        "NPIPB5"        "RP11-166B2.3"  "RP11-416A17.6" "CTD-2014D20.1"
# [6] "LA16c-431H6.6" "RYKP1"         "RP11-49O14.2"  "RP11-192H23.7"

#0.8
# [1] "SGK494"        "AC010976.2"    "ANKRD61"       "GTF3C2-AS1"    "DDX50P1"      
# [6] "RP11-416A17.6" "RP11-29H23.7"  "AL590762.11"  
kwithin11 <- kwithin1[rownames(select1), , drop = FALSE]  # 取df2的第一列
select1$kwithin <- kwithin11[,1] 
write.csv(select1,"0.7&0.2的sub2.csv")
write.csv(select1,"0.8&0.1的sub2.csv")

plotNetworkHeatmap(datExpr_filted,
                   plotGenes = rownames(select1),
                   networkType = "unsigned",
                   useTOM = TRUE,
                   power=powers,
                   main="unsigned correlations")

select2 <- subset(select, geneModuleMembership>=0.8 & abs(geneTraitSignificance)>=0.4)
# [1] "SPARC" "BGN"   "SULF1" "THY1"  "CDH11" "PRRX1" "FAP"   "NOX4" 
kwithin22 <- kwithin2[rownames(select2), , drop = FALSE]  # 取df2的第一列
select2$kwithin <- kwithin22[,1] 
write.csv(select2,"hub3.csv")

plotNetworkHeatmap(datExpr_filted,
                   plotGenes = rownames(select2),
                   networkType = "unsigned",
                   useTOM = TRUE,
                   power=powers,
                   main="unsigned correlations")







dir.create('cytoscape0705', recursive = TRUE)

for (mod in 1:nrow(table(mergedColors))) {
  modules <- names(table(mergedColors))[mod]
  probes <- colnames(datExpr_filted)
  inModule <- (mergedColors == modules)
  modProbes <- probes[inModule]
  modGenes <- modProbes
  modtom_sim <- tom_sim[inModule, inModule]
  dimnames(modtom_sim) <- list(modProbes, modProbes)
  outEdge <- paste('cytoscape0705/', modules, '.edge_list.txt',sep = '')
  outNode <- paste('cytoscape0705/', modules, '.node_list.txt', sep = '')
  exportNetworkToCytoscape(modtom_sim,
                           edgeFile = outEdge,
                           nodeFile = outNode,
                           weighted = TRUE,
                           threshold = 0.3,  #该参数可控制输出的边数量，过滤低权重的边
                           nodeNames = modProbes,
                           altNodeNames = modGenes,
                           nodeAttr = mergedColors[inModule])
}
