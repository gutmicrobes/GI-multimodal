rm(list = ls())
setwd("E:/Graduation_design/WGCNA")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("AnnotationDbi")
# BiocManager::install("impute")
# BiocManager::install("GO.db")
# BiocManager::install("preprocessCore")
# site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
library(dynamicTreeCut)
library(stats)
library(fastcluster)
#biocLite("GO.db")
library(WGCNA)
library(lattice)
library(latticeExtra)
library(ggplot2)

#读入基因表达谱，基因在行，样本在列
rt=read.csv("finaly_mrna.csv",header = T,check.names = F)
rownames(rt) = rt[,1]
rt= rt[,2:ncol(rt)]
#rt = as.numeric(unlist(rt))
exp=as.matrix(rt)

dimname=list(rownames(exp),colnames(exp))
rt=matrix(as.numeric(as.matrix(exp)),nrow = nrow(exp),dimnames = dimname)
rt <- as.data.frame(rt) #转换为样品在行，基因在列的矩阵

#聚类查看是否有离群样品
sampleTree = hclust(dist(rt), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")


#生成包含不同阈值的数列
powers = c(c(1:10), seq(from = 12, to=20, by=2))
power
# 遍历上述给定的阈值，找到符合无标度网络准则的power
sft = pickSoftThreshold(rt, powerVector = powers, verbose = 5)

#将上述遍历计算结果可视化
sizeGrWindow(9, 5)
par(mfrow = c(1,2)); #一页两图
cex1 = 0.8;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
#横轴是Soft threshold (power)，纵轴是无标度网络的评估参数，数值越高，网络越符合无标度特征 (non-scale)
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# 计算基因之间的加权相关系数
softPower <- sft$powerEstimate #得到power的最优值
#softpower = 6
adjacency = adjacency(rt, power = softPower);

# 计算得到拓扑矩阵
TOM = TOMsimilarity(adjacency);

# 计算基因之间的相异度,根据相异程度进行聚类
dissTOM = 1-TOM       
hierTOM = hclust(as.dist(dissTOM),method="average");
#write.csv(dissTOM,"./dissTOM.csv")

# 检验选定的β值下记忆网络是否逼近 scale free
k <- softConnectivity(datE=rt,power=softPower) 
sizeGrWindow(10, 5)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k,main="Check Scale free topology\n")



geneTree = hclust(as.dist(dissTOM), method = "average");
#memory.limit(1000000)
# Plot the resulting clustering tree (dendrogram)
windows()
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.01);

minModuleSize = 30;

dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)


# 鎷撴墤鐑浘锛?
nSelect = 400 
# For reproducibility, we set the random seed 
set.seed(10); 
select = sample(1167, size = nSelect)
selectTOM = dissTOM[select, select]
dynamicColors=labels2colors(dynamicMods)
# There's no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster. 
selectTree = hclust(as.dist(selectTOM), method = "average") 
selectColors = dynamicColors[select]; 
# Open a graphical window 
sizeGrWindow(9,9) 
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing 
# the color palette; setting the diagonal to NA also improves the clarity of the plot 
plotDiss = selectTOM^softPower; 
diag(plotDiss) = NA; 
TOMplot(plotDiss, 
        selectTree, 
        selectColors, 
        main = "Network heatmap plot, selected genes") 

#计算每个模块的特征向量基因，为某一特定模块第一主成分基因E。代表了该模块内基因表达的整体水平
dynamicColors=labels2colors(dynamicMods)
MEList = moduleEigengenes(rt, colors = dynamicColors)
MEs = MEList$eigengenes

#根据模块特征向量基因计算模块相异度：
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result 绘制15个模块特征向量的相关系数热图
plotEigengeneNetworks(MEs, 
                      "Eigengene adjacency heatmap", 
                      marHeatmap = c(3,4,2,2), 
                      plotDendrograms = FALSE, 
                      xLabelsAngle = 90) 

#特征向量基因聚类树状图，红线以下的模块表示相关性>0.8，将被合并
plot(METree, 
     main = "Clustering of module eigengenes",
     xlab = "", 
     sub = "",
     col="purple",
     col.axis = 'green3',
     col.lab = 'purple',
     col.main = 'magenta',
     lty = 1,
     pch = 1,
     #type = 1
     )
# 在聚类图中画出剪切线
MEDissThres = 0.2
abline(h=MEDissThres, col = "red")



# 将相关性系数大于0.8的模块合并掉，即相异性系数小于0.2:最后得到10个模块
merge_modules = mergeCloseModules(rt, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# 合并后的颜色：
mergedColors = merge_modules$colors;
# 新模块的特征向量基因：
mergedMEs = merge_modules$newMEs;
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,cex.lab=1.2,cex.colorLabels = 1.2)

#读入样本性状信息，0代表健康，1代表患病
dataTraits <- read.csv("finally_label.csv",header=TRUE,row.names = 1)
x = c()
for (i in 1:515){
  if (dataTraits[i,] != 2){
    x = c(x,0)
  }
  else {
    x = c(x,1)
  }
}

dataTraits=as.matrix(x)
ncol(dataTraits)
nrow(dataTraits)
ncol(mergedMEs)
nrow(mergedMEs)
#计算性状与特征基因的相关系数以及p值
cor_ADR <- signif(WGCNA::cor(dataTraits,mergedMEs,use="p",method="pearson"),5)
p.values <- corPvalueStudent(cor_ADR,nSamples=515)
#根据基因网络显著性，也就是性状与每个基因表达量相关性在各个模块的均值作为该性状在该模块的显著性，显著性最大的那个模块与该性状最相关：
GS1 <- as.numeric(WGCNA::cor(dataTraits,rt,use="p",method="pearson"))
#显著性是绝对值：
GeneSignificance <- abs(GS1)
write.csv(GeneSignificance,"./third_GS.csv")

# 获得该性状在每个模块中的显著性：
ModuleSignificance <- tapply(GeneSignificance,mergedColors,mean,na.rm=T)
Find_max_ModuleSign <- which.max(ModuleSignificance)
#结果与相关系数和p值所选结果相同



#根据上述计算结果，挑选与性状相关性最高最显著的模块所包含的基因
module = "green";
probes = colnames(rt)
inModule = (mergedColors==module);
modProbes = probes[inModule]; 

modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
write.csv(modTOM,"./third_modTOM.csv")
#将上述基因导出到cytoscape
cyt = exportNetworkToCytoscape(
  modTOM,
  edgeFile = paste("CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
  nodeFile = paste("CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
  weighted = TRUE,
  threshold = 0.02,
  nodeNames = modProbes, 
  nodeAttr = mergedColors[inModule]
);

rtmod <- rt[inModule, inModule];
colorh1 = mergedColors[inModule]


ADJ1=abs(cor(rtmod,use="p"))^softPower 
Alldegrees1=intramodularConnectivity(ADJ1, colorh1) 
write.csv(Alldegrees1,"./third_Alldegrees.csv")

rtKME=signedKME(rt, mergedMEs, outputColumnName="MM.")
head(rtKME)
write.csv(rtKME,"./third_rtKME.csv")

FilterGenes_spe = ((GeneSignificance > 0.1) & (abs(rtKME["MM.green"])>0.5)) 

table(FilterGenes_spe)
trait_hubGenes_spe <- colnames(rt)[FilterGenes_spe] 
write.csv(trait_hubGenes_spe,"./third_genes_spe.csv")
head(trait_hubGenes_spe)
trait_hubGenes_spe
write.csv(trait_hubGenes_spe,"third_cluster_hubGenes_spe.csv")


