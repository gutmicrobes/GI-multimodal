rm(list = ls())
library(dynamicTreeCut)
library(stats)
library(fastcluster)
#biocLite("GO.db")
library(WGCNA)
library(lattice)
library(latticeExtra)
setwd('E:\\Graduation_design\\merge_data\\mRNA')


rt=read.csv('finaly_mrna.csv')

rownames(rt) = rt[,1]
rt = rt[,2:ncol(rt)]
#rt = rt[2:(nrow(rt)-1),]
#rt = rt[,2:ncol(rt)]
rt = t(rt)

rt = data.frame(rt)
rt = as.matrix(rt)
rt = apply(rt,2,as.numeric)

sampleTree = hclust(dist(rt), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")


powers = c(c(1:10), seq(from = 12, to=40, by=2))


sft = pickSoftThreshold(rt, powerVector = powers, verbose = 5)


sizeGrWindow(9, 5)
par(mfrow = c(1,2)); 
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("miRNA Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("miRNA Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

softPower <- sft$powerEstimate
softPower = 36
adjacency = adjacency(rt, power = softPower);
TOM = TOMsimilarity(adjacency);

# 计算基因之间的相异度,根据相异程度进行聚类
dissTOM = 1-TOM 
dissTOM1  = data.frame(dissTOM)
#write.csv(dissTOM1,'pig_disstom.csv')
write.csv(dissTOM1,'new_mRNA_disstom.csv')
