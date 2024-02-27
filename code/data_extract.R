rm(list = ls())
library(limma)
setwd("E:\\Graduation_design\\merge_data\\mRNA")
logFoldChange=1
adjustP=0.05

#rt=read.table("normalize.txt",sep='\t',header = T,check.names = F)

rt = read.csv('colon_stomach_mirna.csv')
coln = colnames(rt)
# for (i in 1:ncol(rt)){
#   if (coln[i] == 'hsa.miR.93.5p'){
#     print(i)
#   }
# }
#2203
# no_need = c(1106,398,2223,201,1557,2227,2203)
# rt = rt[,-no_need]
data = rt   #为后来提取癌症数据作准备
sample = rt$id
rt = subset(rt,select = -c(id))
modType = c()
for(i in 1:length(sample)){
  if(sample[i]==1){
    modType = c(modType,'tumor')
  }
  else{
    modType = c(modType,"normal")
  }
}
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimname=list(rownames(exp),colnames(exp))
rt=matrix(as.numeric(as.matrix(exp)),nrow = nrow(exp),dimnames = dimname)
#differential

#modType=modType<-c(rep("normal",32),rep('tumor',0),rep("tumor",70),rep("normal",12),rep('tumor',35),rep("normal",24),rep("tumor",27),rep("normal",38))
design<-model.matrix(~0+factor(modType))
colnames(design)<-c("con","treat")
rt1 = t(rt)
fit<-lmFit(rt1,design)
cont.matrix<-makeContrasts(treat-con,levels = design)
fit2<-contrasts.fit(fit,cont.matrix)
fit2<-eBayes(fit2)
allDiff=topTable(fit2,adjust.method = 'fdr',number = 200000)
#write.csv(allDiff,file="limmaTab_mirna.csv",quote=F,row.names = T)
#write table
diffSig<-allDiff[with(allDiff,(abs(logFC)>logFoldChange&adj.P.Val<adjustP)),]
#write.csv(diffSig,file='diff_mirna.csv',quote = F,row.names = T)
diffUp<-allDiff[with(allDiff,(logFC>logFoldChange & adj.P.Val<adjustP)),]
#write.csv(diffUp,file='up_diff_mirna.csv',quote = F,row.names = T)
diffDown<-allDiff[with(allDiff,(logFC<(-logFoldChange) & adj.P.Val<adjustP)),]
#write.csv(diffDown,file="down_mirna.csv",quote=F,row.names = T)
#volcano
pdf(file="vol_mirna.pdf")
xMax=25
yMax=max(abs(allDiff$logFC))
plot(-log10(allDiff$adj.P.Val),allDiff$logFC,xlab = "-log10(adj.P.Val)",ylab="logFC",
     main="miRNA_Valcano",xlim = c(0,xMax),ylim = c(-5000,5000),yaxs="i",pch=20,cex=0.8)
diffSub=subset(allDiff,adj.P.Val<adjustP & logFC>logFoldChange)
points(-log10(diffSub$adj.P.Val),diffSub$logFC,pch=20,col="red",cex=0.8)
diffSub=subset(allDiff,adj.P.Val<adjustP& logFC<(-logFoldChange))
points(-log10(diffSub$adj.P.Val),diffSub$logFC,pch=20,col="green",cex=0.8)
abline(h=0,lty=2,lwd=3)
dev.off()


#选出癌症样本
rownames(data) = data[,1]
data = data[,2:ncol(data)]
data_tumor = data[which(data$id == 1),]
mirna_diff = data[,which(colnames(data)%in%rownames(diffSig))]
write.csv(mirna_diff,'mirna_diff.csv')  #保存含有正常样本和差异性基因的数据
mirna_diff_tumor = data_tumor[,which(colnames(data_tumor)%in%rownames(diffSig))]
write.csv(mirna_diff_tumor,'mrna_diff_tumor.csv')    #只保存含有癌症样本的数据
