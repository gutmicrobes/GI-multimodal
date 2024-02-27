rm(list = ls())
library(AnnotationDbi)
library(org.Hs.eg.db)#����ע�Ͱ�
library(clusterProfiler)#������
library(dplyr)
library(ggplot2)#��ͼ��
gene = c('ASCL2','PITX2','LSM7','LEFTY1','CLDN3','HES2','ATP2C2','LINC01558')
#gene = c('UTP14A','GTPBP4','ACTL6A','CBX3','DCAF13','RPP40','TFB2M','PAK1IP1','HSPE1','PNPT1','DKC1','NUDCD1','SSB')
#gene = c('C1orf147','AL031600.1','C2.AS1','NPIPB14P','AC073957.3')

keytypes(org.Hs.eg.db)
gene.df <- bitr(gene,fromType="SYMBOL",toType="ENTREZID", OrgDb = org.Hs.eg.db)#TCGA���ݿ����û�н��л���ע�ͣ���ôfromTypeӦ����Ensembl������ID֮����Ի���ת��,toType������һ���ַ�����Ҳ������һ�����������Լ����� 
          

ego_ALL <- enrichGO(gene = gene.df$ENTREZID,#�������涨����
                    OrgDb=org.Hs.eg.db,
                    keyType = "ENTREZID",
                    ont = "ALL",#������GO����
                    pAdjustMethod = "BH",#������ùܣ�һ�㶼�õ�BH
                    minGSSize = 1,
                    pvalueCutoff = 0.01,#Pֵ����ȡ0.05
                    qvalueCutoff = 0.05,
                    readable = TRUE)

ego_CC <- enrichGO(gene = gene.df$ENTREZID,
                   OrgDb=org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "CC",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05,
                   readable = TRUE)

ego_BP <- enrichGO(gene = gene.df$ENTREZID,
                   OrgDb=org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "BP",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05,
                   readable = TRUE)

ego_MF <- enrichGO(gene = gene.df$ENTREZID,
                   OrgDb=org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "MF",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05,
                   readable = TRUE)

#4����������浽��ǰ·��
ego_ALL <- as.data.frame(ego_ALL)
ego_result_BP <- as.data.frame(ego_BP)
ego_result_CC <- as.data.frame(ego_CC)
ego_result_MF <- as.data.frame(ego_MF)
ego <- rbind(ego_result_BP,ego_result_CC,ego_result_MF)#��������Ҳ�ܵõ�ego_ALLһ���Ľ��
# write.csv(ego_ALL,file = "ego_ALL.csv",row.names = T)
# write.csv(ego_result_BP,file = "ego_result_BP.csv",row.names = T)
# write.csv(ego_result_CC,file = "ego_result_CC.csv",row.names = T)
# write.csv(ego_result_MF,file = "ego_result_MF.csv",row.names = T)
# write.csv(ego,file = "ego.csv",row.names = T)

#5��������ʱ�����Ǹ�������ͨ·̫���ˣ����������ͼƬ������ѡ�񲿷ֻ��ƣ���ʱ��������4����������5��
display_number = c(30,30, 30)#���������ֱַ����ѡȡ��BP��CC��MF��ͨ·����������Լ����þ�����
ego_result_BP <- as.data.frame(ego_BP)[1:display_number[1], ]
ego_result_CC <- as.data.frame(ego_CC)[1:display_number[2], ]
ego_result_MF <- as.data.frame(ego_MF)[1:display_number[3], ]

##����������ժȡ�Ĳ���ͨ·������ϳ����ݿ�
go_enrich_df <- data.frame(
  ID=c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID),                         Description=c(ego_result_BP$Description,ego_result_CC$Description,ego_result_MF$Description),
  GeneNumber=c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count),
  type=factor(c(rep("biological process", display_number[1]), 
                rep("cellular component", display_number[2]),
                rep("molecular function", display_number[3])), 
              levels=c("biological process", "cellular component","molecular function" )))

##ͨ·������̫���ˣ���ѡȡ��ͨ·��ǰ���������Ϊͨ·������
# for(i in 1:nrow(go_enrich_df)){
#   description_splite=strsplit(go_enrich_df$Description[i],split = " ")
#   description_collapse=paste(description_splite[[1]][1:5],collapse = " ") #�����5����ָ5�����ʵ���˼�������Լ�����
#   go_enrich_df$Description[i]=description_collapse
#   go_enrich_df$Description=gsub(pattern = "NA","",go_enrich_df$Description)
# }

##��ʼ����GO��״ͼ
###���ŵ���״ͼ
go_enrich_df$type_order=factor(rev(as.integer(rownames(go_enrich_df))),labels=rev(go_enrich_df$Description))#��һ���Ǳ���ģ�Ϊ�������Ӱ�˳����ʾ�������ں���
COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")#�趨��ɫ

ggplot(data=go_enrich_df, aes(x=type_order,y=GeneNumber, fill=type)) + #������ȡֵ
  geom_bar(stat="identity", width=0.8) + #��״ͼ�Ŀ��ȣ������Լ�����
  scale_fill_manual(values = COLS) + ###��ɫ
  coord_flip() + ##��һ��������״ͼ����������ӵĻ���״ͼ�����ŵ�
  xlab("GO term") + 
  ylab("Gene_Number") + 
  labs(title = "The Most Enriched GO Terms")+
  theme_bw()

###���ŵ���״ͼ 
go_enrich_df$type_order=factor(rev(as.integer(rownames(go_enrich_df))),labels=rev(go_enrich_df$Description))
COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")
ggplot(data=go_enrich_df, aes(x=type_order,y=GeneNumber, fill=type)) + 
  geom_bar(stat="identity", width=0.8) + 
  scale_fill_manual(values = COLS) + 
  theme_bw() + 
  xlab("GO term") + 
  ylab("Num of Genes") + 
  labs(title = "The Most Enriched GO Terms")+ 
  theme(axis.text.x=element_text(face = "bold", color="gray50",angle = 70,vjust = 1, hjust = 1 ))#angle��������������б�ĽǶȣ������Լ�����

#KEGG
#1��KEGG����
kk <- enrichKEGG(gene = gene.df$ENTREZID,keyType = "kegg",organism= "human", qvalueCutoff = 0.1, pvalueCutoff=0.1)

#2�����ӻ�
###��״ͼ
hh <- as.data.frame(kk)#�Լ��ǵñ���������
rownames(hh) <- 1:nrow(hh)
hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Description))
ggplot(hh,aes(y=order,x=Count,fill=p.adjust))+
  geom_bar(stat = "identity",width=0.7)+####���ӿ���
  #coord_flip()+##�ߵ�������
  scale_fill_gradient(low = "red",high ="blue" )+#��ɫ�Լ����Ի�
  labs(title = "KEGG Pathways Enrichment",
       x = "Gene numbers", 
       y = "Pathways")+
  theme(axis.title.x = element_text(face = "bold",size = 16),
        axis.title.y = element_text(face = "bold",size = 16),
        legend.title = element_text(face = "bold",size = 16))+
  theme_bw()

###����ͼ
hh <- as.data.frame(kk)
rownames(hh) <- 1:nrow(hh)
hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Description))
ggplot(hh,aes(y=order,x=Count))+
  geom_point(aes(size=Count,color=-1*p.adjust))+# �޸ĵ�Ĵ�С
  scale_color_gradient(low="green",high = "red")+
  labs(color=expression(p.adjust,size="Count"), 
       x="Gene Number",y="Pathways",title="KEGG Pathway Enrichment")+
  theme_bw()