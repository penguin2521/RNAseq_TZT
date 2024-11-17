

rm(list=ls())
gc()
options(stringsAsFactors = F)

## tpm norm de            ###########################################
library(tidyverse)
load("E:/database/database/GTF_annotation/grch37_annotation_gtf75.RData")
tpm <- read.delim("C:/Users/WIN 10/Desktop/rsem_tpm_exp.txt", stringsAsFactors=FALSE)
table(tpm$X %in% gtf75$gene_id)
colnames(tpm)=gsub("rsem_exp.|_rsem.genes.results","",colnames(tpm))
tpm=merge(gtf75[,-2],tpm,by.x="gene_id",by.y="X")
tpm=limma::avereps(tpm[,3:ncol(tpm)],ID=tpm$gene_name) %>% as.data.frame()
boxplot(tpm)
tpm=log2(tpm+1)
save(tpm,file = "tpm.RData")

library(limma)
library(edgeR)
load("condition.RData")
mm=model.matrix(~0 + condition$Group)
colnames(mm)=gsub("condition[$]Group","",colnames(mm))
table(rowSums(tpm>0.1)>0)
res=lmFit(tpm[rowSums(tpm>0.1)>0,rownames(condition)],mm) %>% contrasts.fit(., makeContrasts(Case - Control, levels = mm)) %>% eBayes(.,trend = T) %>% topTable(., number =Inf,adjust.method = "BH") %>% .[order(.$adj.P.Val),]  %>% na.omit(.)
table(res$adj.P.Val<0.05)

save(res,file = "tpm_limma_res.RData")
write.table(res,"tpm_limma_res.tsv",quote = F,sep = "\t")


################################### KEGG analysis ######################################################################
#引用包
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(ComplexHeatmap)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")library(clusterProfiler)


load("tpm_limma_res.RData")
deg=res[!is.na(res$P.Value) & res$P.Value<0.05,]
write.table(deg,"deg_limma_res.tsv",quote = F,sep = "\t")

de_gene=rownames(res[!is.na(res$P.Value) & res$P.Value<0.05,])
ensem=bitr(de_gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
kegg=enrichKEGG(ensem$ENTREZID, organism = 'hsa', keyType = 'kegg',pAdjustMethod = "none", pvalueCutoff = 0.05,qvalueCutoff = 0.2)
kk=as.data.frame(kegg)

dotplot(kegg, showCategory=10) 

## ploting
#if(T){#
x = kegg
## 
df = data.frame(x)
dd =kk#@result
## 
dd$richFactor =dd$Count / as.numeric(sub("/\\d+", "", dd$BgRatio))

dd <- dd[dd$pvalue < 0.05,]
library(dplyr)
dd=dplyr::arrange(dd,-pvalue)
write.table(dd, file="ddKEGG.txt", sep="\t", quote=F, row.names = F)
dd<-read.table("KEGG.txt", header=T, sep="\t", check.names=F)  #Screening of tumor-related pathways
showNum=20     #show the top 20 pathway
if(nrow(dd)>showNum){
  dd=dd[c(1:showNum),]
}

pdf(file="KEGG_point2.pdf", width=6, height=5)


ggplot(dd,aes(richFactor,forcats::fct_reorder(Description, richFactor))) + 
  ## line
  geom_segment(aes(xend=0, yend = Description)) +
  ## point
  geom_point(aes(color=p.adjust, size = Count)) +
  ## 
  scale_color_viridis_c(begin = 0.1, end = 1) +
  ## buble
  scale_size_continuous(range=c(2, 10)) +
  theme_bw() + 
  xlab("Rich factor") +
  ylab(NULL) + 
  ggtitle("")
dev.off()


## GSEA                   ####################################################
library(msigdf)
library(dplyr)
library(clusterProfiler)
load("tpm_limma_res.RData")
genelist=res$logFC
names(genelist) = rownames(res)
genelist = sort(genelist, decreasing = TRUE)

all_geneset <- msigdf.human %>% as.data.frame()
all_geneset=all_geneset[grepl("^KEGG_",all_geneset$geneset),3:4]

load("E:/proj/LGG/death_pathway3/gene_list.RData")
colnames(gene_list)=c("geneset","symbol")
all_geneset=rbind(all_geneset,gene_list)

res_death_gsea <- GSEA(genelist,TERM2GENE =all_geneset,verbose=FALSE, pvalueCutoff = 0.05,nPerm = 10000,pAdjustMethod = "BH")
#save(res_death_gsea,file = "res_death_gsea.RData")
res_death_gsea1=as.data.frame(res_death_gsea)
#write.table(res_death_gsea1,"res_death_gsea.tsv",quote = F,sep = "\t",row.names = F)
unique(res_death_gsea1$ID)

library(enrichplot)
library(DOSE)
load("res_death_gsea.RData")
pdf("gsea.pdf",16,10)
for (i in unique(res_death_gsea1$ID)) {
    p=gseaplot2(res_death_gsea,i,color = "red",pvalue_table = T)
    print(p)
}
dev.off()


###GSEA分析
#引用包
library(limma)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

load(file = "res_death_gsea.RData")
data=as.data.frame(res_death_gsea)
#termNum=5    
#showTerm="Ferroptosis"#row.names(data)[1:termNum]     
gseaplot=gseaplot2(res_death_gsea,
                   "Ferroptosis", #pathway              
                   base_size=14, 
                   title="Enrichment plot: Ferroptosis ",  #title
                   pvalue_table = T,  #FDR value
                   color = "red",   
                   ES_geom = "line"
)
pdf(file="GSEA_Ferroptosis.pdf", width=6.5, height=5)
print(gseaplot)
dev.off()


#####################15 DEGs in ferroptosis pathway#######################################################################
load(file = "tpm.RData")

colnames(tpm)=c("shEZH2_1_1","shEZH2_1_2","shEZH2_1_3","shEZH2_3_1","sheEZH2_3_2","shEZH2_3_3","shluc_1","shluc_2","shluc_3")
exp=tpm[same,]
write.table(exp,file="Ferroptosisgene.txt",sep="\t",quote=F,row.names = T)

row.names(exp)=sigVec
Type=c(rep("shEZH2",6),rep("shluc",3))
names(Type)=colnames(tpm)
Type=as.data.frame(Type)
#pdf(file="RNAseqheatmap.pdf", width=8, height=4)
pheatmap(exp, 
         annotation=Type, 
         color = colorRampPalette(c(rep("navy",5), "white", rep("firebrick3",5)))(50),
         cluster_cols =F,
         cluster_rows =T,
         scale="row",
         show_colnames = F,
         show_rownames = T,
         fontsize = 8,
         fontsize_row=8,
         fontsize_col=8)
graph2ppt(file="RNAseqheatmap.pptx")
dev.off()