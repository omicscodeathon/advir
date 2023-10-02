library(DESeq2)
library(readxl)
library(readr)
library(ggpubr)
library(tidyr)
library(ggplot2)
library(phyloseq)
library(ComplexHeatmap)
library(EnhancedVolcano)

#cases<-read_xlsx("cases.xlsx")#[,1:26]
posm <- read.csv("malecases.csv",header = T)
# posm <- cases
# posm[sapply(posm, is.character)] <- lapply(posm[sapply(posm, is.character)], as.numeric)
# posm <- as.data.frame(posm)
# posm$Taxon <- as.character(cases$Taxon)
# posm[is.na(posm)] <- 0
colnames(posm)

negm <- read.csv("malecontrols.csv")#[,1:24]
# negm <- ctrls
# negm[sapply(negm, is.character)] <- lapply(negm[sapply(negm, is.character)], as.numeric)
# negm <- as.data.frame(negm)
# negm$Taxon <- as.character(ctrls$Taxon)
# negm[is.na(negm)] <- 0
colnames(negm)

mergeddatam<-merge(x=negm,y=posm, by="Taxon",all=T)
mergeddatam
mergeddatam[is.na(mergeddatam)]<-0
mergeddatam
#mergeddatam <- mergeddatam[c(-5,-7,-18),]
rownames(mergeddatam)<-mergeddatam$Taxon
mergeddatam<-mergeddatam[,-1]
#xx <- mergeddatam[,colSums(mergeddatam[])>0]
#mergeddatam <- xx
#yy <- mergeddatam[rowSums(mergeddatam[])>0,]
#load metadata
metadatam <- read.csv("malemetadata.csv",row.names = 1)
metadatam

all(rownames(metadatam)%in%colnames(mergeddatam))
all(rownames(metadatam)==colnames(mergeddatam))
metadatam<-metadatam[colnames(mergeddatam),drop=F,]
#all(rownames(metadatam)==colnames(mergeddatam))

otu1m<-otu_table(mergeddatam, taxa_are_rows = T)
samplem<-phyloseq::sample_data(metadatam)
ps_objectm<-phyloseq(otu1m, samplem)
#ps_object1 <- ps_objectm
#ps_objectm <- transform_sample_counts(ps_object1, as.integer)

alpha_meas = c("Shannon", "Simpson")
#?plot_richness "Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"
(p <- plot_richness(ps_objectm, "group",measures=alpha_meas))
p$layers <- p$layers[-1]

#p <- ggplot(as.data.frame(ps_objectm))
p + geom_boxplot(data=p$data, aes(x=group, y=value, fill=as.factor(group)), 
                 alpha=0.2, show.legend = F, outlier.shape = NA) +
  theme(axis.text.x = element_text(angle = 0, hjust=0.5)) + theme_bw() + 
  geom_jitter(width = 0.2) + stat_compare_means(method="aov", 
  label.y.npc = "top",label.x.npc = "left") + labs(x="Disease status")
  

# Calculate distances #unifrac, wunifrac, euclidean
DistBC = phyloseq::distance(ps_objectm, method = "jaccard")
ordBC = ordinate(ps_objectm, method = "PCoA", distance = DistBC)
plot_scree(ordBC, "Scree Plot: Jaccard MDS")
plot_ordination(ps_objectm, ordBC, color = "group") +
  geom_point(aes(color=group)) +
  ggtitle("PCoA: Jaccard Distance") + labs(color="Disease Status")

DistB = phyloseq::distance(ps_objectm, method = "bray")
ordB = ordinate(ps_objectm, method = "PCoA", distance = DistB)
plot_scree(ordB, "Scree Plot: Bray-Cutis MDS")
plot_ordination(ps_objectm, ordB, color = "group") +
  geom_point(aes(color=group)) +
  ggtitle("PCoA: Bray-Cutis Distance") + labs(color="Disease Status")
# ordBC$vectors
# ps_objectm
# 
# df<-as.data.frame(ordBC$vectors[,1:2])
# df$condition<-metadatam$group
# df
# 
# ggplot(df, aes(Axis.1, Axis.2, color=condition))+geom_point(size=2) + 
#   labs(color="Disease group", title = "PCoA plot: Euclidean Distance")#PERMANOVA

#############################################################################
diagddsm = phyloseq_to_deseq2(ps_objectm, ~ group)
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

geoMeansm = apply(counts(diagddsm), 1, gm_mean)
#counts(diagddsm) = round(counts(diagddsm)/2)
diagddsm = estimateSizeFactors(diagddsm, geoMeans = geoMeansm)
diagddsm$sizeFactor

#diagddsm <- DESeqDataSetFromMatrix()
######################################################################################
diagddsm = DESeq(diagddsm,
                fitType="local",
                test = "LRT",
                reduced = ~1)
resultsNames(diagddsm)

posm_vs_negm = results(diagddsm,#filterFun = ihw,
                        independentFiltering = FALSE,
                        test = "Wald",
                        name = "group_controls_vs_cases")
summary(posm_vs_negm)

posm_vs_negm<-posm_vs_negm[order(posm_vs_negm$padj),]
posm_vs_negm
posm_vs_negm_sig<-subset(posm_vs_negm, posm_vs_negm$padj<0.05)
posm_vs_negm_sig


rpmm<-fpm(diagddsm)
rpmm

#heatmap(rpm[rownames(posm_vs_negm_sig),])
Heatmap(t(scale(t(log2(rpmm[rownames(posm_vs_negm_sig),]+2)))), column_split = metadatam$group)

#BiocManager::install("ComplexHeatmap")
posm_vs_negm.df <- as.data.frame(posm_vs_negm)
EnhancedVolcano(posm_vs_negm.df,x="log2FoldChange",y="padj",lab = row.names(posm_vs_negm.df),
                FCcutoff = 2,pCutoffCol = "padj", pCutoff = 0.05)
posm_vs_negm.df$sig <- ifelse(posm_vs_negm.df$padj<0.05,"yes","no")
ggplot(posm_vs_negm.df, aes(log2FoldChange, -log10(padj), color=sig)) + 
  geom_point()
