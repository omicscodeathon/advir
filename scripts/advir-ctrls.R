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
malectrl <- read.csv("malecontrols.csv",header = T)
# malectrl <- cases
# malectrl[sapply(malectrl, is.character)] <- lapply(malectrl[sapply(malectrl, is.character)], as.numeric)
# malectrl <- as.data.frame(malectrl)
# malectrl$Taxon <- as.character(cases$Taxon)
# malectrl[is.na(malectrl)] <- 0
colnames(malectrl)

femctrl <- read.csv("trcontroldata.csv")#[,1:24]
# femctrl <- ctrls
# femctrl[sapply(femctrl, is.character)] <- lapply(femctrl[sapply(femctrl, is.character)], as.numeric)
# femctrl <- as.data.frame(femctrl)
# femctrl$Taxon <- as.character(ctrls$Taxon)
# femctrl[is.na(femctrl)] <- 0
colnames(femctrl)

mergeddatam<-merge(x=malectrl,y=femctrl, by="Taxon",all=T)
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
metadatam <- read.csv("ctrlmetadata.csv",row.names = 1)
metadatam

all(rownames(metadatam)%in%colnames(mergeddatam))
all(rownames(metadatam)==colnames(mergeddatam))
#metadatam<-metadatam[colnames(mergeddatam),drop=F,]
#all(rownames(metadatam)==colnames(mergeddatam))

otu1m<-otu_table(mergeddatam, taxa_are_rows = T)
samplem<-phyloseq::sample_data(metadatam)
ps_objectm<-phyloseq(otu1m, samplem)
#ps_object1 <- ps_objectm
#ps_objectm <- transform_sample_counts(ps_object1, as.integer)

alpha_meas = c("Shannon", "Simpson")
#?plot_richness "Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"
(p <- plot_richness(ps_objectm, "gender",measures=alpha_meas))
p$layers <- p$layers[-1]

#p <- ggplot(as.data.frame(ps_objectm))
p + geom_boxplot(data=p$data, aes(x=gender, y=value, fill=as.factor(gender)), 
                 alpha=0.2, show.legend = F, outlier.shape = NA) +
  theme(axis.text.x = element_text(angle = 0, hjust=0.5)) + theme_bw() + 
  geom_jitter(width = 0.2) + stat_compare_means(method="aov", 
  label.y.npc = "top",label.x.npc = "left") + 
  labs(title="Alpha Diversity Controls",x="Gender")
  

# Calculate distances #unifrac, wunifrac, euclidean
DistBC = phyloseq::distance(ps_objectm, method = "jaccard")
ordBC = ordinate(ps_objectm, method = "PCoA", distance = DistBC)
plot_scree(ordBC, "Scree Plot: Jaccard MDS")
plot_ordination(ps_objectm, ordBC, color = "gender") +
  geom_point(aes(color=gender)) +
  ggtitle("PCoA: Jaccard Distance Controls") + labs(color="Gender")

DistB = phyloseq::distance(ps_objectm, method = "bray")
ordB = ordinate(ps_objectm, method = "PCoA", distance = DistB)
plot_scree(ordB, "Scree Plot: Bray-Cutis MDS")
plot_ordination(ps_objectm, ordB, color = "gender") +
  geom_point(aes(color=gender)) +
  ggtitle("PCoA: Bray-Cutis Distance Controls") + labs(color="Gender")
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
diagddsm = phyloseq_to_deseq2(ps_objectm, ~ gender)
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

malectrl_vs_femctrl = results(diagddsm,#filterFun = ihw,
                        independentFiltering = FALSE,
                        test = "Wald",
                        name = "gender_male_vs_female")
summary(malectrl_vs_femctrl)

malectrl_vs_femctrl<-malectrl_vs_femctrl[order(malectrl_vs_femctrl$padj),]
malectrl_vs_femctrl
malectrl_vs_femctrl_sig<-subset(malectrl_vs_femctrl, malectrl_vs_femctrl$padj<0.05)
malectrl_vs_femctrl_sig


rpmm<-fpm(diagddsm)
rpmm

#heatmap(rpm[rownames(malectrl_vs_femctrl_sig),])
Heatmap(t(scale(t(log2(rpmm[rownames(malectrl_vs_femctrl_sig),]+2)))), column_split = metadatam$group)

#BiocManager::install("ComplexHeatmap")
malectrl_vs_femctrl.df <- as.data.frame(malectrl_vs_femctrl)
EnhancedVolcano(malectrl_vs_femctrl.df,x="log2FoldChange",y="padj",lab = row.names(malectrl_vs_femctrl.df),
                FCcutoff = 2,pCutoffCol = "padj", pCutoff = 0.05)
malectrl_vs_femctrl.df$sig <- ifelse(malectrl_vs_femctrl.df$padj<0.05,"yes","no")
ggplot(malectrl_vs_femctrl.df, aes(log2FoldChange, -log10(padj), color=sig)) + 
  geom_point()
