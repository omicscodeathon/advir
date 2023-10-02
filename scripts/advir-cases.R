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
malecase <- read.csv("malecases.csv",header = T)
# malecase <- cases
# malecase[sapply(malecase, is.character)] <- lapply(malecase[sapply(malecase, is.character)], as.numeric)
# malecase <- as.data.frame(malecase)
# malecase$Taxon <- as.character(cases$Taxon)
# malecase[is.na(malecase)] <- 0
colnames(malecase)

femcase <- read.csv("trcasedata.csv")#[,1:24]
# femcase <- cases
# femcase[sapply(femcase, is.character)] <- lapply(femcase[sapply(femcase, is.character)], as.numeric)
# femcase <- as.data.frame(femcase)
# femcase$Taxon <- as.character(cases$Taxon)
# femcase[is.na(femcase)] <- 0
colnames(femcase)

mergeddatam<-merge(x=malecase,y=femcase, by="Taxon",all=T)
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
metadatam <- read.csv("casemetadata.csv",row.names = 1)
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
  labs(title="Alpha Diversity Cases",x="Gender")
  

# Calculate distances #unifrac, wunifrac, euclidean
DistBC = phyloseq::distance(ps_objectm, method = "jaccard")
ordBC = ordinate(ps_objectm, method = "PCoA", distance = DistBC)
plot_scree(ordBC, "Scree Plot: Jaccard MDS")
plot_ordination(ps_objectm, ordBC, color = "gender") +
  geom_point(aes(color=gender)) +
  ggtitle("PCoA: Jaccard Distance Cases") + labs(color="Gender")

DistB = phyloseq::distance(ps_objectm, method = "bray")
ordB = ordinate(ps_objectm, method = "PCoA", distance = DistB)
plot_scree(ordB, "Scree Plot: Bray-Cutis MDS")
plot_ordination(ps_objectm, ordB, color = "gender") +
  geom_point(aes(color=gender)) +
  ggtitle("PCoA: Bray-Cutis Distance Cases") + labs(color="Gender")
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

malecase_vs_femcase = results(diagddsm,#filterFun = ihw,
                        independentFiltering = FALSE,
                        test = "Wald",
                        name = "gender_male_vs_female")
summary(malecase_vs_femcase)

malecase_vs_femcase<-malecase_vs_femcase[order(malecase_vs_femcase$padj),]
malecase_vs_femcase
malecase_vs_femcase_sig<-subset(malecase_vs_femcase, malecase_vs_femcase$padj<0.05)
malecase_vs_femcase_sig


rpmm<-fpm(diagddsm)
rpmm

#heatmap(rpm[rownames(malecase_vs_femcase_sig),])
Heatmap(t(scale(t(log2(rpmm[rownames(malecase_vs_femcase_sig),]+2)))), column_split = metadatam$gender)

#BiocManager::install("ComplexHeatmap")
malecase_vs_femcase.df <- as.data.frame(malecase_vs_femcase)
EnhancedVolcano(malecase_vs_femcase.df,x="log2FoldChange",y="padj",lab = row.names(malecase_vs_femcase.df),
                FCcutoff = 2,pCutoffCol = "padj", pCutoff = 0.05)
malecase_vs_femcase.df$sig <- ifelse(malecase_vs_femcase.df$padj<0.05,"yes","no")
ggplot(malecase_vs_femcase.df, aes(log2FoldChange, -log10(padj), color=sig)) + 
  geom_point()
