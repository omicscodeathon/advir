library(DESeq2)
library(readxl)
library(readr)
library(ggplot2)
#cases<-read_xlsx("cases.xlsx")#[,1:26]
pos <- read.csv("casedata.csv",header = T)
# pos <- cases
# pos[sapply(pos, is.character)] <- lapply(pos[sapply(pos, is.character)], as.numeric)
# pos <- as.data.frame(pos)
# pos$Taxon <- as.character(cases$Taxon)
# pos[is.na(pos)] <- 0
colnames(pos)

neg <- read.csv("controldata.csv")#[,1:24]
# neg <- ctrls
# neg[sapply(neg, is.character)] <- lapply(neg[sapply(neg, is.character)], as.numeric)
# neg <- as.data.frame(neg)
# neg$Taxon <- as.character(ctrls$Taxon)
# neg[is.na(neg)] <- 0
colnames(neg)

mergeddata<-merge(x=neg,y=pos, by="Taxon",all=T)
mergeddata
mergeddata[is.na(mergeddata)]<-0
mergeddata
#mergeddata <- mergeddata[c(-5,-7,-18),]
rownames(mergeddata)<-mergeddata$Taxon
mergeddata<-mergeddata[,-1]
#xx <- mergeddata[,colSums(mergeddata[])>0]
#mergeddata <- xx
#yy <- mergeddata[rowSums(mergeddata[])>0,]
#load metadata
metadata <- read.csv("metadata.csv",row.names = 1)
metadata

all(rownames(metadata)%in%colnames(mergeddata))
all(rownames(metadata)==colnames(mergeddata))
metadata<-metadata[colnames(mergeddata),drop=F,]
all(rownames(metadata)==colnames(mergeddata))

library(phyloseq)
otu1<-otu_table(mergeddata, taxa_are_rows = T)
sample<-phyloseq::sample_data(metadata)
ps_object<-phyloseq(otu1, sample)
#ps_object1 <- ps_object
#ps_object <- transform_sample_counts(ps_object1, as.integer)

alpha_meas = c("Shannon", "Observed")
#?plot_richness "Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"
(p <- plot_richness(ps_object, "status",measures=alpha_meas))
p$layers <- p$layers[-1]
library(ggpubr)
library(tidyr)
library(ggplot2)
#p <- ggplot(as.data.frame(ps_object))
p + geom_boxplot(data=p$data, aes(x=status, y=value, fill=as.factor(status)), 
                 alpha=0.2, show.legend = F, outlier.shape = NA) +
  theme(axis.text.x = element_text(angle = 0, hjust=0.5)) + theme_bw() + 
  geom_jitter(width = 0.2) + stat_compare_means(method="anova", 
  label.y.npc = "top",label.x.npc = "left") + labs(x="Disease status")
  

# Calculate distances #unifrac, wunifrac, euclidean
DistBC = phyloseq::distance(ps_object, method = "jaccard")
ordBC = ordinate(ps_object, method = "PCoA", distance = DistBC)
plot_scree(ordBC, "Scree Plot: Jaccard MDS")
plot_ordination(ps_object, ordBC, color = "status") +
  geom_point(aes(color=status)) +
  ggtitle("PCoA: Jaccard Distance") + labs(color="Disease Status")

DistB = phyloseq::distance(ps_object, method = "bray")
ordB = ordinate(ps_object, method = "PCoA", distance = DistB)
plot_scree(ordB, "Scree Plot: Bray-Cutis MDS")
plot_ordination(ps_object, ordB, color = "status") +
  geom_point(aes(color=status)) +
  ggtitle("PCoA: Bray-Cutis Distance") + labs(color="Disease Status")
# ordBC$vectors
# ps_object
# 
# df<-as.data.frame(ordBC$vectors[,1:2])
# df$condition<-metadata$status
# df
# 
# ggplot(df, aes(Axis.1, Axis.2, color=condition))+geom_point(size=2) + 
#   labs(color="Disease status", title = "PCoA plot: Euclidean Distance")#PERMANOVA

#############################################################################
diagdds = phyloseq_to_deseq2(ps_object, ~ status)
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
library(DESeq2)
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds$sizeFactor



######################################################################################
diagdds = DESeq(diagdds,
                fitType="local",
                test = "LRT",
                reduced = ~1)
resultsNames(diagdds)

pos_vs_neg = results(diagdds,#filterFun = ihw,
                        independentFiltering = FALSE,
                        test = "Wald",
                        name = "status_controls_vs_cases")
summary(pos_vs_neg)

pos_vs_neg<-pos_vs_neg[order(pos_vs_neg$padj),]
pos_vs_neg
pos_vs_neg_sig<-subset(pos_vs_neg, pos_vs_neg$padj<0.05)
pos_vs_neg_sig


rpm<-fpm(diagdds)
rpm

#heatmap(rpm[rownames(pos_vs_neg_sig),])
library(ComplexHeatmap)
Heatmap(t(scale(t(log2(rpm[rownames(pos_vs_neg_sig),]+2)))), column_split = metadata$status)

#BiocManager::install("ComplexHeatmap")
library(EnhancedVolcano)
pos_vs_neg.df <- as.data.frame(pos_vs_neg)
EnhancedVolcano(pos_vs_neg.df,x="log2FoldChange",y="padj",lab = row.names(pos_vs_neg.df),
                FCcutoff = 2,pCutoffCol = "padj", pCutoff = 0.05)
pos_vs_neg.df$sig <- ifelse(pos_vs_neg.df$padj<0.05,"yes","no")
ggplot(pos_vs_neg.df, aes(log2FoldChange, -log10(padj), color=sig)) + 
  geom_point()
