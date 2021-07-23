library(readr)
library(phyloseq)
library(dplyr)
library(ape)
library(ggplot2)
library(vegan)
library(grid)
meta = read.table("sample-metadata-pueria.tsv", header=TRUE, sep="\t")
meta = meta[-c(6,9,11),]
row.names(meta) = meta$Sample
meta = meta[,-which(names(meta) %in% c("Sample"))]
OTU = read.table("ASV-table-306.txt", header=TRUE, sep="\t")
tax = read.table("tree_tax.tsv", header=TRUE, sep="\t")
meta$compartment <- as.factor(meta$compartment)
meta$SampleType <- as.factor(meta$SampleType)
meta$location <- as.factor(meta$location)
row.names(OTU) = OTU$OTU.ID
OTU.clean = OTU[,-which(names(OTU) %in% c("OTU.ID"))]
row.names(tax) = tax$ASV
tax.clean = tax[,-which(names(tax) %in% c("ASV"))]
OTU.UF=otu_table(as.matrix(OTU.clean),taxa_are_rows = T)
tax.UF=tax_table(as.matrix(tax.clean))
meta.UF=sample_data(as.data.frame(meta))
tree <- read.tree("rooted_tree_306.txt")
Ghana <- phyloseq(otu_table(OTU.UF),
                        tax_table(tax.UF),
                        phy_tree(tree),
                        sample_data(meta.UF))

Ghana
Ghana <- prune_samples(sample_sums(Ghana)>1500,Ghana)
Ghana <- prune_taxa(taxa_sums(Ghana) >0,Ghana)
Ghana
# calculate distance matrix 
wUF1_NMDS = ordinate(Ghana, method="NMDS", distance="unifrac", weighted=T)
data.scores = as.data.frame(scores(wUF1_NMDS))
#add columns to data frame 
data.scores$Sample = meta$LabID
data.scores$compartment = meta$compartment
data.scores$SampleType = meta$SampleType
data.scores$location = meta$location
head(data.scores)
p1 = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes(colour = compartment,shape = location))+
  theme(axis.text.y = element_text(colour = "black", size = 12), 
        axis.text.x = element_text(colour = "black", size = 12), 
        legend.text = element_text(size = 12, colour ="black"), 
        legend.position = c(0.89,0.93),
        legend.background = element_blank(),
        legend.spacing = unit(-0.5,"cm"),
        legend.key.size = unit(1,"lines"),
        axis.title.y = element_text( size = 14), 
        axis.title.x = element_text( size = 14, colour = "black"), 
        legend.title = element_blank(),
        #plot.margin =unit(c(1,1,1,1),"cm"),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "compartment", y = "NMDS2")  +
  geom_hline(yintercept =0,linetype="dashed",color="grey")+
  geom_vline(xintercept =0,linetype="dashed",color="grey")+
  scale_colour_manual(values = c("#e41a1c","#2c7fb8"),labels = c("Rhizosphere soil","Root"))
grob <- grobTree(textGrob("stress=0.09", x=0.88,  y=0.02, hjust=0,
                          gp=gpar(col="black", fontsize=10,face="bold")))
p1 <- p1 + annotation_custom(grob)
p1
ggsave("beta_diversity_weighted_compartment_location.tiff",width=8, height=6)
set.seed(1)
wUF.dist1 = UniFrac(Ghana, weighted=TRUE, normalized=T)
adonis(wUF.dist1 ~ compartment*location, data=meta, permutations = 999,p.adjust.methods="fdr")
# here for genus level 
library(microbiome)
library(dplyr)
library(hrbrthemes)
library(gcookbook)
library(tidyverse)
transform <- microbiome::transform
pseq <- transform(Ghana, "compositional")
pseq_Genus <- aggregate_taxa(pseq, level = "Genus")
OTU.genus=as(otu_table(pseq_Genus),"matrix")
OTU.genus=t(OTU.genus)
write.table(OTU.genus,"OTU.genus.telative.abundance.table.csv",row.names = F)
Genus.dist=vegdist(OTU.genus,method = "bray")
Genus.dist
adonis(Genus.dist ~ compartment*location,data = meta,permutations = 1000,p.adjust.methods="fdr")
#here we plot each genus at relative abundance
table(tax_table(Ghana)[,"Genus"])
Ghana_rel = transform_sample_counts(Ghana,function(x){x/sum(x)})
Ghana_rel
Ghana_rel_genus <- tax_glom(Ghana_rel,"Genus")
taxa_names(Ghana_rel_genus) <- tax_table(Ghana_rel_genus)[,"Genus"]
otu_table(Ghana_rel_genus)[1:5,1:5]
psmelt(Ghana_rel_genus) %>%
  ggplot(data= .,aes(x=compartment,y=Abundance))+
  geom_boxplot(outlier.shape = NA)+
  theme_bw()+
  geom_jitter(aes(color=OTU),height =0,width=.2)+
  labs(x="",y="Relative Abundance")+
  facet_wrap(~ OTU, scales = "free_x")
OTU.genus.rel=as(otu_table(Ghana_rel_genus),"matrix")
OTU.genus.rel=t(OTU.genus.rel)
Genus.rel.dist=vegdist(OTU.genus.rel,method = "bray")
adonis(Genus.rel.dist ~ compartment*location,data = meta,permutations = 1000,p.adjust.methods="fdr")
Ghana_genus_rel <- data.frame(otu_table(Ghana_rel_genus))
Ghana_genus_rel_transpose <- as.data.frame(t(as.matrix(Ghana_genus_rel)))
data <- cbind(sample_data(Ghana_rel),Ghana_genus_rel_transpose)
head(data)
shapiro.test(data$Paraglomus)
genus.wilcox <-wilcox.test(Paraglomus ~ compartment,data)
genus.wilcox
shapiro.test(data$Claroideoglomus)
genus.wilcox <-wilcox.test(Claroideoglomus ~ compartment,data)
genus.wilcox
shapiro.test(data$Diversispora)
genus.wilcox <-wilcox.test(Diversispora ~ compartment,data)
genus.wilcox
shapiro.test(data$Acaulospora)
genus.wilcox <-wilcox.test(Acaulospora ~ compartment,data)
genus.wilcox
shapiro.test(data$Septoglomus)
genus.wilcox <-wilcox.test(Septoglomus ~ compartment,data)
genus.wilcox
shapiro.test(data$Funneliformis)
genus.wilcox <-wilcox.test(Funneliformis ~ compartment,data)
genus.wilcox
shapiro.test(data$Dominikia)
genus.wilcox <-wilcox.test(Dominikia ~ compartment,data)
genus.wilcox
shapiro.test(data$Rhizophagus)
genus.wilcox <-wilcox.test(Rhizophagus ~ compartment,data)
genus.wilcox









