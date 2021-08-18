library(readr)
library(phyloseq)
library(dplyr)
library(ape)
library(tidyverse)
library(ggplot2)
library(vegan)
library(ggpubr)
meta = read.table("sample-metadata-pueria.tsv", header=TRUE, sep="\t")
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
# plot rarefaction for each compartment
Ghana_root <- subset_samples(Ghana,compartment=="root")
Ghana_root
otu.r <- otu_table(Ghana_root)
otu.r <- as.data.frame(t(otu.r))
sample_names.r <- rownames(otu.r)
out.r <- rarecurve(otu.r,step=1,label=F) #vegan
rare.r <- lapply(out.r,function(x){
  b <- as.data.frame(x)
  b <- data.frame(OTU=b[,1],raw.read=rownames(b))
  b$raw.read <- as.numeric(gsub("N","",b$raw.read))
  return(b)
})
names(rare.r) <- sample_names.r
rare.r <- map_dfr(rare.r, function(x){
  z <- data.frame(x)
  return(z)
}, .id = "sample")
head(rare.r)
groupings.r <- data.frame(sample=sample_names.r,compartment="Root",stringsAsFactors = F)
groupings.r
rare.r <- map_dfr(groupings.r$sample, function(x){ #loop over samples
  z <- rare.r[rare.r$sample == x,] #subset rare according to sample 
  compartment <- groupings.r$compartment[groupings.r$sample == x] #subset groupings according to sample, if more than one grouping repeat for all
  z <- data.frame(z, compartment) #make a new data frame with the subsets
  return(z)
})
head(rare.r)
p_root_rare <- ggplot(data=rare.r)+
  geom_line(aes(x=raw.read,y=OTU,group=sample,col=compartment))+
  theme_bw()+
  theme(axis.text.y = element_text(colour = "black", size = 12), 
        axis.text.x = element_text(colour = "black", size = 12), 
        legend.position ="none", 
        axis.title.x = element_text( size = 14, colour = "black"), 
        axis.title.y = element_text( size = 14, colour = "black"),
        plot.margin = unit(c(1,1,1,1),"cm"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border=element_blank(),axis.line=element_line(color="black")) + 
  labs(x = "Number of sequences", y = "Number of ASV")+
  scale_x_continuous(expand = c(0, 0),limits = c(0,2500))+
  scale_y_continuous(expand = c(0, 0),limits = c(0,21))+
  geom_vline(xintercept =1501,linetype="dashed",color="black")+
  scale_color_manual(values="#2c7fb8")
p_root_rare
ggsave("rarecurve_root.tiff",width = 5,height = 4)
#rarefaction for rhizosphere
Ghana_rhizosphere <- subset_samples(Ghana,compartment=="rhizosphere")
Ghana_rhizosphere
otu.s <- otu_table(Ghana_rhizosphere)
otu.s <- as.data.frame(t(otu.s))
sample_names.s <- rownames(otu.s)
out.s <- rarecurve(otu.s,step=1,label=F)
rare.s <- lapply(out.s,function(x){
  b <- as.data.frame(x)
  b <- data.frame(OTU=b[,1],raw.read=rownames(b))
  b$raw.read <- as.numeric(gsub("N","",b$raw.read))
  return(b)
})
rare.s
names(rare.s) <- sample_names.s
rare.s <- map_dfr(rare.s, function(x){
  z <- data.frame(x)
  return(z)
}, .id = "sample")
head(rare.s)
groupings.s <- data.frame(sample=sample_names.s,compartment="Rhizosphere soil",stringsAsFactors = F)
groupings.s
rare.s <- map_dfr(groupings.s$sample, function(x){ #loop over samples
  z <- rare.s[rare.s$sample == x,] #subset rare according to sample 
  compartment <- groupings.s$compartment[groupings.s$sample == x] #subset groupings according to sample, if more than one grouping repeat for all
  z <- data.frame(z, compartment) #make a new data frame with the subsets
  return(z)
})
head(rare.s)
p_rhizosphere_rare <- ggplot(data=rare.s)+
  geom_line(aes(x=raw.read,y=OTU,group=sample,col=compartment))+
  theme_bw()+
  theme(axis.text.y = element_text(colour = "black", size = 12), 
        axis.text.x = element_text(colour = "black", size = 12), 
        legend.position ="none", 
        axis.title.x = element_text( size = 14, colour = "black"), 
        axis.title.y = element_text( size = 14, colour = "black"),
        plot.margin = unit(c(1,1,1,1),"cm"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border=element_blank(),axis.line=element_line(color="black")) + 
  labs(x = "Number of sequences", y = "Number of ASV")+
  scale_x_continuous(expand = c(0, 0),limits = c(0,2050))+
  scale_y_continuous(expand = c(0, 0),limits = c(0,31))+
  geom_vline(xintercept =1501,linetype="dashed",color="black")+
  scale_color_manual(values="#e41a1c")
p_rhizosphere_rare 
ggsave("rarecure_rhizosphere.tiff",width = 5,height = 4)
ggarrange(p_rhizosphere_rare,p_root_rare, labels = "AUTO",vjust = 2.5,hjust = -4)
ggsave("rarecurve_rhizosphere_root.tiff",width = 8,height = 6)
