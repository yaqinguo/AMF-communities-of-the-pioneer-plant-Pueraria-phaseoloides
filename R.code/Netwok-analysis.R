library(readr)
library(phyloseq)
library(dplyr)
library(ape)
library(psych)
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
Ghana_root <- subset_samples(Ghana,compartment=="root")
Ghana_root
otu.r <- otu_table(Ghana_root)
otu.r <- as.data.frame(t(otu.r))
Ghana_rhizosphere <- subset_samples(Ghana,compartment=="rhizosphere")
Ghana_rhizosphere
otu.s <- otu_table(Ghana_rhizosphere)
otu.s <- as.data.frame(t(otu.s))
dim(otu.s)
dim(otu.r)
colSums(otu.r)
otu.r <- otu.r[, which(colSums(otu.r) != 0)]
colSums(otu.r)
colSums(otu.s)
otu.s <- otu.s[, which(colSums(otu.s) != 0)]
colSums(otu.s)
dim(otu.s)
dim(otu.r)
#for rhizosphere
s.occor = corr.test(otu.s,use="pairwise",method="spearman",adjust="fdr",alpha=0.05)
s.occor.r = s.occor$r 
s.occor.p = s.occor$p 
head(s.occor.r)
head(s.occor.p)
s.occor.r[s.occor.p>0.05|abs(s.occor.r)<0.6] = 0 
s.occor.p
s.occor.r
write.csv(s.occor.r,"rhizosphere_corror_0.05.csv")
#for root
r.occor = corr.test(otu.r,use="pairwise",method="spearman",adjust="fdr",alpha=.05)
r.occor.r = r.occor$r 
r.occor.p = r.occor$p 
r.occor.r[r.occor.p>0.05|abs(r.occor.r)<0.6] = 0 
write.csv(r.occor.r,"root_corror_0.05.csv")