library(readr)
library(phyloseq)
library(ape)
library(ggplot2)
library(DESeq2)
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
Ghana <- prune_samples(sample_sums(Ghana)>1500,Ghana) #remove low sequence depth
Ghana <- prune_taxa(taxa_sums(Ghana) >0,Ghana)
Ghana
#DESeq2
comp = phyloseq_to_deseq2(Ghana,~ compartment)
comp = DESeq(comp,test = "Wald",fitType ="parametric",sfType = "poscounts")
res = results(comp,cooksCutoff = F)
alpha = 0.05
sigtab = res[which(res$padj < alpha),]
sigtab = cbind(as(sigtab,"data.frame"),as(tax_table(Ghana)[rownames(sigtab),],"matrix"))
head(sigtab)
dim(sigtab)
library(xlsx)
write.xlsx(sigtab,file="DeSeq2_results.xlsx",sheetName="ASV_difference_0.05",append=T)
newcol <- c("ASV170","ASV123","ASV203","ASV274","ASV252","ASV28","ASV133",
            "ASV152","ASV258","ASV114","ASV225","ASV230",
            "ASV290","ASV218","ASV84","ASV32","ASV99","ASV149","ASV238")
data1 <- sigtab
data1$ID <- newcol
data1
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set2", ...) {
  scale_fill_brewer(palette = palname, ...)
}
x = tapply(data1$log2FoldChange, data1$ID, function(x) max(x))
x = sort(x, TRUE)
data1$ID = factor(as.character(data1$ID), levels=names(x))
p <- ggplot(data1, aes(y=ID, x=log2FoldChange,color=Genus))+
  geom_point(aes(size = baseMean))+ 
  labs(x="log2 Fold Change",y="",size="Mean counts")+
  geom_vline(xintercept = 0.0,color="grey",size=0.5,linetype="dashed")+
  theme(legend.spacing = unit(-0.2,"cm"),
        legend.key.size = unit(0.8,"lines"))
p
ggsave("Deseq_results.tiff",width = 5,height = 4)
library(xlsx)
write.xlsx(data1,file="DeSeq2_results_2.xlsx",sheetName="ASV_difference")