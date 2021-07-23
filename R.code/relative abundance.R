library(phyloseq)
library(dplyr)
library(ape)
library(ggplot2)
meta = read.table("sample-metadata-pueria.tsv", header=TRUE,row.names = 1, sep="\t")
OTU = read.table("ASV-table-306.txt", header=TRUE, sep="\t")
tax = read.table("tree_tax.tsv", header=TRUE, sep="\t")
meta$compartment <- as.factor(meta$compartment)
meta$location <- as.factor(meta$location)
row.names(OTU) = OTU$OTU.ID
OTU.clean = OTU[,-which(names(OTU) %in% c("OTU.ID"))]
row.names(tax) = tax$ASV
tax.clean = tax[,-which(names(tax) %in% c("ASV"))]
OTU.UF=otu_table(as.matrix(OTU.clean),taxa_are_rows = TRUE)
tax.UF=tax_table(as.matrix(tax.clean))
meta.UF=sample_data(as.data.frame(meta))
tree <- read.tree("rooted_tree_306.txt")
Ghana <- phyloseq(otu_table(OTU.UF),
                  tax_table(tax.UF),
                  phy_tree(tree),
                  sample_data(meta.UF))
Ghana
Ghana <- prune_samples(sample_sums(Ghana)>1500,Ghana)
Ghana <- prune_taxa(taxa_sums(Ghana)>0,Ghana)
Ghana
###here to plot relative abundance
table(tax_table(Ghana)[,"Genus"])
Ghana_rel = transform_sample_counts(Ghana,function(x){x/sum(x)})
otu_table(Ghana_rel)[1:5,1:5]
new.labs <- c("Rhizosphere soil","Root")
names(new.labs) <- c("rhizosphere","root")
plot_bar(Ghana_rel,fill="Genus")+
  geom_bar(aes(fill=Genus),stat="identity",position =position_stack(reverse=T))+
  labs(x="",y="Relative Abundance (%)")+
  facet_grid(~compartment,scales = "free",labeller = labeller(compartment=new.labs),
             switch = "x")+
  scale_fill_brewer(palette="Set3")+
  scale_y_continuous(expand=c(0,0),labels=function(x) paste0(x*100))+
  guides(fill = guide_legend(reverse = F, keywidth = 0.9, keyheight = 0.9,title = "Genus")) +
  theme(axis.text.y = element_text(color="black",size=10),
        axis.title.y.left = element_text(size=12),
        legend.position = "right",axis.title.y = element_text(size=10),
        legend.title = element_text(color="black",size=12),
        panel.background = element_blank(),axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),strip.background = element_blank(),
        axis.ticks.y = element_line(linetype = "solid",size=0.7),
        axis.line.y = element_line(linetype = "solid",size=1),
        panel.spacing = unit(0.2,"lines"),
        strip.text.x = element_text(color="black",size=12),
        strip.placement = "outside",
        strip.switch.pad.grid = unit(-0.3,"lines"),
        plot.margin = margin(2, 1, 2, 1, "cm"))+
  geom_hline(yintercept = -0.005,linetype="solid",color="black",size=2)
ggsave("relative_abundance_root_rhizosphere.tiff",width = 6,height = 5)
