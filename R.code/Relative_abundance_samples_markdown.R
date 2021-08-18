library(phyloseq)
library(dplyr)
library(ape)
library(ggplot2)
library(RColorBrewer)
library(ggtext)
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
otu <- otu_table(Ghana)
otu <- as.data.frame(t(otu))
library(xlsx)
write.xlsx(otu, "otu_puerira_t.xlsx", sheetName = "otu", 
           col.names = TRUE, row.names = TRUE, append = FALSE)
meta = read.table("sample-metadata-pueria.tsv", header=TRUE, sep="\t")
meta = meta[-c(6,9,11),]
meta = meta %>%
  select(Sample,compartment)
library(readxl)
library(tidyverse)
otu_counts = read_excel("otu_puerira_t.xlsx") 
colnames(otu_counts)[colnames(otu_counts) == "...1"] <- "Sample"
otu_counts = otu_counts %>%
  pivot_longer(-Sample,names_to="ASV",values_to="counts")
taxonomy = read.table("tree_tax.tsv",header = T)
otu_rel_abund <- inner_join(meta, otu_counts,by="Sample") %>%
  inner_join(.,taxonomy,by="ASV") %>%
  group_by(Sample) %>%
  mutate(rel_abund = counts / sum(counts)) %>%
  ungroup() %>%
  select(-counts) %>%
  pivot_longer(c("Phylum","Class","Order","Family","Genus","ASV"),
               names_to="level", values_to="taxon") %>%
  mutate(compartment = factor(compartment,
                              levels = c("rhizosphere","root")))
taxon_rel_abund <- otu_rel_abund %>%
  filter(level=="Genus") %>%
  group_by(compartment,Sample,taxon) %>%
  summarise(rel_abund=100*sum(rel_abund),.groups="drop") %>%
  mutate(taxon=str_replace(taxon,"^(\\S*)$","*\\1*"))
  
taxon_pool <- taxon_rel_abund %>%
  group_by(compartment,taxon) %>%
  summarise(mean=mean(rel_abund),.groups="drop") %>%
  group_by(taxon) %>%
  summarise(pool=max(mean) <10,
            mean=mean(mean),
            .groups="drop")
sample_order <- taxon_rel_abund %>%
  filter(taxon=="*Acaulospora*") %>%
  arrange(desc(rel_abund)) %>%
  mutate(order=1:nrow(.)) %>%
  select(Sample,order)
pretty <- c("rhizosphere"="Rhizosphere soil",
             "root"="Root")
inner_join(taxon_rel_abund,taxon_pool,by="taxon") %>%
   mutate(taxon=if_else(pool,"Others",taxon)) %>%
  group_by(Sample,compartment,taxon) %>%
  summarise(rel_abund=sum(rel_abund),
            mean=min(mean),
            .groups="drop") %>%
  mutate(taxon=factor(taxon),
         taxon=fct_reorder(taxon,mean,.desc = TRUE),
         taxon=fct_shift(taxon,n=1)) %>%
  inner_join(.,sample_order,by="Sample") %>%
  mutate(Sample=factor(Sample),
         Sample=fct_reorder(Sample,order)) %>%
  ggplot(aes(x=Sample,y=rel_abund,fill=taxon)) +
  geom_col(width = 0.8)+
  scale_fill_manual(breaks=c("*Acaulospora*","*Rhizophagus*","*Paraglomus*",
                             "*Claroideoglomus*","*Dominikia*","Others"),
   values =c(brewer.pal(6,"Dark2")))+
  scale_y_continuous(expand = c(0,0))+
  facet_grid(~compartment,scale = "free_x",space = "free",switch = "x",
             labeller = labeller(compartment=pretty))+
  labs(x=NULL, fill= "Genus",
       y="Relative Abundance (%)")+
  theme_classic()+
  theme(axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_markdown(size=8,color = "black"),
        legend.key.size = unit(10,"pt"),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.switch.pad.grid=unit(-0.2,"cm"),
        strip.text.x = element_text(size=8,color = "black"),
        axis.text.y = element_text(size=8,color = "black")
        
  )
ggsave("relative_abundance_samples_markdown.tiff",width = 5,height = 4)
