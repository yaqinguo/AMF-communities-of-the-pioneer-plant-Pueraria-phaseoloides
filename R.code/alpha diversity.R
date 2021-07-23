library(phyloseq)
library(dplyr)
library(vegan)
library(ape)
library(ggplot2)
library(ggpubr)
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
tree
Ghana <- phyloseq(otu_table(OTU.UF),
                  tax_table(tax.UF),
                  phy_tree(tree),
                  sample_data(meta.UF))
Ghana
Ghana <- prune_samples(sample_sums(Ghana)>1500,Ghana) #remove low sequence depth
Ghana <- prune_taxa(taxa_sums(Ghana)>0,Ghana) # remove zero tax
sample_sums(Ghana) #check sample sums
sample_names(Ghana) #check sample names
Ghana
Ghana_rarefied <- rarefy_even_depth(Ghana, sample.size = min(sample_sums(Ghana)) ,
                                    rngseed = 1, replace = F)
Ghana_rarefied
sample_sums(Ghana_rarefied) #check all get rarefied
alpha.diversity_rare <- estimate_richness(Ghana_rarefied,
                                      measures = c("observed","Chao1","Shannon", "InvSimpson"))
head(alpha.diversity_rare)
data_rare <- cbind(sample_data(Ghana), alpha.diversity_rare)
H <- data_rare$Shannon
S1 <- data_rare$Observed
S <- log(S1)
evenness <- H/S
evenness
data_rare$Evenness = evenness
head(data_rare)
data_rare_K <- subset(data_rare, location == "Konongo")
data_rare_B <- subset(data_rare, location == "Bosome-Freho ")
Observed.wilcox <-wilcox.test(Observed ~ compartment,data=data_rare_K)
Observed.wilcox
Shannon.wilcox <-wilcox.test(Shannon ~ compartment,data=data_rare_K)
Shannon.wilcox 
InvSimpson.wilcox <-wilcox.test(InvSimpson ~ compartment,data=data_rare_K)
InvSimpson.wilcox 
Evenness.wilcox <-wilcox.test(Evenness ~ compartment,data=data_rare_K)
Evenness.wilcox
#another location
Observed.wilcox <-wilcox.test(Observed ~ compartment,data=data_rare_B)
Observed.wilcox
Shannon.wilcox <-wilcox.test(Shannon ~ compartment,data=data_rare_B)
Shannon.wilcox 
InvSimpson.wilcox <-wilcox.test(InvSimpson ~ compartment,data=data_rare_B)
InvSimpson.wilcox 
Evenness.wilcox <-wilcox.test(Evenness ~ compartment,data=data_rare_B)
Evenness.wilcox 
#plot
p.observed_box <-ggplot(data_rare,aes(x=compartment,y=Observed,fill=compartment))+
  facet_wrap(~location)+
  facet_grid(~location,scales = "free_x", space="free_x")+
  stat_boxplot(geom = 'errorbar',width=0.1)+
  #geom_violin( trim = F)
  geom_boxplot(outlier.shape = NA,width=0.5)+
  scale_fill_manual(labels = c("Rhizosphere soil","Root"),values = c("#bf812d", "#f6e8c3") )+
  theme_bw()+
  theme(axis.text.y = element_text(colour = "black", size = 12),
        legend.text = element_text(size = 12, colour ="black"), 
        legend.position = "top",
        axis.title.y = element_text( size = 12,colour ="black"),
        axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  labs(y = "ASV richness")+
  scale_y_continuous(expand = c(0, 0),limits = c(0,45))
p.observed_box
ggsave("Observed_compartment_2locations_box.tiff",width=5,height = 4)

p.Evenness_box <-ggplot(data_rare,aes(x=compartment,y=Evenness,fill=compartment))+
  facet_wrap(~location)+
  facet_grid(~location,scales = "free_x", space="free_x")+
  stat_boxplot(geom = 'errorbar',width=0.1)+
  #geom_violin( trim = F)
  geom_boxplot(outlier.shape = NA,width=0.5)+
  scale_fill_manual(values = c("#bf812d", "#f6e8c3"))+
  theme_bw()+
  theme(axis.text.y = element_text(colour = "black", size = 12),
        legend.text = element_text(size = 12, colour ="black"), 
        legend.position = "top",
        legend.background = element_blank(),
        axis.title.y = element_text( size = 12,colour ="black"),
        axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  labs(y = "ASV evenness")+
  scale_y_continuous(expand = c(0, 0),limits = c(0,1),labels = scales::number_format(accuracy = 0.1))
p.Evenness_box
ggsave("Evenness_compartment_2locations_box.tiff",width=5,height = 4)

p.Shannon_box <-ggplot(data_rare,aes(x=compartment,y=Shannon,fill=compartment))+
  facet_wrap(~location)+
  facet_grid(~location,scales = "free_x", space="free_x")+
  stat_boxplot(geom = 'errorbar',width=0.1)+
  #geom_violin( trim = F)
  geom_boxplot(outlier.shape = NA,width=0.5)+
  scale_fill_manual(values = c("#bf812d", "#f6e8c3"))+
  theme_bw()+
  theme(axis.text.y = element_text(colour = "black", size = 12),
        legend.text = element_text(size = 12, colour ="black"), 
        legend.position = "top",
        legend.background = element_blank(),
        axis.title.y = element_text( size = 12,colour ="black"),
        axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  labs(y = "Shannon divsersity")+
  scale_y_continuous(expand = c(0, 0),limits = c(0,4),labels = scales::number_format(accuracy = 0.1))
p.Shannon_box
ggsave("Shannon_compartment_2locations.tiff",width=5,height = 4)

p.InvSimpson_box <-ggplot(data_rare,aes(x=compartment,y=InvSimpson,fill=compartment))+
  facet_wrap(~location)+
  facet_grid(~location,scales = "free_x", space="free_x")+
  stat_boxplot(geom = 'errorbar',width=0.1)+
  #geom_violin( trim = F)
  geom_boxplot(outlier.shape = NA,width=0.5)+
  scale_fill_manual(values = c("#bf812d", "#f6e8c3"))+
  theme_bw()+
  theme(axis.text.y = element_text(colour = "black", size = 12),
        legend.text = element_text(size = 12, colour ="black"), 
        legend.position = "top",
        legend.background = element_blank(),
        axis.title.y = element_text( size = 12,colour ="black"),
        axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  labs(y = "InvSimpson divsersity")+
scale_y_continuous(expand = c(0, 0),limits = c(0,25))
p.InvSimpson_box
ggsave("InvSimpson_compartment_2locations.tiff",width=5,height = 4)
p_box_4figures <- ggarrange(p.observed_box, p.Evenness_box, p.InvSimpson_box,p.Shannon_box, labels = "AUTO",
                               common.legend = TRUE, legend = "top")
p_box_4figures
ggsave("p_box_4figures.tiff",width=8,height = 6)








