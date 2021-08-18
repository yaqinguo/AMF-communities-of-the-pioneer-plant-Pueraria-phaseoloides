library(readr)
library(dplyr)
library(phyloseq)
library(vegan)
library(ape)
library(ggplot2)
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
otu <- otu_table(Ghana)
otu <- as.data.frame(t(otu))
#import environmental factors
EF = read.table("Ghana-environmental-factors.tsv", header=TRUE, sep="\t")
EF = EF[-c(6,9,11),]
rownames(EF) <- EF$X
EF$X <- NULL
head(EF)
# data normolization   (Legendre and Gallagher,2001)
env.log <- log1p(EF) # by log
env.total <- na.omit(env.log) #delete NA
otu.total.hell <- decostand(otu,"hellinger")
#DCA to test axis length to choose rda or cca
sel.total <- decorana(otu.total.hell)
sel.total #DCA1=1
otu.tab.0 <- rda(otu.total.hell ~ 1,env.total)
otu.tab.0
otu.tab.1 <- rda(otu.total.hell ~ .,env.total)
otu.tab.1   
otu.tab.2 <- rda(otu.total.hell ~ AK + Ca + Cu + Fe + K + Mg + Mn + P + Zn + sand + silt + soilMoisture + pH + N + C + C.N,env.total)
vif.cca(otu.tab.2)
max(vif.cca(otu.tab.2),na.rm = T) #remove sand

otu.tab.3 <- rda(otu.total.hell ~ AK + Ca + Cu + Fe + K + Mg + Mn + P + Zn + silt + soilMoisture + pH + N + C + C.N,env.total)
vif.cca(otu.tab.3)
max(vif.cca(otu.tab.3),na.rm = T) #remove K

otu.tab.4 <- rda(otu.total.hell ~ AK + Ca + Cu + Fe + Mg + Mn + P + Zn + silt + soilMoisture + pH + N + C + C.N,env.total)
vif.cca(otu.tab.4)
max(vif.cca(otu.tab.4),na.rm = T) #remove N

otu.tab.5 <- rda(otu.total.hell ~ AK + Ca + Cu + Fe + Mg + Mn + P + Zn + silt + soilMoisture + pH + C + C.N,env.total)
vif.cca(otu.tab.5)
max(vif.cca(otu.tab.5),na.rm = T) #remove Ca

otu.tab.6 <- rda(otu.total.hell ~ AK + Cu + Fe + Mg + Mn + P + Zn + silt + soilMoisture + pH + C + C.N,env.total)
vif.cca(otu.tab.6)
max(vif.cca(otu.tab.6),na.rm = T) #remove Mg

otu.tab.7 <- rda(otu.total.hell ~ AK + Fe +Cu + Mn + P + Zn + silt + soilMoisture + pH + C + C.N,env.total)
vif.cca(otu.tab.7)
max(vif.cca(otu.tab.7),na.rm = T) #remove Cu

otu.tab.8 <- rda(otu.total.hell ~ AK + Fe + Mn + P + Zn + silt + soilMoisture + pH + C + C.N,env.total)
vif.cca(otu.tab.8)
max(vif.cca(otu.tab.8),na.rm = T) #remove Fe

otu.tab.9 <- rda(otu.total.hell ~ AK + Mn + P + Zn + silt + soilMoisture + pH + C + C.N,env.total)
vif.cca(otu.tab.9)
max(vif.cca(otu.tab.9),na.rm = T)
#chaeck AIC for model 
mod.d <- step(otu.tab.0, scope = (list(lower = formula(otu.tab.0), upper = formula(otu.tab.1))))
mod.d
set.seed(1)
anova(otu.tab.9) #overall test
anova(otu.tab.9,by="margin",nperm=49999) #test each parameter
#adjust p.adjust
test_variables <- anova(otu.tab.9,by="margin")
test_variables.adj <- test_variables
test_variables.adj$'P.adjust' <- p.adjust(test_variables$`Pr(>F)`,method = "fdr")
test_variables.adj
#test for rda axis
test_axis <- anova(otu.tab.9, by="axis")
test_axis.adj <- test_axis
test_axis.adj$'P.adjust' <- p.adjust(test_axis$`Pr(>F)`,method = "fdr")
test_axis.adj
wUF.dist1 = UniFrac(Ghana, weighted=TRUE, normalized=F)
set.seed(1)
db.rda <- dbrda(wUF.dist1 ~ AK + Mn + P + Zn + silt + soilMoisture + pH + C + C.N, env.total)
db.rda
db.rda_results <- summary(db.rda)
df1 <- data.frame(db.rda_results$sites[,1:2])
df1 <- cbind(df1,meta)
df2 <- data.frame(db.rda_results$biplot[,1:2])
df3 <- data.frame(db.rda_results$species[,1:2])
#plot biplot
library(ggpubr)
db.rda.biplot <- ggplot()%>%
  ggpar(xlab="dbRDA1 (26.7%)",ylab="dbRDA2(13.7%)",
        font.x=14,font.y=14,font.tickslab=14,ggtheme=theme_bw())+
  theme(plot.title = element_text(hjust=0.5,face="bold"),
        panel.grid = element_blank(),
    aspect.ratio = 0.8,axis.title = element_text(size = 12,face = "bold"),
        axis.ticks.length = unit(-0.25,'cm'),
        axis.title.y.left = element_text(margin = unit(c(0.5,0.5,0.5,0.5),'cm'),size = 10),
        axis.title.x.bottom = element_text(margin = unit(c(0.5,0.5,0.5,0.5),'cm'),size = 10)
        )+
  geom_vline(xintercept = 0,color='black',size=0.6,linetype=2)+
  geom_hline(yintercept = 0,color='black',size=0.6,linetype=2)+
  theme(legend.spacing = unit(0.5,'cm'))+
  #add site score
  geom_point(data = df1,aes(x=dbRDA1,y=dbRDA2,color=compartment),size=3)+
  theme(legend.position = c(0.88,0.95),legend.title = element_blank(),legend.background = element_blank(),
        legend.text = element_text(size=10,face="bold"),legend.spacing =unit(0.5,'cm'))+
  scale_colour_manual(values = c("#e41a1c","#2c7fb8"),labels = c("Rhizosphere soil","Root"))+
  #add environmental arrows
  geom_segment(data = df2,color="black",size=0.6,
               aes(x=0,y=0,xend=dbRDA1,yend=dbRDA2),
               arrow = arrow(length = unit(0.01,"npc")))+
  geom_text(data = df2,aes(x=dbRDA1,y=dbRDA2,label=row.names(df2),
                           hjust=0.4,vjust=0.5*(1-sign(dbRDA2))))
db.rda.biplot
ggsave("db.rda.biplot-remove_Hg.tiff",width = 8,height = 6)




