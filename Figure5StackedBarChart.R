# Figure 5: The relative abundance of the 500 most abundant prokaryotic amplicon sequence variant (ASV) in A) Station 1 free-living, B) Station 1 small particle-associated, C) Station 1 filter large particle-associated, D) Station 2 free-living, E) Station 2 small particle-associated, F) Station 2 large particle-associated, I) Station 3 free-living, J) Station 3 small particle-associated, K) Station 3 large particle-associated in the ETNP OMZ. The large particle-associated sample from 70 m at Station 3 (marked as “NA”) yielded < 20 raw sequencing reads and was excluded from the analysis. 
## Large volume size-fractioned filtration yields novel insights on prokaryotic community composition and interactions in the eastern tropical North Pacific ocean 
## Madeleine A. Thompson, David L. Valentine, Xuefeng Peng
## Code by Madeleine A. Thompson

#set working directory 
setwd("~/OneDrive - University of South Carolina/PengLab/Bioinformatics/etnp18_16s")


### Install Packages 
source("https://raw.githubusercontent.com/joey711/phyloseq/master/inst/scripts/installer.R",
       local = TRUE)
install_phyloseq(branch = "devel")
install.packages("ggplot2", "dplyr", "tidyr", "ape", "grid", "gridExtra", "vegan", "ggpubr", "ecolTest", "BiocManager")

############
## Install Packages
library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library(dplyr)
library(tidyr)
library(vegan)
library(ggpubr)
library(ggrepel)
library("viridis")

############

# Import data 
zotu16s <- read.csv("data/zotutab.csv", header = TRUE, row.names = 1)
metadata <- read.csv("data/metadata.csv", header = TRUE, row.names = 1)
tax16s <- read.csv("data/taxonomy_v1.csv",header = TRUE, row.names = 1)
metadata.phy = phyloseq::sample_data(data.frame(metadata))

#Separate taxa into separate columns and remove beginning string
taxGroup <- tidyr::separate(data = tax16s, col =  Taxon,
                            into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                            sep = ";", remove = TRUE, convert = FALSE, extra = "warn", fill = "warn")
taxGroup$Domain<-gsub("d__","",as.character(taxGroup$Domain))
taxGroup$Phylum<-gsub("p__","",as.character(taxGroup$Phylum))
taxGroup$Class<-gsub("c__","",as.character(taxGroup$Class))
taxGroup$Order <- gsub("o__", "", as.character(taxGroup$Order))
taxGroup$Family<-gsub("f__","",as.character(taxGroup$Family))
taxGroup$Genus<-gsub("g__","",as.character(taxGroup$Genus))
taxGroup$Species<-gsub("s__","",as.character(taxGroup$Species))

############


# Make ZOTU table numeric and a matrix 
# Make Taxonomy a matrix and make the row and column names the same 
mat<- as.matrix(zotu16s)
zotuTab <- matrix(as.numeric(mat),
                  ncol = ncol(mat));
colnames(zotuTab) <- colnames(mat);
rownames(zotuTab) <- rownames(mat);

taxTab <- as.matrix(taxGroup)

setdiff(rownames(taxTab), rownames(zotuTab))
all(rownames(zotuTab) == rownames(taxTab))
setdiff(colnames(zotuTab), rownames(taxTab))

############

# Combine Tax and ZOTU and metadata
OTU <- phyloseq::otu_table(zotuTab, taxa_are_rows = TRUE)
TAX <- phyloseq::tax_table(taxTab) 

physeq = phyloseq::phyloseq(OTU, TAX) # put metadata in this 
phyloseq::taxa_names(physeq)

physeq1 <-phyloseq::merge_phyloseq(physeq, metadata.phy)

#Phyloseq 
physeq_ns <- phyloseq::filter_taxa(physeq1, function(x) sum(x) > 1, prune = TRUE) 
phylor = phyloseq::transform_sample_counts(physeq_ns, function(x) 100 * x / sum(x) )
phylor_nz <- phyloseq::filter_taxa(phylor, function(x) mean(x) > 0, prune = TRUE)
top <- names(sort(phyloseq::taxa_sums(phylor_nz), decreasing=TRUE))[1:1100]
phylor.top <- phyloseq::prune_taxa(top, phylor_nz)

write.csv(phylor_nz@otu_table, "data/relabundzotu.csv")

# Separate by filter size and station
p.02 <- phyloseq::subset_samples(phylor.top, phylor.top@sam_data$Filter == 0.2)
p1.02 <- phyloseq::subset_samples(p.02, p.02@sam_data$Station == 1)
p2.02 <- phyloseq::subset_samples(p.02, p.02@sam_data$Station == 2)
p3.02 <- phyloseq::subset_samples(p.02, p.02@sam_data$Station == 3)

p.2 <- phyloseq::subset_samples(phylor.top, phylor.top@sam_data$Filter == 2.0)
p1.2 <- phyloseq::subset_samples(p.2, p.2@sam_data$Station == 1)
p2.2 <- phyloseq::subset_samples(p.2, p.2@sam_data$Station == 2)
p3.2 <- phyloseq::subset_samples(p.2, p.2@sam_data$Station == 3)

p.22 <- phyloseq::subset_samples(phylor.top, phylor.top@sam_data$Filter == 22)
p1.22 <- phyloseq::subset_samples(p.22, p.22@sam_data$Station == 1)
p2.22 <- phyloseq::subset_samples(p.22, p.22@sam_data$Station == 2)
p3.22 <- phyloseq::subset_samples(p.22, p.22@sam_data$Station == 3)

#Colors 
mycolors <- c("cadetblue", #acido
              "gold", #actino
              "green3", #cyano
              "palegreen", #bacter
              "pink", #gemma
              "deeppink", #marg
              
              "blue", #marini
              "orange", #placto
              "orangered", #verru
              'olivedrab',#chloro 
              "darkcyan", #NB1j
              "cyan", #SAR324
              
              "yellow", #Bdell
              "tomato", #dada
              "mediumpurple", #desulf
              "firebrick", #myxo
              "orchid", #nitrospin
              "navy", #nitrospir
              
              "maroon", #Pate
              "coral", # PAU
              'dodgerblue', # proteo 
              'cadetblue2', # hydro 
              'darkblue', #nano 
              'darkorchid', #thaum 
              
              'deeppink3', #thermo 
              "gray53") #NA

#plot function
my_plot_bar = function (physeq, x = "Sample", y = "Abundance", fill = NULL, title = NULL, 
                        facet_grid = NULL) {
  mdf = psmelt(physeq)
  p = ggplot(mdf, aes_string(x = x, y = y, fill = fill))
  p = p + geom_bar(stat = "identity", position = "stack")
  p = p + theme(axis.text.x = element_text(angle = -90, hjust = 0))
  if (!is.null(facet_grid)) {
    p <- p + facet_grid(facet_grid)
  }
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}

# Bar plots per filter size per station Phylum 
p1.02fig <- my_plot_bar(p1.02,  fill="Phylum", title = expression(paste("Station 1 - Free-living (0.22 – 2 ", mu, "m)")))
stat1.02um <- plot(p1.02fig + 
                     theme_classic()+
                     theme(text=element_text(size=25), 
                           axis.text.x = element_text(size = 25),
                           axis.text.y = element_text(size = 25),
                           legend.title = element_text(size = 40, face = "bold"), 
                           legend.text = element_text(size = 25),
                           legend.key.width = unit(5, "cm")) +
                     coord_flip() + 
                     ylim(0,100) +
                     scale_fill_manual(values = mycolors) +
                     ylab("Relative Abundance (%)") +
                     xlab("Depth (m)") +
                     scale_x_discrete(labels = c("120", "110", "83", "70", "40")))
p1.2fig <- my_plot_bar(p1.2,  fill="Phylum", title = expression(paste("Station 1 - Small Particle-associated (2 – 22 ", mu, "m)")))
stat1.2um <- plot(p1.2fig + 
                    theme_classic()+
                    theme(text=element_text(size=25), 
                          axis.text.x = element_text(size = 25),
                          axis.text.y = element_text(size = 25),
                          legend.title = element_text(size = 40, face = "bold"), 
                          legend.text = element_text(size = 25),
                          legend.key.width = unit(5, "cm")) +
                    coord_flip() + 
                    ylim(0,100) +
                    scale_fill_manual(values = mycolors) +
                    ylab("Relative Abundance (%)") +
                    xlab("Depth (m)") +
                    scale_x_discrete(labels = c("120", "110", "83", "70", "40")))
p1.22fig <- my_plot_bar(p1.22,  fill="Phylum", title = expression(paste("Station 1 - Large Particle-associated (> 22 ", mu, "m)")))
stat1.22um <- plot(p1.22fig + 
                     theme_classic()+
                     theme(text=element_text(size=25), 
                           axis.text.x = element_text(size = 25),
                           axis.text.y = element_text(size = 25),
                           legend.title = element_text(size = 40, face = "bold"), 
                           legend.text = element_text(size = 25),
                           legend.key.width = unit(5, "cm")) +
                     coord_flip() + 
                     ylim(0,100) +
                     scale_fill_manual(values = mycolors) +
                     ylab("Relative Abundance (%)") +
                     xlab("Depth (m)") +
                     scale_x_discrete(labels = c("120", "110", "83", "70", "40")))
p2.02fig <- my_plot_bar(p2.02,  fill="Phylum", title = expression(paste("Station 2 - Free-living (0.22 – 2 ", mu, "m)")))
stat2.02um <- plot(p2.02fig + 
                     theme_classic()+
                     theme(text=element_text(size=25), 
                           axis.text.x = element_text(size = 25),
                           axis.text.y = element_text(size = 25),
                           legend.title = element_text(size = 40, face = "bold"), 
                           legend.text = element_text(size = 25),
                           legend.key.width = unit(5, "cm")) +
                     coord_flip() + 
                     ylim(0,100) +
                     scale_fill_manual(values = mycolors) +
                     ylab("Relative Abundance (%)") +
                     xlab("Depth (m)") +
                     scale_x_discrete(labels = c("900", "120", "113", "90", "30")))
p2.2fig <- my_plot_bar(p2.2,  fill="Phylum", title = expression(paste("Station 2 - Small Particle-associated (2 – 22 ", mu, "m)")))
stat2.2um <- plot(p2.2fig + 
                    theme_classic()+
                    theme(text=element_text(size=25), 
                          axis.text.x = element_text(size = 25),
                          axis.text.y = element_text(size = 25),
                          legend.title = element_text(size = 40, face = "bold"), 
                          legend.text = element_text(size = 25),
                          legend.key.width = unit(5, "cm")) +
                    coord_flip() + 
                    ylim(0,100) +
                    scale_fill_manual(values = mycolors) +
                    ylab("Relative Abundance (%)") +
                    xlab("Depth (m)") +
                    scale_x_discrete(labels = c("900", "120", "113", "90", "30")))
p2.22fig <- my_plot_bar(p2.22,  fill="Phylum", title = expression(paste("Station 2 - Large Particle-associated (> 22 ", mu, "m)")))
stat2.22um <- plot(p2.22fig + 
                     theme_classic()+
                     theme(text=element_text(size=25), 
                           axis.text.x = element_text(size = 25),
                           axis.text.y = element_text(size = 25),
                           legend.title = element_text(size = 40, face = "bold"), 
                           legend.text = element_text(size = 25),
                           legend.key.width = unit(5, "cm")) +                    
                     coord_flip() + 
                     ylim(0,100) +
                     scale_fill_manual(values = mycolors) +
                     ylab("Relative Abundance (%)") +
                     xlab("Depth (m)") +
                     scale_x_discrete(labels = c("900", "120", "113", "90", "30")))
p3.02fig <- my_plot_bar(p3.02,  fill="Phylum", title = expression(paste("Station 3 - Free-living (0.22 – 2 ", mu, "m)")))
stat3.02um <- plot(p3.02fig + 
                     theme_classic()+
                     theme(text=element_text(size=25), 
                           axis.text.x = element_text(size = 25),
                           axis.text.y = element_text(size = 25),
                           legend.title = element_text(size = 40, face = "bold"), 
                           legend.text = element_text(size = 25),
                           legend.key.width = unit(5, "cm")) +
                     coord_flip() + 
                     ylim(0,100) +
                     scale_fill_manual(values = mycolors) +
                     ylab("Relative Abundance (%)") +
                     xlab("Depth (m)") +
                     scale_x_discrete(labels = c("1000", "70", "45", "40", "33", "10")))
p3.2fig <- my_plot_bar(p3.2, fill="Phylum", title = expression(paste("Station 3 - Small Particle-associated (2 – 22 ", mu, "m)")))
stat3.2um <- plot(p3.2fig + 
                    theme_classic()+
                    theme(text=element_text(size=25), 
                          axis.text.x = element_text(size = 25),
                          axis.text.y = element_text(size = 25),
                          legend.title = element_text(size = 40, face = "bold"), 
                          legend.text = element_text(size = 25),
                          legend.key.width = unit(5, "cm")) +
                    coord_flip() + 
                    ylim(0,100) +
                    scale_fill_manual(values = mycolors) +
                    ylab("Relative Abundance (%)") +
                    xlab("Depth (m)") +
                    scale_x_discrete(labels = c("1000", "70", "45", "40", "33", "10")))
p3.22fig <- my_plot_bar(p3.22,  fill="Phylum", title = expression(paste("Station 3 - Large Particle-associated (> 22 ", mu, "m)")))
stat3.22um <- plot(p3.22fig + 
                     theme_classic()+
                     theme(text=element_text(size=25), 
                           axis.text.x = element_text(size = 25),
                           axis.text.y = element_text(size = 25),
                           legend.title = element_text(size = 25, face = "bold"), 
                           legend.text = element_text(size = 25),
                           legend.key.width = unit(5, "cm")) +
                     coord_flip() + 
                     ylim(0,100) +
                     scale_fill_manual(values = mycolors) +
                     ylab("Relative Abundance (%)") +
                     xlab("Depth (m)") +
                     scale_x_discrete(labels = c("1000", "70", "45", "40", "33", "10")))

#put together 
figure <- ggpubr::ggarrange(stat1.02um, stat1.2um, stat1.22um, 
                            stat2.02um, stat2.2um, stat2.22um, 
                            stat3.02um, stat3.2um, stat3.22um, 
                            labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"), 
                            font.label = list(size = 25, color = "black", face = "bold", family = NULL),
                            ncol = 3, nrow = 3,
                            common.legend = TRUE, legend= "bottom")
ggplot2::ggsave("figuresNew/NewPhylumA.png", plot=figure, device="png",
                scale=1, width = 110, height=45, units=c("cm"), dpi=300, limitsize = FALSE)


