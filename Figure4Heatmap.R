#Figure 4: The normalized (by trimmed mean of M values) abundance (log10 of counts per million calculated by EdgeR) of 48 most abundant amplicon sequence variants (rows) that were differentially abundant between surface, upper oxycline and anoxic depths, and lower oxycline depths in all samples (columns) collected from the ETNP OMZ in the free-living (FL), small particle-associated (SPA), and large particle-associated (LPA) size fractions.
## Large volume size-fractioned filtration yields novel insights on prokaryotic community composition and interactions in the eastern tropical North Pacific ocean 
## Madeleine A. Thompson, David L. Valentine, Xuefeng Peng

#Set working directory

setwd("~/Library/CloudStorage/OneDrive-UniversityofSouthCarolina/PengLab/Bioinformatics/etnp18_16s/heatmap")# Below is for Optiplex


#Load the libraries
library(ggplot2); packageVersion("ggplot2")
library(edgeR)
library(statmod)
library(dplyr)
library(tidyr)
library(ggplot2)


#~~~~~~~~~~~~~ First run EdgeR~~~~~~~~~~~~~~~~~~

# Import the ASV table, metadata, and taxonomy
tbl <- read.csv("zotutab_NO_Bacillus.csv", sep = ",", header = TRUE)
metadata <- read.csv("metadata.csv", sep = ",", header = TRUE)
tax <- read.csv("taxonomy.csv", sep = ",", header = TRUE)


# Set up groups. Note that oxcline and anoxic samples are grouped under "Low".
group <- factor(metadata$Category_heatmap)

# Construct edgeR object; Column 1 is ASV ID.
y <- DGEList(tbl[,-1], group=group, genes=tbl[,c(1), drop = FALSE])

# Remove genes with low counts 
# NOTE: this section is not run for ETNP 16S data
#keep <- rowSums(cpm(y) > 0.1) >= 3
#table(keep) # Take a look at how many are left
#y <- y[keep, , keep.lib.sizes=FALSE]

# Normalization for composition bias
y <- calcNormFactors(y, method = "TMM"); y$samples
plotMDS(y) # (Optional visualization step. It should look like NMDS)

# Design matrix
design <- model.matrix(~0+group)
colnames(design) <- levels(group)
design

# Estimate Dispersion
y <- estimateDisp(y, design, robust=TRUE)
plotBCV(y) # Visualize dispersion estimates

# Estimate Quasi-Likelihood (QL) dispersions
fit <- glmQLFit(y, design, robust=TRUE)
head(fit$coefficients)
plotQLDisp(fit) #Visualize Quasi-Likelihood dispersions

# Differential expression analysis
con <- makeContrasts(LvsS = Low - Surface,
                     LvsD = Low - Deep,
                     SvsD = Surface - Deep, levels=design)
res <- glmQLFTest(fit, contrast=con) #QLFTest = Quasi-Likelihood F Test
topTags(res) # view the top DE ASVs

res.table <- topTags(res, n = "Inf")$table
colnames(res.table)[1] <- "ZOTU_ID"

# Count how many ASVs are differentially regulated
is.de <- decideTestsDGE(res)
summary(is.de)

# Visualize differentially regulated ASVs
plotMD(res, status=is.de, values=c(1,-1), col=c("red","blue"), legend="topright")

cpm <- as.data.frame(cpm(y, log=FALSE))
cpm$ZOTU_ID <- y$genes$X.OTU.ID
head(cpm)

#~~~~~~~~Combine CPM table with Taxonomy~~~~~~~~~~~~~~~~~~
# Separate the lineage into different levels and clean up the text
library(stringr)
tax_temp <- as.data.frame(str_split_fixed(tax$Taxon, ';', 7))
colnames(tax_temp) <- c('Domain','Phylum','Class','Order','Family','Genus','Species')
tax_temp$Domain <- gsub("d__","",as.character(tax_temp$Domain))
tax_temp$Phylum <- gsub("p__","",as.character(tax_temp$Phylum))
tax_temp$Class <- gsub("c__","",as.character(tax_temp$Class))
tax_temp$Order <- gsub("o__","",as.character(tax_temp$Order))
tax_temp$Family <- gsub("f__","",as.character(tax_temp$Family))
tax_temp$Genus <- gsub("g__","",as.character(tax_temp$Genus))
tax_temp$Species <- gsub("s__","",as.character(tax_temp$Species))
# Combine into one table
tax_final <- cbind(tax$Feature.ID, tax_temp, tax$Taxon)
colnames(tax_final)[1] <- 'ZOTU_ID'
colnames(tax_final)[9] <- 'Tax_All'

# Merge the EdgeR results with the taxonomy
cpm_tax <- merge(cpm, tax_final, by = "ZOTU_ID", all.x = TRUE)
cpm_tax_res <- merge(cpm_tax, res.table, by = "ZOTU_ID", all.x = TRUE)
# Remove all Chloroplast
cpm_tax_res <- cpm_tax_res[!grepl("Chloroplast", cpm_tax_res$Order),]

# Generate a new column that will be unique for each row
cpm_tax_res$ID_tax <- paste(cpm_tax_res$ZOTU_ID, cpm_tax_res$Tax_All, sep="_")

# Sum the cpm for all depths by groups (surface, Low (oxycline and anoxic), Deep)
all_sample_names <- metadata$Sample
cpm_tax_res$cpm_all <- rowSums(cpm_tax_res[,all_sample_names])

surface_names <- (filter(metadata, metadata$Category_heatmap == "Surface"))$Sample
cpm_tax_res$cpm_surface <- rowSums(cpm_tax_res[,surface_names])

low_names <- (filter(metadata, metadata$Category_heatmap == "Low"))$Sample
cpm_tax_res$cpm_low <- rowSums(cpm_tax_res[,low_names])

deep_names <- (filter(metadata, metadata$Category_heatmap == "Deep"))$Sample
cpm_tax_res$cpm_deep <- rowSums(cpm_tax_res[,deep_names])

#~~~~~~~~~~~~~Generate a Heatmap for most abundant ASVs with GGPLOT~~~~~~~~~~~~~~~~~~

N_ASV = 48; # Number of ASVs to be displayed on the heatmap

ASV_sig <- filter(cpm_tax_res, FDR < 0.05) # Select ASV only, with a threshold of 0.01
ASV_sorted <- ASV_sig[order(-ASV_sig$cpm_all), ] # Sort the table by sum of rpkm across all samples.

top_ASV <- ASV_sorted[1:N_ASV, ] # The number here determines how many ASVs are displayed

# The "small" df selects only the columns needed for pivoting
top_ASV_small <- top_ASV %>% select(starts_with("s", ignore.case = FALSE) | starts_with("ID"))
# Perform pivoting
top_ASV_long <- top_ASV_small %>% pivot_longer(cols=-c('ID_tax'), names_to = 'Sample', values_to = 'CPM')
# Perform log transformation for heat map
top_ASV_long$logCPM <- log10(top_ASV_long$CPM + 0.01)

# Order Samples manually, instead of hierarchical clustering 
Sample_order <- read.delim("Sample_order.txt",header = FALSE)
top_ASV_long$Sample <- factor(top_ASV_long$Sample, levels = as.matrix(Sample_order))

# Order genes by taxonomy
top_ASV_sorted <- top_ASV[with(top_ASV, order(cpm_surface, cpm_low, cpm_deep)), ]
gene_order <- read.delim("gene_order.txt",header = FALSE)
top_ASV_long$ID_tax <- factor(top_ASV_long$ID_tax, levels = as.matrix(gene_order))

write.csv(top_ASV_long, "top_ASV_long.csv")


ylabheat <- c("Woesearchaeales (Nanoarchaeota)",
              "Marine Group II (Thermoplasmatota)",
              "Marine Group III (Thermoplasmatota)", 
              "Crocinitomix", #Flavobacteriales
              "Crocinitomix", #Flavobacteriales
              "Crocinitomix", #Flavobacteriales
              "Flavobacteriaceae", 
              "Pseudofulvibacter",
              "NS5 marine group", 
              "NS5 marine group", 
              "Formosa", 
              "NS5 marine group", 
              "NS2b marine group", 
              "NS4 marine group",
              "NS5 marine group",
              "Formosa",
              "NS4 marine group",
              "NS4 marine group",
              "NS9 marine group",
              "NS9 marine group",
              "Synechococcus (Cyanobacteria)",
              "Desulfobacterota",
              "Marinimicrobia (SAR406 clade)", 
              "NB1-j",
              "OM190", #Planctomycetota
              "OM190", #Planctomycetota
              "Phycisphaerales",
              "Phycisphaeraceae", 
              "Phycisphaeraceae", 
              "Phycisphaeraceae", 
              "Pirellula", 
              "Rhodobacteraceae (Alphaproteobacteria)",
              "Curvibacter", #Gammaproteobacteria
              "Aquabacterium",
              "Comamonas", 
              "Methylomonadaceae", 
              "Woeseia", 
              "Stenotrophomonas", 
              "SAR86 clade",
              "Kiritimatiellaceae", #Verrucomicrobiota
              "Lentimonas", 
              "Lentimonas", 
              "Coraliomargarita", 
              "Lentimonas", 
              "Roseibacillus", 
              "Roseibacillus",
              "Roseibacillus",
              "Roseibacillus")

# Generate the heat map
p_ASV <- ggplot(top_ASV_long, aes(Sample, ID_tax)) +
  geom_tile(aes(fill = logCPM),
            color = "white", lwd = 1.5, linetype = 1)+
  scale_fill_gradient2(low = "#075AFF",
                       mid = "#FFFFCC",
                       high = "#FF0000") +
  scale_y_discrete(labels = ylabheat) + 
  theme(axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size = 14, vjust = 0.5, hjust=1),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.position = "top",
        legend.key.width = unit(2, "cm")
        
        #    legend.key.height = unit(7, "mm"))
  )
p_ASV

ggsave("heatmap_ASV.svg", plot=p_ASV, device="svg",
       scale=1, width = 40, height=30, units=c("cm"), limitsize = FALSE)
