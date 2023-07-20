#Figure 3 - (A – C) Distance triplot of redundancy analysis (RDA) on prokaryotic community composition in free-living, small particle-associated, and large particle-associated communities from the ETNP OMZ, using dissolved oxygen (O2), ammonium concentration (NH4+), nitrite concentration (NO2-), temperature, and chlorophyll a concentration (Chl a, approximated by fluorescence measured by a Seapoint chlorophyll fluorometer) as explanatory variables. The blue arrows are the vectors of the explanatory variables.  Circles represent samples from the coastal OMZ, triangles represent samples from the offshore OMZ, and squares represent samples from the peripheral OMZ. Each symbol was color shaded by the different oxygen regimes (surface/oxic (>150 M of O2), upper oxycline (10 nM – 150 M of O2), anoxic (1 – 10 nM of O2), and lower oxycline (10 nM – 100 M of O2). 
## Large volume size-fractioned filtration yields novel insights on prokaryotic community composition and interactions in the eastern tropical North Pacific ocean 
## Madeleine A. Thompson, David L. Valentine, Xuefeng Peng
## Code by Madeleine A. Thompson

#set working directory 
setwd("~/OneDrive - University of South Carolina/PengLab/Bioinformatics/etnp18_16s")

#load all packages needed
library(ggplot2)
library(vegan)
library(ggrepel)
library(dplyr)
library("viridis")
if (!requireNamespace("viridisLite", quietly = TRUE)) {
  install.packages("viridisLite")
}
library(viridisLite)

#load data 
zotu16s <- read.csv("data/zotutab.csv", header = TRUE, row.names = 1)
metadata <- read.csv("data/metadataETNP2018MT.csv", header = TRUE, row.names = 1)
metaplt <- read.csv("data/metadata.csv", header = TRUE, row.names = 1)

#normalize asv and metadata 

zotu16st <- t(zotu16s)
zotu.chord <- decostand(zotu16st, "normalize")

meta_sqrt <- sqrt(metadata)
## remove NANs
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

meta_sqrt[is.nan(meta_sqrt)] <- 0

meta_sqrt_normal <- decostand(meta_sqrt, "standardize")

#RDA
all.chord <- rda(zotu.chord ~ ., meta_sqrt_normal)
summary(all.chord)

#VIF scores
vif.cca(all.chord)

plot(all.chord)

#ANOVA stat
anova.cca(all.chord, step = 1000, by = "term")

#basic plot
perc <- round(100*(summary(all.chord)$cont$importance[2, 1:2]), 2)

sc_si <- scores(all.chord, display="sites", choices=c(1,2), scaling=1)

sc_bp <- scores(all.chord, display="bp", choices=c(1, 2), scaling=1)

plot(all.chord,
     scaling = 1, # set scaling type 
     type = "none", # this excludes the plotting of any points from the results
     frame = FALSE,
     # set axis limits
     xlim = c(-.3,.3), 
     ylim = c(-.3,.3),
     # label the plot (title, and axes)
     main = "Triplot RDA - scaling 1",
     xlab = paste0("RDA1 (", perc[1], "%)"), 
     ylab = paste0("RDA2 (", perc[2], "%)") 
)
# add points for site scores
points(sc_si, 
       pch = 21, # set shape (here, circle with a fill colour)
       col = "black", # outline colour
       bg = "steelblue", # fill colour
       cex = 1.2) # size
# add arrows for effects of the expanatory variables
arrows(0,0, # start them from (0,0)
       sc_bp[,1], sc_bp[,2], # end them at the score value
       col = "black", 
       lwd = 2)
# add text labels for arrows
text(x = sc_bp[,1] -0.1, # adjust text coordinate to avoid overlap with arrow tip
     y = sc_bp[,2] - 0.03, 
     labels = rownames(sc_bp), 
     col = "black", 
     cex = 1, 
     font = 2)

# plot in ggplot with each size fraction as a different RDA 
sc_bp <- as.data.frame(sc_bp)
zotu_rda = merge(sc_si, metaplt, by = 0)
print(sapply(zotu_rda, class)) 

zotu_rda$Station <- as.factor(zotu_rda$Station)
zotu_rda$OxygenGroup <- factor(zotu_rda$OxygenGroup)

rda.filt.22 <- zotu_rda %>%
  filter(Filter == 22)
rda.filt.2 <- zotu_rda %>%
  filter(Filter == 2)
rda.filt.02 <- zotu_rda %>%
  filter(Filter == 0.2)

zotu.rda22 <- ggplot2::ggplot(data=rda.filt.22, aes(x=RDA1, y=RDA2)) +
  ggplot2::geom_point(aes(x=RDA1, y=RDA2, color=OxygenGroup, shape = Station), size = 5) +
  ggplot2::geom_segment(data = sc_bp, ggplot2::aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
                        size =0.5, alpha = 1, colour = "blue", arrow = arrow(length = unit(0.3,"cm"))) + 
  ggplot2::theme(axis.title = element_text(size = 17, colour = "black"), 
                 panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "black"), 
                 axis.text.x= element_text(color = "black", size = 15), axis.text.y= element_text(color = "black", size = 15), legend.key = element_blank(), 
                 legend.title = element_text(size = 17, colour = "black"), 
                 legend.text = element_text(size = 17, colour = "black")) + 
  geom_text_repel(data = sc_bp, 
                  aes(RDA1*1.06, RDA2 *1.05,
                      label=rownames(sc_bp)), size = 5) +
  geom_hline(yintercept=0, linetype="dashed", color = "grey") +
  geom_vline(xintercept=0, linetype="dashed", color = "grey") +
  xlim(-0.3, 0.35) +
  ylim(-0.3,0.3) +
  labs(colour = "Oxygen Group") + 
  xlab("RDA1 (8.47%)") + ylab("RDA2 (6.72%)") + 
  scale_size_area(max_size = 5) +
  theme(text = element_text(size = 15)) + 
  scale_color_viridis_d(option = "C")+
  ggtitle(expression(paste("Large Particle-associated")))
zotu.rda22

zotu.rda2 <- ggplot2::ggplot(data=rda.filt.2, aes(x=RDA1, y=RDA2)) +
  ggplot2::geom_point(aes(x=RDA1, y=RDA2, color=OxygenGroup, shape = Station), size = 5) +
  ggplot2::geom_segment(data = sc_bp, ggplot2::aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
                        size =0.5, alpha = 1, colour = "blue", arrow = arrow(length = unit(0.3,"cm"))) + 
  ggplot2::theme(axis.title = element_text(size = 17, colour = "black"), 
                 panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "black"), 
                 axis.text.x= element_text(color = "black", size = 15), axis.text.y= element_text(color = "black", size = 15), legend.key = element_blank(), 
                 legend.title = element_text(size = 17, colour = "black"), 
                 legend.text = element_text(size = 17, colour = "black")) + 
  geom_text_repel(data = sc_bp, 
                  aes(RDA1*1.06, RDA2 *1.05,
                      label=rownames(sc_bp)), size = 5) +
  geom_hline(yintercept=0, linetype="dashed", color = "grey") +
  geom_vline(xintercept=0, linetype="dashed", color = "grey") +
  xlim(-0.3, 0.35) +
  ylim(-0.3,0.3) +
  labs(colour = "Oxygen Group") + 
  xlab("RDA1 (8.47%)") + ylab("RDA2 (6.72%)") + 
  scale_size_area(max_size = 5) +
  theme(text = element_text(size = 15)) + 
  scale_color_viridis_d(option = "C")+
  ggtitle(expression(paste("Small Particle-associated")))
zotu.rda2

zotu.rda02 <- ggplot2::ggplot(data=rda.filt.02, aes(x=RDA1, y=RDA2)) +
  ggplot2::geom_point(aes(x=RDA1, y=RDA2, color=OxygenGroup, shape = Station), size = 5) +
  ggplot2::geom_segment(data = sc_bp, ggplot2::aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
                        size =0.5, alpha = 1, colour = "blue", arrow = arrow(length = unit(0.3,"cm"))) + 
  ggplot2::theme(axis.title = element_text(size = 17, colour = "black"), 
                 panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "black"), 
                 axis.text.x= element_text(color = "black", size = 15), axis.text.y= element_text(color = "black", size = 15), legend.key = element_blank(), 
                 legend.title = element_text(size = 17, colour = "black"), 
                 legend.text = element_text(size = 17, colour = "black")) + 
  geom_text_repel(data = sc_bp, 
                  aes(RDA1*1.06, RDA2 *1.05,
                      label=rownames(sc_bp)), size = 5) +
  geom_hline(yintercept=0, linetype="dashed", color = "grey") +
  geom_vline(xintercept=0, linetype="dashed", color = "grey") +
  xlim(-0.3, 0.35) +
  ylim(-0.3,0.3) +
  labs(colour = "Oxygen Group") + 
  xlab("RDA1 (8.47%)") + ylab("RDA2 (6.72%)") + 
  scale_size_area(max_size = 5) +
  theme(text = element_text(size = 15)) + 
  scale_color_viridis_d(option = "C")+
  ggtitle(expression(paste("Free-living"))) 
zotu.rda02


rdafig <- ggpubr::ggarrange(zotu.rda02, zotu.rda2, zotu.rda22, 
                            labels = c("A", "B", "C"), 
                            font.label = list(size = 15, color = "black", face = "bold", family = NULL),
                            ncol = 3, nrow = 1,
                            common.legend = TRUE, legend = c("right"))

ggplot2::ggsave("figuresNew/RDAfilt.png", plot=rdafig, device="png",
                scale=1, width = 35, height=10, units=c("cm"), dpi=300, limitsize = FALSE)

#save as svg and edit in inkscape (i.e. the legend)

