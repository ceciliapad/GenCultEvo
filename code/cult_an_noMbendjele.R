################################
## Analyses cultural networks ##
################################

# DECEMBER 2021 - FULL DATASETS

setwd("/Users/Cecilia/Documents/PhD/African_Pygmies/Spatial analysis culture/My_data/Scripts")
#source utility functions - these now includes ones from Science Advances paper
source("ESC_S2_utilityfunctions.R")
#save.image("SpatialCulture.RData")
require(MetBrewer)

setwd("/Users/Cecilia/Documents/PhD/African_Pygmies/Spatial analysis culture/My_data/Nov_2021")
load("SpatialCulture.RData")

#require(wesanderson)
my_pal <- c(MetPalettes$Isfahan1[[1]], MetPalettes$Isfahan2[[1]][2])

#load data (the excel worksheets should be individually exported as .csv files)

subsistence <- read.csv("P_A_subsistence2.csv", row.names = 1)
music <- read.csv("P_A_music2.csv", row.names = 1)

# Counts per group
subsistence2 <- subsistence
music2 <- music
s <- rowSums(subsistence2, na.rm=T)
m <- rowSums(music2, na.rm=T)
write.csv(s, "results_masked/subsistence_counts.csv")
write.csv(m, "results_masked/music_counts.csv")

#pal <- wes_palette("Cavalcanti1", 11, type = "continuous")
m <- as.data.frame(m)
s <- as.data.frame(s)
m$FID <- rownames(m)
s$FID <- rownames(s)

library(tidyverse)
m <- m %>% arrange(FID)
s <- s %>% arrange(FID)

IDS <-  s$FID
IDS[6] <- "Batwa E"
IDS[7] <- "Batwa W"

# Plot

my_pal_long <- c(MetPalettes$Isfahan1[[1]][1:6],MetPalettes$Isfahan2[[1]][1], MetPalettes$Isfahan1[[1]][7], MetPalettes$Isfahan2[[1]][c(2,3,5)])

png("results_masked/Musical_instruments.png",units="in", width=7, height=5, res=500)
barplot(m$m, col=my_pal_long, names.arg=IDS, ylab="No. of different musical instruments", ylim=c(0,30), cex.names=0.5, cex.axis=0.5 )
dev.off()
png("results_masked/Subsistence_tools.png",units="in", width=7, height=5, res=500)
barplot(s$s, col=my_pal_long, names.arg=IDS, ylab="No. of different subsistence tools", ylim=c(0,30), cex.names=0.5, cex.axis=0.5 )
dev.off()

# Is there a correlation
all_reps <- cbind(m,s)
cor.test(all_reps$m, all_reps$s, method="spearman")

# MASKED HETEROZYGOSITY
# Copy runs_homo.csv into Nov_2021 directory and alter the names to match
# we assume Efe = Sua to be same genetic population and Aka + Mbendjele too
# No genetic data for Batwa (West)

#rohs <- read.csv("runs_homo.csv")

# Make overall diversity df
#div_all <-  merge(rohs,m,by="FID")
#div_all <-  merge(div_all,s,by="FID")

#cor.test(div_all$mean_hom, div_all$m, method="spearman")
#cor.test(div_all$median_hom, div_all$m, method="spearman")
#cor.test(div_all$mean_hom, div_all$s, method="spearman")
#cor.test(div_all$median_hom, div_all$s, method="spearman")

#x <- lm(div_all$m~log(div_all$median_hom))
#t <- lm(div_all$s~log(div_all$median_hom))

# Heterozygosity (total)
het <- read.csv("heterozygosity_masked.csv")
het_all <-  merge(het,m,by="FID")
het_all <-  merge(het_all,s,by="FID")

# No correlation betweem je
cor.test(het_all$Heterozygosity, het_all$s, method="spearman")
cor.test(het_all$Heterozygosity, het_all$m, method="spearman")

# load locations data
locations <- read.csv("location_culture_centroids2.csv", row.names = 1)
locations$FID <- rownames(locations)
area_all <- merge(locations,m,by="FID")
area_all <- merge(area_all,s,by="FID")
cor.test(area_all$Area, area_all$s, method="spearman") 
cor.test(area_all$Area, area_all$m, method="spearman")

# load biome data and /100 to get proportion
biome <- (read.csv("biomes2.csv", row.names = 1))/100

# Maybe is not the % biome that dictates whether you have an object but simply whether
# that biome is present
biome_bin <- read.csv("biomes2.csv", row.names = 1)
biome_bin <- ifelse(biome_bin>0,1,0)

# Read in genetic distances and turn into distance object
Fst_mat<- read.csv("fst_masked.csv", row.names = 1)[-9,-9]
colnames(Fst_mat) <- rownames(Fst_mat)
genDist <- as.dist(Fst_mat)

Fst_mat_IBD<- read.csv("fst_masked_ibd.csv", row.names = 1)[-9,-9]
colnames(Fst_mat_IBD) <- rownames(Fst_mat_IBD)
genDistIBD <- as.dist(Fst_mat_IBD)

Fst_mat_Bantu<- read.csv("fst_bantu.csv", row.names = 1)[-9,-9]
colnames(Fst_mat_Bantu) <- rownames(Fst_mat_Bantu)
genDistBantu <- as.dist(Fst_mat_Bantu)

# Read in sampling_effort file

samples <- read.csv("sampling_effort.csv", row.names = 1)

samples[is.na(samples)] <- 0
samples <- samples %>% mutate_if(is.numeric, ~1 * (. != 0))

samples_group <- as.data.frame(rowSums(samples))
samples_group$FID <- rownames(samples)

write.csv(samples_group, "results_masked/sources_per_group.csv")
write.csv(samples, "results_masked/sample_list.csv")
samples_all_info_het <- merge(het_all,samples_group,by="FID")
samples_all_info_hom <- merge(div_all,samples_group,by="FID")
samples_all_info_area <- merge(area_all,samples_group,by="FID")

reg_hom_music <- lm(samples_all_info_hom$m~log(samples_all_info_hom$`rowSums(samples)`)+samples_all_info_hom$mean_hom)
reg_hom_subsistence <- lm(samples_all_info_hom$s~log(samples_all_info_hom$`rowSums(samples)`)+samples_all_info_hom$mean_hom)

reg_area_subsistence <- lm(samples_all_info_area$s~log(samples_all_info_area$`rowSums(samples)`)+log(samples_all_info_area$Area))
reg_area_music <- lm(samples_all_info_area$m~log(samples_all_info_area$`rowSums(samples)`)+log(samples_all_info_area$Area))

reg_area_subsistence <- lm(samples_all_info_area$s~samples_all_info_area$Area)


#samples <- do.call(cbind,lapply(split(as.list(samples),names(samples)),function(x) Reduce(`+`,x)));
####################
# Data Preparation #
####################

#compute spatial distances
distTable <- data.frame(x=locations$Longitude,y=locations$Latitude)
names_locs <- c("Aka", "Babongo", "Baka", "Bakola", "Bakoya", "Batwa (East)", "Batwa (West)", "Bedzan", "Efe","Sua" )
rownames(distTable) <- names_locs
space<-distMat(data.frame(x=locations$Longitude,y=locations$Latitude))

#library(geosphere)
#spaceSubsistence <- distm(distTable, fun=distGeo)
#rownames(spaceSubsistence) <- names
#colnames(spaceSubsistence) <- names
# here distance is in metres, not kilometres
#spaceSubsistence <- as.dist(spaceSubsistence)

#compute gower dissimilarity between biome composition
library(vegan)
gowBiome<- vegdist(biome,"gower",na.rm=TRUE)
gowBiomeBin<- vegdist(biome_bin,"gower",na.rm=TRUE)

#compute Jaccard distances
#index is bounded between 0 (identical presences in the two sites) and 1 (complete absence of shared traits)
#ignores instances of negative matches (i.e. shared absences) that are common in sparse matrices with many absences (Shennan & Bentley, 2008).
jacSubsistence<-vegdist(subsistence,"jaccard",na.rm=TRUE)
jacMusic<-vegdist(music,"jaccard",na.rm=TRUE)
t.test(jacSubsistence, jacMusic)


# Plot Jaccard distances (higher values = MORE DISSIMILAR!!)

dSub <- as.matrix(jacSubsistence)
dim <- ncol(dSub)

png(file="results_masked/Subsistence_Jaccard.png",
    width=7, height=6, units="in", res=500)
image(1:dim, 1:dim, dSub, axes = FALSE, xlab="", ylab="", col=hcl.colors(12, "Geyser", rev=TRUE), main="Subsistence tools")
axis(1, 1:dim, names_locs, cex.axis = 0.6, las=3)
axis(2, 1:dim, names_locs, cex.axis = 0.6, las=1)
text(expand.grid(1:dim, 1:dim), sprintf("%0.2f", dSub), cex=0.6)
dev.off()

dMus <- as.matrix(jacMusic)
dim <- ncol(dMus)

png(file="results_masked/Music_Jaccard.png",
    width=7, height=6, units="in", res=500)
image(1:dim, 1:dim, dMus, axes = FALSE, xlab="", ylab="", col=hcl.colors(12, "Geyser", rev=TRUE), main="Musical Instruments")
axis(1, 1:dim, names_locs, cex.axis = 0.6, las=3)
axis(2, 1:dim, names_locs, cex.axis = 0.6, las=1)
text(expand.grid(1:dim, 1:dim), sprintf("%0.2f", dMus), cex=0.6)
dev.off()

library(qgraph) # qgraph takes similarity matrices as input
library(wesanderson)
Mus_d <- 1/dMus
png(file="results_masked/Music_Jaccard_net.png",
    width=6, height=6, units="in", res=500)
qgraph(Mus_d, layout='spring', vsize=5, label.cex=1.5, esize=14, edge.color = MetPalettes$Isfahan1[[1]][5])
title("Musical instruments", line = 2.5)
dev.off()

Sub_d <- 1/dSub
png(file="results_masked/Tools_Jaccard_net.png",
    width=6, height=6, units="in", res=500)
qgraph(dSub, layout='spring', vsize=5, label.cex=1.5, edge.color = MetPalettes$Morgenstern[[1]][6])
title("Subsistence tools", line = 2.5)
dev.off()

# USE Plotting_NetInMaps.R script to plot Jaccard in Maps

###################
## Data prep 2.0 ##
###################

# Create reduced distance objects for subsequent analyses without Batwa West
locations_red <- locations[-7,]
biome_red <- biome[-7,]
biome_bin_red <- biome_bin[-7,]
subsistence_red <- subsistence[-7,]
music_red <- music[-7,]

locations_red2 <- locations_red[-c(8:9),]
biome_red2 <- biome_red[-c(8:9),]
biome_bin_red2 <- biome_bin_red[-c(8:9),]
subsistence_red2 <- subsistence_red[-c(8:9),]
music_red2 <- music_red[-c(8:9),]

distTable2 <- data.frame(x=locations_red$Longitude,y=locations_red$Latitude)
rownames(distTable2) <- names_locs[-7]
spaceRed<-distMat(data.frame(x=locations_red$Longitude,y=locations_red$Latitude))

distTable3 <- data.frame(x=locations_red2$Longitude,y=locations_red2$Latitude)
rownames(distTable3) <- rownames(distTable2)[-(8:9)]
spaceRed2<-distMat(data.frame(x=locations_red2$Longitude,y=locations_red2$Latitude))

gowBiomeRed<- vegdist(biome_red,"gower",na.rm=TRUE)
gowBiomeRedBin<- vegdist(biome_bin_red,"gower",na.rm=TRUE)

jacSubsistencered<-vegdist(subsistence_red,"jaccard",na.rm=TRUE)
jacMusicred<-vegdist(music_red,"jaccard",na.rm=TRUE)

gowBiomeRed2<- vegdist(biome_red2,"gower",na.rm=TRUE)
gowBiomeRedBin2<- vegdist(biome_bin_red2,"gower",na.rm=TRUE)

jacSubsistencered2<-vegdist(subsistence_red2,"jaccard",na.rm=TRUE)
jacMusicred2<-vegdist(music_red2,"jaccard",na.rm=TRUE)

# ####################################
# ## Correlations between distances ##
# ####################################
# 
# #make sure reshape package is unloaded
# # plot correlations between predictors
# #
# #get upper triangle of the correlation matrix
# #(source: http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization)
# get_upper_tri <- function(cormat){
#   cormat[lower.tri(cormat)]<- NA
#   return(cormat)
# }
# 
# dr  = d %>%
#   dplyr::select(genDist, spaceRed, gowBiomeRed)
# 
# cormat = dr %>%
#   cor(.) %>%
#   round(., 2) %>%
#   get_upper_tri(.) %>%
#   reshape2::melt(.)
# 
# cor_vars = c("Genetic distance", "Geographical distance", "Biome dissimilarity")
# 
# # calculate correlations and p-values
# cormats = dr %>%
#   as.matrix(.) %>%
#   Hmisc::rcorr(.)
# 
# # format p-values
# cormat.p = cormats$P %>%
#   round(., 2) %>%
#   get_upper_tri(.) %>%
#   reshape2::melt(.)
# 
# # format correlations and merge with p-values
# cormat = cormats$r %>%
#   round(., 2) %>%
#   get_upper_tri(.) %>%
#   reshape2::melt(.) %>%
#   rename(r = value) %>%
#   left_join(cormat.p, by=c("Var1", "Var2")) %>%
#   rename(p = value)
# 
# # rename variables
# cormat$Var1 = factor(cormat$Var1)
# levels(cormat$Var1) = cor_vars
# cormat$Var2 = factor(cormat$Var2)
# levels(cormat$Var2) = cor_vars
# 
# # point out where p < 0.01 (rather than show it as p = 0)
# #cormat$p = as.character(cormat$p)
# #cormat$p = ifelse(cormat$p == "0", "< 0.01", cormat$p)
# 
# cormat$p = as.numeric(cormat$p)
# 
# corrplot <- ggplot(data = cormat, aes(Var2, Var1, fill = p)) +
#   geom_tile(colour = "white") +
#   geom_text(aes(label = r), size=3) +
# 
#   scale_fill_gradient2(low = wes_palette("Darjeeling1")[5], high = wes_palette("Darjeeling1")[1], mid = "white", na.value=wes_palette("Moonrise1")[3],
#                        # midpoint = 0, limit = c(-1,1), space = "Lab",
#                        # name="Pearson\nCorrelation") +
#                        midpoint = 0.05, limit = c(0, 0.15), space = "Lab",
#                        name="P-value") +
# 
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1,
#                                    size = 12, hjust = 1)) +
#   coord_fixed() +
#   xlab("") +
#   ylab("")
# 
# 
# ggsave(corrplot, filename="predictors correlations.png", height=200, width=200, units="mm")

########################
# Partial Mantel Tests #
########################

library(ecodist)
library(geosphere)
#library(vegan)
detach("package:vegan", unload = TRUE)
#partial mantel test between Jaccard Distance and Geography
#pval3 is the one we care about - Null: p = 0
# need to unload vegan
#jacSubsistence2 <- as.data.frame(jacSubsistence)

jacMusic <- as.matrix(jacMusic)
jacMusicred<- as.matrix(jacMusicred)
jacSubsistence<- as.matrix(jacSubsistence)
jacSubsistencered<- as.matrix(jacSubsistencered)
gowBiome<- as.matrix(gowBiome)
gowBiomeBin<- as.matrix(gowBiomeBin)
gowBiomeRed<- as.matrix(gowBiomeRed)
gowBiomeRedBin<- as.matrix(gowBiomeRedBin)
genDist<- as.matrix(genDist)
genDistBantu<- as.matrix(genDistBantu)[-8,-8]
space<- as.matrix(space)
spaceRed<- as.matrix(spaceRed)
spaceRed2<- as.matrix(spaceRed2)

sub_sp_mantel <- mantel(lower(jacSubsistence)~lower(space),nperm=1000, mrank=T)
mus_sp_mantel <-mantel(lower(jacMusic)~lower(space),nperm=1000, mrank=T)

sub_bio_mantel <- mantel(lower(jacSubsistence)~lower(gowBiome),nperm=1000, mrank=T)
mus_bio_mantel <- mantel(lower(jacMusic)~lower(gowBiome),nperm=1000, mrank=T)

sub_bio_mantel_mega <- mantel(lower(jacSubsistence)~lower(gowBiome),nperm=1000, mrank=T)
mus_bio_mantel_mega <- mantel(lower(jacMusic)~lower(gowBiome),nperm=1000, mrank=T)

sub_biobin_mantel <- mantel(lower(jacSubsistence)~lower(gowBiomeBin),nperm=1000, mrank=T)
mus_biobin_mantel <- mantel(lower(jacMusic)~lower(gowBiomeBin),nperm=1000, mrank=T)

sub_gen_mantel <- mantel(lower(jacSubsistencered)~lower(genDist),nperm=1000, mrank=T)
mus_gen_mantel <- mantel(lower(jacMusicred)~lower(genDist),nperm=1000, mrank=T)

sub_gen_ibd_mantel <- mantel(lower(jacSubsistencered)~lower(genDistIBD),nperm=1000, mrank=T)
mus_gen_ibd_mantel <- mantel(lower(jacMusicred)~lower(genDistIBD),nperm=1000, mrank=T)

sub_gen_bantu_mantel <- mantel(lower(jacSubsistencered2)~lower(genDistBantu),nperm=1000, mrank=T)
mus_gen_bantu_mantel <- mantel(lower(jacMusicred2)~lower(genDistBantu),nperm=1000, mrank=T)

# Are biomes and geographical distances related?

sp_bio_mantel <- mantel(lower(space)~lower(gowBiome),nperm=1000, mrank=T)

# Are genetic distances and geographical distances related? 

sp_gen_mantel <- mantel(lower(spaceRed)~lower(genDist),nperm=1000, mrank=T)
sp_gen_ibd_mantel <- mantel(lower(spaceRed)~lower(genDistIBD),nperm=1000, mrank=T)


# Are genetic distances and biome distances related? 
biobin_gen_mantel <- mantel(lower(gowBiomeRedBin)~lower(genDist),nperm=1000, mrank=T)
bio_gen_mantel <- mantel(lower(gowBiomeRed)~lower(genDist),nperm=1000, mrank=T)

# Effect of genetics when controlling for space & biome and other way around
mus_gen_sp_mantel <- mantel(lower(jacMusicred)~lower(genDist)+lower(spaceRed),nperm=1000, mrank=T)
mus_gen_bio_mantel <- mantel(lower(jacMusicred)~lower(genDist)+lower(gowBiomeRed),nperm=1000, mrank=T)

mus_genIBD_sp_mantel <- mantel(lower(jacMusicred)~lower(genDistIBD)+lower(spaceRed),nperm=1000, mrank=T)
mus_genIBD_bio_mantel <- mantel(lower(jacMusicred)~lower(genDistIBD)+lower(gowBiomeRed),nperm=1000, mrank=T)

mus_sp_gen_mantel <- mantel(lower(jacMusicred)~lower(spaceRed)+lower(genDist),nperm=1000, mrank=T)
mus_sp_bio_mantel<- mantel(lower(jacMusic)~lower(space)+lower(gowBiome),nperm=1000, mrank=T)

mus_bio_gen_mantel <- mantel(lower(jacMusicred)~lower(gowBiomeRed)+lower(genDist),nperm=1000, mrank=T)
mus_bio_sp_mantel <-  mantel(lower(jacMusic)~lower(gowBiome)+lower(space),nperm=1000, mrank=T)

sub_gen_sp_mantel <- mantel(lower(jacSubsistencered)~lower(genDist)+lower(spaceRed),nperm=1000, mrank=T)
sub_gen_bio_mantel <- mantel(lower(jacSubsistencered)~lower(genDist)+lower(gowBiomeRed),nperm=1000, mrank=T)

sub_genIBD_sp_mantel <- mantel(lower(jacSubsistencered)~lower(genDistIBD)+lower(spaceRed),nperm=1000, mrank=T)
sub_genIBD_bio_mantel <- mantel(lower(jacSubsistencered)~lower(genDistIBD)+lower(gowBiomeRed),nperm=1000, mrank=T)

sub_sp_gen_mantel <- mantel(lower(jacSubsistencered)~lower(spaceRed)+lower(genDist),nperm=1000, mrank=T)
sub_sp_bio_mantel <- mantel(lower(jacSubsistence)~lower(space)+lower(gowBiome),nperm=1000, mrank=T)

sub_bio_gen_mantel <- mantel(lower(jacSubsistencered)~lower(gowBiomeRed)+lower(genDist),nperm=1000, mrank=T)
sub_bio_sp_mantel <-mantel(lower(jacSubsistence)~lower(gowBiome)+lower(space),nperm=1000, mrank=T)

sub_bio_gen_mantel <- mantel(lower(jacSubsistencered)~lower(gowBiomeRedBin)+lower(genDist),nperm=1000, mrank=T)
sub_bio_sp_mantel <-mantel(lower(jacSubsistence)~lower(gowBiomeBin)+lower(space),nperm=1000, mrank=T)

gen_sp_bio_mantel <-mantel(lower(genDist)~lower(spaceRed)+lower(gowBiomeRed),nperm=1000, mrank=T)
gen_bio_space_mantel <-mantel(lower(genDist)~lower(gowBiomeRed)+lower(spaceRed),nperm=1000, mrank=T)


subsistence_mantel <-  rbind(sub_sp_mantel, sub_bio_mantel, sub_gen_mantel, sub_biobin_mantel, sub_gen_ibd_mantel, sub_gen_bantu_mantel)
music_mantel <-  rbind(mus_sp_mantel, mus_bio_mantel, mus_gen_mantel, mus_biobin_mantel, mus_gen_ibd_mantel, mus_gen_bantu_mantel)                         

pval_adj_subs <- p.adjust(subsistence_mantel[,2], method = "BH")
subsistence_mantel <- cbind(subsistence_mantel, pval_adj_subs)
pval_adj_music <- p.adjust(music_mantel[,2], method = "BH")
music_mantel <- cbind(music_mantel, pval_adj_music)

write.csv(rbind(music_mantel, subsistence_mantel), "final_paper_files/results_final_paper/results_masked/simple_mantel_repertoires_with_bantu.csv")

subsistence_matrix <- rbind(sub_gen_sp_mantel, sub_sp_gen_mantel,sub_bio_gen_mantel,sub_bio_sp_mantel, sub_gen_bio_mantel, sub_sp_bio_mantel)
music_matrix <- rbind(mus_gen_sp_mantel, mus_sp_gen_mantel,mus_bio_gen_mantel,mus_bio_sp_mantel, mus_gen_bio_mantel, mus_sp_bio_mantel)

pval_adj_subs <- p.adjust(subsistence_matrix[,2], method = "BH")
subsistence_matrix <- cbind(subsistence_matrix, pval_adj_subs)
pval_adj_mus <- p.adjust(music_matrix[,2], method = "BH")
music_matrix<- cbind(music_matrix, pval_adj_mus)

write.csv(rbind(music_matrix, subsistence_matrix), "results_masked/m_matrix_repertoires.csv")


#save.image("SpatialCulture.RData")
#table_res[,1:4]

# Plot correlations
reps_df <- read.csv("repertoire_mantel_plot.csv")
c <- ggplot(data = reps_df, aes(response, explanatory, fill = r2_adj))

#####################################################################
## Principal Coordinates Analysis (PCoA) of Tools, Music and Genes ##
#####################################################################

#We perform a principal coordinate analysis (PCoA) on the distance matrices for genetics and music.
#Similar to a PCA, a PCoA produces a set of orthogonal axes whose importance is measured by eigenvalues.
#However, in contrast to the PCA, non-Euclidean distance matrices can be used.
#We correct for negative eigenvalues using the Cailliez procedure.

# Perform PCoA
library(ape)
subsistence_pcoa <- pcoa(jacSubsistence, correction = "cailliez")
genetics_pcoa <- pcoa(genDist, correction = "cailliez")
music_pcoa <- pcoa(jacMusic, correction = "cailliez") 


#We rescale the PCoA components in relation to the explained variance.  

# Genetics does not have negative eigenvalues, no correction needed
for(i in 1:ncol(genetics_pcoa$vectors)) { 
  genetics_pcoa$vectors[,i] <- scale(genetics_pcoa$vectors[,i])*
    genetics_pcoa$values$Rel_corr_eig[i]}

# Subsistence
for(i in 1:ncol(subsistence_pcoa$vectors)) { 
  subsistence_pcoa$vectors[,i] <- scale(subsistence_pcoa$vectors[,i])*
    subsistence_pcoa$values$Relative_eig[i]}

# Music no negative eigenvalues
for(i in 1:ncol(music_pcoa$vectors)) { 
  music_pcoa$vectors[,i] <- scale(music_pcoa$vectors[,i])*
    music_pcoa$values$Relative_eig[i]}

# Extract eigenvalues from PC/PCo results_masked 
genetics_ev <- genetics_pcoa$values$Corr_eig
music_ev <- music_pcoa$values$Eigenvalues
subsistence_ev <- subsistence_pcoa$values$Eigenvalues

# Plot
my_font <- "sans"

png(file="results_masked/pcoaGenetics.png",
    width=6, height=4.5, units="in", res=500)
plot_pc_scree(eigenval=genetics_ev, title="Genetics", type="PCo")
dev.off()
png(file="results_masked/pcoaMusic.png",
    width=6, height=4.5, units="in", res=500)
plot_pc_scree(eigenval=music_ev, title="Musical instruments", type="PCo")
dev.off()
png(file="results_masked/pcoaSubsistence.png",
    width=6, height=4.5, units="in", res=500)
plot_pc_scree(eigenval=subsistence_ev, title="Subsistence tools", type="PCo")
dev.off()

##We extract the principal components/coordinates for all factors: 
#genetics, grammar, music and phonology. 
#We normalize the PCs/PCos to a range from 0 to 1 and then plot a heatmap 
library(reshape2)
# Extract the PCs and the PCos from the PCA and PCoA results_masked
genetics_pco <- genetics_pcoa$vectors
music_pco <- music_pcoa$vectors
subsistence_pco <- subsistence_pcoa$vectors

# Change the column names of all PCs and PCoAs  
colnames(genetics_pco)<- paste("genetics_pco_", seq(1,ncol(genetics_pco)), sep="")
colnames(music_pco) <- paste("music_pco_", seq(1,ncol(music_pco)), sep="")
colnames(subsistence_pco) <- paste("subsistence_pco_", seq(1,ncol(subsistence_pco)), sep="")

png(file="results_masked/HEATpcoaGenetics.png",
    width=10, height=7, units="in", res=500)
heat_map_pc(genetics_pco, type="PCo")
dev.off()

png(file="results_masked/HEATpcoaMusic.png",
    width=10, height=7, units="in", res=500)
heat_map_pc(music_pco, type="PCo")
dev.off()

png(file="results_masked/HEATpcoaSubsistence.png",
    width=10, height=7, units="in", res=500)
heat_map_pc(subsistence_pco, type="PCo")
dev.off()

#########################################
# Distance visualization (NeighborNets) #
#########################################

#We compute Euclidean distance matrices from the dimensionality-reduced data and create NeighborNets. 
# SplitsTree doesnt accept brackets in names
genetics_pco_rn <- genetics_pco
rownames(genetics_pco_rn)[6] <- "BatwaE"

music_pco_rn <- music_pco
rownames(music_pco_rn)[6:7] <- c("BatwaE", "BatwaW")

subsistence_pco_rn <- subsistence_pco
rownames(subsistence_pco_rn)[6:7] <- c("BatwaE", "BatwaW")


dist_list <- lapply(list(genetics=genetics_pco_rn, 
                         music=music_pco_rn, 
                         subsistence=subsistence_pco_rn),
                    function(f) {dist(f)})

# Write nexus files for splitstree
require(phangorn)
writeDist(dist_list$subsistence, file="results_masked/subsistence.nex", format="nexus")
writeDist(dist_list$music, file="results_masked/music.nex", format="nexus")
writeDist(dist_list$genetics, file="results_masked/genes.nex", format="nexus")

# The plotting you do with SplitsTree and Illustrator.

# Single trait networks

# visualise correlations
#Change NAs to blanks for text
PhiMatrix_text <- as.data.frame(PhiMatrix)
PhiMatrix_text <- sapply(PhiMatrix_text, as.character)
PhiMatrix_text[is.na(PhiMatrix_text)] <- " "
PhiMatrix_text <- as.matrix(PhiMatrix_text)

# reps_df <- read.csv("repertoire_simple_plot.csv", row.names = 1)
# MetPalettes$Isfahan2[[1]][2]
# 
# pvals <- as.matrix(reps_df[3:4,])
# cors <- as.matrix(reps_df[1:2,])
# 
# names_y <- c("Geography", "Genetics", "Ecology")
# names_x <- c("Music", "Tools")
# 
# dim_x <- length(names_x)
# dim_y <- length(names_y)
# 
# png(file="Reps_mantel_or.png",
#     width=7, height=6, units="in", res=2500)
# image(1:dim_x, 1:dim_y, cors, axes = FALSE, xlab="", ylab="", col=hcl.colors(12, "Oranges", rev=TRUE),  xlim=c(0,10), ylim=c(0,10))
# axis(3, 1:dim_x, names_x, cex.axis = 0.6, las=1, pos=3.3, col="transparent", col.ticks="transparent")
# axis(2, 1:dim_y, names_y, cex.axis = 0.6, las=1, pos=0.5, col="transparent", col.ticks="transparent")
# text(expand.grid(1:dim_x, 1:dim_y), sprintf("%0.4s", cors), cex=0.6)
# dev.off()

###########
## WORDS ##
###########

subsistence_words <- read.csv("subsistence_words.csv", row.names = 1)
music_words <- read.csv("instrument_words2.csv", row.names = 1)

ncol(subsistence_words)
ncol(music_words)

s_w <- as.data.frame(rowSums(subsistence_words, na.rm=T))
median(s_w$`rowSums(subsistence_words, na.rm = T)`)
range(s_w$`rowSums(subsistence_words, na.rm = T)`)
colnames(s_w)<- "count"
bab <- as.data.frame(0)
bed <- as.data.frame(0)
remaining <- rbind(bab, bed)
colnames(remaining) <- "count"
rownames(remaining) <- c("Babongo", "Bedzan")
rbind(s_w, remaining)
median(s_w$count)
m_w <- rowSums(music_words, na.rm=T)
median(m_w)
write.csv(s_w, "results_masked/subsistence_words_counts.csv")
write.csv(m_w, "results_masked/music_words_counts.csv")
s_w <- read.csv("results_masked/subsistence_words_counts.csv", row.names = 1)
m_w <- read.csv("results_masked/music_words_counts.csv", row.names = 1)

IDS_m <- rownames(m_w)
IDS_s <- rownames(s_w)

png("results_masked/Musical_instruments_words.png",units="in", width=7.5, height=5, res=500)
barplot(m_w$x, col=my_pal_long, names.arg=IDS_m, ylab="Words available musical instruments", ylim=c(0,45), cex.names=0.5, cex.axis=0.5 )
dev.off()
png("results_masked/Subsistence_tools_words.png",units="in", width=7.5, height=5, res=500)
barplot(s_w$count, col=my_pal_long[-c(2,8)], names.arg=IDS_s, ylab="Words available subsistence tools", ylim=c(0,45), cex.names=0.5, cex.axis=0.5 )
dev.off()

library(vegan)
jacSubsistenceWords<-vegdist(subsistence_words,"jaccard",na.rm=TRUE)
jacMusicWords<-vegdist(music_words,"jaccard",na.rm=TRUE)

t.test(jacSubsistenceWords, jacMusicWords)

# Prepare proportion test
# Compute proportions:
music_total <- ncol(music_words)
# create dyads (if you put 3, you can create triads)
# cant do that because for that we would need 1 word per group
# m.dyads <- data.frame(combn(colnames(music_words), 2)) # Creates all possible dyads
music_dyads <- data.frame(ID=rownames(music_words),
                          stringsAsFactors=FALSE)
#with transpose - turn rows into columns
music_dyads <- data.frame(t(combn(music_dyads$ID, 2)))
colnames(music_dyads) <- c("ID1", "ID2")
n_dyads <- nrow(music_dyads)

#replicate lines to match number of shared words per dyad (=word number)
music_dyads <- music_dyads[rep(1:nrow(music_dyads),each=music_total),]

#create plant knowledge column in m.dyads file
#repeat=  number of words; times = number of dyads
music_dyads$word <- rep(1:music_total, times = n_dyads)

#create a file with music words and shared music words
#make a file with ID and words
music_words$ID <- rownames(music_words)
music_words2 <- music_words[rep(1:nrow(music_words), each=music_total),]

#but only want the ID column
music_words2<- subset(music_words2, select=ID)

#change it to ID1
colnames(music_words2)[1] <- "ID1" 

#add plant column
music_words2$word <- rep(1:music_total, times = nrow(music_words))

#create a music word column
#transpose whole music_words file into single column of music words
#but first exclude the columns that are not about music words
music_words3 <- subset(music_words, select = -c(ID))
music_words2$use <- c(t(music_words3)) # make a single vector that contains the transposed columns

# you need to match by ID1 and music words as they are the common columns in both datasets
music.dyads <- merge(music_dyads, music_words2, by = c("ID1", "word")) 

#reorder by ID 
music.dyads <- music.dyads[order(music.dyads$ID1,music.dyads$ID2),]

#Id1 knowledge: change name from know to know1
colnames(music.dyads)[4] <- "use1"

#then ID2 knowledge
#change column names to allow proper matching
colnames(music_words2)[1] <- "ID2" #change the names in the long file with plant knowledge
colnames(music_words2)[3] <- "use2"
music.dyads <- merge(music.dyads, music_words2, by = c("ID2", "word"))

#reorder rows and columns
music.dyads <- music.dyads[order(music.dyads$ID1,music.dyads$ID2),]

#create shared word variable (0,0 is not shared)
music.dyads$shared <- ifelse(music.dyads$use1 == 1 & music.dyads$use2 == 1, 1, 0)

table(music.dyads$shared)
# 6985 possible shared words between pairs of populations; out of which 30 occur

# Now subsistence:
subsistence_total <- ncol(subsistence_words)
# create dyads (if you put 3, you can create triads)
# cant do that because for that we would need 1 word per group
# m.dyads <- data.frame(combn(colnames(subsistence_words), 2)) # Creates all possible dyads
subsistence_dyads <- data.frame(ID=rownames(subsistence_words),
                          stringsAsFactors=FALSE)
#with transpose - turn rows into columns
subsistence_dyads <- data.frame(t(combn(subsistence_dyads$ID, 2)))
colnames(subsistence_dyads) <- c("ID1", "ID2")
n_dyads <- nrow(subsistence_dyads)

#replicate lines to match number of shared words per dyad (=word number)
subsistence_dyads <- subsistence_dyads[rep(1:nrow(subsistence_dyads),each=subsistence_total),]

#create plant knowledge column in m.dyads file
#repeat=  number of words; times = number of dyads
subsistence_dyads$word <- rep(1:subsistence_total, times = n_dyads)

#create a file with subsistence words and shared subsistence words
#make a file with ID and words
subsistence_words$ID <- rownames(subsistence_words)
subsistence_words2 <- subsistence_words[rep(1:nrow(subsistence_words), each=subsistence_total),]

#but only want the ID column
subsistence_words2<- subset(subsistence_words2, select=ID)

#change it to ID1
colnames(subsistence_words2)[1] <- "ID1" 

#add plant column
subsistence_words2$word <- rep(1:subsistence_total, times = nrow(subsistence_words))

#create a subsistence word column
#transpose whole subsistence_words file into single column of subsistence words
#but first exclude the columns that are not about subsistence words
subsistence_words3 <- subset(subsistence_words, select = -c(ID))
subsistence_words2$use <- c(t(subsistence_words3)) # make a single vector that contains the transposed columns

# you need to match by ID1 and subsistence words as they are the common columns in both datasets
subsistence.dyads <- merge(subsistence_dyads, subsistence_words2, by = c("ID1", "word")) 

#reorder by ID 
subsistence.dyads <- subsistence.dyads[order(subsistence.dyads$ID1,subsistence.dyads$ID2),]

#Id1 knowledge: change name from know to know1
colnames(subsistence.dyads)[4] <- "use1"

#then ID2 knowledge
#change column names to allow proper matching
colnames(subsistence_words2)[1] <- "ID2" #change the names in the long file with plant knowledge
colnames(subsistence_words2)[3] <- "use2"
subsistence.dyads <- merge(subsistence.dyads, subsistence_words2, by = c("ID2", "word"))

#reorder rows and columns
subsistence.dyads <- subsistence.dyads[order(subsistence.dyads$ID1,subsistence.dyads$ID2),]

#create shared word variable (0,0 is not shared)
subsistence.dyads$shared <- ifelse(subsistence.dyads$use1 == 1 & subsistence.dyads$use2 == 1, 1, 0)

table(subsistence.dyads$shared)
# 3276 possible subsistence dyads
# create column specifying type of word
subsistence.dyads$type <- "subsistence"
music.dyads$type <- "music"

# words total df
words_df <- rbind(subsistence.dyads, music.dyads)

#proportion test, does cultural domain affect shared?
pt <- table(words_df$type, words_df$shared)
prop.test(pt, correct=F)



###########
## DRUMS ##
###########

# Here you need to reintroduce the Mbendjele in all Matrices
# load locations data
locations <- read.csv("location_culture_centroids.csv", row.names = 1)
locations$FID <- rownames(locations)
# load biome data and /100 to get proportion
biome <- (read.csv("culture_biomes.csv", row.names = 1))/100

# Maybe is not the % biome that dictates whether you have an object but simply whether
# that biome is present
biome_bin <- read.csv("culture_biomes.csv", row.names = 1)
biome_bin <- ifelse(biome_bin>0,1,0)

# Read in genetic distances and turn into distance object
Fst_mat<- read.csv("fst_masked.csv", row.names = 1)
genDist <- as.dist(Fst_mat)
colnames(Fst_mat) <- rownames(Fst_mat)

Fst_mat_IBD<- read.csv("fst_masked_ibd.csv", row.names = 1)
genDistIBD <- as.dist(Fst_mat_IBD)
colnames(Fst_mat_IBD) <- rownames(Fst_mat_IBD)

Fst_mat_Bantu<- read.csv("fst_Bantu.csv", row.names = 1)
genDistBantu<- as.dist(Fst_mat_Bantu)
colnames(Fst_mat_Bantu) <- rownames(Fst_mat_Bantu)


####################
# Data Preparation #
####################

#compute spatial distances
distTable <- data.frame(x=locations$Longitude,y=locations$Latitude)
names_locs <- c("Aka", "Babongo", "Baka", "Bakola", "Bakoya", "Batwa (East)", "Batwa (West)", "Bedzan", "Efe","Mbendjele", "Sua" )
rownames(distTable) <- names_locs
# with ecodist, distMat doesnt work
detach("package:ecodist", unload = TRUE)
space<-distMat(data.frame(x=locations$Longitude,y=locations$Latitude))

#compute gower dissimilarity between biome composition
library(vegan)
gowBiome<- vegdist(biome,"gower",na.rm=TRUE)
gowBiomeBin<- vegdist(biome_bin,"gower",na.rm=TRUE)

locations_red <- locations[-7,]
biome_red <- biome[-7,]
biome_bin_red <- biome_bin[-7,]


distTable2 <- data.frame(x=locations_red$Longitude,y=locations_red$Latitude)
rownames(distTable2) <- names_locs[-7]
spaceRed<-distMat(data.frame(x=locations_red$Longitude,y=locations_red$Latitude))

gowBiomeRed<- vegdist(biome_red,"gower",na.rm=TRUE)
gowBiomeRedBin<- vegdist(biome_bin_red,"gower",na.rm=TRUE)


drums <- read.csv("drums.csv", stringsAsFactors = T,header=T)
ids <- seq(1:nrow(drums))
for (i in ids) {
  ids[i] <- paste("T-", i, sep="")
}

drums$ID <- ids

# For genetics, create another drum dataset with only drums from groups with genes
drums_red <- subset(drums, drums$Culture %in% rownames(Fst_mat))

drums$Trait <- NULL
ids <- drums$ID
ids_red <- drums_red$ID
rownames(drums) <- ids
rownames(drums_red) <- ids_red
drums$ID <- NULL
drums$Culture <- NULL
drums_red$Trait <- NULL
drums_red$ID <- NULL
drums_red$Culture <- NULL

# this file needs to have the coordinates of the culture corresponding to each individual drum
# in the same order as your drums.csv file
locationsdrum <- read.csv("drums_locations.csv")
locationsdrum_red <- subset(locationsdrum, locationsdrum$Culture %in% rownames(Fst_mat))

rownames(locationsdrum) <- rownames(drums)
rownames(locationsdrum_red) <- rownames(drums_red)

biomesdrum <- read.csv("drums_biomes.csv")
biomesdrum_red <- subset(biomesdrum, biomesdrum$Culture %in% rownames(Fst_mat))
biomesdrum$Trait <- NULL
biomesdrum$Culture <- NULL
biomesdrum_red$Trait <- NULL
biomesdrum_red$Culture <- NULL

rownames(biomesdrum) <- rownames(drums)
rownames(biomesdrum_red) <- rownames(drums_red)

# Convert factors into numeric codes (not really needed if stringsAsFactors = F)
indx <- sapply(drums, is.factor)
indx_red <- sapply(drums_red, is.factor)

drums[indx] <- lapply(drums[indx],  function(x) as.numeric(as.factor(x)))
drums_red[indx_red] <- lapply(drums_red[indx_red],  function(x) as.numeric(as.factor(x)))

#compute Gower dissimilarity distance
library(vegan)
gowdrum <- vegdist(drums,"gower",na.rm=TRUE) 
gowdrumRed <- vegdist(drums_red,"gower",na.rm=TRUE) 

gowBiomedrums <- as.matrix(vegdist(biomesdrum,"gower",na.rm=TRUE))
gowBiomedrumsRed <-  as.matrix(vegdist(biomesdrum_red,"gower",na.rm=TRUE))


#compute binary Matrices of cultural affiliation (0:same culture; 1:different culture)
binaryCultDistdrums<-as.matrix(as.dist(BinaryCulture(locationsdrum$Culture)))
binaryCultDistdrumsRed<-as.matrix(as.dist(BinaryCulture(locationsdrum_red$Culture)))

#compute Fst matrix
GenDistdrums<-as.matrix(as.dist(GetFst(locationsdrum_red$Culture)))

#compute spatial distance between cultures and therefore drums
#spacedrums<-distMat(data.frame(x=locationsdrum$Longitude,y=locationsdrum$Latitude))

#compute spatial distances
distTabledrums <- data.frame(x=locationsdrum$Longitude,y=locationsdrum$Latitude)
namesdrum <- as.vector(row.names(locationsdrum))
rownames(distTabledrums) <- namesdrum

distTabledrumsRed <- data.frame(x=locationsdrum_red$Longitude,y=locationsdrum_red$Latitude)
namesdrumRed <- as.vector(row.names(locationsdrum_red))
rownames(distTabledrumsRed) <- namesdrumRed

library(geosphere)
spacedrums <- distm(distTabledrums, fun=distGeo)
rownames(spacedrums) <- namesdrum
colnames(spacedrums) <- namesdrum
# here distance is in metres, not kilometres
#spacedrums <- as.dist(spacedrums)

spacedrumsRed <- distm(distTabledrumsRed, fun=distGeo)
rownames(spacedrumsRed) <- namesdrumRed
colnames(spacedrumsRed) <- namesdrumRed
# here distance is in metres, not kilometres
#spacedrumsRed <- as.dist(spacedrumsRed)


gowdrumRed <- as.matrix(gowdrumRed)
na_rows <- which(is.na(gowdrumRed), arr.ind=TRUE)
na_rows <- unique(na_rows[,1])
gowdrumRed <- gowdrumRed[-na_rows, -na_rows]
spacedrumsRed  <- spacedrumsRed [-na_rows, -na_rows]
gowBiomedrumsRed  <- gowBiomedrumsRed [-na_rows, -na_rows]
binaryCultDistdrumsRed  <- binaryCultDistdrumsRed [-na_rows, -na_rows]
GenDistdrums <- GenDistdrums[-na_rows, -na_rows]
GenDistdrumsIBD <- GenDistdrumsIBD[-na_rows, -na_rows]

gowdrum <- as.matrix(gowdrum)
na_rows <- which(is.na(gowdrum), arr.ind=TRUE)
na_rows <- unique(na_rows[,1])
gowdrum <- gowdrum[-na_rows, -na_rows]
spacedrums <- spacedrums[-na_rows, -na_rows]
gowBiomedrums <- gowBiomedrums[-na_rows, -na_rows]
binaryCultDistdrums <- binaryCultDistdrums[-na_rows, -na_rows]


#partial mantel test between Distance and Geography and between Culture and Geography
library(ecodist)
detach("package:vegan", unload = TRUE)

space_mantel <- mantel(lower(gowdrum)~lower(spacedrums),nperm=1000)
culture_mantel <-mantel(lower(gowdrum)~lower(binaryCultDistdrums),nperm=1000)
genes_mantel <- mantel(lower(gowdrumRed)~lower(GenDistdrums),nperm=1000)
biome_mantel <- mantel(lower(gowBiomedrumsRed)~lower(GenDistdrums),nperm=1000)

#partial mantel test testing for effect of CULTURAL AFFILIATION controlling for rest
cult_space_mantel <- mantel(lower(gowdrum)~lower(binaryCultDistdrums)+lower(spacedrums),nperm=1000)
cult_genes_mantel <- mantel(lower(gowdrumRed)~lower(binaryCultDistdrumsRed)+lower(GenDistdrums),nperm=1000)
cult_biome_mantel <- mantel(lower(gowdrum)~lower(binaryCultDistdrums)+lower(gowBiomedrums),nperm=1000)

#partial mantel test testing for effect of SHARED GENES
genes_space_mantel <- mantel(lower(gowdrumRed)~lower(GenDistdrums)+lower(spacedrumsRed),nperm=1000)
genes_cult_mantel <- mantel(lower(gowdrumRed)~lower(GenDistdrums)+lower(binaryCultDistdrumsRed),nperm=1000)
genes_biome_mantel <- mantel(lower(gowdrumRed)~lower(GenDistdrums)+lower(gowBiomedrumsRed),nperm=1000)

#partial mantel test testing for effect of GEOGRAPHY
space_culture_mantel <- mantel(lower(gowdrum)~lower(spacedrums)+lower(binaryCultDistdrums),nperm=1000)
space_genes_mantel <- mantel(lower(gowdrumRed)~lower(spacedrumsRed)+lower(GenDistdrums),nperm=1000)
space_biome_mantel <- mantel(lower(gowdrum)~lower(spacedrums)+lower(gowBiomedrums),nperm=1000)

#partial mantel test testing for effect of BIOME
biome_culture_mantel <- mantel(lower(gowdrum)~lower(gowBiomedrums)+lower(binaryCultDistdrums),nperm=1000)
biome_genes_mantel <- mantel(lower(gowdrumRed)~lower(gowBiomedrumsRed)+lower(GenDistdrums),nperm=1000)
biome_space_mantel <- mantel(lower(gowdrum)~lower(gowBiomedrums)+lower(spacedrums),nperm=1000)


write.csv(rbind(space_mantel, culture_mantel, genes_mantel,biome_mantel,
                cult_space_mantel,cult_genes_mantel,cult_biome_mantel,
                genes_space_mantel, genes_cult_mantel,genes_biome_mantel,
                space_culture_mantel,space_genes_mantel,space_biome_mantel,
                biome_culture_mantel, biome_genes_mantel, biome_space_mantel), "results_masked/drums_mantel.csv")


#########
# AMOVA #
#########
library(pegas)

#AMOVA
sqgowdrum<-sqrt(gowdrum)
sqgowdrum_mat <- as.data.frame(as.matrix(sqgowdrum))
# get NA
#na_rows <- which(is.na(sqgowdrum_mat), arr.ind=TRUE)

# remove rows with NAs
#na_rows <- unique(na_rows[,1])
#sqgowdrum_mat <- sqgowdrum_mat[-na_rows,-na_rows]
sqgowdrum <- as.dist(sqgowdrum_mat)

groupsdrum<-locationsdrum$Culture
groupsdrum <- as.factor(groupsdrum[-na_rows])
AMOVAdrum=pegas::amova(sqgowdrum~groupsdrum)

# Phi-Stat:
AMOVAdrum$varcomp$sigma2[1]/sum(AMOVAdrum$varcomp$sigma2)
#P-Value:
AMOVAdrum$varcomp$P.value[1]


##################
# pairwise AMOVA #
##################

# Note: This is the important thing to compare against genetic data

#Extract Culture Data
CultureNames<-unique(as.character(groupsdrum))
CultureNames <- sort(CultureNames)
NumberOfCultures<-length(CultureNames)
sqgowdrum<-as.matrix(sqgowdrum)
CultureList<-groupsdrum

#PhiSt and Pval Output Matrices:
PhiMatrix<-matrix(NA,nrow=NumberOfCultures,ncol=NumberOfCultures)
PvalMatrix<-matrix(NA,nrow=NumberOfCultures,ncol=NumberOfCultures)

#Spatial Coordinates Output Data.Frame:
XY<-data.frame(x=numeric(length=NumberOfCultures),y=numeric(length=NumberOfCultures))

for (i in 1:NumberOfCultures) {
  print(paste(i," out of ",NumberOfCultures))
  for (j in 1:NumberOfCultures)
  {
    if (i!=j&i>j)
    { 
      indexI=which(as.character(CultureList)==as.character(CultureNames[i]))
      indexJ=which(as.character(CultureList)==as.character(CultureNames[j]))
      #subsetDistance Matrix
      subDistMat<-sqgowdrum[c(indexI,indexJ),c(indexI,indexJ)]
      subDistMat<-as.dist(subDistMat)
      classes<-as.factor(c(rep("a",length(indexI)),rep("b",length(indexJ))))
      #compute Phi:
      tmp<-pegas::amova(subDistMat~classes)
      N0=tmp$varcoef
      AP=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)
      WP=(tmp$tab$MSD[2])
      PhiMatrix[i,j]=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)/((tmp$tab$MSD[2])+((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef))
      if(PhiMatrix[i,j]<0){PhiMatrix[i,j]=0}
      PvalMatrix[i,j]=tmp$varcomp$P.value[1]
    }
  }
  XY$x[i]<-mean(locationsdrum$Longitude[which(locationsdrum$Culture==as.character(CultureNames[i]))],na.rm=TRUE)
  XY$y[i]<-mean(locationsdrum$Latitude[which(locationsdrum$Culture==as.character(CultureNames[i]))],na.rm=TRUE)
}

dim <- ncol(PhiMatrix)

#Change NAs to blanks for text
PhiMatrix_text <- as.data.frame(PhiMatrix)
PhiMatrix_text <- sapply(PhiMatrix_text, as.character)
PhiMatrix_text[is.na(PhiMatrix_text)] <- " "
PhiMatrix_text <- as.matrix(PhiMatrix_text)

png(file="results_masked/drums_phistats.png",
    width=7, height=6, units="in", res=500)
image(-7:(dim-8), 1:dim, PhiMatrix, axes = FALSE, xlab="", ylab="", col=hcl.colors(12, "Geyser", rev=TRUE), main="drums", xlim=c(-7,5.5), ylim=c(0,11))
axis(1, (-7:(dim-8)), CultureNames, cex.axis = 0.6, las=3, pos=0.5)
axis(4, 1:dim, CultureNames, cex.axis = 0.6, las=1, pos=3.5)
text(expand.grid((-7:(dim-8)), 1:dim), sprintf("%0.4s", PhiMatrix_text), cex=0.6)
dev.off()

# Calculate geographic distance between cultures
#cultSpDist<-distMat(XY)

# Again Partial Mantel Tests using Phi Statistic
library(ecodist)
CultureNames_rn <- CultureNames
CultureNames_rn[CultureNames_rn=="Batwa (East)"] <- "BatwaE"
CultureNames_rn[CultureNames_rn=="Batwa (West)"] <- "BatwaW"
rownames(PhiMatrix) <- CultureNames
colnames(PhiMatrix) <- CultureNames
PhiMatrix_rn <- PhiMatrix
rownames(PhiMatrix_rn) <- CultureNames_rn
colnames(PhiMatrix_rn) <- CultureNames_rn

writeDist(PhiMatrix_rn, file="drums_phi.nex", format="nexus")

PhiMatrixdrums <- PhiMatrix
# Remove Batwa west as we dont have genetics
PhiMatrixdrums <- PhiMatrix[-7, -7]
#PhiMatrixdrums <- as.dist(PhiMatrixdrums)

# Mantel tests between geography, biome, genetics and drums

biome_mantel <- mantel(lower(PhiMatrixdrums)~lower(gowBiomeRed), mrank=T)
space_mantel <- mantel(lower(PhiMatrixdrums)~lower(spaceRed), mrank=T)
genes_mantel <- mantel(lower(PhiMatrixdrums)~lower(genDist), mrank=T)
genes_IBD_mantel <- mantel(lower(PhiMatrixdrums)~lower(genDistIBD), mrank=T)
genes_Bantu_mantel <- mantel(lower(PhiMatrixdrums[-c(8,10), -c(8,10)])~lower(genDistBantu), mrank=T)

#partial mantel test testing for effect of SHARED GENES
genes_space_mantel <- mantel(lower(PhiMatrixdrums)~lower(genDist)+lower(spaceRed), mrank=T)
genes_biome_mantel <- mantel(lower(PhiMatrixdrums)~lower(genDist)+lower(gowBiomeRed), mrank=T)

#partial mantel test testing for effect of GEOGRAPHY
space_genes_mantel <- mantel(lower(PhiMatrixdrums)~lower(spaceRed)+lower(genDist), mrank=T)
space_biome_mantel <- mantel(lower(PhiMatrixdrums)~lower(spaceRed)+lower(gowBiomeRed), mrank=T)

#partial mantel test testing for effect of BIOME
biome_genes_mantel <- mantel(lower(PhiMatrixdrums)~lower(gowBiomeRed)+lower(genDist), mrank=T)
biome_space_mantel <- mantel(lower(PhiMatrixdrums)~lower(gowBiomeRed)+lower(spaceRed), mrank=T)


mantel_drums <- rbind(space_mantel, genes_mantel,
                      biome_mantel,
                      genes_space_mantel, genes_biome_mantel,
                      space_genes_mantel, space_biome_mantel,
                      biome_space_mantel, biome_genes_mantel)

#write.csv(mantel_drums, "results_masked/mantel_drums.csv")


mantel_drums <-  rbind(space_mantel, genes_mantel,
                                      biome_mantel, genes_Bantu_mantel, genes_IBD_mantel)
pval_adj_drums <- p.adjust(mantel_drums[,2], method = "BH")
mantel_drums<- cbind(mantel_drums, pval_adj_drums)

write.csv(mantel_drums, "final_paper_files/results_final_paper/results_masked/mantel_drums.csv")

matrix_drums <- rbind(genes_space_mantel, genes_biome_mantel,
                      space_genes_mantel, space_biome_mantel,
                      biome_space_mantel, biome_genes_mantel)
pval_adj_drums <- p.adjust(matrix_drums[,2], method = "BH")
matrix_drums<- cbind(matrix_drums, pval_adj_drums)

write.csv(matrix_drums, "results_masked/matrix_drums.csv")

# Plot Fst in the same way as drums so it looks pretty
Fst_plot <- as.matrix(Fst_mat_IBD)
Fst_plot[upper.tri(Fst_plot)] <- 0
Fst_plot[Fst_plot==0] <- NA
FstMatrix_text <- as.data.frame(Fst_plot)
FstMatrix_text <- sapply(FstMatrix_text, as.character)
FstMatrix_text[is.na(FstMatrix_text)] <- " "
FstMatrix_text <- as.matrix(FstMatrix_text)

dim <- ncol(Fst_plot)
CultureNames<-unique(as.character(rownames(Fst_plot)))

setwd("/Users/Cecilia/Documents/PhD/African_Pygmies/Spatial analysis culture/My_data/Nov_2021/final_paper_files/results_final_paper")
png(file="results_masked/Fst_plot_IBD.png",
    width=7, height=6, units="in", res=500)
image(-7:(dim-8), 1:dim, Fst_plot, axes = FALSE, xlab="", ylab="", col=hcl.colors(12, "Geyser", rev=TRUE), main="Fst", xlim=c(-7,5.5), ylim=c(0,11))
axis(1, (-7:(dim-8)), CultureNames, cex.axis = 0.6, las=3, pos=0.5)
axis(4, 1:dim, CultureNames, cex.axis = 0.6, las=1, pos=2.5)
text(expand.grid((-7:(dim-8)), 1:dim), sprintf("%0.5s", FstMatrix_text), cex=0.6)
dev.off()

Fst_plot_red <- as.matrix(Fst_mat_IBD[1:8,1:8])
Fst_plot_red[upper.tri(Fst_plot_red)] <- 0
Fst_plot_red[Fst_plot_red==0] <- NA
FstMatrix_text_red <- as.data.frame(Fst_plot_red)
FstMatrix_text_red <- sapply(FstMatrix_text_red, as.character)
FstMatrix_text_red[is.na(FstMatrix_text_red)] <- " "
FstMatrix_text_red <- as.matrix(FstMatrix_text_red)
rownames(Fst_plot_red)[8] <- "Efe & Sua"
colnames(Fst_plot_red)[8] <- rownames(Fst_plot_red)[8] 

dim <- ncol(Fst_plot_red)
CultureNames<-unique(as.character(rownames(Fst_plot_red)))
png(file="results_masked/Fst_plot_red_IBD.png",
    width=7, height=6, units="in", res=500)
image(-7:(dim-8), 1:dim, Fst_plot_red, axes = FALSE, xlab="", ylab="", col=hcl.colors(12, "Geyser", rev=TRUE), main="Masked Fst, no IBD segments", xlim=c(-7,5.5), ylim=c(0,11))
axis(1, (-7:(dim-8)), CultureNames, cex.axis = 0.6, las=3, pos=0.5)
axis(4, 1:dim, CultureNames, cex.axis = 0.6, las=1, pos=0.5)
text(expand.grid((-7:(dim-8)), 1:dim), sprintf("%0.5s", FstMatrix_text_red), cex=0.6)
dev.off()

# bantu fst plot
Fst_plot_red <- as.matrix(Fst_mat_Bantu[1:7,1:7])
Fst_plot_red[upper.tri(Fst_plot_red)] <- 0
Fst_plot_red[Fst_plot_red==0] <- NA
FstMatrix_text_red <- as.data.frame(Fst_plot_red)
FstMatrix_text_red <- sapply(FstMatrix_text_red, as.character)
FstMatrix_text_red[is.na(FstMatrix_text_red)] <- " "
FstMatrix_text_red <- as.matrix(FstMatrix_text_red)
rownames(Fst_plot_red)[1] <- "Aka"
colnames(Fst_plot_red)[1] <- rownames(Fst_plot_red)[1] 

dim <- ncol(Fst_plot_red)
CultureNames<-unique(as.character(rownames(Fst_plot_red)))
png(file="final_paper_files/results_final_paper/results_masked/Fst_plot_red_Bantu.png",
    width=7, height=6, units="in", res=500)
image(-7:(dim-8), 1:dim, Fst_plot_red, axes = FALSE, xlab="", ylab="", col=hcl.colors(12, "Geyser", rev=TRUE), main="Masked Fst, only Bantu segments", xlim=c(-7,5.5), ylim=c(0,11))
axis(1, (-7:(dim-8)), CultureNames, cex.axis = 0.6, las=3, pos=0.5)
axis(4, 1:dim, CultureNames, cex.axis = 0.6, las=1, pos=-0.5)
text(expand.grid((-7:(dim-8)), 1:dim), sprintf("%0.5s", FstMatrix_text_red), cex=0.6)
dev.off()

# Perform PCoA
gowdrum <- as.data.frame(as.matrix(gowdrum))
# get NA
na_rows <- which(is.na(gowdrum), arr.ind=TRUE)

# remove rows with NAs
na_rows <- unique(na_rows[,1])
gowdrum <- gowdrum[-na_rows,-na_rows]
gowdrum <- as.dist(gowdrum)
drums_pcoa <- pcoa(gowdrum, correction = "cailliez")

#We rescale the PCoA components in relation to the explained variance.  

for(i in 1:ncol(drums_pcoa$vectors)) { 
  drums_pcoa$vectors[,i] <- scale(drums_pcoa$vectors[,i])*
    drums_pcoa$values$Rel_corr_eig[i]}

# Extract eigenvalues from PC/PCo results_masked 
drums_ev <- drums_pcoa$values$Corr_eig

# Plot
library(ggplot2)
my_font <- "sans"

png(file="results_masked/pcoadrums.png",
    width=6, height=4.5, units="in", res=500)
plot_pc_scree(eigenval=drums_ev, title="drums", type="PCo")
dev.off()

# Extract the PCs and the PCos from the PCA and PCoA results_masked
drums_pco <- as.data.frame(drums_pcoa$vectors)
drums_pco$Culture <- CultureList

## PCoA with groups

#install.packages('ggfortify')

p <- ggplot(data = drums_pco, aes(x = Axis.1, y = Axis.2, color=Culture)) +
  geom_point(size = 3.5) +
  xlab("PCo1 (13.50%)") +
  ylab("PCo2 (6.56%)") + 
  theme_classic() + 
  scale_color_manual(values = my_pal_long,
                     labels = c("Aka (10)",
                                "Babongo (12)",
                                "Baka (8)",
                                "Bakola (2)",
                                "Bakoya (6)",
                                "East Batwa (23)",
                                "West Batwa (8)",
                                "Bedzan (8)",
                                "Efe (3)",
                                "Mbendjele (2)",
                                "Sua (3)")) +
  geom_hline(yintercept = 0, lty = 2) +
  ggtitle("drums") +
  geom_vline(xintercept = 0, lty = 2) 

p <- p + theme(legend.title = element_blank(),
               legend.text = element_text(size = 12))

ggsave(p, filename="results_masked/PCOA_drums.png", height=150, width=180, units="mm")

#We compute Euclidean distance matrices from the dimensionality-reduced data and create NeighborNets. 
# SplitsTree doesnt accept brackets in names
drums_pco_rn <- drums_pco
drums_pco_rn$Culture <- as.character(drums_pco_rn$Culture)
drums_pco_rn$Culture[drums_pco_rn$Culture=="Batwa (East)"] <- "BatwaE"
drums_pco_rn[drums_pco_rn=="Batwa (West)"] <- "BatwaW"
drums_pco_rn$Culture <- NULL

# Write nexus files for splitstree (Not used so far)
require(phangorn)

dist_list <- lapply(list(drums=drums_pco_rn),
                    function(f) {dist(f)})
writeDist(dist_list$drums, file="drums.nex", format="nexus")

###################################################################
## Now we test for potential mechanisms of cultural transmission ##
###################################################################

# Divide cultures into Western and Eastern Pygmies 

# Question - should we separate between different regions???
drums <- read.csv("drums.csv")

# Remove traits that are only present in one of the two regions (East or West)
drums$Region<-sapply(as.character(drums$Culture),function(x)
  ifelse(x %in% c("Aka", "Babongo", "Baka", "Bakola", "Bakoya",
                  "Bedzan", "Mbendjele", "Batwa (West)"),"West",
         ifelse(x %in% c("Batwa (East)","Efe","Sua"),"East")))
table(drums$Region)

drum_redux<-drums[,c(3:ncol(drums))]
for(i in 1:(ncol(drum_redux)-1)){
  tab<-table(drum_redux[,c("Region",colnames(drum_redux)[i])])  
  check<-apply(tab,2,function(x) ifelse(sum(x!=0)==1,"fix","ok"))  
  drum_redux[,colnames(drum_redux)[i]]<-sapply(drum_redux[,colnames(drum_redux)[i]],
                                               function(x) ifelse(x %in% names(check)[check=="fix"],NA,as.character(x)))
  drum_redux[,colnames(drum_redux)[i]]<-as.factor(drum_redux[,colnames(drum_redux)[i]])}
drum_redux<-drum_redux[,sapply(drum_redux, function(x)!all(is.na(x)))]
drum_redux$Region<-as.factor(drum_redux$Region)

library(party)
library(plyr)
weights_drums<-sapply(drum_redux$Region,
                      function(x) nrow(drum_redux)/nrow(drum_redux[drum_redux$Region==x,]))
rf_redux<-cforest(Region~.,
                  controls=
                    cforest_control(ntree = 500,
                                    mtry=5,
                                    mincriterion=qnorm(0.9),
                                    fraction = 0.632,
                                    testtype ="Teststatistic",
                                    teststat="max",
                                    replace=TRUE,
                                    trace=F,
                                    savesplitstats=F,
                                    minsplit=20,
                                    minbucket=8),
                  data=drum_redux,
                  weights = weights_drums)

# Produce confusion matrix
cmat_redux<-table(drum_redux$Region, predict(rf_redux))
print(cmat_redux)

stats_region<-ldply(c(1:nrow(cmat_redux)),
                    function(x) data.frame(id=rownames(cmat_redux)[x],
                                           N=sum(cmat_redux[,x]),
                                           d=cmat_redux[x,x],
                                           A=sum(cmat_redux[x,]))) %>% transform(p=d/N,
                                                                                 r=d/A)

write.csv(stats_region, "results_masked/stats_region_drums.csv")


drum_redux$predReg<-predict(rf_redux)
missclass<-cbind(drums[rownames(drum_redux[drum_redux$Region!=drum_redux$predReg,]),],
                 data.frame(predReg=drum_redux$predReg[drum_redux$Region!=drum_redux$predReg]))
missclass
write.csv(missclass, "results_masked/drums_missclass.csv")
##########
## FAMD ##
##########

library("FactoMineR")
library("factoextra")
library("missMDA")

drums$Region <- NULL
drums[sapply(drums, is.character)] <- lapply(drums[sapply(drums, is.character)], 
                                       as.factor)
imp.drums <- imputeFAMD(drums[,1:17], ncp = 10)
res.famd <- FAMD(imp.drums$completeObs,   ## Set the target variable "Churn" as a supplementary variable, so it is not included in the analysis for now
                 graph = T, 
                 sup.var=1, # Culture is the supplementary variable
                 ncp=10)

## Inspect principal components
get_eigenvalue(res.famd)
fviz_screeplot(res.famd)
# Plot of variables
fviz_famd_var(res.famd, repel = TRUE, ggtheme=theme_classic())
ggsave("results_masked/drums_FAMD_varcont_dim.png", width=4.5, height=4, units="in", dpi=500)
# Contribution to the first dimension

fviz_contrib(res.famd, "var", axes = 1, ggtheme=theme_classic()) 
ggsave("results_masked/drums_FAMD_varDIM1.png", width=7, height=4.5, units="in", dpi=500)
# Contribution to the second dimension
fviz_contrib(res.famd, "var", axes = 2, ggtheme=theme_classic())
ggsave("results_masked/drums_FAMD_varDIM2.png", width=7, height=4.5, units="in", dpi=500)

fviz_famd_var(res.famd, "quanti.var", col.var = "contrib", 
              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
              repel = TRUE)

famd_plot_drums <- fviz_mfa_ind(res.famd, 
                          habillage = 1, # color by groups 
                          palette = my_pal_long,
                          addEllipses = FALSE, 
                          repel = TRUE, # Avoid text overlapping,
                          label="none",
                          ggtheme = theme_classic(),
) 

famd_plot_drums <- famd_plot_drums + ggtitle("drums - FAMD")
famd_plot_drums
ggsave("results_masked/drums_FAMD_noell.tiff", width=7, height=6, units="in", dpi=500)



############
## ARROWS ##
############

arrows <- read.csv("arrows.csv", stringsAsFactors = T,header=T)
ids <- seq(1:nrow(arrows))
for (i in ids) {
  ids[i] <- paste("T-", i, sep="")
}

arrows$ID <- ids

# For genetics, create another arrow dataset with only arrows from groups with genes
arrows_red <- subset(arrows, arrows$Culture %in% rownames(Fst_mat))

arrows$Trait <- NULL
ids <- arrows$ID
ids_red <- arrows_red$ID
rownames(arrows) <- ids
rownames(arrows_red) <- ids_red
arrows$ID <- NULL
arrows$Culture <- NULL
arrows_red$Trait <- NULL
arrows_red$ID <- NULL
arrows_red$Culture <- NULL

# this file needs to have the coordinates of the culture corresponding to each individual arrow
# in the same order as your arrows.csv file
locationsarrow <- read.csv("arrows_locations.csv")
locationsarrow_red <- subset(locationsarrow, locationsarrow$Culture %in% rownames(Fst_mat))

rownames(locationsarrow) <- rownames(arrows)
rownames(locationsarrow_red) <- rownames(arrows_red)

biomesarrow <- read.csv("arrows_biomes.csv")
biomesarrow_red <- subset(biomesarrow, biomesarrow$Culture %in% rownames(Fst_mat))
biomesarrow$Trait <- NULL
biomesarrow$Culture <- NULL
biomesarrow_red$Trait <- NULL
biomesarrow_red$Culture <- NULL

rownames(biomesarrow) <- rownames(arrows)
rownames(biomesarrow_red) <- rownames(arrows_red)

# Convert factors into numeric codes (not really needed if stringsAsFactors = F)
indx <- sapply(arrows, is.factor)
indx_red <- sapply(arrows_red, is.factor)

arrows[indx] <- lapply(arrows[indx],  function(x) as.numeric(as.factor(x)))
arrows_red[indx_red] <- lapply(arrows_red[indx_red],  function(x) as.numeric(as.factor(x)))

#compute Gower dissimilarity distance
library(vegan)
gowarrow <- vegdist(arrows,"gower",na.rm=TRUE) 
gowarrowRed <- vegdist(arrows_red,"gower",na.rm=TRUE) 

gowBiomearrows <- as.matrix(vegdist(biomesarrow,"gower",na.rm=TRUE))
gowBiomearrowsRed <-  as.matrix(vegdist(biomesarrow_red,"gower",na.rm=TRUE))


#compute binary Matrices of cultural affiliation (0:same culture; 1:different culture)
binaryCultDistarrows<-as.matrix(as.dist(BinaryCulture(locationsarrow$Culture)))
binaryCultDistarrowsRed<-as.matrix(as.dist(BinaryCulture(locationsarrow_red$Culture)))

#compute Fst matrix
GenDistarrows<-as.matrix(as.dist(GetFst(locationsarrow_red$Culture)))

#compute spatial distance between cultures and therefore arrows
#spacearrows<-distMat(data.frame(x=locationsarrow$Longitude,y=locationsarrow$Latitude))

#compute spatial distances
distTablearrows <- data.frame(x=locationsarrow$Longitude,y=locationsarrow$Latitude)
namesarrow <- as.vector(row.names(locationsarrow))
rownames(distTablearrows) <- namesarrow

distTablearrowsRed <- data.frame(x=locationsarrow_red$Longitude,y=locationsarrow_red$Latitude)
namesarrowRed <- as.vector(row.names(locationsarrow_red))
rownames(distTablearrowsRed) <- namesarrowRed

library(geosphere)
spacearrows <- distm(distTablearrows, fun=distGeo)
rownames(spacearrows) <- namesarrow
colnames(spacearrows) <- namesarrow
# here distance is in metres, not kilometres
#spacearrows <- as.dist(spacearrows)

spacearrowsRed <- distm(distTablearrowsRed, fun=distGeo)
rownames(spacearrowsRed) <- namesarrowRed
colnames(spacearrowsRed) <- namesarrowRed
# here distance is in metres, not kilometres
#spacearrowsRed <- as.dist(spacearrowsRed)


gowarrowRed <- as.matrix(gowarrowRed)
na_rows <- which(is.na(gowarrowRed), arr.ind=TRUE)
na_rows <- unique(na_rows[,1])
gowarrowRed <- gowarrowRed[-na_rows, -na_rows]
spacearrowsRed  <- spacearrowsRed [-na_rows, -na_rows]
gowBiomearrowsRed  <- gowBiomearrowsRed [-na_rows, -na_rows]
binaryCultDistarrowsRed  <- binaryCultDistarrowsRed [-na_rows, -na_rows]
GenDistarrows <- GenDistarrows[-na_rows, -na_rows]


gowarrow <- as.matrix(gowarrow)
na_rows <- which(is.na(gowarrow), arr.ind=TRUE)
na_rows <- unique(na_rows[,1])
gowarrow <- gowarrow[-na_rows, -na_rows]
spacearrows <- spacearrows[-na_rows, -na_rows]
gowBiomearrows <- gowBiomearrows[-na_rows, -na_rows]
binaryCultDistarrows <- binaryCultDistarrows[-na_rows, -na_rows]


#partial mantel test between Distance and Geography and between Culture and Geography
library(ecodist)
detach("package:vegan", unload = TRUE)

space_mantel <- mantel(lower(gowarrow)~lower(spacearrows),nperm=1000)
culture_mantel <-mantel(lower(gowarrow)~lower(binaryCultDistarrows),nperm=1000)
genes_mantel <- mantel(lower(gowarrowRed)~lower(GenDistarrows),nperm=1000)
biome_mantel <- mantel(lower(gowBiomearrowsRed)~lower(GenDistarrows),nperm=1000)

#partial mantel test testing for effect of CULTURAL AFFILIATION controlling for rest
cult_space_mantel <- mantel(lower(gowarrow)~lower(binaryCultDistarrows)+lower(spacearrows),nperm=1000)
cult_genes_mantel <- mantel(lower(gowarrowRed)~lower(binaryCultDistarrowsRed)+lower(GenDistarrows),nperm=1000)
cult_biome_mantel <- mantel(lower(gowarrow)~lower(binaryCultDistarrows)+lower(gowBiomearrows),nperm=1000)

#partial mantel test testing for effect of SHARED GENES
genes_space_mantel <- mantel(lower(gowarrowRed)~lower(GenDistarrows)+lower(spacearrowsRed),nperm=1000)
genes_cult_mantel <- mantel(lower(gowarrowRed)~lower(GenDistarrows)+lower(binaryCultDistarrowsRed),nperm=1000)
genes_biome_mantel <- mantel(lower(gowarrowRed)~lower(GenDistarrows)+lower(gowBiomearrowsRed),nperm=1000)

#partial mantel test testing for effect of GEOGRAPHY
space_culture_mantel <- mantel(lower(gowarrow)~lower(spacearrows)+lower(binaryCultDistarrows),nperm=1000)
space_genes_mantel <- mantel(lower(gowarrowRed)~lower(spacearrowsRed)+lower(GenDistarrows),nperm=1000)
space_biome_mantel <- mantel(lower(gowarrow)~lower(spacearrows)+lower(gowBiomearrows),nperm=1000)

#partial mantel test testing for effect of BIOME
biome_culture_mantel <- mantel(lower(gowarrow)~lower(gowBiomearrows)+lower(binaryCultDistarrows),nperm=1000)
biome_genes_mantel <- mantel(lower(gowarrowRed)~lower(gowBiomearrowsRed)+lower(GenDistarrows),nperm=1000)
biome_space_mantel <- mantel(lower(gowarrow)~lower(gowBiomearrows)+lower(spacearrows),nperm=1000)


write.csv(rbind(space_mantel, culture_mantel, genes_mantel,biome_mantel,
                cult_space_mantel,cult_genes_mantel,cult_biome_mantel,
                genes_space_mantel, genes_cult_mantel,genes_biome_mantel,
                space_culture_mantel,space_genes_mantel,space_biome_mantel,
                biome_culture_mantel, biome_genes_mantel, biome_space_mantel), "results_masked/arrows_mantel.csv")


#########
# AMOVA #
#########
library(pegas)

#AMOVA
sqgowarrow<-sqrt(gowarrow)
sqgowarrow_mat <- as.data.frame(as.matrix(sqgowarrow))
# get NA
#na_rows <- which(is.na(sqgowarrow_mat), arr.ind=TRUE)

# remove rows with NAs
#na_rows <- unique(na_rows[,1])
#sqgowarrow_mat <- sqgowarrow_mat[-na_rows,-na_rows]
sqgowarrow <- as.dist(sqgowarrow_mat)

groupsarrow<-locationsarrow$Culture
groupsarrow <- as.factor(groupsarrow[-na_rows])
AMOVAarrow=amova(sqgowarrow~groupsarrow)

# Phi-Stat:
AMOVAarrow$varcomp$sigma2[1]/sum(AMOVAarrow$varcomp$sigma2)
#P-Value:
AMOVAarrow$varcomp$P.value[1]


##################
# pairwise AMOVA #
##################

# Note: This is the important thing to compare against genetic data

# I think you cannot do amova for cultures with only 1 trait!

#Extract Culture Data
CultureNames<-unique(as.character(groupsarrow))
CultureNames <- sort(CultureNames)
NumberOfCultures<-length(CultureNames)
sqgowarrow<-as.matrix(sqgowarrow)
CultureList<-groupsarrow

# Remove Bakoya as only 1 instance
CultureList <- CultureList[-c(7,20)]
CultureNames <- CultureNames[-c(2,5)]
NumberOfCultures<-length(CultureNames)

#PhiSt and Pval Output Matrices:
PhiMatrix<-matrix(NA,nrow=NumberOfCultures,ncol=NumberOfCultures)
PvalMatrix<-matrix(NA,nrow=NumberOfCultures,ncol=NumberOfCultures)

#Spatial Coordinates Output Data.Frame:
XY<-data.frame(x=numeric(length=NumberOfCultures),y=numeric(length=NumberOfCultures))

for (i in 1:NumberOfCultures) {
  print(paste(i," out of ",NumberOfCultures))
  for (j in 1:NumberOfCultures)
  {
    if (i!=j&i>j)
    { 
      indexI=which(as.character(CultureList)==as.character(CultureNames[i]))
      indexJ=which(as.character(CultureList)==as.character(CultureNames[j]))
      #subsetDistance Matrix
      subDistMat<-sqgowarrow[c(indexI,indexJ),c(indexI,indexJ)]
      subDistMat<-as.dist(subDistMat)
      classes<-as.factor(c(rep("a",length(indexI)),rep("b",length(indexJ))))
      #compute Phi:
      tmp<-pegas::amova(subDistMat~classes)
      N0=tmp$varcoef
      AP=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)
      WP=(tmp$tab$MSD[2])
      PhiMatrix[i,j]=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)/((tmp$tab$MSD[2])+((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef))
      if(PhiMatrix[i,j]<0){PhiMatrix[i,j]=0}
      PvalMatrix[i,j]=tmp$varcomp$P.value[1]
    }
  }
  XY$x[i]<-mean(locationsarrow$Longitude[which(locationsarrow$Culture==as.character(CultureNames[i]))],na.rm=TRUE)
  XY$y[i]<-mean(locationsarrow$Latitude[which(locationsarrow$Culture==as.character(CultureNames[i]))],na.rm=TRUE)
}

dim <- ncol(PhiMatrix)

#Change NAs to blanks for text
PhiMatrix_text <- as.data.frame(PhiMatrix)
PhiMatrix_text <- sapply(PhiMatrix_text, as.character)
PhiMatrix_text[is.na(PhiMatrix_text)] <- " "
PhiMatrix_text <- as.matrix(PhiMatrix_text)

png(file="results_masked/arrows_phistats.png",
    width=7, height=6, units="in", res=500)
image(-7:(dim-8), 1:dim, PhiMatrix, axes = FALSE, xlab="", ylab="", col=hcl.colors(12, "Geyser", rev=TRUE), main="Arrows", xlim=c(-7,5.5), ylim=c(0,11))
axis(1, (-7:(dim-8)), CultureNames, cex.axis = 0.6, las=3, pos=0.5)
axis(4, 1:dim, CultureNames, cex.axis = 0.6, las=1, pos=0.5)
text(expand.grid((-7:(dim-8)), 1:dim), sprintf("%0.4s", PhiMatrix_text), cex=0.6)
dev.off()

# Calculate geographic distance between cultures
#cultSpDist<-distMat(XY)

# Again Partial Mantel Tests using Phi Statistic

CultureNames_rn <- CultureNames
CultureNames_rn[CultureNames_rn=="Batwa (East)"] <- "BatwaE"
CultureNames_rn[CultureNames_rn=="Batwa (West)"] <- "BatwaW"
rownames(PhiMatrix) <- CultureNames
colnames(PhiMatrix) <- CultureNames
PhiMatrix_rn <- PhiMatrix
rownames(PhiMatrix_rn) <- CultureNames_rn
colnames(PhiMatrix_rn) <- CultureNames_rn

writeDist(PhiMatrix_rn, file="arrows_phi.nex", format="nexus")

# Mantel tests between geography, biome, genetics and arrows
genDist <- as.matrix(genDist)
genDistIBD <- as.matrix(genDistIBD)
genDistBantu <- as.matrix(genDistBantu)
gowBiomeRed <- as.matrix(gowBiomeRed)
spaceRed <- as.matrix(spaceRed)

genDistArrow <- genDist[-c(2,5),-c(2,5)]
genDistArrowIBD <- genDistIBD[-c(2,5),-c(2,5)]
genDistArrowBantu <- genDistBantu[-c(2,5),-c(2,5)]
BiomeDistArrow <- gowBiomeRed[-c(2,5),-c(2,5)]
spaceArrow <- spaceRed[-c(2,5),-c(2,5)]

PhiMatrixArrows <- as.matrix(PhiMatrix)

biome_mantel <- mantel(lower(PhiMatrixArrows)~lower(BiomeDistArrow), mrank=T)
space_mantel <- mantel(lower(PhiMatrixArrows)~lower(spaceArrow), mrank=T)
genes_mantel <- mantel(lower(PhiMatrixArrows)~lower(genDistArrow), mrank=T)
genes_ibd_mantel <- mantel(lower(PhiMatrixArrows)~lower(genDistArrowIBD), mrank=T)
genes_Bantu_mantel <- mantel(lower(PhiMatrixArrows[-c(7:8), -c(7:8)])~lower(genDistArrowBantu), mrank=T)

#partial mantel test testing for effect of SHARED GENES
genes_space_mantel <- mantel(lower(PhiMatrixArrows)~lower(genDistArrow)+lower(spaceArrow), mrank=T)
genes_biome_mantel <- mantel(lower(PhiMatrixArrows)~lower(genDistArrow)+lower(BiomeDistArrow), mrank=T)

#partial mantel test testing for effect of GEOGRAPHY
space_genes_mantel <- mantel(lower(PhiMatrixArrows)~lower(spaceArrow)+lower(genDistArrow), mrank=T)
space_biome_mantel <- mantel(lower(PhiMatrixArrows)~lower(spaceArrow)+lower(BiomeDistArrow), mrank=T)

#partial mantel test testing for effect of BIOME
biome_genes_mantel <- mantel(lower(PhiMatrixArrows)~lower(BiomeDistArrow)+lower(genDistArrow), mrank=T)
biome_space_mantel <- mantel(lower(PhiMatrixArrows)~lower(BiomeDistArrow)+lower(spaceArrow), mrank=T)


mantel_arrows <- rbind(space_mantel, genes_mantel,
                       biome_mantel,
                       genes_space_mantel, genes_biome_mantel,
                       space_genes_mantel, space_biome_mantel,
                       biome_space_mantel, biome_genes_mantel)

write.csv(mantel_arrows, "results_masked/mantel_arrows.csv")

mantel_arrows <-  rbind(space_mantel, genes_mantel,
                       biome_mantel, genes_ibd_mantel, genes_Bantu_mantel)
pval_adj_arrows <- p.adjust(mantel_arrows[,2], method = "BH")
mantel_arrows<- cbind(mantel_arrows, pval_adj_arrows)

write.csv(mantel_arrows, "final_paper_files/results_final_paper/results_masked/mantel_arrows.csv")

matrix_arrows <- rbind(genes_space_mantel, genes_biome_mantel,
                      space_genes_mantel, space_biome_mantel,
                      biome_space_mantel, biome_genes_mantel)
pval_adj_arrows <- p.adjust(matrix_arrows[,2], method = "BH")
matrix_arrows<- cbind(matrix_arrows, pval_adj_arrows)

write.csv(matrix_arrows, "results_masked/matrix_arrows.csv")

# Perform PCoA
#gowarrow <- as.data.frame(as.matrix(gowarrow))
# get NA
#na_rows <- which(is.na(gowarrow), arr.ind=TRUE)

# remove rows with NAs
#na_rows <- unique(na_rows[,1])
#gowarrow <- gowarrow[-na_rows,-na_rows]
gowarrow <- as.dist(gowarrow)
arrows_pcoa <- pcoa(gowarrow, correction = "cailliez")

#We rescale the PCoA components in relation to the explained variance.  

for(i in 1:ncol(arrows_pcoa$vectors)) { 
  arrows_pcoa$vectors[,i] <- scale(arrows_pcoa$vectors[,i])*
    arrows_pcoa$values$Rel_corr_eig[i]}

# Extract eigenvalues from PC/PCo results_masked  
arrows_ev <- arrows_pcoa$values$Corr_eig

# Plot
library(ggplot2)
my_font <- "sans"

png(file="results_masked/pcoaarrows.png",
    width=6, height=4.5, units="in", res=500)
plot_pc_scree(eigenval=arrows_ev, title="Arrows", type="PCo")
dev.off()

##We extract the principal components/coordinates for all factors: 
#genetics, grammar, music and phonology. 
#We normalize the PCs/PCos to a range from 0 to 1 and then plot a heatmap 

# Extract the PCs and the PCos from the PCA and PCoA results_masked
CultureList<-groupsarrow
arrows_pco <- as.data.frame(arrows_pcoa$vectors)
arrows_pco$Culture <- CultureList

## PCoA with groups

#install.packages('ggfortify')

# remove Mbendjele from palette
my_pal_long_ar <- my_pal_long[-10]

p <- ggplot(data = arrows_pco, aes(x = Axis.1, y = Axis.2, color=Culture)) +
  geom_point(size = 3.5) +
  xlab("PCo1 (22.38%)") +
  ylab("PCo2 (14.32%)") + 
  theme_classic() + 
  scale_color_manual(values = my_pal_long_ar,
                     labels = c("Aka (6)",
                                "Babongo (1)",
                                "Baka (7)",
                                "Bakola (5)",
                                "Bakoya (1)",
                                "East Batwa (16)",
                                "West Batwa (6)",
                                "Bedzan (2)",
                                "Efe (26)",
                                "Sua (21)")) +
  geom_hline(yintercept = 0, lty = 2) +
  ggtitle("Arrows") +
  geom_vline(xintercept = 0, lty = 2) 

p <- p + theme(legend.title = element_blank(),
               legend.text = element_text(size = 12))

ggsave(p, filename="results_masked/PCOA_arrows.png", height=150, width=180, units="mm")

#We compute Euclidean distance matrices from the dimensionality-reduced data and create NeighborNets. 
# SplitsTree doesnt accept brackets in names
arrows_pco_rn <- arrows_pco
arrows_pco_rn$Culture <- as.character(arrows_pco_rn$Culture)
arrows_pco_rn$Culture[arrows_pco_rn$Culture=="Batwa (East)"] <- "BatwaE"
arrows_pco_rn[arrows_pco_rn=="Batwa (West)"] <- "BatwaW"
arrows_pco_rn$Culture <- NULL

# Write nexus files for splitstree (Not used so far)
require(phangorn)

dist_list <- lapply(list(arrows=arrows_pco_rn),
                    function(f) {dist(f)})
writeDist(dist_list$arrows, file="arrows.nex", format="nexus")

########
## RF ##
########

# Divide cultures into Western and Eastern Pygmies 

# Question - should we separate between different regions???
arrows <- read.csv("arrows.csv")

# Remove traits that are only present in one of the two regions (East or West)
arrows$Region<-sapply(as.character(arrows$Culture),function(x)
  ifelse(x %in% c("Aka", "Babongo", "Baka", "Bakola", "Bakoya",
                  "Bedzan", "Mbendjele", "Batwa (West)"),"West",
         ifelse(x %in% c("Batwa (East)","Efe","Sua"),"East")))
table(arrows$Region)

arrow_redux<-arrows[,c(3:ncol(arrows))]
for(i in 1:(ncol(arrow_redux)-1)){
  tab<-table(arrow_redux[,c("Region",colnames(arrow_redux)[i])])  
  check<-apply(tab,2,function(x) ifelse(sum(x!=0)==1,"fix","ok"))  
  arrow_redux[,colnames(arrow_redux)[i]]<-sapply(arrow_redux[,colnames(arrow_redux)[i]],
                                               function(x) ifelse(x %in% names(check)[check=="fix"],NA,as.character(x)))
  arrow_redux[,colnames(arrow_redux)[i]]<-as.factor(arrow_redux[,colnames(arrow_redux)[i]])}
arrow_redux<-arrow_redux[,sapply(arrow_redux, function(x)!all(is.na(x)))]
arrow_redux$Region<-as.factor(arrow_redux$Region)

library(party)
library(plyr)
weights_arrows<-sapply(arrow_redux$Region,
                      function(x) nrow(arrow_redux)/nrow(arrow_redux[arrow_redux$Region==x,]))
rf_redux<-cforest(Region~.,
                  controls=
                    cforest_control(ntree = 500,
                                    mtry=5,
                                    mincriterion=qnorm(0.9),
                                    fraction = 0.632,
                                    testtype ="Teststatistic",
                                    teststat="max",
                                    replace=TRUE,
                                    trace=F,
                                    savesplitstats=F,
                                    minsplit=20,
                                    minbucket=8),
                  data=arrow_redux,
                  weights = weights_arrows)

# Produce confusion matrix
cmat_redux<-table(arrow_redux$Region, predict(rf_redux))
print(cmat_redux)

stats_region<-ldply(c(1:nrow(cmat_redux)),
                    function(x) data.frame(id=rownames(cmat_redux)[x],
                                           N=sum(cmat_redux[,x]),
                                           d=cmat_redux[x,x],
                                           A=sum(cmat_redux[x,]))) %>% transform(p=d/N,
                                                                                 r=d/A)

write.csv(stats_region, "results_masked/stats_region_arrows.csv")


arrow_redux$predReg<-predict(rf_redux)
missclass<-cbind(arrows[rownames(arrow_redux[arrow_redux$Region!=arrow_redux$predReg,]),],
                 data.frame(predReg=arrow_redux$predReg[arrow_redux$Region!=arrow_redux$predReg]))
missclass
write.csv(missclass, "results_masked/arrows_missclass.csv")

##########
## FAMD ##
##########

arrows <- read.csv("arrows.csv", stringsAsFactors = T,header=T)

#arrows$Region <- NULL
arrows[sapply(arrows, is.character)] <- lapply(arrows[sapply(arrows, is.character)], 
                                             as.factor)
imp.arrows <- imputeFAMD(arrows[,2:8], ncp = 10)
res.famd <- FAMD(imp.arrows$completeObs,   ## Set the target variable "Churn" as a supplementary variable, so it is not included in the analysis for now
                 graph = T, 
                 sup.var=1, # Culture is the supplementary variable
                 ncp=10)

## Inspect principal components
get_eigenvalue(res.famd)
fviz_screeplot(res.famd)
# Plot of variables
fviz_famd_var(res.famd, repel = TRUE, ggtheme=theme_classic())
ggsave("results_masked/arrows_FAMD_varcont_dim.png", width=4.5, height=4, units="in", dpi=500)
# Contribution to the first dimension

fviz_contrib(res.famd, "var", axes = 1, ggtheme=theme_classic()) 
ggsave("results_masked/arrows_FAMD_varDIM1.png", width=7, height=4.5, units="in", dpi=500)
# Contribution to the second dimension
fviz_contrib(res.famd, "var", axes = 2, ggtheme=theme_classic())
ggsave("results_masked/arrows_FAMD_varDIM2.png", width=7, height=4.5, units="in", dpi=500)

fviz_famd_var(res.famd, "quanti.var", col.var = "contrib", 
              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
              repel = TRUE)

my_pal_long_ar <- my_pal_long[-10]
famd_plot_arrows <- fviz_mfa_ind(res.famd, 
                          habillage = 1, # color by groups 
                          palette = my_pal_long_ar,
                          addEllipses = FALSE, 
                          repel = TRUE, # Avoid text overlapping,
                          label="none",
                          ggtheme = theme_classic(),
) 

famd_plot_arrows <- famd_plot_arrows + ggtitle("Arrows - FAMD")
famd_plot_arrows
ggsave("results_masked/Arrows_FAMD_noell.tiff", width=7, height=6, units="in", dpi=500)


################
## AEROPHONES ##
################

aerophones <- read.csv("aerophones_5.csv", stringsAsFactors = T,header=T)
ids <- seq(1:nrow(aerophones))
for (i in ids) {
  ids[i] <- paste("T-", i, sep="")
}

aerophones$ID <- ids

# For genetics, create another aerophone dataset with only aerophones from groups with genes
aerophones_red <- subset(aerophones, aerophones$Culture %in% rownames(Fst_mat))

aerophones$Trait <- NULL
ids <- aerophones$ID
ids_red <- aerophones_red$ID
rownames(aerophones) <- ids
rownames(aerophones_red) <- ids_red
aerophones$ID <- NULL
aerophones$Culture <- NULL
aerophones_red$Trait <- NULL
aerophones_red$ID <- NULL
aerophones_red$Culture <- NULL

# this file needs to have the coordinates of the culture corresponding to each individual aerophone
# in the same order as your aerophones.csv file
locationsaerophone <- read.csv("aerophones_locations.csv")
locationsaerophone_red <- subset(locationsaerophone, locationsaerophone$Culture %in% rownames(Fst_mat))

rownames(locationsaerophone) <- rownames(aerophones)
rownames(locationsaerophone_red) <- rownames(aerophones_red)

biomesaerophone <- read.csv("aerophones_biomes.csv")
biomesaerophone_red <- subset(biomesaerophone, biomesaerophone$Culture %in% rownames(Fst_mat))
biomesaerophone$Trait <- NULL
biomesaerophone$Culture <- NULL
biomesaerophone_red$Trait <- NULL
biomesaerophone_red$Culture <- NULL

rownames(biomesaerophone) <- rownames(aerophones)
rownames(biomesaerophone_red) <- rownames(aerophones_red)

# Convert factors into numeric codes (not really needed if stringsAsFactors = F)
indx <- sapply(aerophones, is.factor)
indx_red <- sapply(aerophones_red, is.factor)

aerophones[indx] <- lapply(aerophones[indx],  function(x) as.numeric(as.factor(x)))
aerophones_red[indx_red] <- lapply(aerophones_red[indx_red],  function(x) as.numeric(as.factor(x)))

#compute Gower dissimilarity distance
library(vegan)
gowaerophone <- vegdist(aerophones,"gower",na.rm=TRUE) 
gowaerophoneRed <- vegdist(aerophones_red,"gower",na.rm=TRUE) 

gowBiomeaerophones <- as.matrix(vegdist(biomesaerophone,"gower",na.rm=TRUE))
gowBiomeaerophonesRed <-  as.matrix(vegdist(biomesaerophone_red,"gower",na.rm=TRUE))


#compute binary Matrices of cultural affiliation (0:same culture; 1:different culture)
binaryCultDistaerophones<-as.matrix(as.dist(BinaryCulture(locationsaerophone$Culture)))
binaryCultDistaerophonesRed<-as.matrix(as.dist(BinaryCulture(locationsaerophone_red$Culture)))

#compute Fst matrix
GenDistaerophones<-as.matrix(as.dist(GetFst(locationsaerophone_red$Culture)))

#compute spatial distance between cultures and therefore aerophones
#spaceaerophones<-distMat(data.frame(x=locationsaerophone$Longitude,y=locationsaerophone$Latitude))

#compute spatial distances
distTableaerophones <- data.frame(x=locationsaerophone$Longitude,y=locationsaerophone$Latitude)
namesaerophone <- as.vector(row.names(locationsaerophone))
rownames(distTableaerophones) <- namesaerophone

distTableaerophonesRed <- data.frame(x=locationsaerophone_red$Longitude,y=locationsaerophone_red$Latitude)
namesaerophoneRed <- as.vector(row.names(locationsaerophone_red))
rownames(distTableaerophonesRed) <- namesaerophoneRed

library(geosphere)
spaceaerophones <- distm(distTableaerophones, fun=distGeo)
rownames(spaceaerophones) <- namesaerophone
colnames(spaceaerophones) <- namesaerophone
# here distance is in metres, not kilometres
#spaceaerophones <- as.dist(spaceaerophones)

spaceaerophonesRed <- distm(distTableaerophonesRed, fun=distGeo)
rownames(spaceaerophonesRed) <- namesaerophoneRed
colnames(spaceaerophonesRed) <- namesaerophoneRed
# here distance is in metres, not kilometres
#spaceaerophonesRed <- as.dist(spaceaerophonesRed)

gowaerophone <- as.matrix(gowaerophone)
#na_rows <- which(is.na(gowaerophone), arr.ind=TRUE)
#na_rows <- unique(na_rows[,1])
#gowaerophone <- gowaerophone[-na_rows, -na_rows]
#spaceaerophones <- spaceaerophones[-na_rows, -na_rows]
#gowBiomeaerophones <- gowBiomeaerophones[-na_rows, -na_rows]
#binaryCultDistaerophones <- binaryCultDistaerophones[-na_rows, -na_rows]

gowaerophoneRed <- as.matrix(gowaerophoneRed)
#na_rows <- which(is.na(gowaerophoneRed), arr.ind=TRUE)
#na_rows <- unique(na_rows[,1])
#gowaerophoneRed <- gowaerophoneRed[-na_rows, -na_rows]
#spaceaerophonesRed  <- spaceaerophonesRed [-na_rows, -na_rows]
#gowBiomeaerophonesRed  <- gowBiomeaerophonesRed [-na_rows, -na_rows]
#binaryCultDistaerophonesRed  <- binaryCultDistaerophonesRed [-na_rows, -na_rows]
#GenDistaerophones <- GenDistaerophones[-na_rows, -na_rows]


#partial mantel test between Distance and Geography and between Culture and Geography
library(ecodist)
detach("package:vegan", unload = TRUE)

space_mantel <- mantel(lower(gowaerophone)~lower(spaceaerophones),nperm=1000)
culture_mantel <-mantel(lower(gowaerophone)~lower(binaryCultDistaerophones),nperm=1000)
genes_mantel <- mantel(lower(gowaerophoneRed)~lower(GenDistaerophones),nperm=1000)
biome_mantel <- mantel(lower(gowBiomeaerophonesRed)~lower(GenDistaerophones),nperm=1000)

#partial mantel test testing for effect of CULTURAL AFFILIATION controlling for rest
cult_space_mantel <- mantel(lower(gowaerophone)~lower(binaryCultDistaerophones)+lower(spaceaerophones),nperm=1000)
cult_genes_mantel <- mantel(lower(gowaerophoneRed)~lower(binaryCultDistaerophonesRed)+lower(GenDistaerophones),nperm=1000)
cult_biome_mantel <- mantel(lower(gowaerophone)~lower(binaryCultDistaerophones)+lower(gowBiomeaerophones),nperm=1000)

#partial mantel test testing for effect of SHARED GENES
genes_space_mantel <- mantel(lower(gowaerophoneRed)~lower(GenDistaerophones)+lower(spaceaerophonesRed),nperm=1000)
genes_cult_mantel <- mantel(lower(gowaerophoneRed)~lower(GenDistaerophones)+lower(binaryCultDistaerophonesRed),nperm=1000)
genes_biome_mantel <- mantel(lower(gowaerophoneRed)~lower(GenDistaerophones)+lower(gowBiomeaerophonesRed),nperm=1000)

#partial mantel test testing for effect of GEOGRAPHY
space_culture_mantel <- mantel(lower(gowaerophone)~lower(spaceaerophones)+lower(binaryCultDistaerophones),nperm=1000)
space_genes_mantel <- mantel(lower(gowaerophoneRed)~lower(spaceaerophonesRed)+lower(GenDistaerophones),nperm=1000)
space_biome_mantel <- mantel(lower(gowaerophone)~lower(spaceaerophones)+lower(gowBiomeaerophones),nperm=1000)

#partial mantel test testing for effect of BIOME
biome_culture_mantel <- mantel(lower(gowaerophone)~lower(gowBiomeaerophones)+lower(binaryCultDistaerophones),nperm=1000)
biome_genes_mantel <- mantel(lower(gowaerophoneRed)~lower(gowBiomeaerophonesRed)+lower(GenDistaerophones),nperm=1000)
biome_space_mantel <- mantel(lower(gowaerophone)~lower(gowBiomeaerophones)+lower(spaceaerophones),nperm=1000)

write.csv(rbind(space_mantel, culture_mantel, genes_mantel,biome_mantel,
                cult_space_mantel,cult_genes_mantel,cult_biome_mantel,
                genes_space_mantel, genes_cult_mantel,genes_biome_mantel,
                space_culture_mantel,space_genes_mantel,space_biome_mantel,
                biome_culture_mantel, biome_space_mantel, biome_genes_mantel), "results_masked/aerophones_mantel.csv")

#########
# AMOVA #
#########
library(pegas)

#AMOVA
sqgowaerophone<-sqrt(gowaerophone)
sqgowaerophone_mat <- as.data.frame(as.matrix(sqgowaerophone))
# # get NA
# na_rows <- which(is.na(sqgowaerophone_mat), arr.ind=TRUE)
# 
# # remove rows with NAs
# na_rows <- unique(na_rows[,1])
# sqgowaerophone_mat <- sqgowaerophone_mat[-na_rows,-na_rows]
# sqgowaerophone <- as.dist(sqgowaerophone_mat)

groupsaerophone<-as.factor(locationsaerophone$Culture)
#groupsaerophone <- as.factor(groupsaerophone[-na_rows])
AMOVAaerophone=pegas::amova(sqgowaerophone~groupsaerophone)

# Phi-Stat:
AMOVAaerophone$varcomp$sigma2[1]/sum(AMOVAaerophone$varcomp$sigma2)
#P-Value:
AMOVAaerophone$varcomp$P.value[1]


##################
# pairwise AMOVA #
##################

# Note: This is the important thing to compare against genetic data

# I think you cannot do amova for cultures with only 1 trait!

#Extract Culture Data
CultureNames<-unique(as.character(groupsaerophone))
CultureNames <- sort(CultureNames)
NumberOfCultures<-length(CultureNames)
sqgowaerophone<-as.matrix(sqgowaerophone)
CultureList<-groupsaerophone

#PhiSt and Pval Output Matrices:
PhiMatrix<-matrix(NA,nrow=NumberOfCultures,ncol=NumberOfCultures)
PvalMatrix<-matrix(NA,nrow=NumberOfCultures,ncol=NumberOfCultures)

#Spatial Coordinates Output Data.Frame:
XY<-data.frame(x=numeric(length=NumberOfCultures),y=numeric(length=NumberOfCultures))

for (i in 1:NumberOfCultures) {
  print(paste(i," out of ",NumberOfCultures))
  for (j in 1:NumberOfCultures)
  {
    if (i!=j&i>j)
    { 
      indexI=which(as.character(CultureList)==as.character(CultureNames[i]))
      indexJ=which(as.character(CultureList)==as.character(CultureNames[j]))
      #subsetDistance Matrix
      subDistMat<-sqgowaerophone[c(indexI,indexJ),c(indexI,indexJ)]
      subDistMat<-as.dist(subDistMat)
      classes<-as.factor(c(rep("a",length(indexI)),rep("b",length(indexJ))))
      #compute Phi:
      tmp<-pegas::amova(subDistMat~classes)
      N0=tmp$varcoef
      AP=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)
      WP=(tmp$tab$MSD[2])
      PhiMatrix[i,j]=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)/((tmp$tab$MSD[2])+((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef))
      if(PhiMatrix[i,j]<0){PhiMatrix[i,j]=0}
      PvalMatrix[i,j]=tmp$varcomp$P.value[1]
    }
  }
  XY$x[i]<-mean(locationsaerophone$Longitude[which(locationsaerophone$Culture==as.character(CultureNames[i]))],na.rm=TRUE)
  XY$y[i]<-mean(locationsaerophone$Latitude[which(locationsaerophone$Culture==as.character(CultureNames[i]))],na.rm=TRUE)
}

dim <- ncol(PhiMatrix)

#Change NAs to blanks for text
PhiMatrix_text <- as.data.frame(PhiMatrix)
PhiMatrix_text <- sapply(PhiMatrix_text, as.character)
PhiMatrix_text[is.na(PhiMatrix_text)] <- " "
PhiMatrix_text <- as.matrix(PhiMatrix_text)

png(file="results_masked/aerophones_phistats.png",
    width=7, height=6, units="in", res=500)
image(-7:(dim-8), 1:dim, PhiMatrix, axes = FALSE, xlab="", ylab="", col=hcl.colors(12, "Geyser", rev=TRUE), main="aerophones", xlim=c(-7,3), ylim=c(0,11))
axis(1, (-7:(dim-8)), CultureNames, cex.axis = 0.6, las=3, pos=0.5)
axis(4, 1:dim, CultureNames, cex.axis = 0.6, las=1, pos=-0.5)
text(expand.grid((-7:(dim-8)), 1:dim), sprintf("%0.4s", PhiMatrix_text), cex=0.6)
dev.off()

# Again Partial Mantel Tests using Phi Statistic

CultureNames_rn <- CultureNames
CultureNames_rn[CultureNames_rn=="Batwa (East)"] <- "BatwaE"
CultureNames_rn[CultureNames_rn=="Batwa (West)"] <- "BatwaW"
rownames(PhiMatrix) <- CultureNames
colnames(PhiMatrix) <- CultureNames
PhiMatrix_rn <- PhiMatrix
rownames(PhiMatrix_rn) <- CultureNames_rn
colnames(PhiMatrix_rn) <- CultureNames_rn

library(phangorn)
writeDist(PhiMatrix_rn, file="aerophones_phi.nex", format="nexus")

# Mantel tests between geography, biome, genetics and aerophones
PhiMatrix <- as.matrix(PhiMatrix)
genDist <- as.matrix(genDist)
gowBiomeRed <- as.matrix(gowBiomeRed)
spaceRed <- as.matrix(spaceRed)

# Remove groups with >2 aerophones
genDistaerophone <- genDist[-c(2,4,5,7),-c(2,4,5,7)]
genDistaerophoneIBD <- genDistIBD[-c(2,4,5,7),-c(2,4,5,7)]
genDistaerophoneBantu <- genDistBantu[-c(2,4,5,7),-c(2,4,5,7)]
BiomeDistaerophone <- gowBiomeRed[-c(2,4,5,7),-c(2,4,5,7)]
spaceaerophone <- spaceRed[-c(2,4,5,7),-c(2,4,5,7)]
PhiMatrixaerophones<- PhiMatrix[-4,-4]

detach("package:vegan", unload = TRUE)

biome_mantel <- mantel(lower(PhiMatrixaerophones)~lower(BiomeDistaerophone), mrank=T)
space_mantel <- mantel(lower(PhiMatrixaerophones)~lower(spaceaerophone), mrank=T)
genes_mantel <- mantel(lower(PhiMatrixaerophones)~lower(genDistaerophone), mrank=T)
genes_ibd_mantel <- mantel(lower(PhiMatrixaerophones)~lower(genDistaerophoneIBD), mrank=T)
genes_bantu_mantel <- mantel(lower(PhiMatrixaerophones[-c(4,6),-c(4,6) ])~lower(genDistaerophoneBantu), mrank=T)

#partial mantel test testing for effect of SHARED GENES
genes_space_mantel <- mantel(lower(PhiMatrixaerophones)~lower(genDistaerophone)+lower(spaceaerophone), mrank=T)
genes_biome_mantel <- mantel(lower(PhiMatrixaerophones)~lower(genDistaerophone)+lower(BiomeDistaerophone), mrank=T)

#partial mantel test testing for effect of GEOGRAPHY
space_genes_mantel <- mantel(lower(PhiMatrixaerophones)~lower(spaceaerophone)+lower(genDistaerophone), mrank=T)
space_biome_mantel <- mantel(lower(PhiMatrixaerophones)~lower(spaceaerophone)+lower(BiomeDistaerophone), mrank=T)

#partial mantel test testing for effect of BIOME
biome_genes_mantel <- mantel(lower(PhiMatrixaerophones)~lower(BiomeDistaerophone)+lower(genDistaerophone), mrank=T)
biome_space_mantel <- mantel(lower(PhiMatrixaerophones)~lower(BiomeDistaerophone)+lower(spaceaerophone), mrank=T)


mantel_aerophones <- rbind(space_mantel, genes_mantel,
                           biome_mantel,
                           genes_space_mantel, genes_biome_mantel,
                           space_genes_mantel, space_biome_mantel,
                           biome_space_mantel, biome_genes_mantel)

write.csv(mantel_aerophones, "results_masked/mantel_aerophones.csv")


mantel_aerophones <-  rbind(space_mantel, genes_mantel,
                       biome_mantel,genes_bantu_mantel, genes_ibd_mantel)
pval_adj_aerophones <- p.adjust(mantel_aerophones[,2], method = "BH")
mantel_aerophones<- cbind(mantel_aerophones, pval_adj_aerophones)

write.csv(mantel_aerophones, "final_paper_files/results_final_paper/results_masked/mantel_aerophones.csv")

matrix_aerophones <- rbind(genes_space_mantel, genes_biome_mantel,
                      space_genes_mantel, space_biome_mantel,
                      biome_space_mantel, biome_genes_mantel)
pval_adj_aerophones <- p.adjust(matrix_aerophones[,2], method = "BH")
matrix_aerophones<- cbind(matrix_aerophones, pval_adj_aerophones)

write.csv(matrix_aerophones, "results_masked/matrix_aerophones.csv")


# Perform PCoA
gowaerophone <- as.data.frame(as.matrix(gowaerophone))
# get NA
na_rows <- which(is.na(gowaerophone), arr.ind=TRUE)

# remove rows with NAs
#na_rows <- unique(na_rows[,1])
#gowaerophone <- gowaerophone[-na_rows,-na_rows]
gowaerophone <- as.dist(gowaerophone)
aerophones_pcoa <- pcoa(gowaerophone, correction = "cailliez")

#We rescale the PCoA components in relation to the explained variance.  

for(i in 1:ncol(aerophones_pcoa$vectors)) { 
  aerophones_pcoa$vectors[,i] <- scale(aerophones_pcoa$vectors[,i])*
    aerophones_pcoa$values$Rel_corr_eig[i]}

# Extract eigenvalues from PC/PCo results_masked 
aerophones_ev <- aerophones_pcoa$values$Corr_eig

# Plot
library(ggplot2)
my_font <- "sans"

png(file="results_masked/pcoaaerophones.png",
    width=6, height=4.5, units="in", res=500)
plot_pc_scree(eigenval=aerophones_ev, title="Aerophones", type="PCo")
dev.off()

##We extract the principal components/coordinates for all factors: 
#genetics, grammar, music and phonology. 
#We normalize the PCs/PCos to a range from 0 to 1 and then plot a heatmap 

# Extract the PCs and the PCos from the PCA and PCoA results_masked
CultureList<-groupsaerophone
aerophones_pco <- as.data.frame(aerophones_pcoa$vectors)
aerophones_pco$Culture <- CultureList

## PCoA with groups

#install.packages('ggfortify')

# remove Mbendjele from palette
my_pal_long_ar <- my_pal_long[-c(2,4,5,7)]

p <- ggplot(data = aerophones_pco, aes(x = Axis.1, y = Axis.2, color=Culture)) +
  geom_point(size = 3.5) +
  xlab("PCo1 (7.66%)") +
  ylab("PCo2 (3.77%)") + 
  theme_classic() + 
  scale_color_manual(values = my_pal_long_ar,
                     labels = c("Aka (7)",
                                "Baka (3)",
                                "East Batwa (63)",
                                "West Batwa (4)",
                                "Efe (27)",
                                "Mbendjele (2)",
                                "Sua (16)")) +
  geom_hline(yintercept = 0, lty = 2) +
  ggtitle("Aerophones") +
  geom_vline(xintercept = 0, lty = 2) 

p <- p + theme(legend.title = element_blank(),
               legend.text = element_text(size = 12))

ggsave(p, filename="results_masked/PCOA_aerophones.png", height=150, width=180, units="mm")

#We compute Euclidean distance matrices from the dimensionality-reduced data and create NeighborNets. 
# SplitsTree doesnt accept brackets in names
aerophones_pco_rn <- aerophones_pco
aerophones_pco_rn$Culture <- as.character(aerophones_pco_rn$Culture)
aerophones_pco_rn$Culture[aerophones_pco_rn$Culture=="Batwa (East)"] <- "BatwaE"
aerophones_pco_rn[aerophones_pco_rn=="Batwa (West)"] <- "BatwaW"
aerophones_pco_rn$Culture <- NULL

# Write nexus files for splitstree (Not used so far)
require(phangorn)

dist_list <- lapply(list(aerophones=aerophones_pco_rn),
                    function(f) {dist(f)})
writeDist(dist_list$aerophones, file="aerophones.nex", format="nexus")

###################################################################
## Now we test for potential mechanisms of cultural transmission ##
###################################################################

# Divide cultures into Western and Eastern Pygmies 
aerophones <- read.csv("aerophones_5.csv", stringsAsFactors = T,header=T)

# Question - should we separate between different regions???

# Remove traits that are only present in one of the two regions (East or West)
aerophones$Region<-sapply(as.character(aerophones$Culture),function(x)
  ifelse(x %in% c("Aka", "Babongo", "Baka", "Bakola", "Bakoya",
                  "Bedzan", "Mbendjele", "Batwa (West)"),"West",
         ifelse(x %in% c("Batwa (East)","Efe","Sua"),"East")))
table(aerophones$Region)

aerophone_redux<-aerophones[,c(3:ncol(aerophones))]
for(i in 1:(ncol(aerophone_redux)-1)){
  tab<-table(aerophone_redux[,c("Region",colnames(aerophone_redux)[i])])  
  check<-apply(tab,2,function(x) ifelse(sum(x!=0)==1,"fix","ok"))  
  aerophone_redux[,colnames(aerophone_redux)[i]]<-sapply(aerophone_redux[,colnames(aerophone_redux)[i]],
                                                         function(x) ifelse(x %in% names(check)[check=="fix"],NA,as.character(x)))
  aerophone_redux[,colnames(aerophone_redux)[i]]<-as.factor(aerophone_redux[,colnames(aerophone_redux)[i]])}
aerophone_redux<-aerophone_redux[,sapply(aerophone_redux, function(x)!all(is.na(x)))]
aerophone_redux$Region<-as.factor(aerophone_redux$Region)

library(party)
library(plyr)
weights_aerophones<-sapply(aerophone_redux$Region,
                           function(x) nrow(aerophone_redux)/nrow(aerophone_redux[aerophone_redux$Region==x,]))
rf_redux<-cforest(Region~.,
                  controls=
                    cforest_control(ntree = 500,
                                    mtry=6,
                                    mincriterion=qnorm(0.9),
                                    fraction = 0.632, #the usual,
                                    testtype ="Teststatistic",
                                    teststat="max",
                                    replace=TRUE,
                                    trace=F,
                                    savesplitstats=F,
                                    minsplit=20,
                                    minbucket=8),
                  data=aerophone_redux,
                  weights = weights_aerophones)

# Produce confusion matrix
cmat_redux<-table(aerophone_redux$Region, predict(rf_redux))
print(cmat_redux)

stats_region<-ldply(c(1:nrow(cmat_redux)),
                    function(x) data.frame(id=rownames(cmat_redux)[x],
                                           N=sum(cmat_redux[,x]),
                                           d=cmat_redux[x,x],
                                           A=sum(cmat_redux[x,]))) %>% transform(p=d/N,
                                                                                 r=d/A)
stats_region
write.csv(stats_region, "results_masked/stats_region_aerophones.csv")

aerophone_redux$predReg<-predict(rf_redux)
missclass<-cbind(aerophones[rownames(aerophone_redux[aerophone_redux$Region!=aerophone_redux$predReg,]),],
                 data.frame(predReg=aerophone_redux$predReg[aerophone_redux$Region!=aerophone_redux$predReg]))
missclass

## Import libraries
library(FactoMineR)
library(factoextra)
library(missMDA)

## FAMD

aerophones$Region <- NULL
imp.aerphones <- imputeFAMD(aerophones[,2:15], ncp = 10)
res.famd <- FAMD(imp.aerphones$completeObs,   ## Set the target variable "Churn" as a supplementary variable, so it is not included in the analysis for now
                 graph = T, 
                 sup.var=1, # Culture is the supplementary variable
                 ncp=10)

## Inspect principal components
get_eigenvalue(res.famd)
fviz_screeplot(res.famd)
# Plot of variables
fviz_famd_var(res.famd, repel = TRUE, ggtheme=theme_classic())
ggsave("results_masked/aerophones_FAMD_varcont_dim.png", width=4.5, height=4, units="in", dpi=500)
# Contribution to the first dimension

fviz_contrib(res.famd, "var", axes = 1, ggtheme=theme_classic()) 
ggsave("results_masked/Aerophones_FAMD_varDIM1.png", width=7, height=4.5, units="in", dpi=500)
# Contribution to the second dimension
fviz_contrib(res.famd, "var", axes = 2, ggtheme=theme_classic())
ggsave("results_masked/Aerophones_FAMD_varDIM2.png", width=7, height=4.5, units="in", dpi=500)

fviz_famd_var(res.famd, "quanti.var", col.var = "contrib", 
              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
              repel = TRUE)
my_pal_long_ar <- my_pal_long[-c(2,4,5,7)]
famd_plot_aerophones <- fviz_mfa_ind(res.famd, 
             habillage = 1, # color by groups 
             palette = my_pal_long_ar,
             addEllipses = FALSE, 
             repel = TRUE, # Avoid text overlapping,
             label="none",
             ggtheme = theme_classic(),
) 

famd_plot_aerophones <- famd_plot_aerophones + ggtitle("Aerophones - FAMD")
famd_plot_aerophones
ggsave("results_masked/Aerophones_FAMD_noell.tiff", width=7, height=6, units="in", dpi=500)

############
## FLUTES ##
############

flutes <- read.csv("flutes.csv", stringsAsFactors = T,header=T)
ids <- seq(1:nrow(flutes))
for (i in ids) {
  ids[i] <- paste("T-", i, sep="")
}

flutes$ID <- ids

# For genetics, create another flute dataset with only flutes from groups with genes
flutes_red <- subset(flutes, flutes$Culture %in% rownames(Fst_mat))

flutes$Trait <- NULL
ids <- flutes$ID
ids_red <- flutes_red$ID
groupsflute<-as.factor(flutes$Culture)
groupsflute_red<-as.factor(flutes$Culture)
rownames(flutes) <- ids
rownames(flutes_red) <- ids_red
flutes$ID <- NULL
flutes$Culture <- NULL
flutes_red$Trait <- NULL
flutes_red$ID <- NULL
flutes_red$Culture <- NULL

# Convert factors into numeric codes (not really needed if stringsAsFactors = F)
indx <- sapply(flutes, is.factor)
indx_red <- sapply(flutes_red, is.factor)

flutes[indx] <- lapply(flutes[indx],  function(x) as.numeric(as.factor(x)))
flutes_red[indx_red] <- lapply(flutes_red[indx_red],  function(x) as.numeric(as.factor(x)))

#compute Gower dissimilarity distance
library(vegan)
gowflute <- vegdist(flutes,"gower",na.rm=TRUE) 
gowfluteRed <- vegdist(flutes_red,"gower",na.rm=TRUE) 


#########
# AMOVA #
#########
library(pegas)

#AMOVA
sqgowflute<-sqrt(gowflute)
sqgowflute_mat <- as.data.frame(as.matrix(sqgowflute))
# # get NA
#na_rows <- which(is.na(sqgowflute_mat), arr.ind=TRUE)
# 
# # remove rows with NAs
# na_rows <- unique(na_rows[,1])
# sqgowflute_mat <- sqgowflute_mat[-na_rows,-na_rows]
# sqgowflute <- as.dist(sqgowflute_mat)

#groupsflute <- as.factor(groupsflute[-na_rows])
AMOVAflute=amova(sqgowflute~groupsflute)

# Phi-Stat:
AMOVAflute$varcomp$sigma2[1]/sum(AMOVAflute$varcomp$sigma2)
#P-Value:
AMOVAflute$varcomp$P.value[1]


##################
# pairwise AMOVA #
##################

# Note: This is the important thing to compare against genetic data

# I think you cannot do amova for cultures with only 1 trait!

#Extract Culture Data
CultureNames<-unique(as.character(groupsflute))
CultureNames <- sort(CultureNames)
NumberOfCultures<-length(CultureNames)
sqgowflute<-as.matrix(sqgowflute)
CultureList<-groupsflute

#PhiSt and Pval Output Matrices:
PhiMatrix<-matrix(NA,nrow=NumberOfCultures,ncol=NumberOfCultures)
PvalMatrix<-matrix(NA,nrow=NumberOfCultures,ncol=NumberOfCultures)

#Spatial Coordinates Output Data.Frame:
XY<-data.frame(x=numeric(length=NumberOfCultures),y=numeric(length=NumberOfCultures))

for (i in 1:NumberOfCultures) {
  print(paste(i," out of ",NumberOfCultures))
  for (j in 1:NumberOfCultures)
  {
    if (i!=j&i>j)
    { 
      indexI=which(as.character(CultureList)==as.character(CultureNames[i]))
      indexJ=which(as.character(CultureList)==as.character(CultureNames[j]))
      #subsetDistance Matrix
      subDistMat<-sqgowflute[c(indexI,indexJ),c(indexI,indexJ)]
      subDistMat<-as.dist(subDistMat)
      classes<-as.factor(c(rep("a",length(indexI)),rep("b",length(indexJ))))
      #compute Phi:
      tmp<-pegas::amova(subDistMat~classes)
      N0=tmp$varcoef
      AP=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)
      WP=(tmp$tab$MSD[2])
      PhiMatrix[i,j]=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)/((tmp$tab$MSD[2])+((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef))
      if(PhiMatrix[i,j]<0){PhiMatrix[i,j]=0}
      PvalMatrix[i,j]=tmp$varcomp$P.value[1]
    }
  }
  #XY$x[i]<-mean(locationsflute$Longitude[which(locationsflute$Culture==as.character(CultureNames[i]))],na.rm=TRUE)
  #XY$y[i]<-mean(locationsflute$Latitude[which(locationsflute$Culture==as.character(CultureNames[i]))],na.rm=TRUE)
}

dim <- ncol(PhiMatrix)

#Change NAs to blanks for text
PhiMatrix_text <- as.data.frame(PhiMatrix)
PhiMatrix_text <- sapply(PhiMatrix_text, as.character)
PhiMatrix_text[is.na(PhiMatrix_text)] <- " "
PhiMatrix_text <- as.matrix(PhiMatrix_text)

png(file="results_masked/flutes_phistats.png",
    width=7, height=6, units="in", res=500)
image(-7:(dim-8), 1:dim, PhiMatrix, axes = FALSE, xlab="", ylab="", col=hcl.colors(12, "Geyser", rev=TRUE), main="Flutes", xlim=c(-7,1), ylim=c(0,11))
axis(1, (-7:(dim-8)), CultureNames, cex.axis = 0.6, las=3, pos=0.5)
axis(4, 1:dim, CultureNames, cex.axis = 0.6, las=1, pos=-1.5)
text(expand.grid((-7:(dim-8)), 1:dim), sprintf("%0.4s", PhiMatrix_text), cex=0.6)
dev.off()

# Again Partial Mantel Tests using Phi Statistic

CultureNames_rn <- CultureNames
CultureNames_rn[CultureNames_rn=="Batwa (East)"] <- "BatwaE"
CultureNames_rn[CultureNames_rn=="Batwa (West)"] <- "BatwaW"
rownames(PhiMatrix) <- CultureNames
colnames(PhiMatrix) <- CultureNames
PhiMatrix_rn <- PhiMatrix
rownames(PhiMatrix_rn) <- CultureNames_rn
colnames(PhiMatrix_rn) <- CultureNames_rn

writeDist(PhiMatrix_rn, file="flutes_phi.nex", format="nexus")

# Mantel tests between geography, biome, genetics and flutes
PhiMatrix <- as.matrix(PhiMatrix)
genDist <- as.matrix(genDist)
gowBiomeRed <- as.matrix(gowBiomeRed)
spaceRed <- as.matrix(spaceRed)

# Remove groups with >2 flutes
genDistflute <- genDist[-c(2,4,5,7),-c(2,4,5,7)]
BiomeDistflute <- gowBiomeRed[-c(2,4,5,7),-c(2,4,5,7)]
spaceflute <- spaceRed[-c(2,4,5,7),-c(2,4,5,7)]
#PhiMatrix<- as.dist(PhiMatrix[-4,-4])
PhiMatrixflutes <- PhiMatrix

biome_mantel <- mantel(lower(PhiMatrixflutes)~lower(BiomeDistflute), mrank=T)
space_mantel <- mantel(lower(PhiMatrixflutes)~lower(spaceflute), mrank=T)
genes_mantel <- mantel(lower(PhiMatrixflutes)~lower(genDistflute), mrank=T)

#partial mantel test testing for effect of SHARED GENES
genes_space_mantel <- mantel(lower(PhiMatrixflutes)~lower(genDistflute)+lower(spaceflute), mrank=T)
genes_biome_mantel <- mantel(lower(PhiMatrixflutes)~lower(genDistflute)+lower(BiomeDistflute), mrank=T)

#partial mantel test testing for effect of GEOGRAPHY
space_genes_mantel <- mantel(lower(PhiMatrixflutes)~lower(spaceflute)+lower(genDistflute), mrank=T)
space_biome_mantel <- mantel(lower(PhiMatrixflutes)~lower(spaceflute)+lower(BiomeDistflute), mrank=T)

#partial mantel test testing for effect of BIOME
biome_genes_mantel <- mantel(lower(PhiMatrixflutes)~lower(BiomeDistflute)+lower(genDistflute), mrank=T)
biome_space_mantel <- mantel(lower(PhiMatrixflutes)~lower(BiomeDistflute)+lower(spaceflute), mrank=T)


mantel_flutes <- rbind(space_mantel, genes_mantel,
                       biome_mantel,
                       genes_space_mantel, genes_biome_mantel,
                       space_genes_mantel, space_biome_mantel,
                       biome_space_mantel, biome_genes_mantel)

write.csv(mantel_flutes, "results_masked/mantel_flutes.csv")

# Perform PCoA
gowflute <- as.data.frame(as.matrix(gowflute))
# get NA
#na_rows <- which(is.na(gowflute), arr.ind=TRUE)

# remove rows with NAs
#na_rows <- unique(na_rows[,1])
#gowflute <- gowflute[-na_rows,-na_rows]
#gowflute <- as.dist(gowflute)
flutes_pcoa <- pcoa(gowflute, correction = "cailliez")

#We rescale the PCoA components in relation to the explained variance.  

for(i in 1:ncol(flutes_pcoa$vectors)) { 
  flutes_pcoa$vectors[,i] <- scale(flutes_pcoa$vectors[,i])*
    flutes_pcoa$values$Rel_corr_eig[i]}

# Extract eigenvalues from PC/PCo results_masked 
flutes_ev <- flutes_pcoa$values$Corr_eig

# Plot
library(ggplot2)
my_font <- "sans"

png(file="results_masked/pcoaflutes.png",
    width=6, height=4.5, units="in", res=500)
plot_pc_scree(eigenval=flutes_ev, title="Genetics", type="PCo")
dev.off()

##We extract the principal components/coordinates for all factors: 
#genetics, grammar, music and phonology. 
#We normalize the PCs/PCos to a range from 0 to 1 and then plot a heatmap 

# Extract the PCs and the PCos from the PCA and PCoA results_masked
CultureList<-groupsflute
flutes_pco <- as.data.frame(flutes_pcoa$vectors)
flutes_pco$Culture <- CultureList

## PCoA with groups

#install.packages('ggfortify')

# remove Mbendjele from palette
my_pal_long_ar <- my_pal_long[-c(2,4,5,7,8)]

p <- ggplot(data = flutes_pco, aes(x = Axis.1, y = Axis.2, color=Culture)) +
  geom_point(size = 3.5) +
  xlab("PCo1 (1.26%)") +
  ylab("PCo2 (0.90%)") + 
  theme_classic() + 
  scale_color_manual(values = my_pal_long_ar,
                     labels = c("Aka (5)",
                                "Baka (3)",
                                "East Batwa (4)",
                                "Efe (15)",
                                "Mbendjele (2)",
                                "Sua (4)")) +
  geom_hline(yintercept = 0, lty = 2) +
  ggtitle("flutes") +
  geom_vline(xintercept = 0, lty = 2) 

p <- p + theme(legend.title = element_blank(),
               legend.text = element_text(size = 12))

ggsave(p, filename="results_masked/PCOA_flutes.png", height=150, width=180, units="mm")

#We compute Euclidean distance matrices from the dimensionality-reduced data and create NeighborNets. 
# SplitsTree doesnt accept brackets in names
flutes_pco_rn <- flutes_pco
flutes_pco_rn$Culture <- as.character(flutes_pco_rn$Culture)
flutes_pco_rn$Culture[flutes_pco_rn$Culture=="Batwa (East)"] <- "BatwaE"
flutes_pco_rn[flutes_pco_rn=="Batwa (West)"] <- "BatwaW"
flutes_pco_rn$Culture <- NULL

# Write nexus files for splitstree (Not used so far)
require(phangorn)

dist_list <- lapply(list(flutes=flutes_pco_rn),
                    function(f) {dist(f)})
writeDist(dist_list$flutes, file="flutes.nex", format="nexus")

#############
## BASKETS ##
#############

baskets <- read.csv("baskets.csv", stringsAsFactors = T,header=T)
ids <- seq(1:nrow(baskets))
for (i in ids) {
  ids[i] <- paste("T-", i, sep="")
}

baskets$ID <- ids

# For genetics, create another basket dataset with only baskets from groups with genes
baskets_red <- subset(baskets, baskets$Culture %in% rownames(Fst_mat))

baskets$Trait <- NULL
ids <- baskets$ID
ids_red <- baskets_red$ID
rownames(baskets) <- ids
rownames(baskets_red) <- ids_red
baskets$ID <- NULL
baskets$Culture <- NULL
baskets_red$Trait <- NULL
baskets_red$ID <- NULL
baskets_red$Culture <- NULL

# this file needs to have the coordinates of the culture corresponding to each individual basket
# in the same order as your baskets.csv file
locationsbasket <- read.csv("baskets_locations.csv")
locationsbasket_red <- subset(locationsbasket, locationsbasket$Culture %in% rownames(Fst_mat))

rownames(locationsbasket) <- rownames(baskets)
rownames(locationsbasket_red) <- rownames(baskets_red)

biomesbasket <- read.csv("baskets_biomes.csv")
biomesbasket_red <- subset(biomesbasket, biomesbasket$Culture %in% rownames(Fst_mat))
biomesbasket$Trait <- NULL
biomesbasket$Culture <- NULL
biomesbasket_red$Trait <- NULL
biomesbasket_red$Culture <- NULL

rownames(biomesbasket) <- rownames(baskets)
rownames(biomesbasket_red) <- rownames(baskets_red)

# Convert factors into numeric codes (not really needed if stringsAsFactors = F)
indx <- sapply(baskets, is.factor)
indx_red <- sapply(baskets_red, is.factor)

baskets[indx] <- lapply(baskets[indx],  function(x) as.numeric(as.factor(x)))
baskets_red[indx_red] <- lapply(baskets_red[indx_red],  function(x) as.numeric(as.factor(x)))

#compute Gower dissimilarity distance
library(vegan)
gowbasket <- vegdist(baskets,"gower",na.rm=TRUE) 
gowbasketRed <- vegdist(baskets_red,"gower",na.rm=TRUE) 

gowBiomebaskets <- as.matrix(vegdist(biomesbasket,"gower",na.rm=TRUE))
gowBiomebasketsRed <-  as.matrix(vegdist(biomesbasket_red,"gower",na.rm=TRUE))


#compute binary Matrices of cultural affiliation (0:same culture; 1:different culture)
binaryCultDistbaskets<-as.matrix(as.dist(BinaryCulture(locationsbasket$Culture)))
binaryCultDistbasketsRed<-as.matrix(as.dist(BinaryCulture(locationsbasket_red$Culture)))

#compute Fst matrix
GenDistbaskets<-as.matrix(as.dist(GetFst(locationsbasket_red$Culture)))

#compute spatial distance between cultures and therefore baskets
#spacebaskets<-distMat(data.frame(x=locationsbasket$Longitude,y=locationsbasket$Latitude))

#compute spatial distances
distTablebaskets <- data.frame(x=locationsbasket$Longitude,y=locationsbasket$Latitude)
namesbasket <- as.vector(row.names(locationsbasket))
rownames(distTablebaskets) <- namesbasket

distTablebasketsRed <- data.frame(x=locationsbasket_red$Longitude,y=locationsbasket_red$Latitude)
namesbasketRed <- as.vector(row.names(locationsbasket_red))
rownames(distTablebasketsRed) <- namesbasketRed

library(geosphere)
spacebaskets <- distm(distTablebaskets, fun=distGeo)
rownames(spacebaskets) <- namesbasket
colnames(spacebaskets) <- namesbasket
# here distance is in metres, not kilometres
#spacebaskets <- as.dist(spacebaskets)

spacebasketsRed <- distm(distTablebasketsRed, fun=distGeo)
rownames(spacebasketsRed) <- namesbasketRed
colnames(spacebasketsRed) <- namesbasketRed
# here distance is in metres, not kilometres
#spacebasketsRed <- as.dist(spacebasketsRed)

gowbasket <- as.matrix(gowbasket)
na_rows <- which(is.na(gowbasket), arr.ind=TRUE)
#na_rows <- unique(na_rows[,1])
#gowbasket <- gowbasket[-na_rows, -na_rows]
#spacebaskets <- spacebaskets[-na_rows, -na_rows]
#gowBiomebaskets <- gowBiomebaskets[-na_rows, -na_rows]
#binaryCultDistbaskets <- binaryCultDistbaskets[-na_rows, -na_rows]

gowbasketRed <- as.matrix(gowbasketRed)
na_rows <- which(is.na(gowbasketRed), arr.ind=TRUE)
#na_rows <- unique(na_rows[,1])
#gowbasketRed <- gowbasketRed[-na_rows, -na_rows]
#spacebasketsRed  <- spacebasketsRed [-na_rows, -na_rows]
#gowBiomebasketsRed  <- gowBiomebasketsRed [-na_rows, -na_rows]
#binaryCultDistbasketsRed  <- binaryCultDistbasketsRed [-na_rows, -na_rows]
#GenDistbaskets <- GenDistbaskets[-na_rows, -na_rows]


#partial mantel test between Distance and Geography and between Culture and Geography
library(ecodist)
detach("package:vegan", unload = TRUE)

space_mantel <- mantel(lower(gowbasket)~lower(spacebaskets),nperm=1000)
culture_mantel <-mantel(lower(gowbasket)~lower(binaryCultDistbaskets),nperm=1000)
genes_mantel <- mantel(lower(gowbasketRed)~lower(GenDistbaskets),nperm=1000)
biome_mantel <- mantel(lower(gowBiomebasketsRed)~lower(GenDistbaskets),nperm=1000)

#partial mantel test testing for effect of CULTURAL AFFILIATION controlling for rest
cult_space_mantel <- mantel(lower(gowbasket)~lower(binaryCultDistbaskets)+lower(spacebaskets),nperm=1000)
cult_genes_mantel <- mantel(lower(gowbasketRed)~lower(binaryCultDistbasketsRed)+lower(GenDistbaskets),nperm=1000)
cult_biome_mantel <- mantel(lower(gowbasket)~lower(binaryCultDistbaskets)+lower(gowBiomebaskets),nperm=1000)

#partial mantel test testing for effect of SHARED GENES
genes_space_mantel <- mantel(lower(gowbasketRed)~lower(GenDistbaskets)+lower(spacebasketsRed),nperm=1000)
genes_cult_mantel <- mantel(lower(gowbasketRed)~lower(GenDistbaskets)+lower(binaryCultDistbasketsRed),nperm=1000)
genes_biome_mantel <- mantel(lower(gowbasketRed)~lower(GenDistbaskets)+lower(gowBiomebasketsRed),nperm=1000)

#partial mantel test testing for effect of GEOGRAPHY
space_culture_mantel <- mantel(lower(gowbasket)~lower(spacebaskets)+lower(binaryCultDistbaskets),nperm=1000)
space_genes_mantel <- mantel(lower(gowbasketRed)~lower(spacebasketsRed)+lower(GenDistbaskets),nperm=1000)
space_biome_mantel <- mantel(lower(gowbasket)~lower(spacebaskets)+lower(gowBiomebaskets),nperm=1000)

#partial mantel test testing for effect of BIOME
biome_culture_mantel <- mantel(lower(gowbasket)~lower(gowBiomebaskets)+lower(binaryCultDistbaskets),nperm=1000)
biome_genes_mantel <- mantel(lower(gowbasketRed)~lower(gowBiomebasketsRed)+lower(GenDistbaskets),nperm=1000)
biome_space_mantel <- mantel(lower(gowbasket)~lower(gowBiomebaskets)+lower(spacebaskets),nperm=1000)


write.csv(rbind(space_mantel, culture_mantel, genes_mantel,biome_mantel,
                cult_space_mantel,cult_genes_mantel,cult_biome_mantel,
                genes_space_mantel, genes_cult_mantel,genes_biome_mantel,
                space_culture_mantel,space_genes_mantel,space_biome_mantel,
                biome_culture_mantel, biome_genes_mantel, biome_space_mantel), "results_masked/baskets_mantel.csv")



#########
# AMOVA #
#########
library(pegas)

#AMOVA
sqgowbasket<-sqrt(gowbasket)
sqgowbasket_mat <- as.data.frame(as.matrix(sqgowbasket))
# # get NA
#na_rows <- which(is.na(sqgowbasket_mat), arr.ind=TRUE)
# 
# # remove rows with NAs
# na_rows <- unique(na_rows[,1])
# sqgowbasket_mat <- sqgowbasket_mat[-na_rows,-na_rows]
# sqgowbasket <- as.dist(sqgowbasket_mat)

#groupsbasket <- as.factor(groupsbasket[-na_rows])
groupsbasket <- as.factor(locationsbasket$Culture)
AMOVAbasket=amova(sqgowbasket~groupsbasket)

# Phi-Stat:
AMOVAbasket$varcomp$sigma2[1]/sum(AMOVAbasket$varcomp$sigma2)
#P-Value:
AMOVAbasket$varcomp$P.value[1]


##################
# pairwise AMOVA #
##################

# Note: This is the important thing to compare against genetic data

# I think you cannot do amova for cultures with only 1 trait!

#Extract Culture Data
CultureNames<-unique(as.character(groupsbasket))
CultureNames <- sort(CultureNames)
NumberOfCultures<-length(CultureNames)
sqgowbasket<-as.matrix(sqgowbasket)
CultureList<-groupsbasket

#PhiSt and Pval Output Matrices:
PhiMatrix<-matrix(NA,nrow=NumberOfCultures,ncol=NumberOfCultures)
PvalMatrix<-matrix(NA,nrow=NumberOfCultures,ncol=NumberOfCultures)

#Spatial Coordinates Output Data.Frame:
XY<-data.frame(x=numeric(length=NumberOfCultures),y=numeric(length=NumberOfCultures))

for (i in 1:NumberOfCultures) {
  print(paste(i," out of ",NumberOfCultures))
  for (j in 1:NumberOfCultures)
  {
    if (i!=j&i>j)
    { 
      indexI=which(as.character(CultureList)==as.character(CultureNames[i]))
      indexJ=which(as.character(CultureList)==as.character(CultureNames[j]))
      #subsetDistance Matrix
      subDistMat<-sqgowbasket[c(indexI,indexJ),c(indexI,indexJ)]
      subDistMat<-as.dist(subDistMat)
      classes<-as.factor(c(rep("a",length(indexI)),rep("b",length(indexJ))))
      #compute Phi:
      tmp<-pegas::amova(subDistMat~classes)
      N0=tmp$varcoef
      AP=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)
      WP=(tmp$tab$MSD[2])
      PhiMatrix[i,j]=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)/((tmp$tab$MSD[2])+((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef))
      if(PhiMatrix[i,j]<0){PhiMatrix[i,j]=0}
      PvalMatrix[i,j]=tmp$varcomp$P.value[1]
    }
  }
  #XY$x[i]<-mean(locationsbasket$Longitude[which(locationsbasket$Culture==as.character(CultureNames[i]))],na.rm=TRUE)
  #XY$y[i]<-mean(locationsbasket$Latitude[which(locationsbasket$Culture==as.character(CultureNames[i]))],na.rm=TRUE)
}

dim <- ncol(PhiMatrix)

#Change NAs to blanks for text
PhiMatrix_text <- as.data.frame(PhiMatrix)
PhiMatrix_text <- sapply(PhiMatrix_text, as.character)
PhiMatrix_text[is.na(PhiMatrix_text)] <- " "
PhiMatrix_text <- as.matrix(PhiMatrix_text)

png(file="results_masked/baskets_phistats.png",
    width=7, height=6, units="in", res=500)
image(-7:(dim-8), 1:dim, PhiMatrix, axes = FALSE, xlab="", ylab="", col=hcl.colors(12, "Geyser", rev=TRUE), main="baskets", xlim=c(-7,5.5), ylim=c(0,11))
axis(1, (-7:(dim-8)), CultureNames, cex.axis = 0.6, las=3, pos=0.5)
axis(4, 1:dim, CultureNames, cex.axis = 0.6, las=1, pos=3.5)
text(expand.grid((-7:(dim-8)), 1:dim), sprintf("%0.4s", PhiMatrix_text), cex=0.6)
dev.off()

# Again Partial Mantel Tests using Phi Statistic

CultureNames_rn <- CultureNames
CultureNames_rn[CultureNames_rn=="Batwa (East)"] <- "BatwaE"
CultureNames_rn[CultureNames_rn=="Batwa (West)"] <- "BatwaW"
rownames(PhiMatrix) <- CultureNames
colnames(PhiMatrix) <- CultureNames
PhiMatrix_rn <- PhiMatrix
rownames(PhiMatrix_rn) <- CultureNames_rn
colnames(PhiMatrix_rn) <- CultureNames_rn

writeDist(PhiMatrix_rn, file="baskets_phi.nex", format="nexus")

# Mantel tests between geography, biome, genetics and baskets
PhiMatrix <- as.matrix(PhiMatrix)
#genDist <- as.matrix(genDist)
#gowBiomeRed <- as.matrix(gowBiomeRed)
#spaceRed <- as.matrix(spaceRed)

# Remove groups with <2 baskets
genDistbasket <- genDist
genDistbasketIBD <- genDistIBD
genDistbasketBantu <- genDistBantu
BiomeDistbasket <- gowBiomeRed
spacebasket <- spaceRed
PhiMatrixbaskets<- PhiMatrix[-7,-7]

biome_mantel <- mantel(lower(PhiMatrixbaskets)~lower(BiomeDistbasket), mrank=T)
space_mantel <- mantel(lower(PhiMatrixbaskets)~lower(spacebasket), mrank=T)
genes_mantel <- mantel(lower(PhiMatrixbaskets)~lower(genDistbasket), mrank=T)
genes_ibd_mantel <- mantel(lower(PhiMatrixbaskets)~lower(genDistbasketIBD), mrank=T)
genes_Bantu_mantel <- mantel(lower(PhiMatrixbaskets[-c(8,10), -c(8,10)])~lower(genDistbasketBantu), mrank=T)

#partial mantel test testing for effect of SHARED GENES
genes_space_mantel <- mantel(lower(PhiMatrixbaskets)~lower(genDistbasket)+lower(spacebasket), mrank=T)
genes_biome_mantel <- mantel(lower(PhiMatrixbaskets)~lower(genDistbasket)+lower(BiomeDistbasket), mrank=T)

#partial mantel test testing for effect of GEOGRAPHY
space_genes_mantel <- mantel(lower(PhiMatrixbaskets)~lower(spacebasket)+lower(genDistbasket), mrank=T)
space_biome_mantel <- mantel(lower(PhiMatrixbaskets)~lower(spacebasket)+lower(BiomeDistbasket), mrank=T)

#partial mantel test testing for effect of BIOME
biome_genes_mantel <- mantel(lower(PhiMatrixbaskets)~lower(BiomeDistbasket)+lower(genDistbasket), mrank=T)
biome_space_mantel <- mantel(lower(PhiMatrixbaskets)~lower(BiomeDistbasket)+lower(spacebasket), mrank=T)


mantel_baskets <- rbind(space_mantel, genes_mantel,
                        biome_mantel,
                        genes_space_mantel, genes_biome_mantel,
                        space_genes_mantel, space_biome_mantel,
                        biome_space_mantel, biome_genes_mantel)

write.csv(mantel_baskets, "results_masked/mantel_baskets.csv")


mantel_baskets <-  rbind(space_mantel, genes_mantel,
                       biome_mantel, genes_Bantu_mantel, genes_ibd_mantel)
pval_adj_baskets <- p.adjust(mantel_baskets[,2], method = "BH")
mantel_baskets<- cbind(mantel_baskets, pval_adj_baskets)

write.csv(mantel_baskets, "final_paper_files/results_final_paper/results_masked/mantel_baskets.csv")

matrix_baskets <- rbind(genes_space_mantel, genes_biome_mantel,
                      space_genes_mantel, space_biome_mantel,
                      biome_space_mantel, biome_genes_mantel)
pval_adj_baskets <- p.adjust(matrix_baskets[,2], method = "BH")
matrix_baskets<- cbind(matrix_baskets, pval_adj_baskets)

write.csv(matrix_baskets, "results_masked/matrix_baskets.csv")


# Perform PCoA
gowbasket <- as.data.frame(as.matrix(gowbasket))
# get NA
na_rows <- which(is.na(gowbasket), arr.ind=TRUE)

# remove rows with NAs
#na_rows <- unique(na_rows[,1])
#gowbasket <- gowbasket[-na_rows,-na_rows]
#gowbasket <- as.dist(gowbasket)
baskets_pcoa <- pcoa(gowbasket, correction = "cailliez")

#We rescale the PCoA components in relation to the explained variance.  

for(i in 1:ncol(baskets_pcoa$vectors)) { 
  baskets_pcoa$vectors[,i] <- scale(baskets_pcoa$vectors[,i])*
    baskets_pcoa$values$Rel_corr_eig[i]}

# Extract eigenvalues from PC/PCo results_masked 
baskets_ev <- baskets_pcoa$values$Corr_eig

# Plot
library(ggplot2)
my_font <- "sans"

png(file="results_masked/pcoabaskets.png",
    width=6, height=4.5, units="in", res=500)
plot_pc_scree(eigenval=baskets_ev, title="Genetics", type="PCo")
dev.off()

##We extract the principal components/coordinates for all factors: 
#genetics, grammar, music and phonology. 
#We normalize the PCs/PCos to a range from 0 to 1 and then plot a heatmap 

# Extract the PCs and the PCos from the PCA and PCoA results_masked
CultureList<-groupsbasket
baskets_pco <- as.data.frame(baskets_pcoa$vectors)
baskets_pco$Culture <- CultureList

## PCoA with groups

#install.packages('ggfortify')

# remove Mbendjele from palette
my_pal_long <- my_pal_long

p <- ggplot(data = baskets_pco, aes(x = Axis.1, y = Axis.2, color=Culture)) +
  geom_point(size = 3.5) +
  xlab("PCo1 (17.96%)") +
  ylab("PCo2 (8.04%)") + 
  theme_classic() + 
  scale_color_manual(values = my_pal_long,
                     labels = c("Aka (6)",
                                "Babongo (2)",
                                "Baka (28)",
                                "Bakola (5)",
                                "Bakoya (19)",
                                "East Batwa (6)",
                                "West Batwa (23)",
                                "Bedzan (6)",
                                "Efe (7)",
                                "Mbendjele (13)",
                                "Sua (10)")) +
  geom_hline(yintercept = 0, lty = 2) +
  ggtitle("Baskets") +
  geom_vline(xintercept = 0, lty = 2) 

p <- p + theme(legend.title = element_blank(),
               legend.text = element_text(size = 12))

ggsave(p, filename="results_masked/PCOA_baskets.png", height=150, width=180, units="mm")

#We compute Euclidean distance matrices from the dimensionality-reduced data and create NeighborNets. 
# SplitsTree doesnt accept brackets in names
baskets_pco_rn <- baskets_pco
baskets_pco_rn$Culture <- as.character(baskets_pco_rn$Culture)
baskets_pco_rn$Culture[baskets_pco_rn$Culture=="Batwa (East)"] <- "BatwaE"
baskets_pco_rn[baskets_pco_rn=="Batwa (West)"] <- "BatwaW"
baskets_pco_rn$Culture <- NULL

# Write nexus files for splitstree (Not used so far)
require(phangorn)

dist_list <- lapply(list(baskets=baskets_pco_rn),
                    function(f) {dist(f)})
writeDist(dist_list$baskets, file="baskets.nex", format="nexus")

###################################################################
## Now we test for potential mechanisms of cultural transmission ##
###################################################################

# Divide cultures into Western and Eastern Pygmies 
baskets <- read.csv("baskets.csv", stringsAsFactors = T,header=T)

# Question - should we separate between different regions???

# Remove traits that are only present in one of the two regions (East or West)
baskets$Region<-sapply(as.character(baskets$Culture),function(x)
  ifelse(x %in% c("Aka", "Babongo", "Baka", "Bakola", "Bakoya",
                  "Bedzan", "Mbendjele", "Batwa (West)"),"West",
         ifelse(x %in% c("Batwa (East)","Efe","Sua"),"East")))
table(baskets$Region)

basket_redux<-baskets[,c(3:ncol(baskets))]
for(i in 1:(ncol(basket_redux)-1)){
  tab<-table(basket_redux[,c("Region",colnames(basket_redux)[i])])  
  check<-apply(tab,2,function(x) ifelse(sum(x!=0)==1,"fix","ok"))  
  basket_redux[,colnames(basket_redux)[i]]<-sapply(basket_redux[,colnames(basket_redux)[i]],
                                                   function(x) ifelse(x %in% names(check)[check=="fix"],NA,as.character(x)))
  basket_redux[,colnames(basket_redux)[i]]<-as.factor(basket_redux[,colnames(basket_redux)[i]])}
basket_redux<-basket_redux[,sapply(basket_redux, function(x)!all(is.na(x)))]
basket_redux$Region<-as.factor(basket_redux$Region)

library(party)
library(plyr)
weights_baskets<-sapply(basket_redux$Region,
                        function(x) nrow(basket_redux)/nrow(basket_redux[basket_redux$Region==x,]))
rf_redux<-cforest(Region~.,
                  controls=
                    cforest_control(ntree = 500,
                                    mtry=5,
                                    mincriterion=qnorm(0.9),
                                    fraction = 0.632,
                                    testtype ="Teststatistic",
                                    teststat="max",
                                    replace=TRUE,
                                    trace=F,
                                    savesplitstats=F,
                                    minsplit=20,
                                    minbucket=8),
                  data=basket_redux,
                  weights = weights_baskets)

# Produce confusion matrix
cmat_redux<-table(basket_redux$Region, predict(rf_redux))
print(cmat_redux)

stats_region<-ldply(c(1:nrow(cmat_redux)),
                    function(x) data.frame(id=rownames(cmat_redux)[x],
                                           N=sum(cmat_redux[,x]),
                                           d=cmat_redux[x,x],
                                           A=sum(cmat_redux[x,]))) %>% transform(p=d/N,
                                                                                 r=d/A)
stats_region
write.csv(stats_region, "results_masked/stats_region_baskets.csv")

basket_redux$predReg<-predict(rf_redux)
missclass<-cbind(baskets[rownames(basket_redux[basket_redux$Region!=basket_redux$predReg,]),],
                 data.frame(predReg=basket_redux$predReg[basket_redux$Region!=basket_redux$predReg]))
missclass
write.csv(missclass, "results_masked/missclass_baskets.csv")
##########
## FAMD ##
##########

baskets$Region <- NULL
imp.baskets <- imputeFAMD(baskets[,2:12], ncp = 10)
res.famd <- FAMD(imp.baskets$completeObs,   ## Set the target variable "Churn" as a supplementary variable, so it is not included in the analysis for now
                 graph = T, 
                 sup.var=1, # Culture is the supplementary variable
                 ncp=10)

## Inspect principal components
get_eigenvalue(res.famd)
fviz_screeplot(res.famd)
# Plot of variables
fviz_famd_var(res.famd, repel = TRUE, ggtheme=theme_classic())
ggsave("results_masked/baskets_FAMD_varcont_dim.png", width=4.5, height=4, units="in", dpi=500)
# Contribution to the first dimension

fviz_contrib(res.famd, "var", axes = 1, ggtheme=theme_classic()) 
ggsave("results_masked/baskets_FAMD_varDIM1.png", width=7, height=4.5, units="in", dpi=500)
# Contribution to the second dimension
fviz_contrib(res.famd, "var", axes = 2, ggtheme=theme_classic())
ggsave("results_masked/baskets_FAMD_varDIM2.png", width=7, height=4.5, units="in", dpi=500)

fviz_famd_var(res.famd, "quanti.var", col.var = "contrib", 
              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
              repel = TRUE)

famd_plot_baskets <- fviz_mfa_ind(res.famd, 
                          habillage = 1, # color by groups 
                          palette = my_pal_long,
                          addEllipses = F, 
                          repel = TRUE, # Avoid text overlapping,
                          label="none",
                          ggtheme = theme_classic(),
) 

famd_plot_baskets <- famd_plot_baskets + ggtitle("Baskets - FAMD")
famd_plot_baskets
ggsave("results_masked/baskets_FAMD.tiff", width=7, height=6, units="in", dpi=500)

#############
## RATTLES ##
#############

rattles <- read.csv("rattles.csv", stringsAsFactors = T,header=T)
ids <- seq(1:nrow(rattles))
for (i in ids) {
  ids[i] <- paste("T-", i, sep="")
}

rattles$ID <- ids

# For genetics, create another rattle dataset with only rattles from groups with genes
rattles_red <- subset(rattles, rattles$Culture %in% rownames(Fst_mat))

rattles$Trait <- NULL
ids <- rattles$ID
ids_red <- rattles_red$ID
groupsrattle<-as.factor(rattles$Culture)
groupsrattle_red<-as.factor(rattles$Culture)
rownames(rattles) <- ids
rownames(rattles_red) <- ids_red
rattles$ID <- NULL
rattles$Culture <- NULL
rattles_red$Trait <- NULL
rattles_red$ID <- NULL
rattles_red$Culture <- NULL

# Convert factors into numeric codes (not really needed if stringsAsFactors = F)
indx <- sapply(rattles, is.factor)
indx_red <- sapply(rattles_red, is.factor)

rattles[indx] <- lapply(rattles[indx],  function(x) as.numeric(as.factor(x)))
rattles_red[indx_red] <- lapply(rattles_red[indx_red],  function(x) as.numeric(as.factor(x)))

#compute Gower dissimilarity distance
library(vegan)
gowrattle <- vegdist(rattles,"gower",na.rm=TRUE) 
gowrattleRed <- vegdist(rattles_red,"gower",na.rm=TRUE) 


#########
# AMOVA #
#########
library(pegas)

#AMOVA
sqgowrattle<-sqrt(gowrattle)
sqgowrattle_mat <- as.data.frame(as.matrix(sqgowrattle))
# # get NA
#na_rows <- which(is.na(sqgowrattle_mat), arr.ind=TRUE)
# 
# # remove rows with NAs
# na_rows <- unique(na_rows[,1])
# sqgowrattle_mat <- sqgowrattle_mat[-na_rows,-na_rows]
# sqgowrattle <- as.dist(sqgowrattle_mat)

#groupsrattle <- as.factor(groupsrattle[-na_rows])
AMOVArattle=amova(sqgowrattle~groupsrattle)

# Phi-Stat:
AMOVArattle$varcomp$sigma2[1]/sum(AMOVArattle$varcomp$sigma2)
#P-Value:
AMOVArattle$varcomp$P.value[1]


##################
# pairwise AMOVA #
##################

# Note: This is the important thing to compare against genetic data

# I think you cannot do amova for cultures with only 1 trait!

#Extract Culture Data
CultureNames<-unique(as.character(groupsrattle))
CultureNames <- sort(CultureNames)
NumberOfCultures<-length(CultureNames)
sqgowrattle<-as.matrix(sqgowrattle)
CultureList<-groupsrattle

#PhiSt and Pval Output Matrices:
PhiMatrix<-matrix(NA,nrow=NumberOfCultures,ncol=NumberOfCultures)
PvalMatrix<-matrix(NA,nrow=NumberOfCultures,ncol=NumberOfCultures)

#Spatial Coordinates Output Data.Frame:
XY<-data.frame(x=numeric(length=NumberOfCultures),y=numeric(length=NumberOfCultures))

for (i in 1:NumberOfCultures) {
  print(paste(i," out of ",NumberOfCultures))
  for (j in 1:NumberOfCultures)
  {
    if (i!=j&i>j)
    { 
      indexI=which(as.character(CultureList)==as.character(CultureNames[i]))
      indexJ=which(as.character(CultureList)==as.character(CultureNames[j]))
      #subsetDistance Matrix
      subDistMat<-sqgowrattle[c(indexI,indexJ),c(indexI,indexJ)]
      subDistMat<-as.dist(subDistMat)
      classes<-as.factor(c(rep("a",length(indexI)),rep("b",length(indexJ))))
      #compute Phi:
      tmp<-pegas::amova(subDistMat~classes)
      N0=tmp$varcoef
      AP=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)
      WP=(tmp$tab$MSD[2])
      PhiMatrix[i,j]=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)/((tmp$tab$MSD[2])+((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef))
      if(PhiMatrix[i,j]<0){PhiMatrix[i,j]=0}
      PvalMatrix[i,j]=tmp$varcomp$P.value[1]
    }
  }
  #XY$x[i]<-mean(locationsrattle$Longitude[which(locationsrattle$Culture==as.character(CultureNames[i]))],na.rm=TRUE)
  #XY$y[i]<-mean(locationsrattle$Latitude[which(locationsrattle$Culture==as.character(CultureNames[i]))],na.rm=TRUE)
}

dim <- ncol(PhiMatrix)

#Change NAs to blanks for text
PhiMatrix_text <- as.data.frame(PhiMatrix)
PhiMatrix_text <- sapply(PhiMatrix_text, as.character)
PhiMatrix_text[is.na(PhiMatrix_text)] <- " "
PhiMatrix_text <- as.matrix(PhiMatrix_text)

png(file="results_masked/rattles_phistats.png",
    width=7, height=6, units="in", res=500)
image(-7:(dim-8), 1:dim, PhiMatrix, axes = FALSE, xlab="", ylab="", col=hcl.colors(12, "Geyser", rev=TRUE), main="rattles", xlim=c(-7,1), ylim=c(0,11))
axis(1, (-7:(dim-8)), CultureNames, cex.axis = 0.6, las=3, pos=0.5)
axis(4, 1:dim, CultureNames, cex.axis = 0.6, las=1, pos=0.5)
text(expand.grid((-7:(dim-8)), 1:dim), sprintf("%0.4s", PhiMatrix_text), cex=0.6)
dev.off()

# Again Partial Mantel Tests using Phi Statistic

CultureNames_rn <- CultureNames
CultureNames_rn[CultureNames_rn=="Batwa (East)"] <- "BatwaE"
CultureNames_rn[CultureNames_rn=="Batwa (West)"] <- "BatwaW"
rownames(PhiMatrix) <- CultureNames
colnames(PhiMatrix) <- CultureNames
PhiMatrix_rn <- PhiMatrix
rownames(PhiMatrix_rn) <- CultureNames_rn
colnames(PhiMatrix_rn) <- CultureNames_rn

writeDist(PhiMatrix_rn, file="rattles_phi.nex", format="nexus")

# Mantel tests between geography, biome, genetics and rattles
PhiMatrix <- as.matrix(PhiMatrix)
genDist <- as.matrix(genDist)
gowBiomeRed <- as.matrix(gowBiomeRed)
spaceRed <- as.matrix(spaceRed)

# Remove groups with >2 rattles
genDistrattle <-genDist[-c(5,8,7),-c(5,8,7)]
BiomeDistrattle <- gowBiomeRed[-c(5,8,7),-c(5,8,7)]
spacerattle <- spaceRed[-c(5,8,7),-c(5,8,7)]
PhiMatrixrattles<- PhiMatrix[-6,-6]

biome_mantel <- mantel(lower(PhiMatrixrattles)~lower(BiomeDistrattle), mrank=T)
space_mantel <- mantel(lower(PhiMatrixrattles)~lower(spacerattle), mrank=T)
genes_mantel <- mantel(lower(PhiMatrixrattles)~lower(genDistrattle), mrank=T)

#partial mantel test testing for effect of SHARED GENES
genes_space_mantel <- mantel(lower(PhiMatrixrattles)~lower(genDistrattle)+lower(spacerattle), mrank=T)
genes_biome_mantel <- mantel(lower(PhiMatrixrattles)~lower(genDistrattle)+lower(BiomeDistrattle), mrank=T)

#partial mantel test testing for effect of GEOGRAPHY
space_genes_mantel <- mantel(lower(PhiMatrixrattles)~lower(spacerattle)+lower(genDistrattle), mrank=T)
space_biome_mantel <- mantel(lower(PhiMatrixrattles)~lower(spacerattle)+lower(BiomeDistrattle), mrank=T)

#partial mantel test testing for effect of BIOME
biome_genes_mantel <- mantel(lower(PhiMatrixrattles)~lower(BiomeDistrattle)+lower(genDistrattle), mrank=T)
biome_space_mantel <- mantel(lower(PhiMatrixrattles)~lower(BiomeDistrattle)+lower(spacerattle), mrank=T)


mantel_rattles <- rbind(space_mantel, genes_mantel,
                        biome_mantel,
                        genes_space_mantel, genes_biome_mantel,
                        space_genes_mantel, space_biome_mantel,
                        biome_space_mantel, biome_genes_mantel)

write.csv(mantel_rattles, "results_masked/mantel_rattles.csv")

# Perform PCoA
gowrattle <- as.data.frame(as.matrix(gowrattle))
# get NA
#na_rows <- which(is.na(gowrattle), arr.ind=TRUE)

# remove rows with NAs
#na_rows <- unique(na_rows[,1])
#gowrattle <- gowrattle[-na_rows,-na_rows]
#gowrattle <- as.dist(gowrattle)
rattles_pcoa <- pcoa(gowrattle, correction = "cailliez")

#We rescale the PCoA components in relation to the explained variance.  

for(i in 1:ncol(rattles_pcoa$vectors)) { 
  rattles_pcoa$vectors[,i] <- scale(rattles_pcoa$vectors[,i])*
    rattles_pcoa$values$Rel_corr_eig[i]}

# Extract eigenvalues from PC/PCo results_masked 
rattles_ev <- rattles_pcoa$values$Corr_eig

# Plot
library(ggplot2)
my_font <- "sans"

png(file="results_masked/pcoarattles.png",
    width=6, height=4.5, units="in", res=500)
plot_pc_scree(eigenval=rattles_ev, title="Genetics", type="PCo")
dev.off()

##We extract the principal components/coordinates for all factors: 
#genetics, grammar, music and phonology. 
#We normalize the PCs/PCos to a range from 0 to 1 and then plot a heatmap 

# Extract the PCs and the PCos from the PCA and PCoA results_masked
CultureList<-groupsrattle
rattles_pco <- as.data.frame(rattles_pcoa$vectors)
rattles_pco$Culture <- CultureList

## PCoA with groups

#install.packages('ggfortify')

# remove Mbendjele from palette
my_pal_long_ar <- my_pal_long[-c(5,9,10)]

p <- ggplot(data = rattles_pco, aes(x = Axis.1, y = Axis.2, color=Culture)) +
  geom_point(size = 3.5) +
  xlab("PCo1 (7.98%)") +
  ylab("PCo2 (3.35%)") + 
  theme_classic() + 
  scale_color_manual(values = my_pal_long_ar,
                     labels = c("Aka (3)",
                                "Babongo (4)",
                                "Baka (5)",
                                "Bakola (3)",
                                "East Batwa (3)",
                                "West Batwa (33)",
                                "Bedzan (3)",
                                "Sua (2)")) +
  geom_hline(yintercept = 0, lty = 2) +
  ggtitle("Rattles") +
  geom_vline(xintercept = 0, lty = 2) 

p <- p + theme(legend.title = element_blank(),
               legend.text = element_text(size = 12))

ggsave(p, filename="results_masked/PCOA_rattles.png", height=150, width=180, units="mm")

#We compute Euclidean distance matrices from the dimensionality-reduced data and create NeighborNets. 
# SplitsTree doesnt accept brackets in names
rattles_pco_rn <- rattles_pco
rattles_pco_rn$Culture <- as.character(rattles_pco_rn$Culture)
rattles_pco_rn$Culture[rattles_pco_rn$Culture=="Batwa (East)"] <- "BatwaE"
rattles_pco_rn[rattles_pco_rn=="Batwa (West)"] <- "BatwaW"
rattles_pco_rn$Culture <- NULL

# Write nexus files for splitstree (Not used so far)
require(phangorn)

dist_list <- lapply(list(rattles=rattles_pco_rn),
                    function(f) {dist(f)})
writeDist(dist_list$rattles, file="rattles.nex", format="nexus")

##########
## BOWS ##
##########


bows <- read.csv("bows.csv", stringsAsFactors = T,header=T)
ids <- seq(1:nrow(bows))
for (i in ids) {
  ids[i] <- paste("T-", i, sep="")
}

bows$ID <- ids

# For genetics, create another bow dataset with only bows from groups with genes
bows_red <- subset(bows, bows$Culture %in% rownames(Fst_mat))

bows$Trait <- NULL
ids <- bows$ID
ids_red <- bows_red$ID
rownames(bows) <- ids
rownames(bows_red) <- ids_red
bows$ID <- NULL
bows$Culture <- NULL
bows_red$Trait <- NULL
bows_red$ID <- NULL
bows_red$Culture <- NULL

# this file needs to have the coordinates of the culture corresponding to each individual bow
# in the same order as your bows.csv file
locationsbow <- read.csv("bows_locations.csv")
locationsbow_red <- subset(locationsbow, locationsbow$Culture %in% rownames(Fst_mat))

rownames(locationsbow) <- rownames(bows)
rownames(locationsbow_red) <- rownames(bows_red)

biomesbow <- read.csv("bows_biomes.csv")
biomesbow_red <- subset(biomesbow, biomesbow$Culture %in% rownames(Fst_mat))
biomesbow$Trait <- NULL
biomesbow$Culture <- NULL
biomesbow_red$Trait <- NULL
biomesbow_red$Culture <- NULL

rownames(biomesbow) <- rownames(bows)
rownames(biomesbow_red) <- rownames(bows_red)

# Convert factors into numeric codes (not really needed if stringsAsFactors = F)
indx <- sapply(bows, is.factor)
indx_red <- sapply(bows_red, is.factor)

bows[indx] <- lapply(bows[indx],  function(x) as.numeric(as.factor(x)))
bows_red[indx_red] <- lapply(bows_red[indx_red],  function(x) as.numeric(as.factor(x)))

#compute Gower dissimilarity distance
library(vegan)
gowbow <- vegdist(bows,"gower",na.rm=TRUE) 
gowbowRed <- vegdist(bows_red,"gower",na.rm=TRUE) 

gowBiomebows <- as.matrix(vegdist(biomesbow,"gower",na.rm=TRUE))
gowBiomebowsRed <-  as.matrix(vegdist(biomesbow_red,"gower",na.rm=TRUE))


#compute binary Matrices of cultural affiliation (0:same culture; 1:different culture)
binaryCultDistbows<-as.matrix(as.dist(BinaryCulture(locationsbow$Culture)))
binaryCultDistbowsRed<-as.matrix(as.dist(BinaryCulture(locationsbow_red$Culture)))

#compute Fst matrix
GenDistbows<-as.matrix(as.dist(GetFst(locationsbow_red$Culture)))

#compute spatial distance between cultures and therefore bows
#spacebows<-distMat(data.frame(x=locationsbow$Longitude,y=locationsbow$Latitude))

#compute spatial distances
distTablebows <- data.frame(x=locationsbow$Longitude,y=locationsbow$Latitude)
namesbow <- as.vector(row.names(locationsbow))
rownames(distTablebows) <- namesbow

distTablebowsRed <- data.frame(x=locationsbow_red$Longitude,y=locationsbow_red$Latitude)
namesbowRed <- as.vector(row.names(locationsbow_red))
rownames(distTablebowsRed) <- namesbowRed

library(geosphere)
spacebows <- distm(distTablebows, fun=distGeo)
rownames(spacebows) <- namesbow
colnames(spacebows) <- namesbow
# here distance is in metres, not kilometres
#spacebows <- as.dist(spacebows)

spacebowsRed <- distm(distTablebowsRed, fun=distGeo)
rownames(spacebowsRed) <- namesbowRed
colnames(spacebowsRed) <- namesbowRed
# here distance is in metres, not kilometres
#spacebowsRed <- as.dist(spacebowsRed)


gowbowRed <- as.matrix(gowbowRed)
na_rows <- which(is.na(gowbowRed), arr.ind=TRUE)
#na_rows <- unique(na_rows[,1])
#gowbowRed <- gowbowRed[-na_rows, -na_rows]
#spacebowsRed  <- spacebowsRed [-na_rows, -na_rows]
#gowBiomebowsRed  <- gowBiomebowsRed [-na_rows, -na_rows]
#binaryCultDistbowsRed  <- binaryCultDistbowsRed [-na_rows, -na_rows]
#GenDistbows <- GenDistbows[-na_rows, -na_rows]

gowbow <- as.matrix(gowbow)
na_rows <- which(is.na(gowbow), arr.ind=TRUE)
na_rows <- unique(na_rows[,1])
gowbow <- gowbow[-na_rows, -na_rows]
spacebows <- spacebows[-na_rows, -na_rows]
gowBiomebows <- gowBiomebows[-na_rows, -na_rows]
binaryCultDistbows <- binaryCultDistbows[-na_rows, -na_rows]

#partial mantel test between Distance and Geography and between Culture and Geography
library(ecodist)
detach("package:vegan", unload = TRUE)

space_mantel <- mantel(lower(gowbow)~lower(spacebows),nperm=1000)
culture_mantel <-mantel(lower(gowbow)~lower(binaryCultDistbows),nperm=1000)
genes_mantel <- mantel(lower(gowbowRed)~lower(GenDistbows),nperm=1000)
biome_mantel <- mantel(lower(gowBiomebowsRed)~lower(GenDistbows),nperm=1000)

#partial mantel test testing for effect of CULTURAL AFFILIATION controlling for rest
cult_space_mantel <- mantel(lower(gowbow)~lower(binaryCultDistbows)+lower(spacebows),nperm=1000)
cult_genes_mantel <- mantel(lower(gowbowRed)~lower(binaryCultDistbowsRed)+lower(GenDistbows),nperm=1000)
cult_biome_mantel <- mantel(lower(gowbow)~lower(binaryCultDistbows)+lower(gowBiomebows),nperm=1000)

#partial mantel test testing for effect of SHARED GENES
genes_space_mantel <- mantel(lower(gowbowRed)~lower(GenDistbows)+lower(spacebowsRed),nperm=1000)
genes_cult_mantel <- mantel(lower(gowbowRed)~lower(GenDistbows)+lower(binaryCultDistbowsRed),nperm=1000)
genes_biome_mantel <- mantel(lower(gowbowRed)~lower(GenDistbows)+lower(gowBiomebowsRed),nperm=1000)

#partial mantel test testing for effect of GEOGRAPHY
space_culture_mantel <- mantel(lower(gowbow)~lower(spacebows)+lower(binaryCultDistbows),nperm=1000)
space_genes_mantel <- mantel(lower(gowbowRed)~lower(spacebowsRed)+lower(GenDistbows),nperm=1000)
space_biome_mantel <- mantel(lower(gowbow)~lower(spacebows)+lower(gowBiomebows),nperm=1000)

#partial mantel test testing for effect of BIOME
biome_culture_mantel <- mantel(lower(gowbow)~lower(gowBiomebows)+lower(binaryCultDistbows),nperm=1000)
biome_genes_mantel <- mantel(lower(gowbowRed)~lower(gowBiomebowsRed)+lower(GenDistbows),nperm=1000)
biome_space_mantel <- mantel(lower(gowbow)~lower(gowBiomebows)+lower(spacebows),nperm=1000)


write.csv(rbind(space_mantel, culture_mantel, genes_mantel,biome_mantel,
                cult_space_mantel,cult_genes_mantel,cult_biome_mantel,
                genes_space_mantel, genes_cult_mantel,genes_biome_mantel,
                space_culture_mantel,space_genes_mantel,space_biome_mantel,
                biome_culture_mantel, biome_genes_mantel, biome_space_mantel), "results_masked/bows_mantel.csv")



#########
# AMOVA #
#########
library(pegas)

#AMOVA
sqgowbow<-sqrt(gowbow)
sqgowbow_mat <- as.data.frame(as.matrix(sqgowbow))
# # get NA
#na_rows <- which(is.na(sqgowbow_mat), arr.ind=TRUE)
# 
# # remove rows with NAs
#na_rows <- unique(na_rows[,1])
sqgowbow_mat <- sqgowbow_mat[-na_rows,-na_rows]
sqgowbow <- as.dist(sqgowbow_mat)

groupsbow <- as.factor(groupsbow[-na_rows])
AMOVAbow=amova(sqgowbow~groupsbow)

# Phi-Stat:
AMOVAbow$varcomp$sigma2[1]/sum(AMOVAbow$varcomp$sigma2)
#P-Value:
AMOVAbow$varcomp$P.value[1]


##################
# pairwise AMOVA #
##################

# Note: This is the important thing to compare against genetic data

# I think you cannot do amova for cultures with only 1 trait!

#Extract Culture Data
CultureNames<-unique(as.character(groupsbow))
CultureNames <- sort(CultureNames)
NumberOfCultures<-length(CultureNames)
sqgowbow<-as.matrix(sqgowbow)
CultureList<-groupsbow

#PhiSt and Pval Output Matrices:
PhiMatrix<-matrix(NA,nrow=NumberOfCultures,ncol=NumberOfCultures)
PvalMatrix<-matrix(NA,nrow=NumberOfCultures,ncol=NumberOfCultures)

#Spatial Coordinates Output Data.Frame:
XY<-data.frame(x=numeric(length=NumberOfCultures),y=numeric(length=NumberOfCultures))

for (i in 1:NumberOfCultures) {
  print(paste(i," out of ",NumberOfCultures))
  for (j in 1:NumberOfCultures)
  {
    if (i!=j&i>j)
    { 
      indexI=which(as.character(CultureList)==as.character(CultureNames[i]))
      indexJ=which(as.character(CultureList)==as.character(CultureNames[j]))
      #subsetDistance Matrix
      subDistMat<-sqgowbow[c(indexI,indexJ),c(indexI,indexJ)]
      subDistMat<-as.dist(subDistMat)
      classes<-as.factor(c(rep("a",length(indexI)),rep("b",length(indexJ))))
      #compute Phi:
      tmp<-pegas::amova(subDistMat~classes)
      N0=tmp$varcoef
      AP=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)
      WP=(tmp$tab$MSD[2])
      PhiMatrix[i,j]=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)/((tmp$tab$MSD[2])+((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef))
      if(PhiMatrix[i,j]<0){PhiMatrix[i,j]=0}
      PvalMatrix[i,j]=tmp$varcomp$P.value[1]
    }
  }
  #XY$x[i]<-mean(locationsbow$Longitude[which(locationsbow$Culture==as.character(CultureNames[i]))],na.rm=TRUE)
  #XY$y[i]<-mean(locationsbow$Latitude[which(locationsbow$Culture==as.character(CultureNames[i]))],na.rm=TRUE)
}

dim <- ncol(PhiMatrix)

#Change NAs to blanks for text
PhiMatrix_text <- as.data.frame(PhiMatrix)
PhiMatrix_text <- sapply(PhiMatrix_text, as.character)
PhiMatrix_text[is.na(PhiMatrix_text)] <- " "
PhiMatrix_text <- as.matrix(PhiMatrix_text)

png(file="results_masked/bows_phistats.png",
    width=7, height=6, units="in", res=500)
image(-7:(dim-8), 1:dim, PhiMatrix, axes = FALSE, xlab="", ylab="", col=hcl.colors(12, "Geyser", rev=TRUE), main="Bows", xlim=c(-7,1), ylim=c(0,11))
axis(1, (-7:(dim-8)), CultureNames, cex.axis = 0.6, las=3, pos=0.5)
axis(4, 1:dim, CultureNames, cex.axis = 0.6, las=1, pos=0.5)
text(expand.grid((-7:(dim-8)), 1:dim), sprintf("%0.4s", PhiMatrix_text), cex=0.6)
dev.off()

# Again Partial Mantel Tests using Phi Statistic

CultureNames_rn <- CultureNames
CultureNames_rn[CultureNames_rn=="Batwa (East)"] <- "BatwaE"
CultureNames_rn[CultureNames_rn=="Batwa (West)"] <- "BatwaW"
rownames(PhiMatrix) <- CultureNames
colnames(PhiMatrix) <- CultureNames
PhiMatrix_rn <- PhiMatrix
rownames(PhiMatrix_rn) <- CultureNames_rn
colnames(PhiMatrix_rn) <- CultureNames_rn

writeDist(PhiMatrix_rn, file="bows_phi.nex", format="nexus")

# Mantel tests between geography, biome, genetics and bows
PhiMatrix <- as.matrix(PhiMatrix)
genDist <- as.matrix(genDist)
gowBiomeRed <- as.matrix(gowBiomeRed)
spaceRed <- as.matrix(spaceRed)

# Remove groups with >2 bows
genDistbow <- genDist[-c(4,5,9),-c(4,5,9)]
genDistbowIBD <- genDistIBD[-c(4,5,9),-c(4,5,9)]
genDistbowBantu <- genDistBantu[-c(4,5,8),-c(4,5,8)]
BiomeDistbow <- gowBiomeRed[-c(4,5,9),-c(4,5,9)]
spacebow <- spaceRed[-c(4,5,9),-c(4,5,9)]
PhiMatrixbows<- PhiMatrix[-5,-5]

detach("package:vegan", unload = TRUE)
biome_mantel <- mantel(lower(PhiMatrixbows)~lower(BiomeDistbow), mrank=T)
space_mantel <- mantel(lower(PhiMatrixbows)~lower(spacebow), mrank=T)
genes_mantel <- mantel(lower(PhiMatrixbows)~lower(genDistbow), mrank=T)
genes_ibd_mantel <- mantel(lower(PhiMatrixbows)~lower(genDistbowIBD), mrank=T)
genes_bantu_mantel <- mantel(lower(PhiMatrixbows[-c(6,7),-c(6,7)])~lower(genDistbowBantu), mrank=T)

#partial mantel test testing for effect of SHARED GENES
genes_space_mantel <- mantel(lower(PhiMatrixbows)~lower(genDistbow)+lower(spacebow), mrank=T)
genes_biome_mantel <- mantel(lower(PhiMatrixbows)~lower(genDistbow)+lower(BiomeDistbow), mrank=T)

#partial mantel test testing for effect of GEOGRAPHY
space_genes_mantel <- mantel(lower(PhiMatrixbows)~lower(spacebow)+lower(genDistbow), mrank=T)
space_biome_mantel <- mantel(lower(PhiMatrixbows)~lower(spacebow)+lower(BiomeDistbow), mrank=T)

#partial mantel test testing for effect of BIOME
biome_genes_mantel <- mantel(lower(PhiMatrixbows)~lower(BiomeDistbow)+lower(genDistbow), mrank=T)
biome_space_mantel <- mantel(lower(PhiMatrixbows)~lower(BiomeDistbow)+lower(spacebow), mrank=T)


mantel_bows <- rbind(space_mantel, genes_mantel,
                     biome_mantel,
                     genes_space_mantel, genes_biome_mantel,
                     space_genes_mantel, space_biome_mantel,
                     biome_space_mantel, biome_genes_mantel)

write.csv(mantel_bows, "results_masked/mantel_bows.csv")

mantel_bows <-  rbind(space_mantel, genes_mantel,
                       biome_mantel, genes_ibd_mantel, genes_bantu_mantel)
pval_adj_bows <- p.adjust(mantel_bows[,2], method = "BH")
mantel_bows<- cbind(mantel_bows, pval_adj_bows)

write.csv(mantel_bows, "final_paper_files/results_final_paper/results_masked/mantel_bows.csv")

matrix_bows <- rbind(genes_space_mantel, genes_biome_mantel,
                      space_genes_mantel, space_biome_mantel,
                      biome_space_mantel, biome_genes_mantel)
pval_adj_bows <- p.adjust(matrix_bows[,2], method = "BH")
matrix_bows<- cbind(matrix_bows, pval_adj_bows)

write.csv(matrix_bows, "results_masked/matrix_bows.csv")



# Perform PCoA
gowbow <- as.data.frame(as.matrix(gowbow))
# get NA
na_rows <- which(is.na(gowbow), arr.ind=TRUE)

# remove rows with NAs
na_rows <- unique(na_rows[,1])
gowbow <- gowbow[-na_rows,-na_rows]
gowbow <- as.dist(gowbow)
bows_pcoa <- pcoa(gowbow, correction = "cailliez")

#We rescale the PCoA components in relation to the explained variance.  

for(i in 1:ncol(bows_pcoa$vectors)) { 
  bows_pcoa$vectors[,i] <- scale(bows_pcoa$vectors[,i])*
    bows_pcoa$values$Rel_corr_eig[i]}

# Extract eigenvalues from PC/PCo results_masked 
bows_ev <- bows_pcoa$values$Corr_eig

# Plot
library(ggplot2)
my_font <- "sans"

png(file="results_masked/pcoabows.png",
    width=6, height=4.5, units="in", res=500)
plot_pc_scree(eigenval=bows_ev, title="Genetics", type="PCo")
dev.off()

##We extract the principal components/coordinates for all factors: 
#genetics, grammar, music and phonology. 
#We normalize the PCs/PCos to a range from 0 to 1 and then plot a heatmap 

# Extract the PCs and the PCos from the PCA and PCoA results_masked
CultureList<-groupsbow
bows_pco <- as.data.frame(bows_pcoa$vectors)
bows_pco$Culture <- CultureList

## PCoA with groups

#install.packages('ggfortify')

# remove Mbendjele from palette
my_pal_long_ar <- my_pal_long[-c(4,5,10)]

p <- ggplot(data = bows_pco, aes(x = Axis.1, y = Axis.2, color=Culture)) +
  geom_point(size = 3.5) +
  xlab("PCo1 (17.98%)") +
  ylab("PCo2 (7.99%)") + 
  theme_classic() + 
  scale_color_manual(values = my_pal_long_ar,
                     labels = c("Aka (3)",
                                "Babongo (1)",
                                "Baka (3)",
                                "East Batwa (11)",
                                "West Batwa (9)",
                                "Bedzan (3)",
                                "Efe (22)",
                                "Sua (7)")) +
  geom_hline(yintercept = 0, lty = 2) +
  ggtitle("Bows") +
  geom_vline(xintercept = 0, lty = 2) 

p <- p + theme(legend.title = element_blank(),
               legend.text = element_text(size = 12))

ggsave(p, filename="results_masked/PCOA_bows.png", height=150, width=180, units="mm")

#We compute Euclidean distance matrices from the dimensionality-reduced data and create NeighborNets. 
# SplitsTree doesnt accept brackets in names
bows_pco_rn <- bows_pco
bows_pco_rn$Culture <- as.character(bows_pco_rn$Culture)
bows_pco_rn$Culture[bows_pco_rn$Culture=="Batwa (East)"] <- "BatwaE"
bows_pco_rn[bows_pco_rn=="Batwa (West)"] <- "BatwaW"
bows_pco_rn$Culture <- NULL

# Write nexus files for splitstree (Not used so far)
require(phangorn)

dist_list <- lapply(list(bows=bows_pco_rn),
                    function(f) {dist(f)})
writeDist(dist_list$bows, file="bows.nex", format="nexus")

########
## RF ##
########

# Divide cultures into Western and Eastern Pygmies 

# Question - should we separate between different regions???
bows <- read.csv("bows.csv")

# Remove traits that are only present in one of the two regions (East or West)
bows$Region<-sapply(as.character(bows$Culture),function(x)
  ifelse(x %in% c("Aka", "Babongo", "Baka", "Bakola", "Bakoya",
                  "Bedzan", "Mbendjele", "Batwa (West)"),"West",
         ifelse(x %in% c("Batwa (East)","Efe","Sua"),"East")))
table(bows$Region)

bow_redux<-bows[,c(3:ncol(bows))]
for(i in 1:(ncol(bow_redux)-1)){
  tab<-table(bow_redux[,c("Region",colnames(bow_redux)[i])])  
  check<-apply(tab,2,function(x) ifelse(sum(x!=0)==1,"fix","ok"))  
  bow_redux[,colnames(bow_redux)[i]]<-sapply(bow_redux[,colnames(bow_redux)[i]],
                                                 function(x) ifelse(x %in% names(check)[check=="fix"],NA,as.character(x)))
  bow_redux[,colnames(bow_redux)[i]]<-as.factor(bow_redux[,colnames(bow_redux)[i]])}
bow_redux<-bow_redux[,sapply(bow_redux, function(x)!all(is.na(x)))]
bow_redux$Region<-as.factor(bow_redux$Region)

library(party)
library(plyr)
weights_bows<-sapply(bow_redux$Region,
                       function(x) nrow(bow_redux)/nrow(bow_redux[bow_redux$Region==x,]))
rf_redux<-cforest(Region~.,
                  controls=
                    cforest_control(ntree = 500,
                                    mtry=5,
                                    mincriterion=qnorm(0.9),
                                    fraction = 0.632,
                                    testtype ="Teststatistic",
                                    teststat="max",
                                    replace=TRUE,
                                    trace=F,
                                    savesplitstats=F,
                                    minsplit=20,
                                    minbucket=8),
                  data=bow_redux,
                  weights = weights_bows)

# Produce confusion matrix
cmat_redux<-table(bow_redux$Region, predict(rf_redux))
print(cmat_redux)

stats_region<-ldply(c(1:nrow(cmat_redux)),
                    function(x) data.frame(id=rownames(cmat_redux)[x],
                                           N=sum(cmat_redux[,x]),
                                           d=cmat_redux[x,x],
                                           A=sum(cmat_redux[x,]))) %>% transform(p=d/N,
                                                                                 r=d/A)

write.csv(stats_region, "results_masked/stats_region_bows.csv")


bow_redux$predReg<-predict(rf_redux)
missclass<-cbind(bows[rownames(bow_redux[bow_redux$Region!=bow_redux$predReg,]),],
                 data.frame(predReg=bow_redux$predReg[bow_redux$Region!=bow_redux$predReg]))
missclass
write.csv(missclass, "results_masked/bows_missclass.csv")

##########
## FAMD ##
##########

bows <- read.csv("bows.csv", stringsAsFactors = T,header=T)
#bows$Region <- NULL
imp.bows <- imputeFAMD(bows[,2:9], ncp = 10)
res.famd <- FAMD(imp.bows$completeObs,   ## Set the target variable "Churn" as a supplementary variable, so it is not included in the analysis for now
                 graph = T, 
                 sup.var=1, # Culture is the supplementary variable
                 ncp=10)

## Inspect principal components
get_eigenvalue(res.famd)
fviz_screeplot(res.famd)
# Plot of variables
fviz_famd_var(res.famd, repel = TRUE, ggtheme=theme_classic())
ggsave("results_masked/bows_FAMD_varcont_dim.png", width=4.5, height=4, units="in", dpi=500)
# Contribution to the first dimension

fviz_contrib(res.famd, "var", axes = 1, ggtheme=theme_classic()) 
ggsave("results_masked/bows_FAMD_varDIM1.png", width=7, height=4.5, units="in", dpi=500)
# Contribution to the second dimension
fviz_contrib(res.famd, "var", axes = 2, ggtheme=theme_classic())
ggsave("results_masked/bows_FAMD_varDIM2.png", width=7, height=4.5, units="in", dpi=500)

fviz_famd_var(res.famd, "quanti.var", col.var = "contrib", 
              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
              repel = TRUE)

famd_plot_bows <- fviz_mfa_ind(res.famd, 
                          habillage = 1, # color by groups 
                          palette = my_pal_long,
                          addEllipses = F, ellipse.type = "confidence",
                          repel = TRUE, # Avoid text overlapping,
                          label="none",
                          ggtheme = theme_classic(),
) 

famd_plot_bows <- famd_plot_bows + ggtitle("Bows - FAMD")
famd_plot_bows
ggsave("results_masked/Bows_FAMD_noell.tiff", width=7, height=6, units="in", dpi=500)

############
## SPEARS ##
############

spears <- read.csv("spears.csv", stringsAsFactors = T,header=T)
ids <- seq(1:nrow(spears))
for (i in ids) {
  ids[i] <- paste("T-", i, sep="")
}

spears$ID <- ids

# For genetics, create another spear dataset with only spears from groups with genes
spears_red <- subset(spears, spears$Culture %in% rownames(Fst_mat))

spears$Trait <- NULL
ids <- spears$ID
ids_red <- spears_red$ID
rownames(spears) <- ids
rownames(spears_red) <- ids_red
spears$ID <- NULL
spears$Culture <- NULL
spears_red$Trait <- NULL
spears_red$ID <- NULL
spears_red$Culture <- NULL

# this file needs to have the coordinates of the culture corresponding to each individual spear
# in the same order as your spears.csv file
locationsspear <- read.csv("spears_locations.csv")
locationsspear_red <- subset(locationsspear, locationsspear$Culture %in% rownames(Fst_mat))

rownames(locationsspear) <- rownames(spears)
rownames(locationsspear_red) <- rownames(spears_red)

biomesspear <- read.csv("spears_biomes.csv")
biomesspear_red <- subset(biomesspear, biomesspear$Culture %in% rownames(Fst_mat))
biomesspear$Trait <- NULL
biomesspear$Culture <- NULL
biomesspear_red$Trait <- NULL
biomesspear_red$Culture <- NULL

rownames(biomesspear) <- rownames(spears)
rownames(biomesspear_red) <- rownames(spears_red)

# Convert factors into numeric codes (not really needed if stringsAsFactors = F)
indx <- sapply(spears, is.factor)
indx_red <- sapply(spears_red, is.factor)

spears[indx] <- lapply(spears[indx],  function(x) as.numeric(as.factor(x)))
spears_red[indx_red] <- lapply(spears_red[indx_red],  function(x) as.numeric(as.factor(x)))

#compute Gower dissimilarity distance
library(vegan)
gowspear <- vegdist(spears,"gower",na.rm=TRUE) 
gowspearRed <- vegdist(spears_red,"gower",na.rm=TRUE) 

gowBiomespears <- as.matrix(vegdist(biomesspear,"gower",na.rm=TRUE))
gowBiomespearsRed <-  as.matrix(vegdist(biomesspear_red,"gower",na.rm=TRUE))


#compute binary Matrices of cultural affiliation (0:same culture; 1:different culture)
binaryCultDistspears<-as.matrix(as.dist(BinaryCulture(locationsspear$Culture)))
binaryCultDistspearsRed<-as.matrix(as.dist(BinaryCulture(locationsspear_red$Culture)))

#compute Fst matrix
GenDistspears<-as.matrix(as.dist(GetFst(locationsspear_red$Culture)))

#compute spatial distance between cultures and therefore spears
#spacespears<-distMat(data.frame(x=locationsspear$Longitude,y=locationsspear$Latitude))

#compute spatial distances
distTablespears <- data.frame(x=locationsspear$Longitude,y=locationsspear$Latitude)
namesspear <- as.vector(row.names(locationsspear))
rownames(distTablespears) <- namesspear

distTablespearsRed <- data.frame(x=locationsspear_red$Longitude,y=locationsspear_red$Latitude)
namesspearRed <- as.vector(row.names(locationsspear_red))
rownames(distTablespearsRed) <- namesspearRed

library(geosphere)
spacespears <- distm(distTablespears, fun=distGeo)
rownames(spacespears) <- namesspear
colnames(spacespears) <- namesspear
# here distance is in metres, not kilometres
#spacespears <- as.dist(spacespears)

spacespearsRed <- distm(distTablespearsRed, fun=distGeo)
rownames(spacespearsRed) <- namesspearRed
colnames(spacespearsRed) <- namesspearRed
# here distance is in metres, not kilometres
#spacespearsRed <- as.dist(spacespearsRed)

gowspearRed <- as.matrix(gowspearRed)
na_rows <- which(is.na(gowspearRed), arr.ind=TRUE)
na_rows <- unique(na_rows[,1])
gowspearRed <- gowspearRed[-na_rows, -na_rows]
spacespearsRed  <- spacespearsRed [-na_rows, -na_rows]
gowBiomespearsRed  <- gowBiomespearsRed [-na_rows, -na_rows]
binaryCultDistspearsRed  <- binaryCultDistspearsRed [-na_rows, -na_rows]
GenDistspears <- GenDistspears[-na_rows, -na_rows]

gowspear <- as.matrix(gowspear)
na_rows <- which(is.na(gowspear), arr.ind=TRUE)
na_rows <- unique(na_rows[,1])
gowspear <- gowspear[-na_rows, -na_rows]
spacespears <- spacespears[-na_rows, -na_rows]
gowBiomespears <- gowBiomespears[-na_rows, -na_rows]
binaryCultDistspears <- binaryCultDistspears[-na_rows, -na_rows]




#partial mantel test between Distance and Geography and between Culture and Geography
library(ecodist)
detach("package:vegan", unload = TRUE)

space_mantel <- mantel(lower(gowspear)~lower(spacespears),nperm=1000)
culture_mantel <-mantel(lower(gowspear)~lower(binaryCultDistspears),nperm=1000)
genes_mantel <- mantel(lower(gowspearRed)~lower(GenDistspears),nperm=1000)
biome_mantel <- mantel(lower(gowBiomespearsRed)~lower(GenDistspears),nperm=1000)

#partial mantel test testing for effect of CULTURAL AFFILIATION controlling for rest
cult_space_mantel <- mantel(lower(gowspear)~lower(binaryCultDistspears)+lower(spacespears),nperm=1000)
cult_genes_mantel <- mantel(lower(gowspearRed)~lower(binaryCultDistspearsRed)+lower(GenDistspears),nperm=1000)
cult_biome_mantel <- mantel(lower(gowspear)~lower(binaryCultDistspears)+lower(gowBiomespears),nperm=1000)

#partial mantel test testing for effect of SHARED GENES
genes_space_mantel <- mantel(lower(gowspearRed)~lower(GenDistspears)+lower(spacespearsRed),nperm=1000)
genes_cult_mantel <- mantel(lower(gowspearRed)~lower(GenDistspears)+lower(binaryCultDistspearsRed),nperm=1000)
genes_biome_mantel <- mantel(lower(gowspearRed)~lower(GenDistspears)+lower(gowBiomespearsRed),nperm=1000)

#partial mantel test testing for effect of GEOGRAPHY
space_culture_mantel <- mantel(lower(gowspear)~lower(spacespears)+lower(binaryCultDistspears),nperm=1000)
space_genes_mantel <- mantel(lower(gowspearRed)~lower(spacespearsRed)+lower(GenDistspears),nperm=1000)
space_biome_mantel <- mantel(lower(gowspear)~lower(spacespears)+lower(gowBiomespears),nperm=1000)

#partial mantel test testing for effect of BIOME
biome_culture_mantel <- mantel(lower(gowspear)~lower(gowBiomespears)+lower(binaryCultDistspears),nperm=1000)
biome_genes_mantel <- mantel(lower(gowspearRed)~lower(gowBiomespearsRed)+lower(GenDistspears),nperm=1000)
biome_space_mantel <- mantel(lower(gowspear)~lower(gowBiomespears)+lower(spacespears),nperm=1000)


write.csv(rbind(space_mantel, culture_mantel, genes_mantel,biome_mantel,
                cult_space_mantel,cult_genes_mantel,cult_biome_mantel,
                genes_space_mantel, genes_cult_mantel,genes_biome_mantel,
                space_culture_mantel,space_genes_mantel,space_biome_mantel,
                biome_culture_mantel, biome_genes_mantel, biome_space_mantel), "results_masked/spears_mantel.csv")



#########
# AMOVA #
#########
library(pegas)

#AMOVA
sqgowspear<-sqrt(gowspear)
sqgowspear_mat <- as.data.frame(as.matrix(sqgowspear))
# # get NA
#na_rows <- which(is.na(sqgowspear_mat), arr.ind=TRUE)
# 
# # remove rows with NAs
#na_rows <- unique(na_rows[,1])
#sqgowspear_mat <- sqgowspear_mat[-na_rows,-na_rows]
sqgowspear <- as.dist(sqgowspear_mat)

groupsspear<-as.factor(locationsspear$Culture)
groupsspear <- as.factor(groupsspear[-na_rows])
AMOVAspear=pegas::amova(sqgowspear~groupsspear)

# Phi-Stat:
AMOVAspear$varcomp$sigma2[1]/sum(AMOVAspear$varcomp$sigma2)
#P-Value:
AMOVAspear$varcomp$P.value[1]


##################
# pairwise AMOVA #
##################

# Note: This is the important thing to compare against genetic data

# I think you cannot do amova for cultures with only 1 trait!

#Extract Culture Data
CultureNames<-unique(as.character(groupsspear))
CultureNames <- sort(CultureNames)
NumberOfCultures<-length(CultureNames)
sqgowspear<-as.matrix(sqgowspear)
CultureList<-groupsspear

# Remove the ones with only 1 spear
CultureNames<- CultureNames[-c(4,8)]
NumberOfCultures<-length(CultureNames)

#PhiSt and Pval Output Matrices:
PhiMatrix<-matrix(NA,nrow=NumberOfCultures,ncol=NumberOfCultures)
PvalMatrix<-matrix(NA,nrow=NumberOfCultures,ncol=NumberOfCultures)

#Spatial Coordinates Output Data.Frame:
XY<-data.frame(x=numeric(length=NumberOfCultures),y=numeric(length=NumberOfCultures))

for (i in 1:NumberOfCultures) {
  print(paste(i," out of ",NumberOfCultures))
  for (j in 1:NumberOfCultures)
  {
    if (i!=j&i>j)
    { 
      indexI=which(as.character(CultureList)==as.character(CultureNames[i]))
      indexJ=which(as.character(CultureList)==as.character(CultureNames[j]))
      #subsetDistance Matrix
      subDistMat<-sqgowspear[c(indexI,indexJ),c(indexI,indexJ)]
      subDistMat<-as.dist(subDistMat)
      classes<-as.factor(c(rep("a",length(indexI)),rep("b",length(indexJ))))
      #compute Phi:
      tmp<-pegas::amova(subDistMat~classes)
      N0=tmp$varcoef
      AP=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)
      WP=(tmp$tab$MSD[2])
      PhiMatrix[i,j]=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)/((tmp$tab$MSD[2])+((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef))
      if(PhiMatrix[i,j]<0){PhiMatrix[i,j]=0}
      PvalMatrix[i,j]=tmp$varcomp$P.value[1]
    }
  }
  #XY$x[i]<-mean(locationsspear$Longitude[which(locationsspear$Culture==as.character(CultureNames[i]))],na.rm=TRUE)
  #XY$y[i]<-mean(locationsspear$Latitude[which(locationsspear$Culture==as.character(CultureNames[i]))],na.rm=TRUE)
}

dim <- ncol(PhiMatrix)

#Change NAs to blanks for text
PhiMatrix_text <- as.data.frame(PhiMatrix)
PhiMatrix_text <- sapply(PhiMatrix_text, as.character)
PhiMatrix_text[is.na(PhiMatrix_text)] <- " "
PhiMatrix_text <- as.matrix(PhiMatrix_text)

png(file="results_masked/spears_phistats.png",
    width=7, height=6, units="in", res=500)
image(-7:(dim-8), 1:dim, PhiMatrix, axes = FALSE, xlab="", ylab="", col=hcl.colors(12, "Geyser", rev=TRUE), main="spears", xlim=c(-7,1), ylim=c(0,11))
axis(1, (-7:(dim-8)), CultureNames, cex.axis = 0.6, las=3, pos=0.5)
axis(4, 1:dim, CultureNames, cex.axis = 0.6, las=1, pos=-0.5)
text(expand.grid((-7:(dim-8)), 1:dim), sprintf("%0.4s", PhiMatrix_text), cex=0.6)
dev.off()

# Again Partial Mantel Tests using Phi Statistic

CultureNames_rn <- CultureNames
CultureNames_rn[CultureNames_rn=="Batwa (East)"] <- "BatwaE"
CultureNames_rn[CultureNames_rn=="Batwa (West)"] <- "BatwaW"
rownames(PhiMatrix) <- CultureNames
colnames(PhiMatrix) <- CultureNames
PhiMatrix_rn <- PhiMatrix
rownames(PhiMatrix_rn) <- CultureNames_rn
colnames(PhiMatrix_rn) <- CultureNames_rn

writeDist(PhiMatrix_rn, file="spears_phi.nex", format="nexus")

# Mantel tests between geography, biome, genetics and spears
PhiMatrix <- as.matrix(PhiMatrix)
genDist <- as.matrix(genDist)
gowBiomeRed <- as.matrix(gowBiomeRed)
spaceRed <- as.matrix(spaceRed)

# Remove groups with >2 spears
genDistspear <- genDist[-c(2,5,7),-c(2,5,7)]
BiomeDistspear <- gowBiomeRed[-c(2,5,7),-c(2,5,7)]
spacespear <- spaceRed[-c(2,5,7),-c(2,5,7)]
PhiMatrixspears <- PhiMatrix[-5,-5]

biome_mantel <- mantel(lower(PhiMatrixspears)~lower(BiomeDistspear), mrank=T)
space_mantel <- mantel(lower(PhiMatrixspears)~lower(spacespear), mrank=T)
genes_mantel <- mantel(lower(PhiMatrixspears)~lower(genDistspear), mrank=T)

#partial mantel test testing for effect of SHARED GENES
genes_space_mantel <- mantel(lower(PhiMatrixspears)~lower(genDistspear)+lower(spacespear), mrank=T)
genes_biome_mantel <- mantel(lower(PhiMatrixspears)~lower(genDistspear)+lower(BiomeDistspear), mrank=T)

#partial mantel test testing for effect of GEOGRAPHY
space_genes_mantel <- mantel(lower(PhiMatrixspears)~lower(spacespear)+lower(genDistspear), mrank=T)
space_biome_mantel <- mantel(lower(PhiMatrixspears)~lower(spacespear)+lower(BiomeDistspear), mrank=T)

#partial mantel test testing for effect of BIOME
biome_genes_mantel <- mantel(lower(PhiMatrixspears)~lower(BiomeDistspear)+lower(genDistspear), mrank=T)
biome_space_mantel <- mantel(lower(PhiMatrixspears)~lower(BiomeDistspear)+lower(spacespear), mrank=T)


mantel_spears <- rbind(space_mantel, genes_mantel,
                            biome_mantel,
                            genes_space_mantel, genes_biome_mantel,
                            space_genes_mantel, space_biome_mantel,
                            biome_space_mantel, biome_genes_mantel)

write.csv(mantel_spears, "results_masked/mantel_spears.csv")



# Perform PCoA
gowspear <- as.data.frame(as.matrix(gowspear))
# get NA
na_rows <- which(is.na(gowspear), arr.ind=TRUE)

# remove rows with NAs
#na_rows <- unique(na_rows[,1])
#gowspear <- gowspear[-na_rows,-na_rows]
#gowspear <- as.dist(gowspear)
spears_pcoa <- pcoa(gowspear, correction = "cailliez")

#We rescale the PCoA components in relation to the explained variance.  

for(i in 1:ncol(spears_pcoa$vectors)) { 
  spears_pcoa$vectors[,i] <- scale(spears_pcoa$vectors[,i])*
    spears_pcoa$values$Rel_corr_eig[i]}

# Extract eigenvalues from PC/PCo results_masked 
spears_ev <- spears_pcoa$values$Corr_eig

# Plot
library(ggplot2)
my_font <- "sans"

png(file="results_masked/pcoaspears.png",
    width=6, height=4.5, units="in", res=500)
plot_pc_scree(eigenval=spears_ev, title="Genetics", type="PCo")
dev.off()

##We extract the principal components/coordinates for all factors: 
#genetics, grammar, music and phonology. 
#We normalize the PCs/PCos to a range from 0 to 1 and then plot a heatmap 

# Extract the PCs and the PCos from the PCA and PCoA results_masked
CultureList<-groupsspear
spears_pco <- as.data.frame(spears_pcoa$vectors)
spears_pco$Culture <- CultureList

## PCoA with groups

#install.packages('ggfortify')

# remove Mbendjele from palette
my_pal_long_ar <- my_pal_long[-c(2,8)]

p <- ggplot(data = spears_pco, aes(x = Axis.1, y = Axis.2, color=Culture)) +
  geom_point(size = 3.5) +
  xlab("PCo1 (13.56%)") +
  ylab("PCo2 (9.15%)") + 
  theme_classic() + 
  scale_color_manual(values = my_pal_long_ar,
                     labels = c("Aka (8)",
                                "Baka (2)",
                                "Bakola (8)",
                                "Bakoya (1)",
                                "East Batwa (8)",
                                "West Batwa (14)",
                                "Efe (22)",
                                "Mbendjele (1)",
                                "Sua (7)")) +
  geom_hline(yintercept = 0, lty = 2) +
  ggtitle("Spears") +
  geom_vline(xintercept = 0, lty = 2) 

p <- p + theme(legend.title = element_blank(),
               legend.text = element_text(size = 12))

ggsave(p, filename="results_masked/PCOA_spears.png", height=150, width=180, units="mm")

#We compute Euclidean distance matrices from the dimensionality-reduced data and create NeighborNets. 
# SplitsTree doesnt accept brackets in names
spears_pco_rn <- spears_pco
spears_pco_rn$Culture <- as.character(spears_pco_rn$Culture)
spears_pco_rn$Culture[spears_pco_rn$Culture=="Batwa (East)"] <- "BatwaE"
spears_pco_rn[spears_pco_rn=="Batwa (West)"] <- "BatwaW"
spears_pco_rn$Culture <- NULL

# Write nexus files for splitstree (Not used so far)
require(phangorn)

dist_list <- lapply(list(spears=spears_pco_rn),
                    function(f) {dist(f)})
writeDist(dist_list$spears, file="results_masked/spears.nex", format="nexus")




#############
## ZITHERS ##
#############

zitherharps <- read.csv("zitherharps.csv", stringsAsFactors = T,header=T)
ids <- seq(1:nrow(zitherharps))
for (i in ids) {
  ids[i] <- paste("T-", i, sep="")
}

zitherharps$ID <- ids

# For genetics, create another zitherharp dataset with only zitherharps from groups with genes
zitherharps_red <- subset(zitherharps, zitherharps$Culture %in% rownames(Fst_mat))

zitherharps$Trait <- NULL
ids <- zitherharps$ID
ids_red <- zitherharps_red$ID
rownames(zitherharps) <- ids
rownames(zitherharps_red) <- ids_red
zitherharps$ID <- NULL
zitherharps$Culture <- NULL
zitherharps_red$Trait <- NULL
zitherharps_red$ID <- NULL
zitherharps_red$Culture <- NULL

# this file needs to have the coordinates of the culture corresponding to each individual zitherharp
# in the same order as your zitherharps.csv file
locationszitherharp <- read.csv("zitherharps_locations.csv")
locationszitherharp_red <- subset(locationszitherharp, locationszitherharp$Culture %in% rownames(Fst_mat))

rownames(locationszitherharp) <- rownames(zitherharps)
rownames(locationszitherharp_red) <- rownames(zitherharps_red)

biomeszitherharp <- read.csv("zitherharps_biomes.csv")
biomeszitherharp_red <- subset(biomeszitherharp, biomeszitherharp$Culture %in% rownames(Fst_mat))
biomeszitherharp$Trait <- NULL
biomeszitherharp$Culture <- NULL
biomeszitherharp_red$Trait <- NULL
biomeszitherharp_red$Culture <- NULL

rownames(biomeszitherharp) <- rownames(zitherharps)
rownames(biomeszitherharp_red) <- rownames(zitherharps_red)

# Convert factors into numeric codes (not really needed if stringsAsFactors = F)
indx <- sapply(zitherharps, is.factor)
indx_red <- sapply(zitherharps_red, is.factor)

zitherharps[indx] <- lapply(zitherharps[indx],  function(x) as.numeric(as.factor(x)))
zitherharps_red[indx_red] <- lapply(zitherharps_red[indx_red],  function(x) as.numeric(as.factor(x)))

#compute Gower dissimilarity distance
library(vegan)
gowzitherharp <- vegdist(zitherharps,"gower",na.rm=TRUE) 
gowzitherharpRed <- vegdist(zitherharps_red,"gower",na.rm=TRUE) 

gowBiomezitherharps <- as.matrix(vegdist(biomeszitherharp,"gower",na.rm=TRUE))
gowBiomezitherharpsRed <-  as.matrix(vegdist(biomeszitherharp_red,"gower",na.rm=TRUE))


#compute binary Matrices of cultural affiliation (0:same culture; 1:different culture)
binaryCultDistzitherharps<-as.matrix(as.dist(BinaryCulture(locationszitherharp$Culture)))
binaryCultDistzitherharpsRed<-as.matrix(as.dist(BinaryCulture(locationszitherharp_red$Culture)))

#compute Fst matrix
GenDistzitherharps<-as.matrix(as.dist(GetFst(locationszitherharp_red$Culture)))

#compute spatial distance between cultures and therefore zitherharps
#spacezitherharps<-distMat(data.frame(x=locationszitherharp$Longitude,y=locationszitherharp$Latitude))

#compute spatial distances
distTablezitherharps <- data.frame(x=locationszitherharp$Longitude,y=locationszitherharp$Latitude)
nameszitherharp <- as.vector(row.names(locationszitherharp))
rownames(distTablezitherharps) <- nameszitherharp

distTablezitherharpsRed <- data.frame(x=locationszitherharp_red$Longitude,y=locationszitherharp_red$Latitude)
nameszitherharpRed <- as.vector(row.names(locationszitherharp_red))
rownames(distTablezitherharpsRed) <- nameszitherharpRed

library(geosphere)
spacezitherharps <- distm(distTablezitherharps, fun=distGeo)
rownames(spacezitherharps) <- nameszitherharp
colnames(spacezitherharps) <- nameszitherharp
# here distance is in metres, not kilometres
#spacezitherharps <- as.dist(spacezitherharps)

spacezitherharpsRed <- distm(distTablezitherharpsRed, fun=distGeo)
rownames(spacezitherharpsRed) <- nameszitherharpRed
colnames(spacezitherharpsRed) <- nameszitherharpRed
# here distance is in metres, not kilometres
#spacezitherharpsRed <- as.dist(spacezitherharpsRed)

gowzitherharp <- as.matrix(gowzitherharp)
na_rows <- which(is.na(gowzitherharp), arr.ind=TRUE)
#na_rows <- unique(na_rows[,1])
#gowzitherharp <- gowzitherharp[-na_rows, -na_rows]
#spacezitherharps <- spacezitherharps[-na_rows, -na_rows]
#gowBiomezitherharps <- gowBiomezitherharps[-na_rows, -na_rows]
#binaryCultDistzitherharps <- binaryCultDistzitherharps[-na_rows, -na_rows]

gowzitherharpRed <- as.matrix(gowzitherharpRed)
na_rows <- which(is.na(gowzitherharpRed), arr.ind=TRUE)
#na_rows <- unique(na_rows[,1])
#gowzitherharpRed <- gowzitherharpRed[-na_rows, -na_rows]
#spacezitherharpsRed  <- spacezitherharpsRed [-na_rows, -na_rows]
#gowBiomezitherharpsRed  <- gowBiomezitherharpsRed [-na_rows, -na_rows]
#binaryCultDistzitherharpsRed  <- binaryCultDistzitherharpsRed [-na_rows, -na_rows]
#GenDistzitherharps <- GenDistzitherharps[-na_rows, -na_rows]


#partial mantel test between Distance and Geography and between Culture and Geography
library(ecodist)
detach("package:vegan", unload = TRUE)

space_mantel <- mantel(lower(gowzitherharp)~lower(spacezitherharps),nperm=1000)
culture_mantel <-mantel(lower(gowzitherharp)~lower(binaryCultDistzitherharps),nperm=1000)
genes_mantel <- mantel(lower(gowzitherharpRed)~lower(GenDistzitherharps),nperm=1000)
biome_mantel <- mantel(lower(gowBiomezitherharpsRed)~lower(GenDistzitherharps),nperm=1000)

#partial mantel test testing for effect of CULTURAL AFFILIATION controlling for rest
cult_space_mantel <- mantel(lower(gowzitherharp)~lower(binaryCultDistzitherharps)+lower(spacezitherharps),nperm=1000)
cult_genes_mantel <- mantel(lower(gowzitherharpRed)~lower(binaryCultDistzitherharpsRed)+lower(GenDistzitherharps),nperm=1000)
cult_biome_mantel <- mantel(lower(gowzitherharp)~lower(binaryCultDistzitherharps)+lower(gowBiomezitherharps),nperm=1000)

#partial mantel test testing for effect of SHARED GENES
genes_space_mantel <- mantel(lower(gowzitherharpRed)~lower(GenDistzitherharps)+lower(spacezitherharpsRed),nperm=1000)
genes_cult_mantel <- mantel(lower(gowzitherharpRed)~lower(GenDistzitherharps)+lower(binaryCultDistzitherharpsRed),nperm=1000)
genes_biome_mantel <- mantel(lower(gowzitherharpRed)~lower(GenDistzitherharps)+lower(gowBiomezitherharpsRed),nperm=1000)

#partial mantel test testing for effect of GEOGRAPHY
space_culture_mantel <- mantel(lower(gowzitherharp)~lower(spacezitherharps)+lower(binaryCultDistzitherharps),nperm=1000)
space_genes_mantel <- mantel(lower(gowzitherharpRed)~lower(spacezitherharpsRed)+lower(GenDistzitherharps),nperm=1000)
space_biome_mantel <- mantel(lower(gowzitherharp)~lower(spacezitherharps)+lower(gowBiomezitherharps),nperm=1000)

#partial mantel test testing for effect of BIOME
biome_culture_mantel <- mantel(lower(gowzitherharp)~lower(gowBiomezitherharps)+lower(binaryCultDistzitherharps),nperm=1000)
biome_genes_mantel <- mantel(lower(gowzitherharpRed)~lower(gowBiomezitherharpsRed)+lower(GenDistzitherharps),nperm=1000)
biome_space_mantel <- mantel(lower(gowzitherharp)~lower(gowBiomezitherharps)+lower(spacezitherharps),nperm=1000)


write.csv(rbind(space_mantel, culture_mantel, genes_mantel,biome_mantel,
                cult_space_mantel,cult_genes_mantel,cult_biome_mantel,
                genes_space_mantel, genes_cult_mantel,genes_biome_mantel,
                space_culture_mantel,space_genes_mantel,space_biome_mantel,
                biome_culture_mantel, biome_genes_mantel, biome_space_mantel), "results_masked/zitherharps_mantel.csv")


#########
# AMOVA #
#########
library(pegas)

#AMOVA
sqgowzitherharp<-sqrt(gowzitherharp)
sqgowzitherharp_mat <- as.data.frame(as.matrix(sqgowzitherharp))
# # get NA
# na_rows <- which(is.na(sqgowzitherharp_mat), arr.ind=TRUE)
# # 
# # # remove rows with NAs
# na_rows <- unique(na_rows[,1])
# sqgowzitherharp_mat <- sqgowzitherharp_mat[-na_rows,-na_rows]
# sqgowzitherharp <- as.dist(sqgowzitherharp_mat)

groupszitherharp<- as.factor(locationszitherharp$Culture)
# groupszitherharp <- as.factor(groupszitherharp[-na_rows])
AMOVAzitherharp=amova(sqgowzitherharp~groupszitherharp)

# Phi-Stat:
AMOVAzitherharp$varcomp$sigma2[1]/sum(AMOVAzitherharp$varcomp$sigma2)
#P-Value:
AMOVAzitherharp$varcomp$P.value[1]


##################
# pairwise AMOVA #
##################

# Note: This is the important thing to compare against genetic data

# I think you cannot do amova for cultures with only 1 trait!

#Extract Culture Data
CultureNames<-unique(as.character(groupszitherharp))
CultureNames <- sort(CultureNames)
NumberOfCultures<-length(CultureNames)
sqgowzitherharp<-as.matrix(sqgowzitherharp)
CultureList<-groupszitherharp

# Remove the ones with only 1 zitherharp
CultureNames<- CultureNames[-c(7)]
NumberOfCultures<-length(CultureNames)

#PhiSt and Pval Output Matrices:
PhiMatrix<-matrix(NA,nrow=NumberOfCultures,ncol=NumberOfCultures)
PvalMatrix<-matrix(NA,nrow=NumberOfCultures,ncol=NumberOfCultures)

#Spatial Coordinates Output Data.Frame:
XY<-data.frame(x=numeric(length=NumberOfCultures),y=numeric(length=NumberOfCultures))

for (i in 1:NumberOfCultures) {
  print(paste(i," out of ",NumberOfCultures))
  for (j in 1:NumberOfCultures)
  {
    if (i!=j&i>j)
    { 
      indexI=which(as.character(CultureList)==as.character(CultureNames[i]))
      indexJ=which(as.character(CultureList)==as.character(CultureNames[j]))
      #subsetDistance Matrix
      subDistMat<-sqgowzitherharp[c(indexI,indexJ),c(indexI,indexJ)]
      subDistMat<-as.dist(subDistMat)
      classes<-as.factor(c(rep("a",length(indexI)),rep("b",length(indexJ))))
      #compute Phi:
      tmp<-pegas::amova(subDistMat~classes)
      N0=tmp$varcoef
      AP=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)
      WP=(tmp$tab$MSD[2])
      PhiMatrix[i,j]=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)/((tmp$tab$MSD[2])+((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef))
      if(PhiMatrix[i,j]<0){PhiMatrix[i,j]=0}
      PvalMatrix[i,j]=tmp$varcomp$P.value[1]
    }
  }
  #XY$x[i]<-mean(locationszitherharp$Longitude[which(locationszitherharp$Culture==as.character(CultureNames[i]))],na.rm=TRUE)
  #XY$y[i]<-mean(locationszitherharp$Latitude[which(locationszitherharp$Culture==as.character(CultureNames[i]))],na.rm=TRUE)
}

dim <- ncol(PhiMatrix)

#Change NAs to blanks for text
PhiMatrix_text <- as.data.frame(PhiMatrix)
PhiMatrix_text <- sapply(PhiMatrix_text, as.character)
PhiMatrix_text[is.na(PhiMatrix_text)] <- " "
PhiMatrix_text <- as.matrix(PhiMatrix_text)

png(file="results_masked/zitherharps_phistats.png",
    width=7, height=6, units="in", res=500)
image(-7:(dim-8), 1:dim, PhiMatrix, axes = FALSE, xlab="", ylab="", col=hcl.colors(12, "Geyser", rev=TRUE), main="Zither harps", xlim=c(-7,1), ylim=c(0,11))
axis(1, (-7:(dim-8)), CultureNames, cex.axis = 0.6, las=3, pos=0.5)
axis(4, 1:dim, CultureNames, cex.axis = 0.6, las=1, pos=-0.5)
text(expand.grid((-7:(dim-8)), 1:dim), sprintf("%0.4s", PhiMatrix_text), cex=0.6)
dev.off()

# Again Partial Mantel Tests using Phi Statistic

CultureNames_rn <- CultureNames
CultureNames_rn[CultureNames_rn=="Batwa (East)"] <- "BatwaE"
CultureNames_rn[CultureNames_rn=="Batwa (West)"] <- "BatwaW"
rownames(PhiMatrix) <- CultureNames
colnames(PhiMatrix) <- CultureNames
PhiMatrix_rn <- PhiMatrix
rownames(PhiMatrix_rn) <- CultureNames_rn
colnames(PhiMatrix_rn) <- CultureNames_rn

writeDist(PhiMatrix_rn, file="zitherharps_phi.nex", format="nexus")

# Mantel tests between geography, biome, genetics and zitherharps
PhiMatrix <- as.matrix(PhiMatrix)
genDist <- as.matrix(genDist)
gowBiomeRed <- as.matrix(gowBiomeRed)
spaceRed <- as.matrix(spaceRed)

# Remove groups with <2 zitherharps
genDistzitherharp <- genDist[-c(2,8,9),-c(2,8,9)]
genDistzitherharpIBD <- genDistIBD[-c(2,8,9),-c(2,8,9)]
genDistzitherharpBantu <- genDistBantu[-c(2,8),-c(2,8)]
BiomeDistzitherharp <- gowBiomeRed[-c(2,8,9),-c(2,8,9)]
spacezitherharp <- spaceRed[-c(2,8,9),-c(2,8,9)]
PhiMatrixzitherharps <- PhiMatrix
#PhiMatrix<- as.dist(PhiMatrix[-5,-5])

biome_mantel <- mantel(lower(PhiMatrixzitherharps)~lower(BiomeDistzitherharp), mrank=T)
space_mantel <- mantel(lower(PhiMatrixzitherharps)~lower(spacezitherharp), mrank=T)
genes_mantel <- mantel(lower(PhiMatrixzitherharps)~lower(genDistzitherharp), mrank=T)
genes_ibd_mantel <- mantel(lower(PhiMatrixzitherharps)~lower(genDistzitherharpIBD), mrank=T)
genes_bantu_mantel <- mantel(lower(PhiMatrixzitherharps[-7,-7])~lower(genDistzitherharpBantu), mrank=T)

#partial mantel test testing for effect of SHARED GENES
genes_space_mantel <- mantel(lower(PhiMatrixzitherharps)~lower(genDistzitherharp)+lower(spacezitherharp), mrank=T)
genes_biome_mantel <- mantel(lower(PhiMatrixzitherharps)~lower(genDistzitherharp)+lower(BiomeDistzitherharp), mrank=T)

#partial mantel test testing for effect of GEOGRAPHY
space_genes_mantel <- mantel(lower(PhiMatrixzitherharps)~lower(spacezitherharp)+lower(genDistzitherharp), mrank=T)
space_biome_mantel <- mantel(lower(PhiMatrixzitherharps)~lower(spacezitherharp)+lower(BiomeDistzitherharp), mrank=T)

#partial mantel test testing for effect of BIOME
biome_genes_mantel <- mantel(lower(PhiMatrixzitherharps)~lower(BiomeDistzitherharp)+lower(genDistzitherharp), mrank=T)
biome_space_mantel <- mantel(lower(PhiMatrixzitherharps)~lower(BiomeDistzitherharp)+lower(spacezitherharp), mrank=T)


mantel_zitherharps <- rbind(space_mantel, genes_mantel,
                            biome_mantel,
                            genes_space_mantel, genes_biome_mantel,
                            space_genes_mantel, space_biome_mantel,
                            biome_space_mantel, biome_genes_mantel)

write.csv(mantel_zitherharps, "results_masked/mantel_zitherharps.csv")

mantel_zitherharps <-  rbind(space_mantel, genes_mantel,
                       biome_mantel, genes_ibd_mantel, genes_bantu_mantel)
pval_adj_zitherharps <- p.adjust(mantel_zitherharps[,2], method = "BH")
mantel_zitherharps<- cbind(mantel_zitherharps, pval_adj_zitherharps)

write.csv(mantel_zitherharps, "final_paper_files/results_final_paper/results_masked/mantel_zitherharps.csv")

matrix_zitherharps <- rbind(genes_space_mantel, genes_biome_mantel,
                      space_genes_mantel, space_biome_mantel,
                      biome_space_mantel, biome_genes_mantel)
pval_adj_zitherharps <- p.adjust(matrix_zitherharps[,2], method = "BH")
matrix_zitherharps<- cbind(matrix_zitherharps, pval_adj_zitherharps)

write.csv(matrix_zitherharps, "results_masked/matrix_zitherharps.csv")


# Perform PCoA
gowzitherharp <- as.data.frame(as.matrix(gowzitherharp))
# get NA
na_rows <- which(is.na(gowzitherharp), arr.ind=TRUE)

# remove rows with NAs
#na_rows <- unique(na_rows[,1])
#gowzitherharp <- gowzitherharp[-na_rows,-na_rows]
#gowzitherharp <- as.dist(gowzitherharp)
zitherharps_pcoa <- pcoa(gowzitherharp, correction = "cailliez")

#We rescale the PCoA components in relation to the explained variance.  

for(i in 1:ncol(zitherharps_pcoa$vectors)) { 
  zitherharps_pcoa$vectors[,i] <- scale(zitherharps_pcoa$vectors[,i])*
    zitherharps_pcoa$values$Rel_corr_eig[i]}

# Extract eigenvalues from PC/PCo results_masked 
zitherharps_ev <- zitherharps_pcoa$values$Corr_eig

# Plot
library(ggplot2)
my_font <- "sans"

png(file="results_masked/pcoazitherharps.png",
    width=6, height=4.5, units="in", res=500)
plot_pc_scree(eigenval=zitherharps_ev, title="Genetics", type="PCo")
dev.off()

##We extract the principal components/coordinates for all factors: 
#genetics, grammar, music and phonology. 
#We normalize the PCs/PCos to a range from 0 to 1 and then plot a heatmap 

# Extract the PCs and the PCos from the PCA and PCoA results_masked
CultureList<-groupszitherharp
zitherharps_pco <- as.data.frame(zitherharps_pcoa$vectors)
zitherharps_pco$Culture <- CultureList

## PCoA with groups

#install.packages('ggfortify')

# remove Mbendjele from palette
my_pal_long_ar <- my_pal_long[-c(2,7,9)]

p <- ggplot(data = zitherharps_pco, aes(x = Axis.1, y = Axis.2, color=Culture)) +
  geom_point(size = 3.5) +
  xlab("PCo1 (1.91%)") +
  ylab("PCo2 (0.92%)") + 
  theme_classic() + 
  scale_color_manual(values = my_pal_long_ar,
                     labels = c("Aka (4)",
                                "Baka (7)",
                                "Bakola (2)",
                                "Bakoya (3)",
                                "East Batwa (2)",
                                "Bedzan (4)",
                                "Mbendjele (1)",
                                "Sua (3)")) +
  geom_hline(yintercept = 0, lty = 2) +
  ggtitle("Zither harps") +
  geom_vline(xintercept = 0, lty = 2) 

p <- p + theme(legend.title = element_blank(),
               legend.text = element_text(size = 12))

ggsave(p, filename="results_masked/PCOA_zitherharps.png", height=150, width=180, units="mm")

#We compute Euclidean distance matrices from the dimensionality-reduced data and create NeighborNets. 
# SplitsTree doesnt accept brackets in names
zitherharps_pco_rn <- zitherharps_pco
zitherharps_pco_rn$Culture <- as.character(zitherharps_pco_rn$Culture)
zitherharps_pco_rn$Culture[zitherharps_pco_rn$Culture=="Batwa (East)"] <- "BatwaE"
zitherharps_pco_rn[zitherharps_pco_rn=="Batwa (West)"] <- "BatwaW"
zitherharps_pco_rn$Culture <- NULL

# Write nexus files for splitstree (Not used so far)
require(phangorn)

dist_list <- lapply(list(zitherharps=zitherharps_pco_rn),
                    function(f) {dist(f)})
writeDist(dist_list$zitherharps, file="zitherharps.nex", format="nexus")

# Divide cultures into Western and Eastern Pygmies 

# Question - should we separate between different regions???
zitherharps <- read.csv("zitherharps.csv")

# Remove traits that are only present in one of the two regions (East or West)
zitherharps$Region<-sapply(as.character(zitherharps$Culture),function(x)
  ifelse(x %in% c("Aka", "Babongo", "Baka", "Bakola", "Bakoya",
                  "Bedzan", "Mbendjele", "Batwa (West)"),"West",
         ifelse(x %in% c("Batwa (East)","Efe","Sua"),"East")))
table(zitherharps$Region)

zitherharp_redux<-zitherharps[,c(3:ncol(zitherharps))]
for(i in 1:(ncol(zitherharp_redux)-1)){
  tab<-table(zitherharp_redux[,c("Region",colnames(zitherharp_redux)[i])])  
  check<-apply(tab,2,function(x) ifelse(sum(x!=0)==1,"fix","ok"))  
  zitherharp_redux[,colnames(zitherharp_redux)[i]]<-sapply(zitherharp_redux[,colnames(zitherharp_redux)[i]],
                                             function(x) ifelse(x %in% names(check)[check=="fix"],NA,as.character(x)))
  zitherharp_redux[,colnames(zitherharp_redux)[i]]<-as.factor(zitherharp_redux[,colnames(zitherharp_redux)[i]])}
zitherharp_redux<-zitherharp_redux[,sapply(zitherharp_redux, function(x)!all(is.na(x)))]
zitherharp_redux$Region<-as.factor(zitherharp_redux$Region)

library(party)
library(plyr)
weights_zitherharps<-sapply(zitherharp_redux$Region,
                     function(x) nrow(zitherharp_redux)/nrow(zitherharp_redux[zitherharp_redux$Region==x,]))
rf_redux<-cforest(Region~.,
                  controls=
                    cforest_control(ntree = 500,
                                    mtry=5,
                                    mincriterion=qnorm(0.9),
                                    fraction = 0.632,
                                    testtype ="Teststatistic",
                                    teststat="max",
                                    replace=TRUE,
                                    trace=F,
                                    savesplitstats=F,
                                    minsplit=20,
                                    minbucket=8),
                  data=zitherharp_redux,
                  weights = weights_zitherharps)

# Produce confusion matrix
cmat_redux<-table(zitherharp_redux$Region, predict(rf_redux))
print(cmat_redux)

stats_region<-ldply(c(1:nrow(cmat_redux)),
                    function(x) data.frame(id=rownames(cmat_redux)[x],
                                           N=sum(cmat_redux[,x]),
                                           d=cmat_redux[x,x],
                                           A=sum(cmat_redux[x,]))) %>% transform(p=d/N,
                                                                                 r=d/A)

write.csv(stats_region, "results_masked/stats_region_zitherharps.csv")


zitherharp_redux$predReg<-predict(rf_redux)
missclass<-cbind(zitherharps[rownames(zitherharp_redux[zitherharp_redux$Region!=zitherharp_redux$predReg,]),],
                 data.frame(predReg=zitherharp_redux$predReg[zitherharp_redux$Region!=zitherharp_redux$predReg]))
missclass
write.csv(missclass, "results_masked/zitherharps_missclass.csv")

##########
## FAMD ##
##########

zitherharps <- read.csv("zitherharps.csv", stringsAsFactors = T,header=T)
#zitherharps$Region <- NULL
imp.zitherharps <- imputeFAMD(zitherharps[,2:11], ncp = 10)
res.famd <- FAMD(imp.zitherharps$completeObs,   ## Set the target variable "Churn" as a supplementary variable, so it is not included in the analysis for now
                 graph = T, 
                 sup.var=1, # Culture is the supplementary variable
                 ncp=10)

## Inspect principal components
get_eigenvalue(res.famd)
fviz_screeplot(res.famd)
# Plot of variables
fviz_famd_var(res.famd, repel = TRUE, ggtheme=theme_classic())
ggsave("results_masked/zitherharps_FAMD_varcont_dim.png", width=4.5, height=4, units="in", dpi=500)
# Contribution to the first dimension

fviz_contrib(res.famd, "var", axes = 1, ggtheme=theme_classic()) 
ggsave("results_masked/zitherharps_FAMD_varDIM1.png", width=7, height=4.5, units="in", dpi=500)
# Contribution to the second dimension
fviz_contrib(res.famd, "var", axes = 2, ggtheme=theme_classic())
ggsave("results_masked/zitherharps_FAMD_varDIM2.png", width=7, height=4.5, units="in", dpi=500)

fviz_famd_var(res.famd, "quanti.var", col.var = "contrib", 
              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
              repel = TRUE)

my_pal_long_ar <- my_pal_long[-c(2,7,9)]

famd_plot_zitherharps <- fviz_mfa_ind(res.famd, 
                          habillage = 1, # color by groups 
                          palette = my_pal_long_ar,
                          addEllipses = F,
                          repel = TRUE, # Avoid text overlapping,
                          label="none",
                          ggtheme = theme_classic(),
) 

famd_plot_zitherharps <- famd_plot_zitherharps + ggtitle("Zither harps - FAMD")
famd_plot_zitherharps 
ggsave("results_masked/Zitherharps_FAMD.tiff", width=7, height=6, units="in", dpi=500)

library(ggpubr)
bb <- ggpubr::ggarrange(famd_plot_bows, famd_plot_aerophones, famd_plot_baskets,
                   famd_plot_baskets, famd_plot_arrows, famd_plot_zitherharps, ncol=2, nrow=3)
bb
ggsave("results_masked/all_FAMD.tiff", bb, width=8, height=9, units="in", dpi=1000, device = "tiff")

