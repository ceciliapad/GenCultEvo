################################
## Analyses cultural networks ##
################################

source("utility_functions.R")
require(MetBrewer)

load("SpatialCulture.RData")

#require(wesanderson)
my_pal <- c(MetPalettes$Isfahan1[[1]], MetPalettes$Isfahan2[[1]][2])

#load data (the excel worksheets should be individually exported as .csv files)

subsistence <- read.csv("P_A_subsistence.csv", row.names = 1)
music <- read.csv("P_A_music.csv", row.names = 1)

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

write.csv(rbind(music_mantel, subsistence_mantel), "simple_mantel_repertoires_with_bantu.csv")

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


###########
## WORDS ##
###########

# Do the same procedure first for UNIQUE words (unique to CAHG) and then for all shared
# words. Change filenames for "subsistence_words_all.csv" and "music_words_all.csv".

subsistence_words <- read.csv("subsistence_words_unique.csv", row.names = 1)[-10,]
music_words <- read.csv("instrument_words_unique.csv", row.names = 1)[-10,]

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
m_w <- as.data.frame(rowSums(music_words, na.rm=T))
median(m_w)
write.csv(s_w, "../results_sep_23/subsistence_words_counts_unique.csv")
write.csv(m_w, "../results_sep_23/music_words_counts_unique.csv")
s_w <- read.csv("../results_sep_23/subsistence_words_counts.csv", row.names = 1)
m_w <- read.csv("../results_sep_23/music_words_counts.csv", row.names = 1)

IDS_m <- rownames(m_w)
IDS_s <- rownames(s_w)

png("../results_sep_23/Musical_instruments_words_unique.png",units="in", width=7.5, height=5, res=500)
barplot(m_w$`rowSums(music_words, na.rm = T)`, col=my_pal_long, names.arg=IDS_m, ylab="Words available musical instruments", ylim=c(0,20), cex.names=0.5, cex.axis=0.5 )
dev.off()
png("../results_sep_23/Subsistence_tools_words_unique.png",units="in", width=7.5, height=5, res=500)
barplot(s_w$`rowSums(subsistence_words, na.rm = T)`, col=my_pal_long[-c(2,8)], names.arg=IDS_s, ylab="Words available subsistence tools", ylim=c(0,20), cex.names=0.5, cex.axis=0.5 )
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
# Now read in

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


require(pscl)
require(MASS)
require(boot)

# name dyads

words_df[is.na(words_df)]  <- 0

words_shared <- words_df %>%
  group_by(ID2, ID1, type) %>%
  summarise_at(vars(shared),
               list(total_shared=sum))

write.csv(words_shared, "../shared_words_counts_unique.csv")

# Read in again space data and genetic distance data
genes <- read.csv("fst_masked.csv", row.names = 1)
colnames(genes) <- rownames(genes)
genes_list <- melt(as.matrix(genes))
colnames(genes_list) <- c("ID2", "ID1", "gen_dist")

genes_words <- merge(words_shared, genes_list, by=c("ID2", "ID1"))
space <- read.csv("location_culture_centroids.csv")
distTable <- data.frame(x=space$Longitude,y=space$Latitude)
names_locs <- c("Aka", "Babongo", "Baka", "Bakola", "Bakoya", "Batwa (East)", "Batwa (West)", "Bedzan", "Efe", "Mbendjele", "Sua" )
rownames(distTable) <- names_locs
space<-as.matrix(distMat(distTable))
rownames(space) <- names_locs
colnames(space) <- names_locs
write.csv(space, "geographical_distances.csv")

# Read in distances DF
library(pscl)
word_dists <- read.csv("shared_words_counts.csv")
word_dists_sub <- subset(word_dists, word_dists$linguistic_distance > 0)

png("PMI_MUSIC_WORDS.png", width = 5, height = 5, units = "in", res = 500)
plot(word_dists_sub$linguistic_distance, word_dists_sub$total_shared, xlab="Linguistic distance (PMI)",
     ylab="No. of shared musical instrument words", pch=16)
dev.off()

png("PMI_MUSIC_WORDS_WITH0.png", width = 5, height = 5, units = "in", res = 500)
plot(word_dists$linguistic_distance, word_dists$total_shared, xlab="Linguistic distance (PMI)",
     ylab="No. of shared musical instrument words", pch=16)
dev.off()

words_music <- subset(word_dists, word_dists$type=="music")
words_subsistence <- subset(word_dists, word_dists$type=="subsistence")


m2 <- zeroinfl(total_shared ~ type + spatial_dist + gen_dist_masked + linguistic_distance, dist = "poisson", data=word_dists)
summary(m2)
m2 <- zeroinfl(total_shared ~ spatial_dist + gen_dist_masked + linguistic_distance, dist = "poisson", data=words_music)
summary(m2)
m5 <- zeroinfl(shared_unique ~ spatial_dist + gen_dist_masked + linguistic_distance, dist = "poisson", data=words_music)
summary(m5)

unwords_complete <- words_music[complete.cases(words_music),]
m2l <- zeroinfl(total_shared ~ linguistic_distance, dist = "poisson", data=words_music)
m2g <- zeroinfl(total_shared ~ gen_dist_masked, dist = "poisson", data=words_music)
m2s <- zeroinfl(total_shared ~ spatial_dist, dist = "poisson", data=words_music)
lmtest::lrtest(m2, m2l)

m3l <- zeroinfl(shared_unique ~ linguistic_distance, dist = "poisson", data=words_music)
m3g <- zeroinfl(shared_unique ~ gen_dist_masked, dist = "poisson", data=words_music)
m3s <- zeroinfl(shared_unique ~ spatial_dist, dist = "poisson", data=words_music)

m2 <- zeroinfl(total_shared ~ spatial_dist + gen_dist_masked + linguistic_distance, dist = "poisson", data=words_subsistence)
summary(m2)


# DO THE SAME WITH PATRISTIC DISTANCES

bantu_music <- subset(words_music, !(words_music$ID2 %in% c("Efe", "Bedzan", "Baka")))
bantu_music <- subset(bantu_music, !(bantu_music$ID1 %in% c("Efe", "Bedzan", "Baka")))

m2 <- zeroinfl(total_shared ~ type + spatial_dist + gen_dist_masked + koile_patristic, dist = "poisson", data=word_dists)
summary(m2)
m3 <- zeroinfl(total_shared ~ spatial_dist + gen_dist_masked + koile_patristic, dist = "poisson", data=words_music)
summary(m3)
m4 <- zeroinfl(total_shared ~ spatial_dist + gen_dist_masked + koile_patristic, dist = "poisson", data=bantu_music)
summary(m4)

m3l <- zeroinfl(shared_unique ~ koile_patristic, dist = "poisson", data=words_music)
m3g <- zeroinfl(shared_unique ~ gen_dist_masked, dist = "poisson", data=words_music)
m3s <- zeroinfl(shared_unique ~ spatial_dist, dist = "poisson", data=words_music)

bantu_complete <- bantu_music[complete.cases(bantu_music),]
words_complete <- words_music[complete.cases(words_music),]
m2l <- zeroinfl(total_shared ~ koile_patristic, dist = "poisson", data=words_music)
m2g <- zeroinfl(total_shared ~ gen_dist_masked, dist = "poisson", data=words_music)
m2s <- zeroinfl(total_shared ~ spatial_dist, dist = "poisson", data=words_music)
lmtest::lrtest(m2, m2l)

cor.test(words_music$total_shared, words_music$koile_patristic, method="pearson")
cor.test(words_complete$total_shared, words_complete$linguistic_distance, method="pearson")
cor.test(words_complete$total_shared, words_complete$gen_dist_masked, method="pearson")

words_non0 <- subset(words_music, words_music$koile_patristic > 0)
cor.test(words_music$koile_patristic, words_music$linguistic_distance)
plot(words_non0$koile_patristic, words_non0$linguistic_distance)
m5 <- lm(koile_patristic ~ linguistic_distance, data=words_music)
summary(m5)
m6 <- lm(koile_patristic ~ linguistic_distance, data=words_non0)
summary(m6)
m7 <- lm(linguistic_distance ~ koile_patristic, data=words_non0)
summary(m7)
m8 <- lm(linguistic_distance ~ koile_patristic, data=words_music)
summary(m8)

scatter_plot <- ggplot(words_non0, aes(y = linguistic_distance, x= koile_patristic)) +
  geom_point() +
  ylab("PMI") +
  xlab("Patristic distance")

# Add the regression line to the scatter plot
scatter_plot_with_line <- scatter_plot +
  geom_abline(intercept = coef(m7)[1], slope = coef(m7)[2], color = "turquoise", lwd=1.2)

# Display the scatter plot with the regression line
print(scatter_plot_with_line)




