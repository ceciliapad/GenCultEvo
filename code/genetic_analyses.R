## This will produce genetic descriptives we need for Gene-Culture coevolution paper

pops <- read.csv("pop_desc.csv")
pops[pops=="Biaka"] <- "Aka"
pops[pops=="Mbuti"] <- "Efe & Sua"
pops[pops=="Bezan"] <- "Bedzan"
pops[pops=="BakaG"] <- "Baka"
pops[pops=="BabongoS"] <- "Babongo"
pops[pops=="BabongoE"] <- "Babongo"
coord_pop<- data.frame(long=pops$Longitude, lat=pops$Latitude)

# Change the data set to match 

# read .fam file from prunned dataset (needed for PCA & FST)
fam_prunned <- read.table("Patin_Jarvis_cleaned_KING_pruned.fam")

# change names to match groups
fam_prunned[fam_prunned=="BakaG"] <- "Baka"
fam_prunned[fam_prunned=="Bezan"] <- "Bedzan"
fam_prunned[fam_prunned=="BabongoS"] <- "Babongo"
fam_prunned[fam_prunned=="BabongoE"] <- "Babongo"

# write fam file again
write.table(fam_prunned, "Patin_Jarvis_cleaned_KING_pruned2.fam", sep=" ", row.names = F,  col.names = F,quote=F)

# write list of cahg to keep in filtering process
cahg <- c("Biaka", "Baka", "Bakola", "Bakoya", "Babongo", "Bedzan", "Mbuti",
          "Batwa")

# write file to subset
fam_prunned_cahg <- subset(fam_prunned, fam_prunned$V1 %in% cahg)
write.table(fam_prunned_cahg[,1:2], "cahg_keep.txt", row.names = F,  sep=" ",col.names = F,quote=F)

# read .fam file from non prunned dataset (needed for HOMOZYG)
fam_non <- read.table("Patin_Jarvis_cleaned_KING.fam")

# change names to match groups
fam_non[fam_non=="BakaG"] <- "Baka"
fam_non[fam_non=="Bezan"] <- "Bedzan"
fam_non[fam_non=="BabongoS"] <- "Babongo"
fam_non[fam_non=="BabongoE"] <- "Babongo"

# write fam file again
write.table(fam_non, "Patin_Jarvis_cleaned_KING_grouped.fam", sep=" ", row.names = F,  col.names = F,quote=F)

# write list of cahg to keep in filtering process
cahg <- c("Biaka", "Baka", "Bakola", "Bakoya", "Babongo", "Bedzan", "Mbuti",
          "Batwa")

# Plot genetic populations within areas sampled for material culture

# load shapefile
setwd("/Users/Cecilia/Documents/PhD/African_Pygmies/Spatial analysis culture/My_data/shapefile_cultures")
pygmies <- rgdal::readOGR ("cultural.groupsfixed.shp")
setwd("/Users/Cecilia/Documents/PhD/African_Pygmies/Spatial analysis culture/My_data/shapefile_borders")
borders <- rgdal::readOGR ("borders.congo.shp")

setwd("/Users/Cecilia/Documents/PhD/African_Pygmies/Spatial analysis culture/My_data/Nov_2021")
# Load natural earth data
library(raster)
g <- list.files(pattern="NE1")
g <- raster(g)

# extent of cultures map
my_extent <- extent(8,31,-7,7) #extent of the pygmies shapefile
g <- crop(g, my_extent)

library(reshape2)
library(geosphere)
library(wesanderson)

my_pal <- c(MetPalettes$Isfahan1[[1]], MetPalettes$Isfahan2[[1]][2])

tiff ("map_culture_genes.tiff", units="in", width=6, height=4, res=500) 
par (mar=c(3,4,3,0))
#par (mar=c(0,4,0,0))
plot(g, col=rev(hcl.colors(100, palette="Earth")), legend=F, box=F, axes=F)
plot(borders, add=T)
plot(pygmies, add=T, col = MetPalettes$Isfahan2[[1]][1])
points(coord_pop, pch=17, col=MetPalettes$Isfahan1[[1]][7], cex=1, lwd=1.5)
text(coord_pop, labels = pops$Population, pos =3, offset = 0.6, cex =0.7)
dev.off()

# All the results including PCA, HOM and FST (Lane's program) are here - only comprising CAHG:
# For results with full dataset go to "/Users/Cecilia/Documents/PhD/African_Pygmies/Genetic_studies_pygmies/Final_Data/New genetics march 2021/SNP chip data/Patin_Jarvis"
# see folders there

setwd("/Users/Cecilia/Documents/PhD/African_Pygmies/Genetic_studies_pygmies/Final_Data/New genetics march 2021/SNP chip data/Patin_Jarvis/ibdne_trial/results_summary")
het=as.matrix(read.table("homozygosity_CAHG_PJ.het",head=T))
het_df <- as.data.frame(het)

# Absolute Heterozygosity (observed/total genotypes)
hetf=1-as.numeric(het[,3])/as.numeric(het[,5])
het_df$H_abs <- hetf

hist(hetf,xlab="Heterozygosity") # Do it by population

library(dplyr)
het_pops <- het_df

het_pops <- het_pops %>%                                    # Specify data frame
  group_by(FID) %>%                             # Specify group indicator
  summarise_at(vars(H_abs),                     # Specify column
               list(name = mean))               # Specify function

write.csv(het_pops, "inputs_paper/het_obs_total.csv")


# Now observed - expected
het_df$O.HOM. <- as.numeric(het_df$O.HOM.)
het_df$E.HOM. <- as.numeric(het_df$E.HOM.)
het_df$H_OE <- ((het_df$O.HOM.)/(het_df$E.HOM.))

het_pops_2 <- het_df

het_pops_2 <- het_pops_2 %>%                    # Specify data frame
  group_by(FID) %>%                             # Specify group indicator
  summarise_at(vars(H_OE),                     # Specify column
               list(name = mean))               # Specify function

write.csv(het_pops_2, "inputs_paper/het_obs_exp.csv")

het=as.matrix(read.table("merged_data_flipped_clean.het",head=T))
het_df <- as.data.frame(het)
# Absolute Heterozygosity (observed/total genotypes)
hetf=1-as.numeric(het[,3])/as.numeric(het[,5])
het_df$H_abs <- hetf


neighbors <- c("Badwee", "Nzime", "Akele", "Benga", "Duma", "Eshira", "Eviya", "Fang", "Galoa", "Makina", "Ndumu", "Obamba","Shake", "Bakiga")
east <- c("Bakiga")
cameroon <- c("Badwee", "Nzime")
nb_het <- subset(het_df, het_df$FID %in% neighbors)
nb_het$H_abs  <- round(nb_het$H_abs,3)

pal_nb <- c(MetPalettes$Isfahan2[[1]][3], MetPalettes$Isfahan2[[1]][5], MetPalettes$Isfahan2[[1]][4])
#pal <- wes_palette("Darjeeling1", 11, type = "continuous")
#pal <- hcl.colors(11, "Geyser", rev=TRUE)

nb_het$meta <- "Gabon"
for (i in 1:nrow(nb_het)) {
  if (nb_het[i,1] %in% east) {
    nb_het[i,8] <- "East"
  } else if (nb_het[i,1] %in% cameroon) {
    nb_het[i,8] <- "Cameroon"
  }
}

nb_het_mean <- nb_het  %>%
  group_by(FID) %>%
  summarise(mean=mean(H_abs), sd=sd(H_abs))

nb_het_mean$meta <- "Gabon"
for (i in 1:nrow(nb_het_mean)) {
  if (nb_het_mean[i,1] %in% east) {
    nb_het_mean[i,4] <- "East"
  } else if (nb_het_mean[i,1] %in% cameroon) {
    nb_het_mean[i,4] <- "Cameroon"
  }
}


p <- ggplot(nb_het_mean, aes(x=FID, y=mean, colour=meta)) + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, size=.7) +
  geom_point(cex=2) +
  theme_classic() + 
  ylab("Heterozygosity") +
  xlab("") + 
  scale_colour_manual(values=pal_nb) + 
  theme(legend.position="none")

p <- p + theme(legend.position = "none", text = element_text(size = 14))
ggsave(p, filename="het_masked_new_nb.tiff", height=120, width=220, units="mm")



# Now FST
fst <- read.csv("output_CAHG/pairwise_fst.csv", row.names = 1)
fst <- as.matrix(fst)
pygmy_names <- row.names(fst)
# edit name directly in .csv file
#pygmy_names[pygmy_names=="Biaka"] <- "Aka"
#fst[fst=="Biaka"] <- "Aka"

colnames(fst) <- rownames(fst)
dim <- ncol(fst)

heatmap(fst, scale = "none", col=hcl.colors(10, "Geyser", rev=TRUE))
text(expand.grid(1:dim, 1:dim), sprintf("%0.8f", fst), cex=0.6)
notes_fst <- round(fst,3)

library("gplots")
png(file="fst_3dp.png",
    width=7, height=7, units="in", res=500)
heatmap.2(fst, # to build the heat map you need all digits as otherwise you'll get blank cells
          cellnote = notes_fst,  # same data set for cell labels
          main = "Fst", # heat map title
          notecol="black", 
          key= FALSE,# change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          margins=c(7,7),
          trace="none",  
          dendrogram="none",
          Rowv=FALSE,
          Colv=FALSE,# turns off trace lines inside the heat map
          col=hcl.colors(12, "Geyser", rev=TRUE)) # use on color palette defined earlier
dev.off()   

## PCA with groups

pca1 <- read.table("pca_Patin_Jarvis_CAHG_groups.eigenval",sep=" ",header=F)
pca2 <- read.table("pca_Patin_Jarvis_CAHG_groups.eigenvec",sep=" ",header=F)

#pca1 <- read.table("/Users/Cecilia/Documents/PhD/African_Pygmies/Genetic_studies_pygmies/Final_Data/New genetics march 2021/SNP chip data/Patin_Jarvis/pca_groups.eigenval",sep=" ",header=F)
#pca2 <- read.table("/Users/Cecilia/Documents/PhD/African_Pygmies/Genetic_studies_pygmies/Final_Data/New genetics march 2021/SNP chip data/Patin_Jarvis/pca_groups.eigenvec",sep=" ",header=F)

pca2[pca2=="Biaka"]<-"Aka"
population=as.data.frame(table(pca2$V1))$Var1
sample.id = pca2$V1

tab <- data.frame(sample.id = pca2$V1,
                  EV1 = pca2$V3,    # the first eigenvector
                  EV2 = pca2$V4,    # the second eigenvector
                  stringsAsFactors = FALSE)
head(tab)

#install.packages('ggfortify')
library(ggfortify)
library(wesanderson)

p <- ggplot(data = pca2, aes(x = V3, y = V4, color=V1)) +
  geom_point(alpha = 0.8, size = 1.5) +
  xlab("PCA1 (20.21%)") +
  ylab("PCA2 (7.36%)") + 
  theme_classic() + 
  scale_color_manual(values = my_pal,
                     labels = c("Aka (20)",
                                "Babongo (65)",
                                "Baka (229)",
                                "Bakola (29)",
                                "Bakoya (25)",
                                "Batwa (176)",
                                "Bedzan (60)",
                                "Mbuti (15)")) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) 

p <- p + theme(legend.title = element_blank(),
               legend.text = element_text(size = 12))

ggsave(p, filename="inputs_paper/PCA_groups.png", height=150, width=180, units="mm")

# Runs of homozygosity
# sum_sample <- read.csv("Summary_all.csv")
# subsample <- subset(sum_sample, summ_sample$Subsample==1)
# rn_subssample <- subsample$FID

roh <- read.table("merged_data_flipped_clean.hom",sep="", header = T)
roh_ind <- read.table("merged_data_flipped_clean.hom.indiv",sep="", header = T)

roh_pops_2 <- roh_ind <-  as.data.frame(roh_ind)
pygmies_roh <- roh_ind

mean(pygmies_roh$KB)

roh_pops_2 <- roh_pops_2 %>%                 # Specify data frame
  group_by(FID) %>%                          # Population
  summarise_at(vars(KB),                     # Total Genome in Runs of homozygosity
               list(mean_hom = mean,
                    sd_hom = sd,
                    median_hom= median))            # Specify function
write.csv(roh_pops_2, "runs_homo.csv")


cahg <- c("Biaka", "Baka", "Bakola", "Bakoya", "Babongo", "Bedzan", "Mbuti",
          "Batwa")

pygmies_roh <- subset(pygmies_roh, pygmies_roh$FID %in% cahg)

pygmies_roh$KB  <- as.integer(round(pygmies_roh$KB,0))
# Violin plot for pygmies
pygmies_roh$FID[pygmies_roh$FID=="Biaka"]<-"Aka"
pygmies_roh$FID[pygmies_roh$FID=="Batwa"]<-"Batwa (East)"
pygmies_roh$FID[pygmies_roh$FID=="Mbuti"]<-"Efe & Sua"

library(MetBrewer)

pal <- c(MetPalettes$Isfahan1[[1]], MetPalettes$Isfahan2[[1]][2])
#pal <- wes_palette("Darjeeling1", 11, type = "continuous")
#pal <- hcl.colors(11, "Geyser", rev=TRUE)

p <- ggplot(pygmies_roh, aes(x=FID, y=KB, fill=FID)) + 
  ggtitle("CAHG") + 
  geom_violin(trim=TRUE) + 
  xlab("") +
  ylab("Total Genome in RoH (KB)") + 
  geom_boxplot(width=0.1, outlier.shape = NA, fill="white")+
  theme_classic() 
p <- p + scale_fill_manual(values=pal) + theme(legend.position = "none", text = element_text(size = 12))
ggsave(p, filename="runs_homo_pygmies_new.png", height=150, width=250, units="mm")

write.csv(pygmies_roh, "pygmies_roh_new.csv")


other_hg <- c("Hadza", "Ogiek", "Wata", "Boni", "Sandawe")
all_hg <- c(cahg, other_hg)
all_hg_roh <- subset(roh_ind, roh_ind$FID %in% all_hg)
all_hg_roh$KB  <- as.integer(round(all_hg_roh$KB,0))

pal_hg <- c(MetPalettes$Isfahan1[[1]][3], MetPalettes$Isfahan1[[1]][4])
#pal <- wes_palette("Darjeeling1", 11, type = "continuous")
#pal <- hcl.colors(11, "Geyser", rev=TRUE)

all_hg_roh$meta <- "CAHG"
for (i in 1:nrow(all_hg_roh)) {
  if (all_hg_roh[i,1] %in% other_hg)
  all_hg_roh[i,7] <- "Other HG"
}

all_hg_roh$FID[all_hg_roh$FID=="Biaka"]<-"Aka"
all_hg_roh$FID[all_hg_roh$FID=="Batwa"]<-"Batwa (East)"
all_hg_roh$FID[all_hg_roh$FID=="Mbuti"]<-"Efe & Sua"


p <- ggplot(all_hg_roh, aes(x=FID, y=KB, fill=meta)) + 
  ggtitle("African HG") + 
  geom_violin(trim=TRUE) + 
  xlab("") +
  ylab("Total Genome in RoH (KB)") + 
  geom_boxplot(width=0.1, outlier.shape = NA, fill="white")+
  theme_classic() 
p <- p + scale_fill_manual(values=pal_hg) + theme(legend.position = "none", text = element_text(size = 12))
ggsave(p, filename="runs_homo_all_hg_new.png", height=150, width=300, units="mm")


neighbors <- c("Badwee", "Nzime", "Akele", "Benga", "Duma", "Eshira", "Eviya", "Fang", "Galoa", "Makina", "Ndumu", "Obamba","Shake", "Bakiga")
east <- c("Bakiga")
cameroon <- c("Badwee", "Nzime")
nb_roh <- subset(roh_ind, roh_ind$FID %in% neighbors)
nb_roh$KB  <- as.integer(round(nb_roh$KB,0))

pal_nb <- c(MetPalettes$Isfahan2[[1]][3], MetPalettes$Isfahan2[[1]][5], MetPalettes$Isfahan2[[1]][4])
#pal <- wes_palette("Darjeeling1", 11, type = "continuous")
#pal <- hcl.colors(11, "Geyser", rev=TRUE)

nb_roh$meta <- "Gabon"
for (i in 1:nrow(nb_roh)) {
  if (nb_roh[i,1] %in% east) {
    nb_roh[i,7] <- "East"
  } else if (nb_roh[i,1] %in% cameroon) {
    nb_roh[i,7] <- "Cameroon"
  }
}

p <- ggplot(nb_roh, aes(x=FID, y=KB, fill=meta)) + 
  ggtitle("Neighbours") + 
  geom_violin(trim=TRUE) + 
  xlab("") +
  ylim(0, 400000) +
  ylab("Total Genome in RoH (KB)") + 
  geom_boxplot(width=0.1, outlier.shape = NA, fill="white")+
  theme_classic() 
p <- p + scale_fill_manual(values=pal_nb) + theme(legend.position = "none", text = element_text(size = 12))
ggsave(p, filename="runs_homo_all_nb.png", height=150, width=300, units="mm")

write.csv(pygmies_roh, "pygmies_all_hg_new.csv")


