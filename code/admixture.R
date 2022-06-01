#################################################
## ADMIXTURE PLOT PATIN JARVIS DATA (GROUPED) ##
################################################

setwd("/Users/Cecilia/ncbi/public/downloads/dbGaP-25678/86588/PhenoGenotypeFiles/RootStudyConsentSet_phs001780.AfricanDemographicHistory.v1.p1.c1.GRU-PUB-NPU/GenotypeFiles/Genotype_files/ADMIXTURE_REDUCED")
my_dir <- "/Users/Cecilia/ncbi/public/downloads/dbGaP-25678/86588/PhenoGenotypeFiles/RootStudyConsentSet_phs001780.AfricanDemographicHistory.v1.p1.c1.GRU-PUB-NPU/GenotypeFiles/Genotype_files/ADMIXTURE_REDUCED"
directory <- my_dir


#results <- plot.admixture(my_dir)

west_hg <- c("Babongo", "Baka", "Bedzan",  "Bakoya", "Bakola",
             "Biaka")
east_hg <- c("Batwa", "Mbuti")

#tbl=read.table("/Users/Cecilia/Documents/PhD/African_Pygmies/Genetic_studies_pygmies/Final_Data/New genetics march 2021/SNP chip data/Patin_Jarvis/Patin_Jarvis_cleaned_KING_pruned2.3.Q")
tbl=read.table("merged_data_reduced.3.Q")
fam<-dir(directory, pattern="merged_data_reduced.fam")
temp.name<-read.csv(paste(directory, fam, sep=""), sep=" ", header=FALSE)
dd <- cbind(temp.name, tbl)
dd$PopGroup <- "AGR"

for (i in 1:nrow(dd)) {
  if (dd[i,1] %in% west_hg) {
    dd[i,10] <- "West_HG"
  } else if (dd[i,1] %in% east_hg) {
    dd[i,10] <- "East_HG"
  }
}

barNaming <- function(vec) {
  retVec <- vec
  for(k in 2:length(vec)) {
    if(vec[k-1] == vec[k])
      retVec[k] <- ""
  }
  return(retVec)
}

bant <- subset(dd, dd$PopGroup=="AGR")
hg <- subset(dd, dd$PopGroup!="AGR")
 
mergedAdmWithPopGroups <- bant

ordered = mergedAdmWithPopGroups[order(mergedAdmWithPopGroups$Pop),]
colnames(ordered) <- c("Pop", "Ind", "x1", "x2","x3","x4","V1","V2","V3", "PopGroup")
barplot(t(as.matrix(subset(ordered, select=V1:V3))), col=rainbow(3), border=NA)
names <- barNaming(ordered$Pop)

cameroon_here <- c("Baka", "Biaka","Badwee", "Nzime")
cam <- subset(bant, bant$V1 %in% cameroon_here)
aggregate(cam[, 8], list(cam$V1), sd)

bantu_gabon <- c("Akele", "Benga", "Duma", "Eshira", "Eviya", "Fang", "Galoa",
                 "Makina", "Ndumu", "Obamba",
                 "Orungu", "Shake")
gab <- subset(bant, bant$V1 %in% bantu_gabon)
aggregate(gab[, 8], list(gab$V1), sd)

eastern_here <- c("Batwa", "Mbuti","Bakiga", "Luhya")
east <- subset(bant, bant$V1 %in% eastern_here)
aggregate(east[, 8], list(east$V1), mean)

require(MetBrewer)
my_pal <- MetPalettes$Isfahan1[[1]][3:5]
            
#png("admix.png",units="in", width=10, height=3, res=500)
barplot(t(as.matrix(ordered[,7:9])), col=my_pal, border=NA,
        names.arg=barNaming(ordered$Pop), las=2,cex.names = 0.8, cex.axis = 0.8)
#dev.off()

png("admix_bantu.png",units="in", width=10, height=3, res=500)
## Plot, but suppress the labels
midpts <- barplot(t(as.matrix(ordered[,7:9])), col=my_pal, border=NA,
                  cex.axis = 0.5, names.arg=rep("", nrow(ordered)))

## Use grid to add the labels    
vps <- baseViewports()
pushViewport(vps$inner, vps$figure, vps$plot)

grid.text(names,
          x = unit(midpts, "native"), y=unit(-1, "lines"),
          just="right", rot=90, gp=gpar(fontsize=6))

popViewport(3)
dev.off()

mergedAdmWithPopGroups <- hg

ordered = mergedAdmWithPopGroups[order(mergedAdmWithPopGroups$Pop),]
colnames(ordered) <- c("Pop", "Ind", "x1", "x2","x3","x4","V1","V2","V3", "PopGroup")
names <- barNaming(ordered$Pop)

#png("admix.png",units="in", width=10, height=3, res=500)
barplot(t(as.matrix(ordered[,7:9])), col=my_pal, border=NA,
        names.arg=barNaming(ordered$Pop), las=2,cex.names = 0.8, cex.axis = 0.8)
#dev.off()

png("admix_hg.png",units="in", width=10, height=3, res=500)
## Plot, but suppress the labels
midpts <- barplot(t(as.matrix(ordered[,7:9])), col=my_pal, border=NA,
                  cex.axis = 0.5, names.arg=rep("", nrow(ordered)))

## Use grid to add the labels    
vps <- baseViewports()
pushViewport(vps$inner, vps$figure, vps$plot)

grid.text(names,
          x = unit(midpts, "native"), y=unit(-1, "lines"),
          just="right", rot=90, gp=gpar(fontsize=6))

popViewport(3)
dev.off()

ordered_mean_admix <- ordered %>% group_by(Pop)
ordered_mean_admix <- ordered_mean_admix %>% summarise(
  CAHG_1 = mean(V1),
  CAHG_2 = mean(V2),
  AGR = mean(V3)
)

write.csv(ordered_mean_admix, "mean_admix_pops.csv")

##############################
## RUN WITH ALL POPULATIONS ##
##############################

setwd("/Users/Cecilia/ncbi/public/downloads/dbGaP-25678/86588/PhenoGenotypeFiles/RootStudyConsentSet_phs001780.AfricanDemographicHistory.v1.p1.c1.GRU-PUB-NPU/GenotypeFiles/Genotype_files")
directory <- getwd()

#tbl=read.table("/Users/Cecilia/Documents/PhD/African_Pygmies/Genetic_studies_pygmies/Final_Data/New genetics march 2021/SNP chip data/Patin_Jarvis/Patin_Jarvis_cleaned_KING_pruned2.3.Q")
tbl=read.table("merged_data_flipped_clean_prunned.5.Q")
temp.name<-read.csv("merged_data_flipped_clean_prunned.fam", sep=" ", header=FALSE)
dd <- cbind(temp.name, tbl)

dd <- subset(dd, dd$V1 != "Twa")

ordered = dd[order(dd$V1),]
colnames(ordered) <- c("Pop", "Ind", "x1", "x2","x3","x4","V1","V2","V3", "V4", "V5")
barplot(t(as.matrix(subset(ordered, select=V1:V5))), col=rainbow(3), border=NA)
names <- barNaming(ordered$Pop)

#png("admix.png",units="in", width=10, height=3, res=500)
barplot(t(as.matrix(ordered[,7:9])), col=my_pal, border=NA,
        names.arg=barNaming(ordered$Pop), las=2,cex.names = 0.8, cex.axis = 0.8)
#dev.off()

png("admix_all_pops_k5.png",units="in", width=15, height=3, res=500)
## Plot, but suppress the labels
pal <- my_pal <- MetPalettes$Isfahan1[[1]]
midpts <- barplot(t(as.matrix(ordered[,7:11])), col=pal, border=NA,
                  cex.axis = 0.5, names.arg=rep("", nrow(ordered)))

## Use grid to add the labels    
vps <- baseViewports()
pushViewport(vps$inner, vps$figure, vps$plot)
grid.text(names,
          x = unit(midpts, "native"), y=unit(-5, "lines"),
          just="right", rot=270, gp=gpar(fontsize=6))

popViewport(3)
dev.off()

#### NOW K=4

#tbl=read.table("/Users/Cecilia/Documents/PhD/African_Pygmies/Genetic_studies_pygmies/Final_Data/New genetics march 2021/SNP chip data/Patin_Jarvis/Patin_Jarvis_cleaned_KING_pruned2.3.Q")
tbl=read.table("merged_data_flipped_clean_prunned.4.Q")
temp.name<-read.csv("merged_data_flipped_clean_prunned.fam", sep=" ", header=FALSE)
dd <- cbind(temp.name, tbl)

dd <- subset(dd, dd$V1 != "Twa")

ordered = dd[order(dd$V1),]
colnames(ordered) <- c("Pop", "Ind", "x1", "x2","x3","x4","V1","V2","V3", "V4")
barplot(t(as.matrix(subset(ordered, select=V1:V4))), col=rainbow(3), border=NA)
names <- barNaming(ordered$Pop)

#png("admix.png",units="in", width=10, height=3, res=500)
barplot(t(as.matrix(ordered[,7:9])), col=my_pal[2:5], border=NA,
        names.arg=barNaming(ordered$Pop), las=2,cex.names = 0.8, cex.axis = 0.8)
#dev.off()

png("admix_all_pops_k4.png",units="in", width=15, height=3, res=500)
## Plot, but suppress the labels
pal <- my_pal <- MetPalettes$Isfahan1[[1]]
midpts <- barplot(t(as.matrix(ordered[,7:10])), col=pal[2:5], border=NA,
                  cex.axis = 0.5, names.arg=rep("", nrow(ordered)))

## Use grid to add the labels    
vps <- baseViewports()
pushViewport(vps$inner, vps$figure, vps$plot)
grid.text(names,
          x = unit(midpts, "native"), y=unit(-5, "lines"),
          just="right", rot=270, gp=gpar(fontsize=6))

popViewport(3)
dev.off()


east <- subset(dd, dd$PopGroup=="East_HG")
west <- subset(dd, dd$PopGroup=="West_HG")

## Select Hunter-Gatherers with >98% HG-related ancestry (Eastern and Western components together)

ref_subset <- subset(hg, hg$V2.1 < 0.03)
ref_subset_bantu <- subset(ordered, ordered$V2 > 0.98)
setwd("/Users/Cecilia/Documents/PhD/African_Pygmies/Genetic_studies_pygmies/Final_Data/New genetics march 2021/SNP chip data/Patin_Jarvis/gnomix_snp/Reference panels")
write.csv(ordered, "all_hg_and_bantu.csv")
write.csv(hg, "all_hg.csv")
write.csv(ref_subset, "all_hg_refs.csv")
write.csv(ref_subset_bantu, "bantu_refs.csv")

all_inds <- read.table("inds_keep_gnomix.txt")
ref_inds <- subset(all_inds, all_inds$V1 %in% ref_subset$V2 |  all_inds$V1 %in% ref_subset_bantu$Ind)
admix_inds <- subset(all_inds, !(all_inds$V1 %in% ref_inds$V1))

write.table(ref_inds, "ref_inds.txt", row.names = F, col.names = F, quote = FALSE)
write.table(admix_inds, "admix_ind.txt", row.names = F, col.names = F, quote = FALSE)

colnames(ref_subset) <- colnames(ref_subset_bantu)
ref_all <- rbind(ref_subset, ref_subset_bantu)
ref_all <- ref_all[,c(2,10)]

for (i in 1:nrow(ref_all)){
  ifelse(ref_all[i,2]=="Bantu", ref_all[i,2] <- "AGR", ref_all[i,2] <- "CAHG")
}

write.table(ref_all, "references.smap", row.names = F, col.names = F, quote = FALSE, sep="\t")
