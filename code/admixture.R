#################################################
## ADMIXTURE PLOT PATIN JARVIS DATA (GROUPED) ##
################################################

# Code to make plots with results from ADMIXTURE software
# Example provided for plot at K=3. The same code can be
# Used to plot results at any K value.

setwd("ADMIXTURE_REDUCED")
my_dir <- "ADMIXTURE_REDUCED" # complete with directories
directory <- my_dir


#results <- plot.admixture(my_dir)

west_hg <- c("Babongo", "Baka", "Bedzan",  "Bakoya", "Bakola",
             "Aka")
east_hg <- c("Batwa", "Mbuti")
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

