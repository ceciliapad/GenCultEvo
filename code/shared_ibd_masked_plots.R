###################################################
## Script for Shared IBD segments on masked data ##
###################################################

# Processed datasets not available for public distribution
# Contact data administrators for acess (see Patin et al. 2017 and Jarvis et al. 2012)
# Code to build the figures below

ped=read.table("admix_combination_masked_bin.ped")
ped_snps <- ped[,-c(1:6)]
#ped_snps[ped_snps!=0]<-1
#ped_snps<-as.data.frame(lapply(ped_snps,as.numeric))
ind_ancestries <- as.data.frame(rowSums(ped_snps))
ped_ind <- ped[,1:6]

masked_per_ind <- cbind(ped_ind, ind_ancestries)
write.csv(masked_per_ind, "masked_per_ind.csv")

# Average mask per group
require(dplyr)
masked_by_fam <- masked_per_ind %>% group_by(V1) %>% 
  summarise(mean_hg_ancestry = mean(`rowSums(ped_snps)`))

masked_by_fam$mean_prop_hg <- masked_by_fam$mean_hg_ancestry/1111260

# Shared IBD segments
### R
ibd<-read.table("all.refinedIBD.filled", as.is=T)
colnames(ibd)<-c("firstID","firstHapIndex","secondID","secondHapIndex", "chromosome", "start","end","LOD","length")
infocomplete<-read.csv("pop_desc.csv", header=T, as.is=T , comment.char = "", fill=T) #info file each line one population
infoID<-read.table("../admix_combination_masked_filtered.fam",header=F, as.is=T)  #info file each line one individual

minimuminfo<-infoID[,c(2,1)]   # select the columns which have the sample ID name and the corresponding population
colnames(minimuminfo)[1]<-"firstID"
ibdmatch<-merge(ibd,minimuminfo,all.x=TRUE)   # associate the population source for the first sample ID of the couple
colnames(ibdmatch)[10]<-"source1"
colnames(minimuminfo)[1]<-"secondID"
ibdmatch2<-merge(ibdmatch,minimuminfo,all.x=TRUE) # associate the population source for the second sample ID of the couple
colnames(ibdmatch2)[11]<-"source2"

write.table(ibdmatch2, "refinedIBD_merged_withinfo.txt", row.names = F, quote=F)
ibd<-ibdmatch2

### table sharing per population

poporder<-unique(infoID$V1)
pops<-table(infoID$V1)

perpop<-matrix(NA,length(poporder),11)
colnames(perpop)<-c("population","samplesize","numberSharingTot","numbersharingWithin","numberSharingOut","FreqSharingTot","FreqsharingWithin","FreqSharingOut","Mean_lengthsharingWithin","totallenghtsharing","howmanypops")
perpop[,1]<-poporder
perpop[,2]<-pops[poporder]

for (i in 1:nrow(perpop)){
  popp<-poporder[i]
  within<-which(ibd$source1%in%popp & ibd$source2%in%popp)
  tempWithin<-ibd[within,]
  tempTOT<-ibd[union(which(ibd$source1%in%popp),which(ibd$source2%in%popp)),]
  tempOut<-tempTOT[-which(tempTOT$source1 == tempTOT$source2),]
  perpop[i,3]<-nrow(tempTOT)
  perpop[i,4]<-nrow(tempWithin)
  perpop[i,5]<-nrow(tempOut)
  perpop[i,6]<-nrow(tempTOT)/as.numeric(perpop[i,2])
  perpop[i,7]<-nrow(tempWithin)/as.numeric(perpop[i,2])
  perpop[i,8]<-nrow(tempOut)/as.numeric(perpop[i,2])
  perpop[i,9]<-mean(tempTOT$length)
  popvarie<-c(tempTOT$source1,tempTOT$source2)
  perpop[i,10]<-sum(tempTOT$length)
  perpop[i,11]<-length(unique(popvarie))
}
#perpop2<-merge(perpop, infoID[,c(1,8,11,12)], by.x=perpop[,1], by.y=infoID[,1])

write.table(perpop,"popInfoIBDsharing.txt",sep="\t", row.names = F, quote=F)
perpop<-read.table("popInfoIBDsharing.txt",sep="\t",header=T, as.is=T,comment.char = "", fill=T, quote="")

# matrix with the total number of shared blocks
matrixIBD<-matrix(NA, length(poporder),length(poporder), dimnames=list(poporder, poporder))

for (i in 1:length(poporder)){
  for (k in 1:length(poporder)){
    pop.i=poporder[i]
    pop.k=poporder[k]
    temp<-ibd[union(which(ibd$source1==pop.i),which(ibd$source2==pop.i)),]
    if (pop.i==pop.k){
      matrixIBD[i,k]<- length(which(temp$source1==temp$source2))
    } else {
      tempp<-rbind(temp[which(temp$source1==pop.k),],temp[which(temp$source2==pop.k),])
      matrixIBD[i,k]<-nrow(tempp)
    }
  }
}
write.table(matrixIBD,"matrix_refinedIBD_merge_sharing.txt", sep="\t")

# make a matrix with the average length of blocks
matrixIBDAverageLength<-matrix(NA, length(poporder),length(poporder), dimnames=list(poporder, poporder))

for (i in 1:length(poporder)){
  for (k in 1:length(poporder)){
    pop.i=poporder[i]
    pop.k=poporder[k]
    temp<-ibd[union(which(ibd$source1==pop.i),which(ibd$source2==pop.i)),]
    if (pop.i==pop.k){
      temp2<-temp[(which(temp$source1==temp$source2)),]
      matrixIBDAverageLength[i,k]<- mean(temp2$length)
    } else {
      tempp<-rbind(temp[which(temp$source1==pop.k),],temp[which(temp$source2==pop.k),])
      matrixIBDAverageLength[i,k]<-mean(tempp$length)
    }
  }
}
write.table(matrixIBDAverageLength,"matrix_IBD_averageLength.txt", sep="\t")


# make a matrix with the TOTAL length of blocks
matrixIBDTotLength<-matrix(NA, length(poporder),length(poporder), dimnames=list(poporder, poporder))

for (i in 1:length(poporder)){
  for (k in 1:length(poporder)){
    pop.i=poporder[i]
    pop.k=poporder[k]
    temp<-ibd[union(which(ibd$source1==pop.i),which(ibd$source2==pop.i)),]
    if (pop.i==pop.k){
      temp2<-temp[(which(temp$source1==temp$source2)),]
      matrixIBDTotLength[i,k]<- sum(temp2$length)
    } else {
      tempp<-rbind(temp[which(temp$source1==pop.k),],temp[which(temp$source2==pop.k),])
      matrixIBDTotLength[i,k]<-sum(tempp$length)
    }
  }
}
write.table(matrixIBDTotLength,"matrix_IBD_totalLength.txt", sep="\t")

pops<-table(infoID$V1)
#pops<-pops[which(pops>0)]
pops<-pops[rownames(matrixIBD)]

#adjust for population size
matrixIBDadjustpopsize<-matrixIBD
for (i in 1:length(poporder)){
  for (k in 1:length(poporder)){
    pop.i=poporder[i]
    pop.k=poporder[k]
    matrixIBDadjustpopsize[i,k]<- matrixIBD[i,k]/(pops[i]*pops[k])
  }
}
write.table(matrixIBDadjustpopsize,"matrix_IBDsharingAdjustPopSize.txt", sep="\t")

######################
# Make plots in maps #
######################

# Plot below uses all IBD segments <1cM. The same code can be used to make plots subsetting
# IBD segments of the desired length (see SM for examples)

library(reshape)
library(reshape)
library(ggplot2)
melted_df<-melt(matrixIBD)
colnames(melted_df)<-c("source1", "source2", "n_sharing")

melted_dfaverage<-melt(matrixIBDAverageLength)
melted_df$averageLength<-melted_dfaverage$value
melted_dflength<-melt(matrixIBDTotLength)
melted_df$totalLength<-melted_dflength$value
melted_dfadjuxt<-melt(matrixIBDadjustpopsize)
melted_df$sharingadjust<-melted_dfadjuxt$value
melted_dflengthadjuxt<-melt(matrixIBDadjustpopsizelength)
melted_df$lengthadjust<-melted_dflengthadjuxt$value


melted_df2<-melted_df[-(which(melted_df$source1==melted_df$source2)),] #exclude same pop sharing

infoID<-infoID[which(infoID$V1%in%poporder),]
pops<-table(infoID$V1)
pops<-pops[poporder]
infopops <- infocomplete[,c(1,3,4)]
#colnames(infopops) <- infopops[1,]
#infopops <- infopops[-1,]
infoID <- infoID[,c(1,2)]
colnames(infoID)[1] <- "FID"
infoID <- merge(infoID, infopops, by="FID", all.x = T)

library(fields)
lista<-(cbind(as.numeric(infoID$longitude_1), as.numeric(infoID$latitude_1)))  #with longitude and latitude coordinates
MatrixGeo<-rdist.earth (lista, miles=FALSE)  #matrix of distance in km between locations
rownames(MatrixGeo)<-infoID$FID
colnames(MatrixGeo)<-infoID$FID

geomelt<-melt(MatrixGeo)

melted_df3<-melted_df2[which(melted_df2$source1%in%geomelt$X1),]
melted_df3<-melted_df3[which(melted_df3$source2%in%geomelt$X1),]
melted_df3<-melted_df3[which(melted_df3$source2!=melted_df3$source1),]

geomelt<-geomelt[which(geomelt$X1!=geomelt$X2),]
geomelt=geomelt[!duplicated(geomelt[c(1,2)]),]

melted_df3$geodist<-as.character(geomelt$value)
melted_df3<-melted_df3[which(melted_df3$n_sharing!=0),]

melted_df3$lengthGeo<-melted_df3$lengthadjust*as.numeric(melted_df3$geodist)
melted_df3$sharingGeo<-melted_df3$sharingadjust*as.numeric(melted_df3$geodist)

library(maps)
library('geosphere')

cahg <- c("Baka", "Babongo",  "Batwa", "Mbuti", "Bakoya", "Bakola",
          "Aka", "Bedzan")

## ALL LENGTHS > 1

png ("map_ibd_all_length.png", units="in", width=6, height=4, res=500) 
par (mar=c(3,4,3,0))
#par (mar=c(0,4,0,0))
plot(g, col=rev(hcl.colors(100, palette="Earth")), legend=F, main="Shared IBD segments >1cM", box=F, axes=F)
map(add=T)
points(x=infoID$longitude_1, y=infoID$latitude_1, pch=19,  cex=0.5, col=MetPalettes$Isfahan1[[1]][7])
for(i in 1:nrow(melted_df3))  {
  node1 <- infoID[infoID$FID == as.character(melted_df3[i,]$source1),]
  node2 <- infoID[infoID$FID == as.character(melted_df3[i,]$source2),]
  
  arc <- gcIntermediate(as.numeric(c(node1[1,]$longitude_1, node1[1,]$latitude_1)), 
                        as.numeric(c(node2[1,]$longitude_1, node2[1,]$latitude_1)), 
                        n=1, addStartEnd=TRUE )
  edge.ind <- round(round((20*melted_df3[i,]$sharingadjust) / max(melted_df3$sharingadjust)))
  
  lines(arc, col=edge.col[edge.ind], lwd=edge.ind)
}
text(x=as.numeric(infoID$longitude_1), y=as.numeric(infoID$latitude_1), labels=infoID_more_10$FID,  cex=0.5, col="black")
dev.off()

