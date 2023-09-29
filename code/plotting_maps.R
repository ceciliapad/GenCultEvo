## SCRIPT FOR PLOTTING NETWORKS OF JACCARD DISTANCES ON MAP OF CENTRAL AFRICA


#######################################
# Plot networks on a map

library(raster)
library(rgdal)
library(maps)

# load shapefile
cahg <- rgdal::readOGR ("cultural.groupsfixed.shp")
borders <- rgdal::readOGR ("borders.congo.shp")

# Load natural earth data
g <- list.files(pattern="NE1_HR")
g <- raster(g)

# extent of cultures map
my_extent <- extent(8,31,-7.5,9) #extent of the cahg shapefile
g <- crop(g, my_extent)

# Start with musical instruments

dMus <- as.matrix(jacMusic)
Mus_d <- 1/dMus

locations$names <- rownames(locations)
locations$names[2] <- "Babongo"

library(reshape2)
library(geosphere)
library(wesanderson)
df_mus <- melt(Mus_d, varnames = c("Cult_1", "Cult_2"))
df_mus <- subset(df_mus, df_mus$Cult_1 != df_mus$Cult_2)

col.1 <- adjustcolor(MetPalettes$Isfahan1[[1]][6], alpha=0.4)
col.2 <-  adjustcolor(MetPalettes$Isfahan1[[1]][4], alpha=0.4)

#col.1 <- 
#col.2 <- 
edge.pal <- colorRampPalette(c(col.1, col.2), alpha = TRUE) 
edge.col <- edge.pal(100)


png ("results/map_music.png", units="in", width=6, height=4, res=500) 
par (mar=c(3,4,3,0))
#par (mar=c(0,4,0,0))
plot(g, col=rev(hcl.colors(100, palette="Earth")), legend=F, main="Musical instruments", box=F, axes=F)
plot(borders, add=T)
plot(cahg, add=T, col = MetPalettes$Isfahan1[[1]][4])
points(locations$Longitude, locations$Latitude, cex=0.5, pch=19, col=MetPalettes$Isfahan1[[1]][7])
for(i in 1:nrow(df_mus)) {
  node1 <- locations[locations$names == df_mus[i,]$Cult_1,] 
  node2 <- locations[locations$names == df_mus[i,]$Cult_2,]
  arc <- geosphere::gcIntermediate( c(node1[1,]$Longitude, node1[1,]$Latitude), 
                                    c(node2[1,]$Longitude, node2[1,]$Latitude),
                                    n=0, addStartEnd=T)
  edge.ind <- df_mus[i,]$value^2
  if (edge.ind > 1) {
    lines(arc, col=edge.col[edge.ind], lwd=edge.ind) 
  }
}
dev.off()

# Now subsistence tools

dSub <- as.matrix(jacSubsistence)
Sub_d <- 1/dSub

df_sub <- melt(Sub_d, varnames = c("Cult_1", "Cult_2"))
df_sub <- subset(df_sub, df_sub$Cult_1 != df_sub$Cult_2)

col.1 <- adjustcolor(MetPalettes$Morgenstern[[1]][7], alpha=0.4)
col.2 <- adjustcolor(MetPalettes$Morgenstern[[1]][6], alpha=0.4)
edge.pal <- colorRampPalette(c(col.1, col.2), alpha = TRUE) 
edge.col <- edge.pal(100)

locations$names <- rownames(locations)
locations$names[2] <- "Babongo"


png ("results/map_subs.png", units="in", width=6, height=4, res=500) 
par (mar=c(3,4,3,0))
#par (mar=c(0,4,0,0))
plot(g, col=rev(hcl.colors(100, palette="Earth")), legend=F, main="Subsistence tools", box=F, axes=F)
plot(borders, add=T)
plot(cahg, add=T, col = MetPalettes$Isfahan2[[1]][1])
points(locations$Longitude, locations$Latitude, cex=0.5, pch=19, col=MetPalettes$Isfahan1[[1]][2])
for(i in 1:nrow(df_sub)) {
  node1 <- locations[locations$names == df_sub[i,]$Cult_1,] 
  node2 <- locations[locations$names == df_sub[i,]$Cult_2,]
  arc <- geosphere::gcIntermediate( c(node1[1,]$Longitude, node1[1,]$Latitude), 
                                    c(node2[1,]$Longitude, node2[1,]$Latitude),
                                    n=0, addStartEnd=T)
  edge.ind <- df_sub[i,]$value^2
  if (edge.ind > 1) {
    lines(arc, col=edge.col[edge.ind], lwd=edge.ind) 
  }
}
dev.off()


# Start with musical instruments WORDS

locations <- read.csv("location_culture_centroids.csv", row.names = 1)
locations$FID <- rownames(locations)

dMus <- as.matrix(jacMusicWords)
Mus_d <- 1/dMus

locations$names <- rownames(locations)
locations$names[2] <- "Babongo"

library(reshape2)
library(geosphere)
library(wesanderson)
df_mus <- melt(Mus_d, varnames = c("Cult_1", "Cult_2"))
df_mus <- subset(df_mus, df_mus$Cult_1 != df_mus$Cult_2)

col.1 <- adjustcolor(MetPalettes$Isfahan1[[1]][6], alpha=0.4)
col.2 <-  adjustcolor(MetPalettes$Isfahan1[[1]][4], alpha=0.4)

#col.1 <- 
#col.2 <- 
edge.pal <- colorRampPalette(c(col.1, col.2), alpha = TRUE) 
edge.col <- edge.pal(100)


png ("results/map_music_words.png", units="in", width=6, height=4, res=500) 
par (mar=c(3,4,3,0))
#par (mar=c(0,4,0,0))
plot(g, col=rev(hcl.colors(100, palette="Earth")), legend=F, main="Words for musical instruments", box=F, axes=F)
plot(borders, add=T)
plot(cahg, add=T, col = MetPalettes$Isfahan1[[1]][4])
points(locations$Longitude, locations$Latitude, cex=0.5, pch=19, col=MetPalettes$Isfahan1[[1]][7])
for(i in 1:nrow(df_mus)) {
  node1 <- locations[locations$names == df_mus[i,]$Cult_1,] 
  node2 <- locations[locations$names == df_mus[i,]$Cult_2,]
  arc <- geosphere::gcIntermediate( c(node1[1,]$Longitude, node1[1,]$Latitude), 
                                    c(node2[1,]$Longitude, node2[1,]$Latitude),
                                    n=0, addStartEnd=T)
  edge.ind <- df_mus[i,]$value^10
  if (edge.ind > 1) {
    lines(arc, col=edge.col[edge.ind], lwd=edge.ind) 
  }
}
dev.off()

# Start with musical instruments WORDS

locations <- read.csv("location_culture_centroids.csv", row.names = 1)
locations$FID <- rownames(locations)

locations$names <- rownames(locations)
locations$names[2] <- "Babongo"

dSub <- as.matrix(jacSubsistenceWords)
Sub_d <- 1/dSub

df_sub <- melt(Sub_d, varnames = c("Cult_1", "Cult_2"))
df_sub <- subset(df_sub, df_sub$Cult_1 != df_sub$Cult_2)

col.1 <- adjustcolor(MetPalettes$Morgenstern[[1]][7], alpha=0.4)
col.2 <- adjustcolor(MetPalettes$Morgenstern[[1]][6], alpha=0.4)
edge.pal <- colorRampPalette(c(col.1, col.2), alpha = TRUE) 
edge.col <- edge.pal(100)


png ("results/map_subs_words.png", units="in", width=6, height=4, res=500) 
par (mar=c(3,4,3,0))
#par (mar=c(0,4,0,0))
plot(g, col=rev(hcl.colors(100, palette="Earth")), legend=F, main="Words for subsistence tools", box=F, axes=F)
plot(borders, add=T)
plot(cahg, add=T, col = MetPalettes$Isfahan2[[1]][1])
points(locations$Longitude, locations$Latitude, cex=0.5, pch=19, col=MetPalettes$Isfahan1[[1]][2])
for(i in 1:nrow(df_sub)) {
  node1 <- locations[locations$names == df_sub[i,]$Cult_1,] 
  node2 <- locations[locations$names == df_sub[i,]$Cult_2,]
  arc <- geosphere::gcIntermediate( c(node1[1,]$Longitude, node1[1,]$Latitude), 
                                    c(node2[1,]$Longitude, node2[1,]$Latitude),
                                    n=0, addStartEnd=T)
  edge.ind <- df_sub[i,]$value^10
  if (edge.ind > 1) {
    lines(arc, col=edge.col[edge.ind], lwd=edge.ind) 
  }
}
dev.off()
