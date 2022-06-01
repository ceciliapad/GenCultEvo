## SCRIPT FOR PLOTTING NETWORKS OF JACCARD DISTANCES ON MAP OF CENTRAL AFRICA
# Author: Cecilia Padilla-Iglesias

setwd("/Users/Cecilia/Documents/PhD/African_Pygmies/Spatial analysis culture/My_data/Nov_2021")
load("SpatialCulture.RData")

#######################################
# Plot networks on a map

library(raster)
library(rgdal)
library(maps)

# load shapefile
setwd("/Users/Cecilia/Documents/PhD/African_Pygmies/Spatial analysis culture/My_data/shapefile_cultures")
pygmies <- rgdal::readOGR ("cultural.groupsfixed.shp")
setwd("/Users/Cecilia/Documents/PhD/African_Pygmies/Spatial analysis culture/My_data/shapefile_borders")
borders <- rgdal::readOGR ("borders.congo.shp")
setwd("/Users/Cecilia/Documents/PhD/African_Pygmies/Spatial analysis culture/My_data/Nov_2021")

# Load natural earth data
g <- list.files(pattern="NE1_HR")
g <- raster(g)

# extent of cultures map
my_extent <- extent(8,31,-7,7) #extent of the pygmies shapefile
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
plot(pygmies, add=T, col = MetPalettes$Isfahan1[[1]][4])
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
plot(pygmies, add=T, col = MetPalettes$Isfahan2[[1]][1])
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
plot(pygmies, add=T, col = MetPalettes$Isfahan1[[1]][4])
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
plot(pygmies, add=T, col = MetPalettes$Isfahan2[[1]][1])
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



#################################
# Sampling of spatial locations #
#################################
require(sf)
require(rgdal)

polys_sf <- st_as_sf(pygmies)
# Polygons are really small, so let's take 500 samples per polygon
# 9500 samples in total
# First sample 500 random points from each culture polygon
all_samples <- st_sample(polys_sf, size=9500, type="random", by_polygon=T)
nc_df2 <- as.data.frame(all_samples)



# create point sample within polygons
## followed http://casoilresource.lawr.ucdavis.edu/drupal/node/644 

polyPts = sapply(slot(pygmies, 'polygons'), function(i) spsample(i, n=500,type='random', offset=c(0,0)))  # works but doesn't vary sample by
#poly_df <- as.data.frame(polyPts)

polyPts.merged <- do.call('rbind', polyPts) # this works when n above is constant but not when varied
poly_df <- as.data.frame(polyPts.merged)
# add an id, based on source polygon:
ids <- sapply(slot(pygmies, 'polygons'), function(i) slot(i, 'ID')) #

pt_id=c()
for (i in ids) {
  ids_here <- rep(i, 500)
  pt_id <- c(pt_id, ids_here)
}

# promote to SpatialPointsDataFrame
#polyPts.final <- SpatialPointsDataFrame(polyPts.merged,data=data.frame(poly_id=pt_id))
poly_df$poly_id <- pt_id

# rename with actual names of cultural groups
#df_points <- as.data.frame(polyPts.final)
df_points <- as.data.frame(poly_df)
df_points$poly_id[df_points$poly_id=="0"]<-"Bakola"
df_points$poly_id[df_points$poly_id=="1"]<-"Baka"
df_points$poly_id[df_points$poly_id=="2"]<-"Baka"
df_points$poly_id[df_points$poly_id=="3"]<-"Baka"
df_points$poly_id[df_points$poly_id=="4"]<-"Bakoya"
df_points$poly_id[df_points$poly_id=="5"]<-"Aka"
df_points$poly_id[df_points$poly_id=="6"]<-"Mbendjele"
df_points$poly_id[df_points$poly_id=="7"]<-"Batwa (West)"
df_points$poly_id[df_points$poly_id=="8"]<-"Bedzan"
df_points$poly_id[df_points$poly_id=="9"]<-"Babongo"
df_points$poly_id[df_points$poly_id=="10"]<-"Babongo"
df_points$poly_id[df_points$poly_id=="11"]<-"Babongo"
df_points$poly_id[df_points$poly_id=="12"]<-"Babongo"
df_points$poly_id[df_points$poly_id=="13"]<-"Babongo"
df_points$poly_id[df_points$poly_id=="14"]<-"Babongo"
df_points$poly_id[df_points$poly_id=="15"]<-"Babongo"
df_points$poly_id[df_points$poly_id=="16"]<-"Batwa (East)"
df_points$poly_id[df_points$poly_id=="17"]<-"Sua"
df_points$poly_id[df_points$poly_id=="18"]<-"Efe"

# sample 500 points by ID
df_points_new <- df_points %>% group_by(poly_id) %>% slice_sample(n=500)
write.csv(df_points_new, "geo_random_points.csv", row.names = F)


