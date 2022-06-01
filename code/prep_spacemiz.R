# Prepare spacemix input

# If you have use is_weighted option you'll already have freqs adjusted for sample size

setwd("/Users/Cecilia/Documents/PhD/African_Pygmies/Genetic_studies_pygmies/Final_Data/maasMDS/data_modify")
pops <- c("Baka", "Bakoya", "Bakola", "Mbuti", "Batwa", "Biaka", "Bedzan", "Babongo")
i <- "Babongo"
freq <- read.table(paste(i,".freq", sep=""), sep=",", header=1)
loci_list <- as.vector(freq$CHROM_IDX)
j <- 1
for (i in pops) {
  freq <- read.table(paste(i,".freq", sep=""), sep=",", header=1)
  freqs <- freq$FREQ
  loci_list <- as.data.frame(cbind(loci_list, freqs))
  colnames(loci_list)[(j+1)] <- i
  j <- j + 1
}

library(tibble)
colnames(loci_list)[5] <- "Efe"
Sua <- loci_list[,5]
add_column(loci_list, Sua = Sua, .after = "Efe")
loci_list <- add_column(loci_list, Sua = Sua, .after = "Efe")
loci_list_t <- as.data.frame(t(loci_list))
colnames(loci_list_t) <- loci_list_t[1,]
loci_list_t[1,] <- NULL
loci_list_t <- loci_list_t[-1,]
write.csv(loci_list_t, "loci_list.csv", col.names = T, row.names = T)

i <- "Babongo"
freq <- read.table(paste(i,".freq", sep=""), sep=",", header=1)
count_list <- as.vector(freq$CHROM_IDX)
j <- 1
for (i in pops) {
  freq <- read.table(paste(i,".freq", sep=""), sep=",", header=1)
  count <- freq$CT
  count_list <- as.data.frame(cbind(count_list, count))
  colnames(count_list)[(j+1)] <- i
  j <- j + 1
}
library(tibble)
colnames(count_list)[5] <- "Efe"
Sua <- count_list[,5]
count_list <- add_column(loci_list, Sua = Sua, .after = "Efe")
count_list_t <- as.data.frame(t(count_list))
colnames(count_list_t) <- count_list_t[1,]
count_list_t <- count_list_t[-1,]
write.csv(count_list_t, "count_list.csv", row.names = T, col.names = T)

## Plot genetic groups
setwd("/Users/Cecilia/Documents/PhD/African_Pygmies/Genetic_studies_pygmies/Final_Data/New genetics march 2021/SNP chip data/Patin_Jarvis")
pops <- read.csv("pop_desc.csv")
pops[pops=="Biaka"] <- "Aka"
coord_pop<- data.frame(long=pops$Longitude, lat=pops$Latitude)
coord_pop
pops
pops[pops=="Mbuti"] <- "Efe & Sua"
coord_pop<- data.frame(long=pops$Longitude, lat=pops$Latitude)
pops[pops=="Bezan"] <- "Bedzan"
coord_pop<- data.frame(long=pops$Longitude, lat=pops$Latitude)
# load shapefile
setwd("/Users/Cecilia/Documents/PhD/African_Pygmies/Spatial analysis culture/My_data/shapefile_cultures")
pygmies <- rgdal::readOGR ("cultural.groupsfixed.shp")
setwd("/Users/Cecilia/Documents/PhD/African_Pygmies/Spatial analysis culture/My_data/shapefile_borders")
borders <- rgdal::readOGR ("borders.congo.shp")
setwd("/Users/Cecilia/Documents/PhD/African_Pygmies/Spatial analysis culture/My_data/Nov_2021")

my_pal <- c(MetPalettes$Isfahan1[[1]], MetPalettes$Isfahan2[[1]][2])

pops[pops=="BakaG"] <- "Baka"
pops[pops=="BabongoS"] <- "Babongo"
pops[pops=="BabongoE"] <- "Babongo"

tiff ("results_masked/map_culture_genes.tiff", units="in", width=6, height=4, res=500)
par (mar=c(3,4,3,0))
#par (mar=c(0,4,0,0))
plot(g, col=rev(hcl.colors(100, palette="Earth")), legend=F, box=F, axes=F)
plot(borders, add=T)
plot(pygmies, add=T, col = MetPalettes$Isfahan2[[1]][1])
points(coord_pop, pch=17, col=MetPalettes$Isfahan1[[1]][7], cex=1, lwd=1.5)
text(coord_pop, labels = pops$Population, pos =3, offset = 0.6, cex =0.7)
dev.off()

#####################
locations <- read.csv("/Users/Cecilia/Documents/PhD/African_Pygmies/Spatial analysis culture/My_data/Nov_2021/location_culture_centroids2.csv", row.names = 1)
setwd("/Users/Cecilia/Documents/PhD/African_Pygmies/Spatial analysis culture/My_data/Nov_2021/spacemix")

freq_mat <- read.table("loci_list.csv", sep=",")
count_mat <- read.table("count_list.csv", sep=",")
rownames(freq_mat) <- freq_mat[,1]
colnames(freq_mat) <- freq_mat[1,]
freq_mat <- freq_mat[-1,-1]
freq_mat <- as.matrix(freq_mat)
tmp.sample.sizes <- count_mat
sample.colors <- rainbow(n=nrow(locations),start=4/6,end=6/6)[as.numeric(cut(locations[,2],nrow(locations)))]

# Spacemix for genes
require(MetBrewer)
my_pal_long <- c(MetPalettes$Isfahan1[[1]][1:6],MetPalettes$Isfahan2[[1]][1], MetPalettes$Isfahan1[[1]][7], MetPalettes$Isfahan2[[1]][c(2,3,5)])
freq_mat <- loci_list_t
count_mat <- count_list_t
sample_df <- as.data.frame(matrix(nrow=nrow(freq_mat), ncol=ncol(freq_mat)))
freq_mat <- as.matrix(freq_mat)
count_mat <- as.matrix(count_mat)
tmp.sample.sizes <- count_mat
#tmp.sample.sizes[which(ex.sample.sizes == 0)] <- NA
mean.sample.sizes <- rowMeans(tmp.sample.sizes[-1,-1],na.rm=TRUE)
#mean.sample.sizes <- rowMeans(tmp.sample.sizes,na.rm=TRUE)
# Do sample means with actual number of individuals with allele (non-masked) per population
run.spacemix.analysis(n.fast.reps=10, fast.MCMC.ngen=1e5, fast.model.option="target", long.model.option="source_and_target",
                                      data.type="sample.frequencies", sample.frequencies = freq_mat, mean.sample.sizes = mean.sample.sizes,
                                      counts = NULL, sample.sizes = NULL, sample.covariance = NULL, target.spatial.prior.scale = NULL,
                                      source.spatial.prior.scale = NULL, spatial.prior.X.coordinates=locations$Longitude, spatial.prior.Y.coordinates=locations$Latitude, round.earth=FALSE,
                                      long.run.initial.parameters = NULL, k=nrow(freq_mat), loci=ncol(freq_mat), ngen=1e6, printfreq=1e2, samplefreq=1e3,
                                      mixing.diagn.freq = 50, savefreq=1e5, directory = NULL, prefix = "cahg_spacemix_all")


# Spacemix with music
count_music <- read.csv("P_A_music2.csv", row.names = 1)
count_music <- as.matrix(count_music)
total_music <- read.csv("music_sample_size.csv", row.names = 1)
total_music <- as.matrix(total_music)
run.spacemix.analysis(n.fast.reps=10, fast.MCMC.ngen=1e5, fast.model.option="target", long.model.option="source_and_target",
                                    data.type="counts", sample.frequencies = NULL, mean.sample.sizes = NULL,
                                    counts = count_music, sample.sizes = total_music, sample.covariance = NULL, target.spatial.prior.scale = NULL,
                                    source.spatial.prior.scale = NULL, spatial.prior.X.coordinates=locations$Longitude, spatial.prior.Y.coordinates=locations$Latitude, round.earth=FALSE,
                                    long.run.initial.parameters = NULL, k=nrow(count_music), loci=ncol(count_music), ngen=1e6, printfreq=1e2, samplefreq=1e3,
                                    mixing.diagn.freq = 50, savefreq=1e5, directory = NULL, prefix = "cahg_spacemix_music")


# Spacemix with subsistence
count_subsistence <- read.csv("P_A_subsistence2.csv", row.names = 1)
count_subsistence <- as.matrix(count_subsistence)
total_subsistence <- read.csv("subsistence_sample_size.csv", row.names = 1)
total_subsistence <- as.matrix(total_subsistence)
run.spacemix.analysis(n.fast.reps=10, fast.MCMC.ngen=1e5, fast.model.option="target", long.model.option="source_and_target",
                                    data.type="counts", sample.frequencies = NULL, mean.sample.sizes = NULL,
                                    counts = count_subsistence, sample.sizes = total_subsistence, sample.covariance = NULL, target.spatial.prior.scale = NULL,
                                    source.spatial.prior.scale = NULL, spatial.prior.X.coordinates=locations$Longitude, spatial.prior.Y.coordinates=locations$Latitude, round.earth=FALSE,
                                    long.run.initial.parameters = NULL, k=nrow(count_subsistence), loci=ncol(count_subsistence), ngen=1e6, printfreq=1e2, samplefreq=1e3,
                                    mixing.diagn.freq = 50, savefreq=1e5, directory = NULL, prefix = "cahg_spacemix_subsistence")

# Make maps

subsistence.map.list <- make.spacemix.map.list(MCMC.output.file = "/Users/Cecilia/Documents/PhD/African_Pygmies/Spatial analysis culture/My_data/Nov_2021/spacemix/run_611922/cahg_spacemix_subsistence_LongRun/cahg_spacemix_subsistence_space_MCMC_output1.Robj",
                                               geographic.locations = as.matrix(locations[,2:3]), name.vector = rownames(freq_mat),color.vector = sample.colors, quantile=0.95, burnin=0)
music.map.list <- make.spacemix.map.list(MCMC.output.file = "/Users/Cecilia/Documents/PhD/African_Pygmies/Spatial analysis culture/My_data/Nov_2021/spacemix/run_299078/cahg_spacemix_music_LongRun/cahg_spacemix_music_space_MCMC_output1.Robj",
                                              geographic.locations = as.matrix(locations[,2:3]), name.vector = rownames(freq_mat),color.vector = sample.colors, quantile=0.95, burnin=0)

cahg.map.list <- make.spacemix.map.list(MCMC.output.file = "/Users/Cecilia/Documents/PhD/African_Pygmies/Spatial analysis culture/My_data/Nov_2021/spacemix/run_793742/cahg_spacemix_all_LongRun/cahg_spacemix_all_space_MCMC_output1.Robj",
                                        geographic.locations = as.matrix(locations[,2:3]), name.vector = rownames(freq_mat),color.vector = sample.colors, quantile=0.95, burnin=0)

make.spacemix.map(spacemix.map.list = subsistence.map.list,
                                text=TRUE,
                                ellipses=T,
                                source.option=FALSE)
make.spacemix.map(spacemix.map.list = music.map.list,
                                text=TRUE,
                                ellipses=T,
                                source.option=FALSE)

make.spacemix.map(spacemix.map.list = cahg.map.list,
                                text=TRUE,
                                ellipses=T,
                                source.option=FALSE)

# Plot in actual map
# Make map a bit bigger
range_x <- c((range(cahg.map.list$geographic.locations[,1])[1]-1), (range(cahg.map.list$geographic.locations[,1])[2]+1))
range_y <- c((range(cahg.map.list$geographic.locations[,2])[1]-1), (range(cahg.map.list$geographic.locations[,2])[2]+1))

library(rworldmap)
newmap <- getMap(resolution = "high")

png(file="map_music.png", units="in", width=6, height=5, res=500)
plot(newmap, xlim = range_x, ylim = range_y, asp = 1, col="lightgray", main="Musical instruments and genetics")
points(cahg.map.list$geographic.locations[,1], cahg.map.list$geographic.locations[,2],col=locations$Col, pch=15)
points(music.map.list$MAPP.geogen.coords[,1], music.map.list$MAPP.geogen.coords[,2], col=locations$Col, pch=16)
points(cahg.map.list$MAPP.geogen.coords[,1], cahg.map.list$MAPP.geogen.coords[,2], col=locations$Col, pch=17)
segments(cahg.map.list$geographic.locations[,1], cahg.map.list$geographic.locations[,2], music.map.list$MAPP.geogen.coords[,1], music.map.list$MAPP.geogen.coords[,2], col=locations$Col, lwd=1.5)
segments(cahg.map.list$geographic.locations[,1], cahg.map.list$geographic.locations[,2], cahg.map.list$MAPP.geogen.coords[,1], cahg.map.list$MAPP.geogen.coords[,2], col=locations$Col, lwd=1.5)
segments(cahg.map.list$MAPP.geogen.coords[,1], cahg.map.list$MAPP.geogen.coords[,2], music.map.list$MAPP.geogen.coords[,1], music.map.list$MAPP.geogen.coords[,2], col=locations$Col, lwd=1.5)
text(cahg.map.list$geographic.locations[,1], (cahg.map.list$geographic.locations[,2])+0.6, cahg.map.list$name.vector, col=locations$Col, cex=0.8, font=2)
dev.off()

newmap <- getMap(resolution = "high")
png(file="map_subsistence.png", units="in", width=6, height=5, res=500)
plot(newmap, xlim = range_x, ylim = range_y, asp = 1, col="lightgray", main="Subsistence tools and genetics")
points(cahg.map.list$geographic.locations[,1], cahg.map.list$geographic.locations[,2],col=locations$Col, pch=15)
points(subsistence.map.list$MAPP.geogen.coords[,1], subsistence.map.list$MAPP.geogen.coords[,2], col=locations$Col, pch=16)
points(cahg.map.list$MAPP.geogen.coords[,1], cahg.map.list$MAPP.geogen.coords[,2], col=locations$Col, pch=17)
segments(cahg.map.list$geographic.locations[,1], cahg.map.list$geographic.locations[,2], subsistence.map.list$MAPP.geogen.coords[,1], subsistence.map.list$MAPP.geogen.coords[,2], col=locations$Col, lwd=1.5)
segments(cahg.map.list$geographic.locations[,1], cahg.map.list$geographic.locations[,2], cahg.map.list$MAPP.geogen.coords[,1], cahg.map.list$MAPP.geogen.coords[,2], col=locations$Col, lwd=1.5)
segments(cahg.map.list$MAPP.geogen.coords[,1], cahg.map.list$MAPP.geogen.coords[,2], subsistence.map.list$MAPP.geogen.coords[,1], subsistence.map.list$MAPP.geogen.coords[,2], col=locations$Col, lwd=1.5)
text(cahg.map.list$geographic.locations[,1], (cahg.map.list$geographic.locations[,2])+0.6, cahg.map.list$name.vector, col=locations$Col, cex=0.8, font=2)
dev.off()

setwd("/Users/Cecilia/Documents/PhD/African_Pygmies/Spatial analysis culture/My_data/Nov_2021")
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
my_extent <- extent(8,31,-7.5,8) #extent of the pygmies shapefile
g <- crop(g, my_extent)

# Color for Baka not very clear
locations$Col[1] <- "#d7aca1"


tiff(file="map_music_real.tiff", units="in", width=7, height=6, res=500)
plot(g, col=rev(hcl.colors(100, palette="Earth")), legend=F, main="Musical instruments", box=F, axes=F)
map(add=T)
points(cahg.map.list$geographic.locations[,1], cahg.map.list$geographic.locations[,2],col=locations$Col, pch=15, cex=1.5)
points(music.map.list$MAPP.geogen.coords[,1], music.map.list$MAPP.geogen.coords[,2], col=locations$Col, pch=16, cex=1.5)
points(cahg.map.list$MAPP.geogen.coords[,1], cahg.map.list$MAPP.geogen.coords[,2], col=locations$Col, pch=17, cex=1.5)
segments(cahg.map.list$geographic.locations[,1], cahg.map.list$geographic.locations[,2], music.map.list$MAPP.geogen.coords[,1], music.map.list$MAPP.geogen.coords[,2], col=locations$Col, lwd=2)
segments(cahg.map.list$geographic.locations[,1], cahg.map.list$geographic.locations[,2], cahg.map.list$MAPP.geogen.coords[,1], cahg.map.list$MAPP.geogen.coords[,2], col=locations$Col, lwd=2)
segments(cahg.map.list$MAPP.geogen.coords[,1], cahg.map.list$MAPP.geogen.coords[,2], music.map.list$MAPP.geogen.coords[,1], music.map.list$MAPP.geogen.coords[,2], col=locations$Col, lwd=2)
text(cahg.map.list$geographic.locations[,1], (cahg.map.list$geographic.locations[,2])+0.8, cahg.map.list$name.vector, col=locations$Col, cex=0.8, font=2)
dev.off()

tiff(file="map_subsistence_real.tiff", units="in", width=7, height=6, res=500)
plot(g, col=rev(hcl.colors(100, palette="Earth")), legend=F, main="Subsistence tools", box=F, axes=F)
map(add=T)
points(cahg.map.list$geographic.locations[,1], cahg.map.list$geographic.locations[,2],col=locations$Col, pch=15)
points(subsistence.map.list$MAPP.geogen.coords[,1], subsistence.map.list$MAPP.geogen.coords[,2], col=locations$Col, pch=16)
points(cahg.map.list$MAPP.geogen.coords[,1], cahg.map.list$MAPP.geogen.coords[,2], col=locations$Col, pch=17)
segments(cahg.map.list$geographic.locations[,1], cahg.map.list$geographic.locations[,2], subsistence.map.list$MAPP.geogen.coords[,1], subsistence.map.list$MAPP.geogen.coords[,2], col=locations$Col, lwd=1.5)
segments(cahg.map.list$geographic.locations[,1], cahg.map.list$geographic.locations[,2], cahg.map.list$MAPP.geogen.coords[,1], cahg.map.list$MAPP.geogen.coords[,2], col=locations$Col, lwd=1.5)
segments(cahg.map.list$MAPP.geogen.coords[,1], cahg.map.list$MAPP.geogen.coords[,2], subsistence.map.list$MAPP.geogen.coords[,1], subsistence.map.list$MAPP.geogen.coords[,2], col=locations$Col, lwd=1.5)
text(cahg.map.list$geographic.locations[,1], (cahg.map.list$geographic.locations[,2])+0.8, cahg.map.list$name.vector, col=locations$Col, cex=0.8, font=2)
dev.off()

# Make admixture map

fade.admixture.source.points <- function(pop.cols,admix.proportions){
  faded.colors <- numeric(length(pop.cols))
  for(i in 1:length(pop.cols)){
    faded.colors[i] <- adjustcolor(pop.cols[i],admix.proportions[i], alpha.f = 0.4, alpha.f = 0.4, alpha.f = 0.4)
  }
  return(faded.colors)
}


make.spacemix.map.list <- function(MCMC.output.file,geographic.locations,name.vector,color.vector,quantile=0.95,burnin=0){
  MCMC.output <- load_MCMC_output(MCMC.output.file)
  x <- seq((burnin+1),length(MCMC.output$Prob),by=1)
  best.iter <- which.max(MCMC.output$Prob[x])
  if(max(MCMC.output$Prob) > max(MCMC.output$Prob[x])){
    message("warning: the MCMC iteration with the\n 
				highest posterior probability is being\n 
				discarded as burn-in")
  }
  k <- MCMC.output$last.params$k
  if(is.null(color.vector)){
    color.vector <- rep(1,k)
  }
  admix.source.color.vector <- fade.admixture.source.points(color.vector,rowMeans(MCMC.output$admix.proportions[,x]))
  MAPP.geogen.coords <- spacemix.procrustes(X = geographic.locations,
                                            Y = MCMC.output$population.coordinates[[best.iter]][1:k,],
                                            k = k,option = 1)
  MAPP.admix.source.coords <- spacemix.procrustes(X = geographic.locations,
                                                  Y = MCMC.output$population.coordinates[[best.iter]][1:k,],
                                                  k = k,
                                                  admix.source.locs = MCMC.output$population.coordinates[[best.iter]][(k+1):(2*k),],
                                                  option=2)
  procrustes.coord.posterior.lists <- get.procrustes.locations.posterior.list(geographic.locations = geographic.locations,
                                                                              population.coordinates.posterior = MCMC.output$population.coordinates[x])
  pp.geogen.location.matrices <- lapply(1:k,get.posterior.location.matrix.from.list,posterior.list=procrustes.coord.posterior.lists$geogen.coords.list)
  pp.admix.source.location.matrices <- lapply(1:k,get.posterior.location.matrix.from.list,posterior.list=procrustes.coord.posterior.lists$admix.source.coords.list)
  pp.geogen.ellipses <- lapply(pp.geogen.location.matrices,get.credible.ellipse,quantile)
  pp.admix.source.ellipses <- lapply(pp.admix.source.location.matrices,get.credible.ellipse,quantile)
  spacemix.map.list <- c(list(MCMC.output=MCMC.output),
                         list(geographic.locations=geographic.locations),
                         list(name.vector=name.vector),list(color.vector=color.vector),
                         list(quantile=quantile),list(best.iter = best.iter),
                         list(admix.source.color.vector = admix.source.color.vector),
                         list(k = k),list(MAPP.geogen.coords = MAPP.geogen.coords),
                         list(MAPP.admix.source.coords = MAPP.admix.source.coords),
                         list(procrustes.coord.posterior.lists = procrustes.coord.posterior.lists),
                         list(pp.geogen.location.matrices = pp.geogen.location.matrices),
                         list(pp.admix.source.location.matrices = pp.admix.source.location.matrices),
                         list(pp.geogen.ellipses = pp.geogen.ellipses),
                         list(pp.admix.source.ellipses = pp.admix.source.ellipses))
  return(spacemix.map.list)
}



cahg.map.list <- make.spacemix.map.list(MCMC.output.file = "/Users/Cecilia/Documents/PhD/African_Pygmies/Spatial analysis culture/My_data/Nov_2021/spacemix/run_793742/cahg_spacemix_all_LongRun/cahg_spacemix_all_space_MCMC_output1.Robj",
                                        geographic.locations = as.matrix(locations[,2:3]), name.vector = locations$Culture,color.vector = locations$Col, quantile=0.95, burnin=0)


make.spacemix.map <- function(spacemix.map.list,text=FALSE,ellipses=TRUE,source.option=TRUE,xlim=NULL,ylim=NULL){
  with(spacemix.map.list,{ 
    plot(MAPP.geogen.coords,type='n',xlim=xlim,ylim=ylim,xlab="",ylab="")
    if(ellipses){
      lapply(1:k,FUN=function(i){plot.credible.ellipse(pp.geogen.ellipses[[i]],color.vector[i])})
    }
    if(text){
      text(MAPP.geogen.coords+0.7,col="black",font=2,labels=name.vector,cex=1)
    }
    if(source.option){
      if(ellipses){
        lapply(1:k,FUN=function(i){plot.credible.ellipse(pp.admix.source.ellipses[[i]],admix.source.color.vector[i],fading=3,lty=2)})
      }
      text(MAPP.admix.source.coords,col= admix.source.color.vector,font=3,labels=name.vector,cex=0.7)
      plot.admix.arrows(MAPP.admix.source.coords, MAPP.geogen.coords,
                        admix.proportions=MCMC.output$admix.proportions[,best.iter],
                        colors=admix.source.color.vector,length=0.1)
      print(MCMC.output$admix.proportions[,best.iter])
    }
    box(lwd=1)
  })
  return(invisible("spacemix map!"))
}


# Manually input admixture props
admix_proportions <- c (0.1099353843, 0.0540065167, 0.0023568149, 0.0006655948, 0.0041017206,
                        0.0630641436, 0.0608503573, 0.1282241622, 0.0100971584)
# Modify Spacemix function
query.spacemix.map <- function(focal.pops,spacemix.map.list,ellipses=TRUE,source.option=TRUE){
  with(spacemix.map.list,{
    # browser()
    focal.indices <- match(focal.pops,name.vector)
    if(ellipses){
      for(i in 1:length(focal.indices)){
        plot.credible.ellipse(pp.geogen.ellipses[[focal.indices[i]]],color.vector[focal.indices[i]],)
      }
    }
    if(source.option){
      if(ellipses){
        for(i in 1:length(focal.indices)){
          plot.credible.ellipse(pp.admix.source.ellipses[[focal.indices[i]]], color.vector[focal.indices[i]],lty=2, fading=0.2)
        }
      }
      edge.ind <- admix_proportions[focal.indices]^-0.4
      print(edge.ind)
      arrows(	x0 = MAPP.admix.source.coords[focal.indices,1],
              y0 = MAPP.admix.source.coords[focal.indices,2],
              x1 = MAPP.geogen.coords[focal.indices,1],
              y1 = MAPP.geogen.coords[focal.indices,2],
              col= color.vector[focal.indices],
              lwd=edge.ind,
              length=0.1)
      text(MAPP.admix.source.coords[focal.indices,,drop=FALSE],col=1,font=3,labels=name.vector[focal.indices])
    }
    text(MAPP.geogen.coords[focal.indices,,drop=FALSE]+0.7,col=1,font=2,labels=name.vector[focal.indices],cex=1)
    box(lwd=1)
  })
  return(invisible("highlighted samples!"))
}


# Edited this function so I can see the admixture proportions from the besst run

tiff(file="admixture_sources.tiff", units="in", width=7.5, height=6, res=500)
make.spacemix.map(spacemix.map.list = cahg.map.list,
                  text=T,
                  ellipses=T,
                  source.option=F,xlim=c(5,35),ylim=c(-10,15))

# Select populations with >1% admixture from their source location
# all but Efe, Sua and Bakola
query.spacemix.map(focal.pops = locations$Culture[-c(3,4,5)],
                                 spacemix.map.list = cahg.map.list,
                                 ellipses=T,
                                 source.option=T)
dev.off()


