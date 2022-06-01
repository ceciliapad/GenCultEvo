# Script to extract environmental variables from cultural group ranges

setwd("/Users/Cecilia/Documents/PhD/African Pygmies/Spatial analysis culture/My_data")

# IMPORT SHAPEFILE 

library(sp)
library(rgdal)
library(raster)
library(RColorBrewer)
library(tidyverse)
require(sf) # Simple Features for R
library(stars) # Spatiotemporal Arrays, Raster and Vector Data Cubes


cultures <- readOGR("../My_data/shapefile_cultures/cultural.groupsfixed.shp", stringsAsFactors = TRUE)
plot(cultures)

## remember that we had worldclim data already downloaded
setwd("/Users/Cecilia/Documents/PhD/Courses_Tutorial/ENM_15_06/Data/worldclim")
var<- list.files( pattern=".tif") # do this in the path you have the worldclim variables
wc<- stack (var) 
setwd("/Users/Cecilia/Documents/PhD/African Pygmies/Spatial analysis culture/My_data")

## we can extract the climatic data for the polygons with cultures
bio_cult <- extract(wc, cultures, fun=mean, na.rm=TRUE, sp = T)

environment_cult <- as.data.frame(bio_cult@data)
write.csv(environment_cult, "environmennt_cultures.csv")

# Analyse relationship between environment and heterozygosity

env_gen <- read.csv("environmennt_cultures_gen.csv")
env_gen <- env_gen[1:10,] # Eliminate Twa (looks very weird)
z <- lm(env_gen$median_hom ~ env_gen$prec_season)

# but it seems like minimum annual temperature and maximum annual temperature
# related negatively to homozygosity. That is, populations living in places
# with higher minimum temperatures and lower maximum ones (less extreme)
# have less RoH

# What about cultural richess?
env_cult <- read.csv("environmennt_cultures_traits.csv")
z <- lm(env_cult$music_counts ~ env_cult$prec_season)

library(lme4) ## Run the same model in non-Bayesian framework

g <- lm (formula=music_counts ~ prec_season, data=env_cult)

drop1(g, test="Chisq")

