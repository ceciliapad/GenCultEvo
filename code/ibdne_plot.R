## PLOT IBDNE PATIN AND JARVIS

setwd("/Users/Cecilia/ncbi/public/downloads/dbGaP-25678/86588/PhenoGenotypeFiles/RootStudyConsentSet_phs001780.AfricanDemographicHistory.v1.p1.c1.GRU-PUB-NPU/GenotypeFiles/Genotype_files/ibdne")

t_col <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color
  
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  
  ## Save the color
  invisible(t.col)
}

require(Jodorowsky)
require(MetBrewer)

# it used to be 70% but like that the visibility is crap
#bantu_red <- t_col(holy_mountain(2)[5], percent=50)
bantu_red <-  t_col(MetPalettes$Isfahan2[[1]][1], percent=50)

# CAMEROON AND CAR

cameroon_here <- c("Baka", "Aka", "Bakola", "Badwee", "Nzime")
colors_cam <-  c(MetPalettes$Isfahan1[[1]][3],MetPalettes$Isfahan1[[1]][1],MetPalettes$Isfahan1[[1]][4], bantu_red, bantu_red)

tiff("cameroon_ibdne_pres.tiff",units="in", width=7, height=5, res=500)
par(xpd = T, mar = par()$mar + c(0,0,0,7))
plot(x[,1],x[,2],type="n",log="y",xlim=c(-150, -20),ylim=c(1e2,1e8),xlab="g (generations before present)",ylab="Ne (effective population size)",xaxs="i",yaxt="n",main=paste("Northwestern region"), cex.lab=0.8, cex.axis=0.7)
axis(2,at=c(1e3,1e4,1e5,1e6,1e7,1e8),labels=c(expression(10^3),expression(10^4),expression(10^5),expression(10^6),expression(10^7), expression(10^8)),las=1, cex.axis=0.7)
i <- 1
for(anc in cameroon_here){
  filename=paste(anc,".ne",sep="")
  # only read rows from g=4 to g=150
  x=read.table(filename,header=T)[21:151,]
  x$GEN <- -x$GEN
  lines(x[,1],x[,2],las=1,lwd=3,col=colors_cam[i],)
  i <- i + 1
}
legend("right", legend=c("Baka", "Aka", "Bakola", "Bantu"), col=colors_cam,lty=c(1), lwd=c(2),cex=0.8,inset=c(-0.3,0),
       bty = "n", xpd=TRUE)
dev.off()

#############
### GABON ###
#############

gabon_here <- c("Babongo", "Bakoya",
                "Bedzan", "Akele", "Benga", "Duma", "Eshira", "Eviya", "Fang", "Galoa",
                "Makina", "Ndumu", "Obamba",
                 "Shake")
bantu_gabon <- c("Akele", "Benga", "Duma", "Eshira", "Eviya", "Fang", "Galoa",
                 "Makina", "Ndumu", "Obamba",
                  "Shake")

colors_gabon <-  c(MetPalettes$Isfahan1[[1]][2],MetPalettes$Isfahan1[[1]][5],
                   MetPalettes$Isfahan1[[1]][7],
                   bantu_red)

tiff("gabon_ibdne_pres.tiff",units="in", width=7, height=5, res=500)
par(xpd = T, mar = par()$mar + c(0,0,0,7))
plot(x[,1],x[,2],type="n",log="y",xlim=c(-150, -20),ylim=c(1e2,1e8),xlab="g (generations before present)",ylab="Ne (effective population size)",xaxs="i",yaxt="n",main=paste("Western region (Gabon)"), cex.lab=0.8, cex.axis=0.7)
axis(2,at=c(1e2,1e3,1e4,1e5,1e6,1e7, 1e8),labels=c(expression(10^2),expression(10^3),expression(10^4),expression(10^5),expression(10^6),expression(10^7), expression(10^8)),las=1, cex.axis=0.7)
i <- 1
col_this <- 0
for(anc in gabon_here){
  filename=paste(anc,".ne",sep="")
  # only read rows from g=4 to g=150
  x=read.table(filename,header=T)[21:151,]
  x$GEN <- -x$GEN
  if (anc =="Babongo") {
        col_this <- 1 } else if (anc == "Bakoya") {
          col_this <- 2 } else if (anc == "Bedzan") {
            col_this <- 3 } else if (anc %in% bantu_gabon) {
              col_this <- 4 }
  lines(x[,1],x[,2],las=1,lwd=3,col=colors_gabon[col_this])
  i <- i + 1
}
legend("right", legend=c("Babongo", "Bakoya", "Bedzan", "Bantu"), col=colors_gabon,lty=c(1), lwd=c(2),cex=0.8,inset=c(-0.3,0),
       bty = "n", xpd=TRUE)
dev.off()

#############
### EAST ####
#############

eastern_here <- c("Batwa", "Mbuti","Bakiga")
colors_east <-  c(MetPalettes$Isfahan1[[1]][6],MetPalettes$Isfahan2[[1]][2], bantu_red, bantu_red)

png("eastern_ibdne.tiff",units="in", width=7, height=5, res=500)
par(xpd = T, mar = par()$mar + c(0,0,0,7))
plot(x[,1],x[,2],type="n",log="y",xlim=c(-150, -20),ylim=c(1e2,1e8),xlab="g (generations before present)",ylab="Ne (effective population size)",xaxs="i",yaxt="n",main=paste("Eastern region"), cex.lab=0.8, cex.axis=0.7)
axis(2,at=c(1e3,1e4,1e5,1e6,1e7, 1e8),labels=c(expression(10^3),expression(10^4),expression(10^5),expression(10^6),expression(10^7), expression(10^8)),las=1, cex.axis=0.7)
i <- 1
for(anc in eastern_here){
  filename=paste(anc,".ne",sep="")
  # only read rows from g=4 to g=150
  x=read.table(filename,header=T)[21:151,]
  x$GEN <- -x$GEN
  lines(x[,1],x[,2],las=1,lwd=3,col=colors_east[i],)
  i <- i + 1
}
legend("right", legend=c("Batwa (East)", "Efe & Sua", "Bantu"), col=colors_east,lty=c(1), lwd=c(2),cex=0.8,inset=c(-0.3,0),
       bty = "n", xpd=TRUE)
dev.off()

#################
### OTHER HG ####
#################

hg_here <- c("Hadza", "Ogiek","Sandawe", "Boni")
colors_hg <-  c(MetPalettes$Isfahan2[[1]])

png("hg_ibdne.tiff",units="in", width=7, height=5, res=500)
par(xpd = T, mar = par()$mar + c(0,0,0,7))
plot(x[,1],x[,2],type="n",log="y",xlim=c(-150, -20),ylim=c(1e2,1e8),xlab="g (generations before present)",ylab="Ne (effective population size)",xaxs="i",yaxt="n",main=paste("Other African hunter-gatherers"), cex.lab=0.8, cex.axis=0.7)
axis(2,at=c(1e3,1e4,1e5,1e6,1e7, 1e8),labels=c(expression(10^3),expression(10^4),expression(10^5),expression(10^6),expression(10^7), expression(10^8)),las=1, cex.axis=0.7)
i <- 1
for(anc in hg_here){
  filename=paste(anc,".ne",sep="")
  # only read rows from g=4 to g=150
  x=read.table(filename,header=T)[21:151,]
  x$GEN <- -x$GEN
  lines(x[,1],x[,2],las=1,lwd=3,col=colors_hg[i],)
  i <- i + 1
}
legend("right", legend=c("Hadza", "Ogiek", "Sandawe", "Boni"), col=colors_hg,lty=c(1), lwd=c(2),cex=0.8,inset=c(-0.3,0),
       bty = "n", xpd=TRUE)
dev.off()
