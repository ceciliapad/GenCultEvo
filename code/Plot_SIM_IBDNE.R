####################################
### Plot simulated IBDNE results ###
####################################

setwd("/Users/Cecilia/Documents/PhD/African_Pygmies/Genetic_studies_pygmies/Final_Data/New genetics march 2021/Revision_october_Lane/sims/constant_OUT")

# plot the results in R

pops <- c(15,20,50,150)
colors_cam <-  c(holy_mountain(2)[2],holy_mountain(5)[5],holy_mountain(2)[4], holy_mountain(5)[6])

png("constant_04_550000.png",units="in", width=7, height=5, res=500)
par(xpd = T, mar = par()$mar + c(0,0,0,7))
plot(x[,1],x[,2],type="n",log="y",xlim=c(-150, -4),ylim=c(1e2,1e20),xlab="g (generations before present)",ylab="Ne (effective population size)",xaxs="i",yaxt="n",main=paste("Constant population MAF=0.4, 550,000SNP"), cex.lab=0.8, cex.axis=0.7)
axis(2,at=c(1e3,1e4,1e5,1e6,1e7, 1e8, 1e9,1e10,1e11,1e12,1e13,1e14,1e15, 1e16, 1e17, 1e18, 1e19),labels=c(expression(10^3),expression(10^4),expression(10^5),expression(10^6),expression(10^7), expression(10^8), expression(10^9), expression(10^10),
                                                                                                          expression(10^11), expression(10^12), expression(10^13), expression(10^14), expression(10^15), expression(10^16),
                                                                                                          expression(10^17), expression(10^18), expression(10^19)),las=1, cex.axis=0.7)
i <- 1
for(pop in pops){
  filename=paste("constant.",pop,".maf0.4.ALL_sorted.550000_OUT.ne",sep="")
  # only read rows from g=4 to g=150
  x=read.table(filename,header=T)[5:151,]
  x$GEN <- -x$GEN
  lines(x[,1],x[,2],las=1,lwd=3,col=colors_cam[i],)
  col2="#00000080"
  #col2="azure3"
  polygon(c(x[,1],x[nrow(x):1,1]),c(x[,3],x[nrow(x):1,4]),col=col2,density=-1,border=col2)
  i <- i + 1
}
legend("right", legend=c("15","20", "50", "150"), col=colors_cam,lty=c(1), lwd=c(2),cex=0.8,inset=c(-0.3,0),
       bty = "n", xpd=TRUE)
dev.off()

png("constant_04_all.png",units="in", width=7, height=5, res=500)
par(xpd = T, mar = par()$mar + c(0,0,0,7))
plot(x[,1],x[,2],type="n",log="y",xlim=c(-150, -4),ylim=c(1e2,1e20),xlab="g (generations before present)",ylab="Ne (effective population size)",xaxs="i",yaxt="n",main=paste("Constant population MAF=0.4"), cex.lab=0.8, cex.axis=0.7)
axis(2,at=c(1e3,1e4,1e5,1e6,1e7, 1e8, 1e9,1e10,1e11,1e12,1e13,1e14,1e15, 1e16, 1e17, 1e18, 1e19),labels=c(expression(10^3),expression(10^4),expression(10^5),expression(10^6),expression(10^7), expression(10^8), expression(10^9), expression(10^10),
                                                                                                          expression(10^11), expression(10^12), expression(10^13), expression(10^14), expression(10^15), expression(10^16),
                                                                                                          expression(10^17), expression(10^18), expression(10^19)),las=1, cex.axis=0.7)
i <- 1
for(pop in pops){
  filename=paste("constant.",pop,".fixed.maf0.4_OUT.ne",sep="")
  # only read rows from g=4 to g=150
  x=read.table(filename,header=T)[5:151,]
  x$GEN <- -x$GEN
  lines(x[,1],x[,2],las=1,lwd=3,col=colors_cam[i],)
  col2="#00000080"
  #col2="azure3"
  polygon(c(x[,1],x[nrow(x):1,1]),c(x[,3],x[nrow(x):1,4]),col=col2,density=-1,border=col2)
  i <- i + 1
}
legend("right", legend=c("15","20", "50", "150"), col=colors_cam,lty=c(1), lwd=c(2),cex=0.8,inset=c(-0.3,0),
       bty = "n", xpd=TRUE)
dev.off()

png("constant_005_550000.png",units="in", width=7, height=5, res=500)
par(xpd = T, mar = par()$mar + c(0,0,0,7))
plot(x[,1],x[,2],type="n",log="y",xlim=c(-150, -4),ylim=c(1e2,1e20),xlab="g (generations before present)",ylab="Ne (effective population size)",xaxs="i",yaxt="n",main=paste("Constant population MAF=0.05, 550,000SNP"), cex.lab=0.8, cex.axis=0.7)
axis(2,at=c(1e3,1e4,1e5,1e6,1e7, 1e8, 1e9,1e10,1e11,1e12,1e13,1e14,1e15, 1e16, 1e17, 1e18, 1e19),labels=c(expression(10^3),expression(10^4),expression(10^5),expression(10^6),expression(10^7), expression(10^8), expression(10^9), expression(10^10),
                                                                                                          expression(10^11), expression(10^12), expression(10^13), expression(10^14), expression(10^15), expression(10^16),
                                                                                                          expression(10^17), expression(10^18), expression(10^19)),las=1, cex.axis=0.7)
i <- 1
for(pop in pops){
  filename=paste("constant.",pop,".maf0.05.ALL_sorted.550000_OUT.ne",sep="")
  # only read rows from g=4 to g=150
  x=read.table(filename,header=T)[5:151,]
  x$GEN <- -x$GEN
  lines(x[,1],x[,2],las=1,lwd=3,col=colors_cam[i],)
  col2="#00000080"
  #col2="azure3"
  polygon(c(x[,1],x[nrow(x):1,1]),c(x[,3],x[nrow(x):1,4]),col=col2,density=-1,border=col2)
  i <- i + 1
}
legend("right", legend=c("15","20",  "50", "150"), col=colors_cam,lty=c(1), lwd=c(2),cex=0.8,inset=c(-0.3,0),
       bty = "n", xpd=TRUE)
dev.off()

png("constant_005_all.png",units="in", width=7, height=5, res=500)
par(xpd = T, mar = par()$mar + c(0,0,0,7))
plot(x[,1],x[,2],type="n",log="y",xlim=c(-150, -4),ylim=c(1e2,1e20),xlab="g (generations before present)",ylab="Ne (effective population size)",xaxs="i",yaxt="n",main=paste("Constant population MAF=0.05"), cex.lab=0.8, cex.axis=0.7)
axis(2,at=c(1e3,1e4,1e5,1e6,1e7, 1e8, 1e9,1e10,1e11,1e12,1e13,1e14,1e15, 1e16, 1e17, 1e18, 1e19),labels=c(expression(10^3),expression(10^4),expression(10^5),expression(10^6),expression(10^7), expression(10^8), expression(10^9), expression(10^10),
                                                                                                          expression(10^11), expression(10^12), expression(10^13), expression(10^14), expression(10^15), expression(10^16),
                                                                                                          expression(10^17), expression(10^18), expression(10^19)),las=1, cex.axis=0.7)
i <- 1
for(pop in pops){
  filename=paste("constant.",pop,".fixed.maf0.05_OUT.ne",sep="")
  # only read rows from g=4 to g=150
  x=read.table(filename,header=T)[5:151,]
  x$GEN <- -x$GEN
  lines(x[,1],x[,2],las=1,lwd=3,col=colors_cam[i],)
  col2="#00000080"
  #col2="azure3"
  polygon(c(x[,1],x[nrow(x):1,1]),c(x[,3],x[nrow(x):1,4]),col=col2,density=-1,border=col2)
  i <- i + 1
}
legend("right", legend=c("15","20", "50", "150"), col=colors_cam,lty=c(1), lwd=c(2),cex=0.8,inset=c(-0.3,0),
       bty = "n", xpd=TRUE)
dev.off()

###########

setwd("/Users/Cecilia/Documents/PhD/African_Pygmies/Genetic_studies_pygmies/Final_Data/New genetics march 2021/Revision_october_Lane/sims/expansion_OUT")

png("expansion_04_550000.png",units="in", width=7, height=5, res=500)
par(xpd = T, mar = par()$mar + c(0,0,0,7))
plot(x[,1],x[,2],type="n",log="y",xlim=c(-150, -4),ylim=c(1e2,1e20),xlab="g (generations before present)",ylab="Ne (effective population size)",xaxs="i",yaxt="n",main=paste("Expansion MAF=0.4, 550,000SNP"), cex.lab=0.8, cex.axis=0.7)
axis(2,at=c(1e3,1e4,1e5,1e6,1e7, 1e8, 1e9,1e10,1e11,1e12,1e13,1e14,1e15, 1e16, 1e17, 1e18, 1e19),labels=c(expression(10^3),expression(10^4),expression(10^5),expression(10^6),expression(10^7), expression(10^8), expression(10^9), expression(10^10),
                                                                                                          expression(10^11), expression(10^12), expression(10^13), expression(10^14), expression(10^15), expression(10^16),
                                                                                                          expression(10^17), expression(10^18), expression(10^19)),las=1, cex.axis=0.7)
i <- 1
for(pop in pops){
  filename=paste("expansion.",pop,".maf0.4_sorted.550000_OUT.ne",sep="")
  # only read rows from g=4 to g=150
  x=read.table(filename,header=T)[5:151,]
  x$GEN <- -x$GEN
  lines(x[,1],x[,2],las=1,lwd=3,col=colors_cam[i],)
  col2="#00000080"
  #col2="azure3"
  polygon(c(x[,1],x[nrow(x):1,1]),c(x[,3],x[nrow(x):1,4]),col=col2,density=-1,border=col2)
  i <- i + 1
}
legend("right", legend=c("15", "20", "50", "150"), col=colors_cam,lty=c(1), lwd=c(2),cex=0.8,inset=c(-0.3,0),
       bty = "n", xpd=TRUE)
dev.off()

png("expansion_04_all.png",units="in", width=7, height=5, res=500)
par(xpd = T, mar = par()$mar + c(0,0,0,7))
plot(x[,1],x[,2],type="n",log="y",xlim=c(-150, -4),ylim=c(1e2,1e20),xlab="g (generations before present)",ylab="Ne (effective population size)",xaxs="i",yaxt="n",main=paste("Expansion MAF=0.4"), cex.lab=0.8, cex.axis=0.7)
axis(2,at=c(1e3,1e4,1e5,1e6,1e7, 1e8, 1e9,1e10,1e11,1e12,1e13,1e14,1e15, 1e16, 1e17, 1e18, 1e19),labels=c(expression(10^3),expression(10^4),expression(10^5),expression(10^6),expression(10^7), expression(10^8), expression(10^9), expression(10^10),
                                                                                                          expression(10^11), expression(10^12), expression(10^13), expression(10^14), expression(10^15), expression(10^16),
                                                                                                          expression(10^17), expression(10^18), expression(10^19)),las=1, cex.axis=0.7)
i <- 1
for(pop in pops){
  filename=paste("expansion.",pop,".fixed.maf0.4_OUT.ne",sep="")
  # only read rows from g=4 to g=150
  x=read.table(filename,header=T)[5:151,]
  x$GEN <- -x$GEN
  lines(x[,1],x[,2],las=1,lwd=3,col=colors_cam[i],)
  col2="#00000080"
  #col2="azure3"
  polygon(c(x[,1],x[nrow(x):1,1]),c(x[,3],x[nrow(x):1,4]),col=col2,density=-1,border=col2)
  i <- i + 1
}
legend("right", legend=c("15", "20","50", "150"), col=colors_cam,lty=c(1), lwd=c(2),cex=0.8,inset=c(-0.3,0),
       bty = "n", xpd=TRUE)
dev.off()

png("expansion_005_550000.png",units="in", width=7, height=5, res=500)
par(xpd = T, mar = par()$mar + c(0,0,0,7))
plot(x[,1],x[,2],type="n",log="y",xlim=c(-150, -4),ylim=c(1e2,1e20),xlab="g (generations before present)",ylab="Ne (effective population size)",xaxs="i",yaxt="n",main=paste("Expansion MAF=0.05, 550,000SNP"), cex.lab=0.8, cex.axis=0.7)
axis(2,at=c(1e3,1e4,1e5,1e6,1e7, 1e8, 1e9,1e10,1e11,1e12,1e13,1e14,1e15, 1e16, 1e17, 1e18, 1e19),labels=c(expression(10^3),expression(10^4),expression(10^5),expression(10^6),expression(10^7), expression(10^8), expression(10^9), expression(10^10),
                                                                                                          expression(10^11), expression(10^12), expression(10^13), expression(10^14), expression(10^15), expression(10^16),
                                                                                                          expression(10^17), expression(10^18), expression(10^19)),las=1, cex.axis=0.7)
i <- 1
for(pop in pops){
  filename=paste("expansion.",pop,".maf0.05_sorted.550000_OUT.ne",sep="")
  # only read rows from g=4 to g=150
  x=read.table(filename,header=T)[5:151,]
  x$GEN <- -x$GEN
  lines(x[,1],x[,2],las=1,lwd=3,col=colors_cam[i],)
  col2="#00000080"
  #col2="azure3"
  polygon(c(x[,1],x[nrow(x):1,1]),c(x[,3],x[nrow(x):1,4]),col=col2,density=-1,border=col2)
  i <- i + 1
}
legend("right", legend=c("15", "20", "50", "150"), col=colors_cam,lty=c(1), lwd=c(2),cex=0.8,inset=c(-0.3,0),
       bty = "n", xpd=TRUE)
dev.off()

png("expansion_005_all.png",units="in", width=7, height=5, res=500)
par(xpd = T, mar = par()$mar + c(0,0,0,7))
plot(x[,1],x[,2],type="n",log="y",xlim=c(-150, -4),ylim=c(1e2,1e20),xlab="g (generations before present)",ylab="Ne (effective population size)",xaxs="i",yaxt="n",main=paste("Expansion MAF=0.05"), cex.lab=0.8, cex.axis=0.7)
axis(2,at=c(1e3,1e4,1e5,1e6,1e7, 1e8, 1e9,1e10,1e11,1e12,1e13,1e14,1e15, 1e16, 1e17, 1e18, 1e19),labels=c(expression(10^3),expression(10^4),expression(10^5),expression(10^6),expression(10^7), expression(10^8), expression(10^9), expression(10^10),
                                                                                                          expression(10^11), expression(10^12), expression(10^13), expression(10^14), expression(10^15), expression(10^16),
                                                                                                          expression(10^17), expression(10^18), expression(10^19)),las=1, cex.axis=0.7)
i <- 1
for(pop in pops){
  filename=paste("expansion.",pop,".fixed.maf0.05_OUT.ne",sep="")
  # only read rows from g=4 to g=150
  x=read.table(filename,header=T)[5:151,]
  x$GEN <- -x$GEN
  lines(x[,1],x[,2],las=1,lwd=3,col=colors_cam[i],)
  col2="#00000080"
  #col2="azure3"
  polygon(c(x[,1],x[nrow(x):1,1]),c(x[,3],x[nrow(x):1,4]),col=col2,density=-1,border=col2)
  i <- i + 1
}
legend("right", legend=c("15","20", "50", "150"), col=colors_cam,lty=c(1), lwd=c(2),cex=0.8,inset=c(-0.3,0),
       bty = "n", xpd=TRUE)
dev.off()

###########

setwd("/Users/Cecilia/Documents/PhD/African_Pygmies/Genetic_studies_pygmies/Final_Data/New genetics march 2021/Revision_october_Lane/sims/collapse_OUT")

png("collapse_04_550000.png",units="in", width=7, height=5, res=500)
par(xpd = T, mar = par()$mar + c(0,0,0,7))
plot(x[,1],x[,2],type="n",log="y",xlim=c(-150, -4),ylim=c(1e2,1e20),xlab="g (generations before present)",ylab="Ne (effective population size)",xaxs="i",yaxt="n",main=paste("Collapse MAF=0.4, 550,000SNP"), cex.lab=0.8, cex.axis=0.7)
axis(2,at=c(1e3,1e4,1e5,1e6,1e7, 1e8, 1e9,1e10,1e11,1e12,1e13,1e14,1e15, 1e16, 1e17, 1e18, 1e19),labels=c(expression(10^3),expression(10^4),expression(10^5),expression(10^6),expression(10^7), expression(10^8), expression(10^9), expression(10^10),
                                                                                                          expression(10^11), expression(10^12), expression(10^13), expression(10^14), expression(10^15), expression(10^16),
                                                                                                          expression(10^17), expression(10^18), expression(10^19)),las=1, cex.axis=0.7)
i <- 1
for(pop in pops){
  filename=paste("collapse.",pop,".maf0.4_sorted.550000_OUT.ne",sep="")
  # only read rows from g=4 to g=150
  x=read.table(filename,header=T)[5:151,]
  x$GEN <- -x$GEN
  lines(x[,1],x[,2],las=1,lwd=3,col=colors_cam[i],)
  col2="#00000080"
  #col2="azure3"
  polygon(c(x[,1],x[nrow(x):1,1]),c(x[,3],x[nrow(x):1,4]),col=col2,density=-1,border=col2)
  i <- i + 1
}
legend("right", legend=c("15", "20", "50", "150"), col=colors_cam,lty=c(1), lwd=c(2),cex=0.8,inset=c(-0.3,0),
       bty = "n", xpd=TRUE)
dev.off()

png("collapse_04_all.png",units="in", width=7, height=5, res=500)
par(xpd = T, mar = par()$mar + c(0,0,0,7))
plot(x[,1],x[,2],type="n",log="y",xlim=c(-150, -4),ylim=c(1e2,1e20),xlab="g (generations before present)",ylab="Ne (effective population size)",xaxs="i",yaxt="n",main=paste("Collapse MAF=0.4"), cex.lab=0.8, cex.axis=0.7)
axis(2,at=c(1e3,1e4,1e5,1e6,1e7, 1e8, 1e9,1e10,1e11,1e12,1e13,1e14,1e15, 1e16, 1e17, 1e18, 1e19),labels=c(expression(10^3),expression(10^4),expression(10^5),expression(10^6),expression(10^7), expression(10^8), expression(10^9), expression(10^10),
                                                                                                          expression(10^11), expression(10^12), expression(10^13), expression(10^14), expression(10^15), expression(10^16),
                                                                                                          expression(10^17), expression(10^18), expression(10^19)),las=1, cex.axis=0.7)
i <- 1
for(pop in pops){
  filename=paste("collapse.",pop,".fixed.maf0.4_OUT.ne",sep="")
  # only read rows from g=4 to g=150
  x=read.table(filename,header=T)[5:151,]
  x$GEN <- -x$GEN
  lines(x[,1],x[,2],las=1,lwd=3,col=colors_cam[i],)
  col2="#00000080"
  #col2="azure3"
  polygon(c(x[,1],x[nrow(x):1,1]),c(x[,3],x[nrow(x):1,4]),col=col2,density=-1,border=col2)
  i <- i + 1
}
legend("right", legend=c("15", "20","50", "150"), col=colors_cam,lty=c(1), lwd=c(2),cex=0.8,inset=c(-0.3,0),
       bty = "n", xpd=TRUE)
dev.off()

png("collapse_005_550000.png",units="in", width=7, height=5, res=500)
par(xpd = T, mar = par()$mar + c(0,0,0,7))
plot(x[,1],x[,2],type="n",log="y",xlim=c(-150, -4),ylim=c(1e2,1e20),xlab="g (generations before present)",ylab="Ne (effective population size)",xaxs="i",yaxt="n",main=paste("Collapse MAF=0.05, 550,000SNP"), cex.lab=0.8, cex.axis=0.7)
axis(2,at=c(1e3,1e4,1e5,1e6,1e7, 1e8, 1e9,1e10,1e11,1e12,1e13,1e14,1e15, 1e16, 1e17, 1e18, 1e19),labels=c(expression(10^3),expression(10^4),expression(10^5),expression(10^6),expression(10^7), expression(10^8), expression(10^9), expression(10^10),
                                                                                                          expression(10^11), expression(10^12), expression(10^13), expression(10^14), expression(10^15), expression(10^16),
                                                                                                          expression(10^17), expression(10^18), expression(10^19)),las=1, cex.axis=0.7)
i <- 1
for(pop in pops){
  filename=paste("collapse.",pop,".maf0.05_sorted.550000_OUT.ne",sep="")
  # only read rows from g=4 to g=150
  x=read.table(filename,header=T)[5:151,]
  x$GEN <- -x$GEN
  lines(x[,1],x[,2],las=1,lwd=3,col=colors_cam[i],)
  col2="#00000080"
  #col2="azure3"
  polygon(c(x[,1],x[nrow(x):1,1]),c(x[,3],x[nrow(x):1,4]),col=col2,density=-1,border=col2)
  i <- i + 1
}
legend("right", legend=c("15", "20", "50", "150"), col=colors_cam,lty=c(1), lwd=c(2),cex=0.8,inset=c(-0.3,0),
       bty = "n", xpd=TRUE)
dev.off()

png("collapse_005_all.png",units="in", width=7, height=5, res=500)
par(xpd = T, mar = par()$mar + c(0,0,0,7))
plot(x[,1],x[,2],type="n",log="y",xlim=c(-150, -4),ylim=c(1e2,1e20),xlab="g (generations before present)",ylab="Ne (effective population size)",xaxs="i",yaxt="n",main=paste("Collapse MAF=0.05"), cex.lab=0.8, cex.axis=0.7)
axis(2,at=c(1e3,1e4,1e5,1e6,1e7, 1e8, 1e9,1e10,1e11,1e12,1e13,1e14,1e15, 1e16, 1e17, 1e18, 1e19),labels=c(expression(10^3),expression(10^4),expression(10^5),expression(10^6),expression(10^7), expression(10^8), expression(10^9), expression(10^10),
                                                                                                          expression(10^11), expression(10^12), expression(10^13), expression(10^14), expression(10^15), expression(10^16),
                                                                                                          expression(10^17), expression(10^18), expression(10^19)),las=1, cex.axis=0.7)
i <- 2
for(pop in pops[-1]){
  filename=paste("collapse.",pop,".fixed.maf0.05_OUT.ne",sep="")
  # only read rows from g=4 to g=150
  x=read.table(filename,header=T)[5:151,]
  x$GEN <- -x$GEN
  lines(x[,1],x[,2],las=1,lwd=3,col=colors_cam[i],)
  col2="#00000080"
  #col2="azure3"
  polygon(c(x[,1],x[nrow(x):1,1]),c(x[,3],x[nrow(x):1,4]),col=col2,density=-1,border=col2)
  i <- i + 1
}
legend("right", legend=c("20", "50", "150"), col=colors_cam[-1],lty=c(1), lwd=c(2),cex=0.8,inset=c(-0.3,0),
       bty = "n", xpd=TRUE)
dev.off()

####################################
### Plot simulated IBDNE results ###
####################################

setwd("/Users/Cecilia/Documents/PhD/African_Pygmies/Genetic_studies_pygmies/Final_Data/New genetics march 2021/Revision_october_Lane/sims/constant_OUT")

# plot the results in R

pops <- c(15,20,50,150)
colors_cam <-  c(holy_mountain(2)[2],holy_mountain(5)[5],holy_mountain(2)[4], holy_mountain(5)[6])

png("constant_04_550000.png",units="in", width=7, height=5, res=500)
par(xpd = T, mar = par()$mar + c(0,0,0,7))
plot(x[,1],x[,2],type="n",log="y",xlim=c(-150, -20),ylim=c(1e2,1e12),xlab="g (generations before present)",ylab="Ne (effective population size)",xaxs="i",yaxt="n",main=paste("Constant population MAF=0.4, 550,000SNP"), cex.lab=0.8, cex.axis=0.7)
axis(2,at=c(1e3,1e4,1e5,1e6,1e7, 1e8, 1e9,1e10,1e11,1e12),labels=c(expression(10^3),expression(10^4),expression(10^5),expression(10^6),expression(10^7), expression(10^8), expression(10^9), expression(10^10),
                                                                                                          expression(10^11), expression(10^12)),las=1, cex.axis=0.7)
i <- 1
for(pop in pops){
  filename=paste("constant.",pop,".maf0.4.ALL_sorted.550000_OUT.ne",sep="")
  # only read rows from g=4 to g=150
  x=read.table(filename,header=T)[21:151,]
  x$GEN <- -x$GEN
  lines(x[,1],x[,2],las=1,lwd=3,col=colors_cam[i],)
  x$real <- 10000
  lines(x[,1],x[,5],las=1,lty=2,lwd=3,col="black")
  col2="#00000080"
  #col2="azure3"
  polygon(c(x[,1],x[nrow(x):1,1]),c(x[,3],x[nrow(x):1,4]),col=col2,density=-1,border=col2)
  i <- i + 1
}
legend("right", legend=c("15","20", "50", "150"), col=colors_cam,lty=c(1), lwd=c(2),cex=0.8,inset=c(-0.3,0),
       bty = "n", xpd=TRUE)
dev.off()

png("constant_04_all.png",units="in", width=7, height=5, res=500)
par(xpd = T, mar = par()$mar + c(0,0,0,7))
plot(x[,1],x[,2],type="n",log="y",xlim=c(-150, -20),ylim=c(1e2,1e12),xlab="g (generations before present)",ylab="Ne (effective population size)",xaxs="i",yaxt="n",main=paste("Constant population MAF=0.4"), cex.lab=0.8, cex.axis=0.7)
axis(2,at=c(1e3,1e4,1e5,1e6,1e7, 1e8, 1e9,1e10,1e11,1e12),labels=c(expression(10^3),expression(10^4),expression(10^5),expression(10^6),expression(10^7), expression(10^8), expression(10^9), expression(10^10),
                                                                                                          expression(10^11), expression(10^12)),las=1, cex.axis=0.7)
i <- 1
for(pop in pops){
  filename=paste("constant.",pop,".fixed.maf0.4_OUT.ne",sep="")
  # only read rows from g=4 to g=150
  x=read.table(filename,header=T)[21:151,]
  x$GEN <- -x$GEN
  x$real <- 10000
  lines(x[,1],x[,2],las=1,lwd=3,col=colors_cam[i],)
  x$real <- 10000
  lines(x[,1],x[,5],las=1,lty=2,lwd=3,col="black")
  col2="#00000080"
  #col2="azure3"
  polygon(c(x[,1],x[nrow(x):1,1]),c(x[,3],x[nrow(x):1,4]),col=col2,density=-1,border=col2)
  i <- i + 1
}
legend("right", legend=c("15","20", "50", "150"), col=colors_cam,lty=c(1), lwd=c(2),cex=0.8,inset=c(-0.3,0),
       bty = "n", xpd=TRUE)
dev.off()

png("constant_005_550000.png",units="in", width=7, height=5, res=500)
par(xpd = T, mar = par()$mar + c(0,0,0,7))
plot(x[,1],x[,2],type="n",log="y",xlim=c(-150, -20),ylim=c(1e2,1e12),xlab="g (generations before present)",ylab="Ne (effective population size)",xaxs="i",yaxt="n",main=paste("Constant population MAF=0.05, 550,000SNP"), cex.lab=0.8, cex.axis=0.7)
axis(2,at=c(1e3,1e4,1e5,1e6,1e7, 1e8, 1e9,1e10,1e11,1e12),labels=c(expression(10^3),expression(10^4),expression(10^5),expression(10^6),expression(10^7), expression(10^8), expression(10^9), expression(10^10),
                                                                                                          expression(10^11), expression(10^12)),las=1, cex.axis=0.7)
i <- 1
for(pop in pops){
  filename=paste("constant.",pop,".maf0.05.ALL_sorted.550000_OUT.ne",sep="")
  # only read rows from g=4 to g=150
  x=read.table(filename,header=T)[21:151,]
  x$GEN <- -x$GEN
  x$real <- 10000
  lines(x[,1],x[,2],las=1,lwd=3,col=colors_cam[i],)
  x$real <- 10000
  lines(x[,1],x[,5],las=1,lty=2,lwd=3,col="black")
  col2="#00000080"
  #col2="azure3"
  polygon(c(x[,1],x[nrow(x):1,1]),c(x[,3],x[nrow(x):1,4]),col=col2,density=-1,border=col2)
  i <- i + 1
}
legend("right", legend=c("15","20",  "50", "150"), col=colors_cam,lty=c(1), lwd=c(2),cex=0.8,inset=c(-0.3,0),
       bty = "n", xpd=TRUE)
dev.off()

png("constant_005_all.png",units="in", width=7, height=5, res=500)
par(xpd = T, mar = par()$mar + c(0,0,0,7))
plot(x[,1],x[,2],type="n",log="y",xlim=c(-150, -20),ylim=c(1e2,1e12),xlab="g (generations before present)",ylab="Ne (effective population size)",xaxs="i",yaxt="n",main=paste("Constant population MAF=0.05"), cex.lab=0.8, cex.axis=0.7)
axis(2,at=c(1e3,1e4,1e5,1e6,1e7, 1e8, 1e9,1e10,1e11,1e12),labels=c(expression(10^3),expression(10^4),expression(10^5),expression(10^6),expression(10^7), expression(10^8), expression(10^9), expression(10^10),
                                                                                                          expression(10^11), expression(10^12)),las=1, cex.axis=0.7)
i <- 1
for(pop in pops){
  filename=paste("constant.",pop,".fixed.maf0.05_OUT.ne",sep="")
  # only read rows from g=4 to g=150
  x=read.table(filename,header=T)[21:151,]
  x$GEN <- -x$GEN
  x$real <- 10000
  lines(x[,1],x[,2],las=1,lwd=3,col=colors_cam[i],)
  x$real <- 10000
  lines(x[,1],x[,5],las=1,lty=2,lwd=3,col="black")
  col2="#00000080"
  #col2="azure3"
  polygon(c(x[,1],x[nrow(x):1,1]),c(x[,3],x[nrow(x):1,4]),col=col2,density=-1,border=col2)
  i <- i + 1
}
legend("right", legend=c("15","20", "50", "150"), col=colors_cam,lty=c(1), lwd=c(2),cex=0.8,inset=c(-0.3,0),
       bty = "n", xpd=TRUE)
dev.off()

png("constant_01_550000.png",units="in", width=7, height=5, res=500)
par(xpd = T, mar = par()$mar + c(0,0,0,7))
plot(x[,1],x[,2],type="n",log="y",xlim=c(-150, -20),ylim=c(1e2,1e12),xlab="g (generations before present)",ylab="Ne (effective population size)",xaxs="i",yaxt="n",main=paste("Constant population MAF=0.1, 550,000SNP"), cex.lab=0.8, cex.axis=0.7)
axis(2,at=c(1e3,1e4,1e5,1e6,1e7, 1e8, 1e9,1e10,1e11,1e12),labels=c(expression(10^3),expression(10^4),expression(10^5),expression(10^6),expression(10^7), expression(10^8), expression(10^9), expression(10^10),
                                                                   expression(10^11), expression(10^12)),las=1, cex.axis=0.7)
i <- 1
for(pop in pops){
  filename=paste("constant.",pop,".maf0.1.ALL_sorted.550000_OUT.ne",sep="")
  # only read rows from g=4 to g=150
  x=read.table(filename,header=T)[21:151,]
  x$GEN <- -x$GEN
  lines(x[,1],x[,2],las=1,lwd=3,col=colors_cam[i],)
  x$real <- 10000
  lines(x[,1],x[,5],las=1,lty=2,lwd=3,col="black")
  col2="#00000080"
  #col2="azure3"
  polygon(c(x[,1],x[nrow(x):1,1]),c(x[,3],x[nrow(x):1,4]),col=col2,density=-1,border=col2)
  i <- i + 1
}
legend("right", legend=c("15","20", "50", "150"), col=colors_cam,lty=c(1), lwd=c(2),cex=0.8,inset=c(-0.3,0),
       bty = "n", xpd=TRUE)
dev.off()

png("constant_01_all.png",units="in", width=7, height=5, res=500)
par(xpd = T, mar = par()$mar + c(0,0,0,7))
plot(x[,1],x[,2],type="n",log="y",xlim=c(-150, -20),ylim=c(1e2,1e12),xlab="g (generations before present)",ylab="Ne (effective population size)",xaxs="i",yaxt="n",main=paste("Constant population MAF=0.1"), cex.lab=0.8, cex.axis=0.7)
axis(2,at=c(1e3,1e4,1e5,1e6,1e7, 1e8, 1e9,1e10,1e11,1e12),labels=c(expression(10^3),expression(10^4),expression(10^5),expression(10^6),expression(10^7), expression(10^8), expression(10^9), expression(10^10),
                                                                   expression(10^11), expression(10^12)),las=1, cex.axis=0.7)
i <- 1
for(pop in pops){
  filename=paste("constant.",pop,".fixed.maf0.1_OUT.ne",sep="")
  # only read rows from g=4 to g=150
  x=read.table(filename,header=T)[21:151,]
  x$GEN <- -x$GEN
  x$real <- 10000
  lines(x[,1],x[,2],las=1,lwd=3,col=colors_cam[i],)
  x$real <- 10000
  lines(x[,1],x[,5],las=1,lty=2,lwd=3,col="black")
  col2="#00000080"
  #col2="azure3"
  polygon(c(x[,1],x[nrow(x):1,1]),c(x[,3],x[nrow(x):1,4]),col=col2,density=-1,border=col2)
  i <- i + 1
}
legend("right", legend=c("15","20", "50", "150"), col=colors_cam,lty=c(1), lwd=c(2),cex=0.8,inset=c(-0.3,0),
       bty = "n", xpd=TRUE)
dev.off()


###########

setwd("/Users/Cecilia/Documents/PhD/African_Pygmies/Genetic_studies_pygmies/Final_Data/New genetics march 2021/Revision_october_Lane/sims/expansion_OUT")

r <- 0.03     #Assign the value of 0.1 to the object "r", or per-capita growth rate (discrete)
lambda <- 1+r
N0 <- 11109   #Assign "N0", or initial population size
nyears <- 150 #Assign  the number of time steps to simulate
years <- seq(from=0, to=nyears, by=1)
N <- numeric(nyears+1)
N[1] <- N0
for (i in 2:(nyears+1)){  # This for-loop will run through the line of code between the curly brackets {}. "i" is simply the name of a variable (you can use "j", or "k", instead -- any variable name will do). "i" changes each time the loop iterates; basically, it will increase by 1 each time the loop is run, starting at "2" up until the specified maximum number of loops "50+1".
  N[i] <- lambda*N[i-1]   # This takes the [i - 1] element of "N", multiplies that element by the value of lambda, then assigns that calculated result to the [i] element of "N".
}
# These are the log files for 10,20,50,100,200
#exp_plot <- c(740818.22068172,548811.6360940, 223130.16014843,  49787.06836786,11108.99653824)
#nyears <- 150 #Assign  the number of time steps to simulate
#years <- c(10,20,50,100,150) 
#lo <- loess(exp_plot~years, x.is.log=T)

#exp_plot<-predict(lo, seq(from=0, to=150, by=1), x.is.log=T)
#exp_plot <- exp_plot[-c(1:20)]
exp_plot <- N[-c(1:20)]


png("expansion_04_550000.png",units="in", width=7, height=5, res=500)
par(xpd = T, mar = par()$mar + c(0,0,0,7))
plot(x[,1],x[,2],type="n",log="y",xlim=c(-150, -20),ylim=c(1e2,1e12),xlab="g (generations before present)",ylab="Ne (effective population size)",xaxs="i",yaxt="n",main=paste("Expansion MAF=0.4, 550,000SNP"), cex.lab=0.8, cex.axis=0.7)
axis(2,at=c(1e3,1e4,1e5,1e6,1e7, 1e8, 1e9,1e10,1e11,1e12),labels=c(expression(10^3),expression(10^4),expression(10^5),expression(10^6),expression(10^7), expression(10^8), expression(10^9), expression(10^10),
                                                                                                          expression(10^11), expression(10^12)),las=1, cex.axis=0.7)
i <- 1
for(pop in pops){
  filename=paste("expansion.",pop,".maf0.4_sorted.550000_OUT.ne",sep="")
  # only read rows from g=4 to g=150
  x=read.table(filename,header=T)[21:151,]
  x$real <- rev(exp_plot)
  x$GEN <- -x$GEN
  lines(x[,1],x[,2],las=1,lwd=3,col=colors_cam[i],)
  lines(x$GEN, x$real,las=1,lty=2,lwd=3,col="black")
  col2="#00000080"
  #col2="azure3"
  polygon(c(x[,1],x[nrow(x):1,1]),c(x[,3],x[nrow(x):1,4]),col=col2,density=-1,border=col2)
  i <- i + 1
}
legend("right", legend=c("15", "20", "50", "150"), col=colors_cam,lty=c(1), lwd=c(2),cex=0.8,inset=c(-0.3,0),
       bty = "n", xpd=TRUE)
dev.off()

png("expansion_04_all.png",units="in", width=7, height=5, res=500)
par(xpd = T, mar = par()$mar + c(0,0,0,7))
plot(x[,1],x[,2],type="n",log="y",xlim=c(-150, -20),ylim=c(1e2,1e12),xlab="g (generations before present)",ylab="Ne (effective population size)",xaxs="i",yaxt="n",main=paste("Expansion MAF=0.4"), cex.lab=0.8, cex.axis=0.7)
axis(2,at=c(1e3,1e4,1e5,1e6,1e7, 1e8, 1e9,1e10,1e11,1e12),labels=c(expression(10^3),expression(10^4),expression(10^5),expression(10^6),expression(10^7), expression(10^8), expression(10^9), expression(10^10),
                                                                                                          expression(10^11), expression(10^12)),las=1, cex.axis=0.7)
i <- 1
for(pop in pops){
  filename=paste("expansion.",pop,".fixed.maf0.4_OUT.ne",sep="")
  # only read rows from g=4 to g=150
  x=read.table(filename,header=T)[21:151,]
  x$GEN <- -x$GEN
  lines(x[,1],x[,2],las=1,lwd=3,col=colors_cam[i],)
  x$real <- rev(exp_plot)
  lines(x$GEN, x$real,las=1,lty=2,lwd=3,col="black")
  col2="#00000080"
  #col2="azure3"
  polygon(c(x[,1],x[nrow(x):1,1]),c(x[,3],x[nrow(x):1,4]),col=col2,density=-1,border=col2)
  i <- i + 1
}
legend("right", legend=c("15", "20","50", "150"), col=colors_cam,lty=c(1), lwd=c(2),cex=0.8,inset=c(-0.3,0),
       bty = "n", xpd=TRUE)
dev.off()

png("expansion_005_550000.png",units="in", width=7, height=5, res=500)
par(xpd = T, mar = par()$mar + c(0,0,0,7))
plot(x[,1],x[,2],type="n",log="y",xlim=c(-150, -20),ylim=c(1e2,1e12),xlab="g (generations before present)",ylab="Ne (effective population size)",xaxs="i",yaxt="n",main=paste("Expansion MAF=0.05, 550,000SNP"), cex.lab=0.8, cex.axis=0.7)
axis(2,at=c(1e3,1e4,1e5,1e6,1e7, 1e8, 1e9,1e10,1e11,1e12),labels=c(expression(10^3),expression(10^4),expression(10^5),expression(10^6),expression(10^7), expression(10^8), expression(10^9), expression(10^10),
                                                                                                          expression(10^11), expression(10^12)),las=1, cex.axis=0.7)
i <- 1
for(pop in pops){
  filename=paste("expansion.",pop,".maf0.05_sorted.550000_OUT.ne",sep="")
  # only read rows from g=4 to g=150
  x=read.table(filename,header=T)[21:151,]
  x$GEN <- -x$GEN
  x$real <- rev(exp_plot)
  lines(x$GEN, x$real,las=1,lty=2,lwd=3,col="black")
  lines(x[,1],x[,2],las=1,lwd=3,col=colors_cam[i],)
  col2="#00000080"
  #col2="azure3"
  polygon(c(x[,1],x[nrow(x):1,1]),c(x[,3],x[nrow(x):1,4]),col=col2,density=-1,border=col2)
  i <- i + 1
}
legend("right", legend=c("15", "20", "50", "150"), col=colors_cam,lty=c(1), lwd=c(2),cex=0.8,inset=c(-0.3,0),
       bty = "n", xpd=TRUE)
dev.off()

png("expansion_005_all.png",units="in", width=7, height=5, res=500)
par(xpd = T, mar = par()$mar + c(0,0,0,7))
plot(x[,1],x[,2],type="n",log="y",xlim=c(-150, -20),ylim=c(1e2,1e12),xlab="g (generations before present)",ylab="Ne (effective population size)",xaxs="i",yaxt="n",main=paste("Expansion MAF=0.05"), cex.lab=0.8, cex.axis=0.7)
axis(2,at=c(1e3,1e4,1e5,1e6,1e7, 1e8, 1e9,1e10,1e11,1e12),labels=c(expression(10^3),expression(10^4),expression(10^5),expression(10^6),expression(10^7), expression(10^8), expression(10^9), expression(10^10),
                                                                                                          expression(10^11), expression(10^12)),las=1, cex.axis=0.7)
i <- 1
for(pop in pops){
  filename=paste("expansion.",pop,".fixed.maf0.05_OUT.ne",sep="")
  # only read rows from g=4 to g=150
  x=read.table(filename,header=T)[21:151,]
  x$GEN <- -x$GEN
  lines(x[,1],x[,2],las=1,lwd=3,col=colors_cam[i],)
  x$real <- rev(exp_plot)
  lines(x$GEN, x$real,las=1,lty=2,lwd=3,col="black")
  col2="#00000080"
  #col2="azure3"
  polygon(c(x[,1],x[nrow(x):1,1]),c(x[,3],x[nrow(x):1,4]),col=col2,density=-1,border=col2)
  i <- i + 1
}
legend("right", legend=c("15","20", "50", "150"), col=colors_cam,lty=c(1), lwd=c(2),cex=0.8,inset=c(-0.3,0),
       bty = "n", xpd=TRUE)
dev.off()

png("expansion_01_550000.png",units="in", width=7, height=5, res=500)
par(xpd = T, mar = par()$mar + c(0,0,0,7))
plot(x[,1],x[,2],type="n",log="y",xlim=c(-150, -20),ylim=c(1e2,1e12),xlab="g (generations before present)",ylab="Ne (effective population size)",xaxs="i",yaxt="n",main=paste("Expansion MAF=0.1, 550,000SNP"), cex.lab=0.8, cex.axis=0.7)
axis(2,at=c(1e3,1e4,1e5,1e6,1e7, 1e8, 1e9,1e10,1e11,1e12),labels=c(expression(10^3),expression(10^4),expression(10^5),expression(10^6),expression(10^7), expression(10^8), expression(10^9), expression(10^10),
                                                                   expression(10^11), expression(10^12)),las=1, cex.axis=0.7)
i <- 1
for(pop in pops){
  filename=paste("expansion.",pop,".maf0.1_sorted.550000_OUT.ne",sep="")
  # only read rows from g=4 to g=150
  x=read.table(filename,header=T)[21:151,]
  x$real <- rev(exp_plot)
  x$GEN <- -x$GEN
  lines(x[,1],x[,2],las=1,lwd=3,col=colors_cam[i],)
  lines(x$GEN, x$real,las=1,lty=2,lwd=3,col="black")
  col2="#00000080"
  #col2="azure3"
  polygon(c(x[,1],x[nrow(x):1,1]),c(x[,3],x[nrow(x):1,4]),col=col2,density=-1,border=col2)
  i <- i + 1
}
legend("right", legend=c("15", "20", "50", "150"), col=colors_cam,lty=c(1), lwd=c(2),cex=0.8,inset=c(-0.3,0),
       bty = "n", xpd=TRUE)
dev.off()

png("expansion_01_all.png",units="in", width=7, height=5, res=500)
par(xpd = T, mar = par()$mar + c(0,0,0,7))
plot(x[,1],x[,2],type="n",log="y",xlim=c(-150, -20),ylim=c(1e2,1e12),xlab="g (generations before present)",ylab="Ne (effective population size)",xaxs="i",yaxt="n",main=paste("Expansion MAF=0.1"), cex.lab=0.8, cex.axis=0.7)
axis(2,at=c(1e3,1e4,1e5,1e6,1e7, 1e8, 1e9,1e10,1e11,1e12),labels=c(expression(10^3),expression(10^4),expression(10^5),expression(10^6),expression(10^7), expression(10^8), expression(10^9), expression(10^10),
                                                                   expression(10^11), expression(10^12)),las=1, cex.axis=0.7)
i <- 1
for(pop in pops){
  filename=paste("expansion.",pop,".fixed.maf0.1_OUT.ne",sep="")
  # only read rows from g=4 to g=150
  x=read.table(filename,header=T)[21:151,]
  x$GEN <- -x$GEN
  lines(x[,1],x[,2],las=1,lwd=3,col=colors_cam[i],)
  x$real <- rev(exp_plot)
  lines(x$GEN, x$real,las=1,lty=2,lwd=3,col="black")
  col2="#00000080"
  #col2="azure3"
  polygon(c(x[,1],x[nrow(x):1,1]),c(x[,3],x[nrow(x):1,4]),col=col2,density=-1,border=col2)
  i <- i + 1
}
legend("right", legend=c("15", "20","50", "150"), col=colors_cam,lty=c(1), lwd=c(2),cex=0.8,inset=c(-0.3,0),
       bty = "n", xpd=TRUE)
dev.off()

###########

setwd("/Users/Cecilia/Documents/PhD/African_Pygmies/Genetic_studies_pygmies/Final_Data/New genetics march 2021/Revision_october_Lane/sims/collapse_OUT")


r <- 0.03     #Assign the value of 0.1 to the object "r", or per-capita growth rate (discrete)
lambda <- 1+r
N0 <- 500  #Assign "N0", or initial population size
nyears <- 150 #Assign  the number of time steps to simulate
years <- seq(from=0, to=nyears, by=1)
N <- numeric(nyears+1)
N[1] <- N0
for (i in 2:(nyears+1)){  # This for-loop will run through the line of code between the curly brackets {}. "i" is simply the name of a variable (you can use "j", or "k", instead -- any variable name will do). "i" changes each time the loop iterates; basically, it will increase by 1 each time the loop is run, starting at "2" up until the specified maximum number of loops "50+1".
  N[i] <- lambda*N[i-1]   # This takes the [i - 1] element of "N", multiplies that element by the value of lambda, then assigns that calculated result to the [i] element of "N".
}
coll_plot <- N[-c(1:20)]

png("collapse_04_550000.png",units="in", width=7, height=5, res=500)
par(xpd = T, mar = par()$mar + c(0,0,0,7))
plot(x[,1],x[,2],type="n",log="y",xlim=c(-150, -20),ylim=c(1e2,1e12),xlab="g (generations before present)",ylab="Ne (effective population size)",xaxs="i",yaxt="n",main=paste("Collapse MAF=0.4, 550,000SNP"), cex.lab=0.8, cex.axis=0.7)
axis(2,at=c(1e3,1e4,1e5,1e6,1e7, 1e8, 1e9,1e10,1e11,1e12),labels=c(expression(10^3),expression(10^4),expression(10^5),expression(10^6),expression(10^7), expression(10^8), expression(10^9), expression(10^10),
                                                                                                          expression(10^11), expression(10^12)),las=1, cex.axis=0.7)
i <- 1
for(pop in pops){
  filename=paste("collapse.",pop,".maf0.4_sorted.550000_OUT.ne",sep="")
  # only read rows from g=4 to g=150
  x=read.table(filename,header=T)[21:151,]
  x$GEN <- -x$GEN
  lines(x[,1],x[,2],las=1,lwd=3,col=colors_cam[i],)
  x$real <- coll_plot
  lines(x$GEN, x$real,las=1,lty=2,lwd=3,col="black")
  col2="#00000080"
  #col2="azure3"
  polygon(c(x[,1],x[nrow(x):1,1]),c(x[,3],x[nrow(x):1,4]),col=col2,density=-1,border=col2)
  i <- i + 1
}
legend("right", legend=c("15", "20", "50", "150"), col=colors_cam,lty=c(1), lwd=c(2),cex=0.8,inset=c(-0.3,0),
       bty = "n", xpd=TRUE)
dev.off()

png("collapse_04_all.png",units="in", width=7, height=5, res=500)
par(xpd = T, mar = par()$mar + c(0,0,0,7))
plot(x[,1],x[,2],type="n",log="y",xlim=c(-150, -20),ylim=c(1e2,1e12),xlab="g (generations before present)",ylab="Ne (effective population size)",xaxs="i",yaxt="n",main=paste("Collapse MAF=0.4"), cex.lab=0.8, cex.axis=0.7)
axis(2,at=c(1e3,1e4,1e5,1e6,1e7, 1e8, 1e9,1e10,1e11,1e12),labels=c(expression(10^3),expression(10^4),expression(10^5),expression(10^6),expression(10^7), expression(10^8), expression(10^9), expression(10^10),
                                                                                                          expression(10^11), expression(10^12)),las=1, cex.axis=0.7)
i <- 1
for(pop in pops){
  filename=paste("collapse.",pop,".fixed.maf0.4_OUT.ne",sep="")
  # only read rows from g=4 to g=150
  x=read.table(filename,header=T)[21:151,]
  x$GEN <- -x$GEN
  lines(x[,1],x[,2],las=1,lwd=3,col=colors_cam[i],)
  x$real <- coll_plot
  lines(x$GEN, x$real,las=1,lty=2,lwd=3,col="black")
  col2="#00000080"
  #col2="azure3"
  polygon(c(x[,1],x[nrow(x):1,1]),c(x[,3],x[nrow(x):1,4]),col=col2,density=-1,border=col2)
  i <- i + 1
}
legend("right", legend=c("15", "20","50", "150"), col=colors_cam,lty=c(1), lwd=c(2),cex=0.8,inset=c(-0.3,0),
       bty = "n", xpd=TRUE)
dev.off()

png("collapse_005_550000.png",units="in", width=7, height=5, res=500)
par(xpd = T, mar = par()$mar + c(0,0,0,7))
plot(x[,1],x[,2],type="n",log="y",xlim=c(-150, -20),ylim=c(1e2,1e12),xlab="g (generations before present)",ylab="Ne (effective population size)",xaxs="i",yaxt="n",main=paste("Collapse MAF=0.05, 550,000SNP"), cex.lab=0.8, cex.axis=0.7)
axis(2,at=c(1e3,1e4,1e5,1e6,1e7, 1e8, 1e9,1e10,1e11,1e12),labels=c(expression(10^3),expression(10^4),expression(10^5),expression(10^6),expression(10^7), expression(10^8), expression(10^9), expression(10^10),
                                                                                                          expression(10^11), expression(10^12)),las=1, cex.axis=0.7)
i <- 1
for(pop in pops){
  filename=paste("collapse.",pop,".maf0.05_sorted.550000_OUT.ne",sep="")
  # only read rows from g=4 to g=150
  x=read.table(filename,header=T)[21:151,]
  x$GEN <- -x$GEN
  lines(x[,1],x[,2],las=1,lwd=3,col=colors_cam[i],)
  x$real <- coll_plot
  lines(x$GEN, x$real,las=1,lty=2,lwd=3,col="black")
  col2="#00000080"
  #col2="azure3"
  polygon(c(x[,1],x[nrow(x):1,1]),c(x[,3],x[nrow(x):1,4]),col=col2,density=-1,border=col2)
  i <- i + 1
}
legend("right", legend=c("15", "20", "50", "150"), col=colors_cam,lty=c(1), lwd=c(2),cex=0.8,inset=c(-0.3,0),
       bty = "n", xpd=TRUE)
dev.off()

png("collapse_005_all.png",units="in", width=7, height=5, res=500)
par(xpd = T, mar = par()$mar + c(0,0,0,7))
plot(x[,1],x[,2],type="n",log="y",xlim=c(-150, -20),ylim=c(1e2,1e20),xlab="g (generations before present)",ylab="Ne (effective population size)",xaxs="i",yaxt="n",main=paste("Collapse MAF=0.05"), cex.lab=0.8, cex.axis=0.7)
axis(2,at=c(1e3,1e4,1e5,1e6,1e7, 1e8, 1e9,1e10,1e11,1e12),labels=c(expression(10^3),expression(10^4),expression(10^5),expression(10^6),expression(10^7), expression(10^8), expression(10^9), expression(10^10),
                                                                                                          expression(10^11), expression(10^12)),las=1, cex.axis=0.7)
i <- 2
for(pop in pops[-1]){
  filename=paste("collapse.",pop,".fixed.maf0.05_OUT.ne",sep="")
  # only read rows from g=4 to g=150
  x=read.table(filename,header=T)[21:151,]
  x$GEN <- -x$GEN
  lines(x[,1],x[,2],las=1,lwd=3,col=colors_cam[i],)
  x$real <- coll_plot
  lines(x$GEN, x$real,las=1,lty=2,lwd=3,col="black")
  col2="#00000080"
  #col2="azure3"
  polygon(c(x[,1],x[nrow(x):1,1]),c(x[,3],x[nrow(x):1,4]),col=col2,density=-1,border=col2)
  i <- i + 1
}
legend("right", legend=c("20", "50", "150"), col=colors_cam[-1],lty=c(1), lwd=c(2),cex=0.8,inset=c(-0.3,0),
       bty = "n", xpd=TRUE)
dev.off()

png("collapse_01_550000.png",units="in", width=7, height=5, res=500)
par(xpd = T, mar = par()$mar + c(0,0,0,7))
plot(x[,1],x[,2],type="n",log="y",xlim=c(-150, -20),ylim=c(1e2,1e12),xlab="g (generations before present)",ylab="Ne (effective population size)",xaxs="i",yaxt="n",main=paste("Collapse MAF=0.1, 550,000SNP"), cex.lab=0.8, cex.axis=0.7)
axis(2,at=c(1e3,1e4,1e5,1e6,1e7, 1e8, 1e9,1e10,1e11,1e12),labels=c(expression(10^3),expression(10^4),expression(10^5),expression(10^6),expression(10^7), expression(10^8), expression(10^9), expression(10^10),
                                                                   expression(10^11), expression(10^12)),las=1, cex.axis=0.7)
i <- 1
for(pop in pops){
  filename=paste("collapse.",pop,".maf0.1_sorted.550000_OUT.ne",sep="")
  # only read rows from g=4 to g=150
  x=read.table(filename,header=T)[21:151,]
  x$GEN <- -x$GEN
  lines(x[,1],x[,2],las=1,lwd=3,col=colors_cam[i],)
  x$real <- coll_plot
  lines(x$GEN, x$real,las=1,lty=2,lwd=3,col="black")
  col2="#00000080"
  #col2="azure3"
  polygon(c(x[,1],x[nrow(x):1,1]),c(x[,3],x[nrow(x):1,4]),col=col2,density=-1,border=col2)
  i <- i + 1
}
legend("right", legend=c("15", "20", "50", "150"), col=colors_cam,lty=c(1), lwd=c(2),cex=0.8,inset=c(-0.3,0),
       bty = "n", xpd=TRUE)
dev.off()

png("collapse_01_all.png",units="in", width=7, height=5, res=500)
par(xpd = T, mar = par()$mar + c(0,0,0,7))
plot(x[,1],x[,2],type="n",log="y",xlim=c(-150, -20),ylim=c(1e2,1e12),xlab="g (generations before present)",ylab="Ne (effective population size)",xaxs="i",yaxt="n",main=paste("Collapse MAF=0.1"), cex.lab=0.8, cex.axis=0.7)
axis(2,at=c(1e3,1e4,1e5,1e6,1e7, 1e8, 1e9,1e10,1e11,1e12),labels=c(expression(10^3),expression(10^4),expression(10^5),expression(10^6),expression(10^7), expression(10^8), expression(10^9), expression(10^10),
                                                                   expression(10^11), expression(10^12)),las=1, cex.axis=0.7)
i <- 1
for(pop in pops){
  filename=paste("collapse.",pop,".fixed.maf0.1_OUT.ne",sep="")
  # only read rows from g=4 to g=150
  x=read.table(filename,header=T)[21:151,]
  x$GEN <- -x$GEN
  lines(x[,1],x[,2],las=1,lwd=3,col=colors_cam[i],)
  x$real <- coll_plot
  lines(x$GEN, x$real,las=1,lty=2,lwd=3,col="black")
  col2="#00000080"
  #col2="azure3"
  polygon(c(x[,1],x[nrow(x):1,1]),c(x[,3],x[nrow(x):1,4]),col=col2,density=-1,border=col2)
  i <- i + 1
}
legend("right", legend=c("15", "20","50", "150"), col=colors_cam,lty=c(1), lwd=c(2),cex=0.8,inset=c(-0.3,0),
       bty = "n", xpd=TRUE)
dev.off()





















png("constant_nopoly.png",units="in", width=7, height=5, res=500)
par(xpd = T, mar = par()$mar + c(0,0,0,7))
plot(x[,1],x[,2],type="n",log="y",xlim=c(-150, -4),ylim=c(1e2,1e20),xlab="g (generations before present)",ylab="Ne (effective population size)",xaxs="i",yaxt="n",main=paste("collapse"), cex.lab=0.8, cex.axis=0.7)
axis(2,at=c(1e3,1e4,1e5,1e6,1e7, 1e8, 1e9,1e10,1e11,1e12,1e13,1e14,1e15, 1e16, 1e17, 1e18, 1e19),labels=c(expression(10^3),expression(10^4),expression(10^5),expression(10^6),expression(10^7), expression(10^8), expression(10^9), expression(10^10),
                                                                                                          expression(10^11), expression(10^12), expression(10^13), expression(10^14), expression(10^15), expression(10^16),
                                                                                                          expression(10^17), expression(10^18), expression(10^19)),las=1, cex.axis=0.7)
i <- 1
for(pop in pops){
  filename=paste("constant.",pop,".fixed.maf0.05_OUT.ne",sep="")
  # only read rows from g=4 to g=150
  x=read.table(filename,header=T)[5:151,]
  x$GEN <- -x$GEN
  lines(x[,1],x[,2],las=1,lwd=3,col=colors_cam[i],)
  col2="#00000080"
  #col2="azure3"
  #polygon(c(x[,1],x[nrow(x):1,1]),c(x[,3],x[nrow(x):1,4]),col=col2,density=-1,border=col2)
  i <- i + 1
}
legend("right", legend=c("15", "50", "150"), col=colors_cam,lty=c(1), lwd=c(2),cex=0.8,inset=c(-0.3,0),
       bty = "n", xpd=TRUE)
dev.off()

png("expansion_nopoly.png",units="in", width=7, height=5, res=500)
par(xpd = T, mar = par()$mar + c(0,0,0,7))
plot(x[,1],x[,2],type="n",log="y",xlim=c(-150, -4),ylim=c(1e2,1e20),xlab="g (generations before present)",ylab="Ne (effective population size)",xaxs="i",yaxt="n",main=paste("Expansion"), cex.lab=0.8, cex.axis=0.7)
axis(2,at=c(1e3,1e4,1e5,1e6,1e7, 1e8, 1e9,1e10,1e11,1e12,1e13,1e14,1e15, 1e16, 1e17, 1e18, 1e19),labels=c(expression(10^3),expression(10^4),expression(10^5),expression(10^6),expression(10^7), expression(10^8), expression(10^9), expression(10^10),
                                                                                                          expression(10^11), expression(10^12), expression(10^13), expression(10^14), expression(10^15), expression(10^16),
                                                                                                          expression(10^17), expression(10^18), expression(10^19)),las=1, cex.axis=0.7)
i <- 1
for(pop in pops){
  filename=paste("expansion.",pop,"_OUT.ne",sep="")
  # only read rows from g=4 to g=150
  x=read.table(filename,header=T)[5:151,]
  x$GEN <- rev(x$GEN) # expansion is reversed
  x$GEN <- -x$GEN
  lines(x[,1],x[,2],las=1,lwd=3,col=colors_cam[i],)
  col2="#00000080"
  #col2="azure3"
  #polygon(c(x[,1],x[nrow(x):1,1]),c(x[,3],x[nrow(x):1,4]),col=col2,density=-1,border=col2)
  i <- i + 1
}
legend("right", legend=c("15", "50", "150"), col=colors_cam,lty=c(1), lwd=c(2),cex=0.8,inset=c(-0.3,0),
       bty = "n", xpd=TRUE)
dev.off()

##########################################
## Plot the downsampled BIAKA from HGDP ##
##########################################

setwd("/Users/Cecilia/Documents/PhD/African_Pygmies/Genetic_studies_pygmies/Final_Data/New genetics march 2021/Revision_october_Lane/biaka_hgdp")

biaka <- c("Biaka_hgdp_sorted_maf0.1_550k_OUT.ne", "Biaka_hgdp_sorted_maf0.1_750k_OUT.ne",
           "Biaka_hgdp_sorted_maf0.1_1mil_OUT.ne")

colors_biaka <-  c(holy_mountain(2)[3], holy_mountain(2)[4],
                   holy_mountain(2)[5])

png("hgdp_biaka_ds.png",units="in", width=7, height=5, res=500)
par(xpd = T, mar = par()$mar + c(0,0,0,7))
plot(x[,1],x[,2],type="n",log="y",xlim=c(-150, -4),ylim=c(1e2,1e10),xlab="g (generations before present)",ylab="Ne (effective population size)",xaxs="i",yaxt="n",main=paste("Biaka HGDP"), cex.lab=0.8, cex.axis=0.7)
axis(2,at=c(1e3,1e4,1e5,1e6,1e7, 1e8, 1e9, 1e10),labels=c(expression(10^3),expression(10^4),expression(10^5),expression(10^6),expression(10^7), expression(10^8), expression(10^9), expression(10^10)),las=1, cex.axis=0.7)
i <- 1
for(anc in biaka){
  filename=anc
  # only read rows from g=4 to g=150
  x=read.table(filename,header=T)[5:151,]
  x$GEN <- -x$GEN
  lines(x[,1],x[,2],las=1,lwd=2,col=colors_biaka[i],)
  i <- i + 1
}
legend("right", legend=c("550,000 SNPs", "750,000 SNPs", "1,000,000 SNPs"), col=colors_biaka,lty=c(1), lwd=c(2),cex=0.8,inset=c(-0.3,0),
       bty = "n", xpd=TRUE)
dev.off()

## MAF 0.4

biaka <- c("Biaka_hgdp_sorted_maf0.4_550k_OUT.ne", "Biaka_hgdp_sorted_maf0.4_750k_OUT.ne",
           "Biaka_hgdp_sorted_maf0.4_1mil_OUT.ne")

colors_biaka <-  c(holy_mountain(2)[3], holy_mountain(2)[4],
                   holy_mountain(2)[5])

png("hgdp_biaka_ds_maf04.png",units="in", width=7, height=5, res=500)
par(xpd = T, mar = par()$mar + c(0,0,0,7))
plot(x[,1],x[,2],type="n",log="y",xlim=c(-150, -4),ylim=c(1e2,1e10),xlab="g (generations before present)",ylab="Ne (effective population size)",xaxs="i",yaxt="n",main=paste("Biaka HGDP"), cex.lab=0.8, cex.axis=0.7)
axis(2,at=c(1e3,1e4,1e5,1e6,1e7, 1e8, 1e9, 1e10),labels=c(expression(10^3),expression(10^4),expression(10^5),expression(10^6),expression(10^7), expression(10^8), expression(10^9), expression(10^10)),las=1, cex.axis=0.7)
i <- 1
for(anc in biaka){
  filename=anc
  # only read rows from g=4 to g=150
  x=read.table(filename,header=T)[5:151,]
  x$GEN <- -x$GEN
  lines(x[,1],x[,2],las=1,lwd=2,col=colors_biaka[i],)
  i <- i + 1
}
legend("right", legend=c("550,000 SNPs", "750,000 SNPs", "1,000,000 SNPs"), col=colors_biaka,lty=c(1), lwd=c(2),cex=0.8,inset=c(-0.3,0),
       bty = "n", xpd=TRUE)
dev.off()

#############################
## Plot the regional IBDNE ##
#############################

setwd("/Users/Cecilia/Documents/PhD/African_Pygmies/Genetic_studies_pygmies/Final_Data/New genetics march 2021/Revision_october_Lane/regional_out")

regions <- c("Northwestern_HG_OUT.ne", "Western_HG_OUT.ne",
           "Eastern_HG_OUT.ne")

colors_biaka <-  c(holy_mountain(2)[3], holy_mountain(2)[4],
                   holy_mountain(2)[5])

png("regions.png",units="in", width=7, height=5, res=500)
par(xpd = T, mar = par()$mar + c(0,0,0,7))
plot(x[,1],x[,2],type="n",log="y",xlim=c(-150, -4),ylim=c(1e2,1e10),xlab="g (generations before present)",ylab="Ne (effective population size)",xaxs="i",yaxt="n",main=paste("Regional trends"), cex.lab=0.8, cex.axis=0.7)
axis(2,at=c(1e3,1e4,1e5,1e6,1e7, 1e8, 1e9, 1e10),labels=c(expression(10^3),expression(10^4),expression(10^5),expression(10^6),expression(10^7), expression(10^8), expression(10^9), expression(10^10)),las=1, cex.axis=0.7)
i <- 1
for(anc in regions){
  filename=anc
  # only read rows from g=4 to g=150
  x=read.table(filename,header=T)[5:151,]
  x$GEN <- -x$GEN
  lines(x[,1],x[,2],las=1,lwd=2,col=colors_biaka[i],)
  i <- i + 1
}
legend("right", legend=c("Northwestern HG", "Western HG", "Eastern HG"), col=colors_biaka,lty=c(1), lwd=c(2),cex=0.8,inset=c(-0.4,0),
       bty = "n", xpd=TRUE)
dev.off()


