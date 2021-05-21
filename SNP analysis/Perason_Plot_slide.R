setwd("~/Dropbox/popbuild/Mark_Scripts/AFS_Het_Ne")
library(reshape2)
library(stringr)
library(tidyverse)
library(GenWin)
par(mfrow=c(6,1),mar=c(2,5,0.5,2),mgp=c(3,1,0), oma=c(3,4,0,0))

#candidate genes cubillios
candidate <- read.table("Candidate_regs_cubillos.txt", header= TRUE)
candidate$start <- candidate$start/1e6
candidate$end <- candidate$end/1e6

#candidate genes Burke 2014
candidate_2 <- read.table("candidate_burke2014.txt", header= TRUE)
candidate_2$start <- candidate_2$start/1e6
candidate_2$end <- candidate_2$end/1e6

#4 Founders
B4 <- read.table("B4_Pearson.txt", header = TRUE)

B4$chr <- as.character(B4$chr)
replacement <- c('C01','C02','C03','C04','C05','C06','C07','C08','C09','C10','C11','C12','C13','C14','C15','C16')
old <- as.vector(unique(B4$chr))
for (i in 1:16){
  
  B4[B4==old[i]] <- replacement[i]
}

#S4
B4_S <- na.omit(B4[,1:3])
chr_index <- unique(B4_S$chr)

B4_S_final <- data.frame()
for (i in 1:16){
  
  temp <- subset(B4_S, chr == chr_index[i])
  
  temp2 <- splineAnalyze(-log10(temp$S4_p),temp$pos,smoothness=2000,
                        plotRaw=FALSE,plotWindows=FALSE,method=3 )
  
  p <- temp2[["windowData"]][["MeanY"]]
  pos <- temp2[["windowData"]][["WindowStop"]]
  chr <- rep( chr_index[i], length(p))
  temp3 <- data.frame(chr,pos,p )
  B4_S_final <- rbind(B4_S_final, temp3)
}



#Make a genome-wide scale so the entire genome can go on a single plot
#Add the max pos of C01 to C02, C01+C02 to C03, etc
#Asssuming you are starting with a data frame called xx with one colname = pos, and another colname = var1 (some variable you want to plot)
Gaxis <- numeric(length=0)
chrs<-c('C01','C02','C03','C04','C05','C06','C07','C08','C09','C10','C11','C12','C13','C14','C15','C16')

B4_S_final <- arrange(transform(B4_S_final,
                              chr=factor(chr,levels=chrs)),chr)

for (k in 1:length(chrs)){
  data.samp<-subset(B4_S_final,chr==chrs[k])
  if (chrs[k]=='C01'){Gaxis.samp<-data.samp$pos}
  if (chrs[k]=='C02'){Gaxis.samp<-data.samp$pos+230218}
  if (chrs[k]=='C03'){Gaxis.samp<-data.samp$pos+230218+813184}
  if (chrs[k]=='C04'){Gaxis.samp<-data.samp$pos+230218+813184+316620}
  if (chrs[k]=='C05'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933}
  if (chrs[k]=='C06'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874}
  if (chrs[k]=='C07'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161}
  if (chrs[k]=='C08'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940}
  if (chrs[k]=='C09'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643}
  if (chrs[k]=='C10'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888}
  if (chrs[k]=='C11'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751}
  if (chrs[k]=='C12'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816} 
  if (chrs[k]=='C13'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816+1078177}
  if (chrs[k]=='C14'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816+1078177+924431}
  if (chrs[k]=='C15'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816+1078177+924431+784333}
  if (chrs[k]=='C16'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816+1078177+924431+784333+1091291}
  Gaxis<-c(Gaxis,Gaxis.samp) 
}
MB=Gaxis/1e6

plot(MB,B4_S_final$p,
     xlab="",
     ylab="",
     main="",
     ylim=c(0,30),
     cex.main=3, #size of text for main title
     cex.lab=2, #size of text for axis labels
     type="n",	#type = none: don't make plot yet, need to make rectangles first
     axes=FALSE  #we will do the axes by hand next
)

#Here draw are gray rectangles that delineate chromosomes
rect(0.23,-10,1.04,80,col="grey80",lty=0)  
rect(1.36,-10,2.89,80,col="grey80",lty=0)  
rect(3.46,-10,3.73,80,col="grey80",lty=0)
rect(4.82,-10,5.39,80,col="grey80",lty=0)  
rect(5.83,-10,6.57,80,col="grey80",lty=0) 
rect(7.24,-10,8.32,80,col="grey80",lty=0)
rect(9.24,-10,10.03,80,col="grey80",lty=0)
rect(11.12,-10,12.07,80,col="grey80",lty=0) 

#Here draw points or lines.  Often it is useful to play around with point type, font size, color...
points(MB,B4_S_final$p, col="black", pch = 20, cex = 0.60, type = "l")
box() 
#Now draw axes back in (you have more flexibility this way)
axis(2,cex.axis=1.5, las = 1)

####Candidate regions from Cubillos et al. 
rect(candidate$start,25,candidate$end,27,col="red",lty=0) #red

####Candidate regions from Burke 2014 
rect(candidate_2$start,25,candidate_2$end,27,col="#009E73",lty=0) #green

#99.9 quantile
#abline(h=quantile(B4_S_final$p,probs=.999, na.rm=TRUE),v=NA, col="red")

legend("right",inset=-.05, "S4",bty="n", cex=2)

legend("left",inset=-.080, "(A)",bty="n", cex=2)

#K4

B4_K <- na.omit(B4[,c(1,2,4)])
chr_index <- unique(B4_K$chr)

B4_K_final <- data.frame()
for (i in 1:16){
  
  temp <- subset(B4_K, chr == chr_index[i])
  
  temp2 <- splineAnalyze(-log10(temp$K4_p),temp$pos,smoothness=2000,
                         plotRaw=FALSE,plotWindows=FALSE,method=3 )
  
  p <- temp2[["windowData"]][["MeanY"]]
  pos <- temp2[["windowData"]][["WindowStop"]]
  chr <- rep( chr_index[i], length(p))
  temp3 <- data.frame(chr,pos,p )
  B4_K_final <- rbind(B4_K_final, temp3)
}

#Make a genome-wide scale so the entire genome can go on a single plot
#Add the max pos of C01 to C02, C01+C02 to C03, etc
#Asssuming you are starting with a data frame called xx with one colname = pos, and another colname = var1 (some variable you want to plot)
Gaxis <- numeric(length=0)
chrs<-c('C01','C02','C03','C04','C05','C06','C07','C08','C09','C10','C11','C12','C13','C14','C15','C16')

B4_K_final <- arrange(transform(B4_K_final,
                                chr=factor(chr,levels=chrs)),chr)

for (k in 1:length(chrs)){
  data.samp<-subset(B4_K_final,chr==chrs[k])
  if (chrs[k]=='C01'){Gaxis.samp<-data.samp$pos}
  if (chrs[k]=='C02'){Gaxis.samp<-data.samp$pos+230218}
  if (chrs[k]=='C03'){Gaxis.samp<-data.samp$pos+230218+813184}
  if (chrs[k]=='C04'){Gaxis.samp<-data.samp$pos+230218+813184+316620}
  if (chrs[k]=='C05'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933}
  if (chrs[k]=='C06'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874}
  if (chrs[k]=='C07'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161}
  if (chrs[k]=='C08'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940}
  if (chrs[k]=='C09'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643}
  if (chrs[k]=='C10'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888}
  if (chrs[k]=='C11'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751}
  if (chrs[k]=='C12'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816} 
  if (chrs[k]=='C13'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816+1078177}
  if (chrs[k]=='C14'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816+1078177+924431}
  if (chrs[k]=='C15'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816+1078177+924431+784333}
  if (chrs[k]=='C16'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816+1078177+924431+784333+1091291}
  Gaxis<-c(Gaxis,Gaxis.samp) 
}
MB=Gaxis/1e6

plot(MB,B4_K_final$p,
     xlab="",
     ylab="",
     main="",
     ylim=c(0,30),
     cex.main=3, #size of text for main title
     cex.lab=2, #size of text for axis labels
     type="n",	#type = none: don't make plot yet, need to make rectangles first
     axes=FALSE  #we will do the axes by hand next
)

#Here draw are gray rectangles that delineate chromosomes
rect(0.23,-10,1.04,80,col="grey80",lty=0)  
rect(1.36,-10,2.89,80,col="grey80",lty=0)  
rect(3.46,-10,3.73,80,col="grey80",lty=0)
rect(4.82,-10,5.39,80,col="grey80",lty=0)  
rect(5.83,-10,6.57,80,col="grey80",lty=0) 
rect(7.24,-10,8.32,80,col="grey80",lty=0)
rect(9.24,-10,10.03,80,col="grey80",lty=0)
rect(11.12,-10,12.07,80,col="grey80",lty=0) 

#Here draw points or lines.  Often it is useful to play around with point type, font size, color...
points(MB,B4_K_final$p, col="black",pch = 20, cex = 0.70, type="l")

####Candidate regions from Cubillos et al. 
rect(candidate$start,25,candidate$end,27,col="red",lty=0)

####Candidate regions from Burke 2014 
rect(candidate_2$start,25,candidate_2$end,27,col="#009E73",lty=0)

#99.9 quantile
#abline(h=quantile(B4_K_final$p,probs=.999, na.rm=TRUE),v=NA, col="red")


box() 
#Now draw axes back in (you have more flexibility this way)
axis(2,cex.axis=1.5, las = 1)

legend("right",inset=-.05, "K4",bty="n", cex=2)

legend("left",inset=-.08, "(B)",bty="n", cex=2)

#8 founders

B8 <- read.table("B8_Pearson.txt", header = TRUE)

B8$chr <- as.character(B8$chr)
replacement <- c('C01','C02','C03','C04','C05','C06','C07','C08','C09','C10','C11','C12','C13','C14','C15','C16')
old <- as.vector(unique(B8$chr))
for (i in 1:16){
  
  B8[B8==old[i]] <- replacement[i]
}

#S8
B8_S <- na.omit(B8[,1:3])
chr_index <- unique(B8_S$chr)

B8_S_final <- data.frame()
for (i in 1:16){
  
  temp <- subset(B8_S, chr == chr_index[i])
  
  temp2 <- splineAnalyze(-log10(temp$S8_p),temp$pos,smoothness=2000,
                         plotRaw=FALSE,plotWindows=FALSE,method=3 )
  
  p <- temp2[["windowData"]][["MeanY"]]
  pos <- temp2[["windowData"]][["WindowStop"]]
  chr <- rep( chr_index[i], length(p))
  temp3 <- data.frame(chr,pos,p )
  B8_S_final <- rbind(B8_S_final, temp3)
}


#Make a genome-wide scale so the entire genome can go on a single plot
#Add the max pos of C01 to C02, C01+C02 to C03, etc
#Asssuming you are starting with a data frame called xx with one colname = pos, and another colname = var1 (some variable you want to plot)
Gaxis <- numeric(length=0)
chrs<-c('C01','C02','C03','C04','C05','C06','C07','C08','C09','C10','C11','C12','C13','C14','C15','C16')

B8_S_final <- arrange(transform(B8_S_final,
                        chr=factor(chr,levels=chrs)),chr)

for (k in 1:length(chrs)){
  data.samp<-subset(B8_S_final,chr==chrs[k])
  if (chrs[k]=='C01'){Gaxis.samp<-data.samp$pos}
  if (chrs[k]=='C02'){Gaxis.samp<-data.samp$pos+230218}
  if (chrs[k]=='C03'){Gaxis.samp<-data.samp$pos+230218+813184}
  if (chrs[k]=='C04'){Gaxis.samp<-data.samp$pos+230218+813184+316620}
  if (chrs[k]=='C05'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933}
  if (chrs[k]=='C06'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874}
  if (chrs[k]=='C07'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161}
  if (chrs[k]=='C08'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940}
  if (chrs[k]=='C09'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643}
  if (chrs[k]=='C10'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888}
  if (chrs[k]=='C11'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751}
  if (chrs[k]=='C12'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816} 
  if (chrs[k]=='C13'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816+1078177}
  if (chrs[k]=='C14'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816+1078177+924431}
  if (chrs[k]=='C15'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816+1078177+924431+784333}
  if (chrs[k]=='C16'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816+1078177+924431+784333+1091291}
  Gaxis<-c(Gaxis,Gaxis.samp) 
}
MB=Gaxis/1e6

#S8
plot(MB,B8_S_final$p,
     xlab="",
     ylab="",
     main="",
     ylim=c(0,30),
     cex.main=3, #size of text for main title
     cex.lab=2, #size of text for axis labels
     type="n",	#type = none: don't make plot yet, need to make rectangles first
     axes=FALSE  #we will do the axes by hand next
)

#Here draw are gray rectangles that delineate chromosomes
rect(0.23,-10,1.04,80,col="grey80",lty=0)  
rect(1.36,-10,2.89,80,col="grey80",lty=0)  
rect(3.46,-10,3.73,80,col="grey80",lty=0)
rect(4.82,-10,5.39,80,col="grey80",lty=0)  
rect(5.83,-10,6.57,80,col="grey80",lty=0) 
rect(7.24,-10,8.32,80,col="grey80",lty=0)
rect(9.24,-10,10.03,80,col="grey80",lty=0)
rect(11.12,-10,12.07,80,col="grey80",lty=0) 

#Here draw points or lines.  Often it is useful to play around with point type, font size, color...
points(MB,B8_S_final$p, col="black", pch = 20, cex = 0.60,type="l")

####Candidate regions from Cubillos et al. 
rect(candidate$start,25,candidate$end,27,col="red",lty=0)

####Candidate regions from Burke 2014 
rect(candidate_2$start,25,candidate_2$end,27,col="#009E73",lty=0)

#99.9 quantile
#abline(h=quantile(B8_S_final$p,probs=.999, na.rm=TRUE),v=NA, col="red")

box() 
#Now draw axes back in (you have more flexibility this way)
axis(2,cex.axis=1.5, las = 1)

legend("right",inset=-.05, "S8",bty="n", cex=2)

legend("left",inset=-.08, "(C)",bty="n", cex=2)

#K8
B8_K <- na.omit(B8[,c(1,2,4)])
chr_index <- unique(B8_K$chr)

B8_K_final <- data.frame()
for (i in 1:16){
  
  temp <- subset(B8_K, chr == chr_index[i])
  
  temp2 <- splineAnalyze(-log10(temp$K8_p),temp$pos,smoothness=2000,
                         plotRaw=FALSE,plotWindows=FALSE,method=3 )
  
  p <- temp2[["windowData"]][["MeanY"]]
  pos <- temp2[["windowData"]][["WindowStop"]]
  chr <- rep( chr_index[i], length(p))
  temp3 <- data.frame(chr,pos,p )
  B8_K_final <- rbind(B8_K_final, temp3)
}


#Make a genome-wide scale so the entire genome can go on a single plot
#Add the max pos of C01 to C02, C01+C02 to C03, etc
#Asssuming you are starting with a data frame called xx with one colname = pos, and another colname = var1 (some variable you want to plot)
Gaxis <- numeric(length=0)
chrs<-c('C01','C02','C03','C04','C05','C06','C07','C08','C09','C10','C11','C12','C13','C14','C15','C16')

B8_K_final <- arrange(transform(B8_K_final,
                                chr=factor(chr,levels=chrs)),chr)

for (k in 1:length(chrs)){
  data.samp<-subset(B8_K_final,chr==chrs[k])
  if (chrs[k]=='C01'){Gaxis.samp<-data.samp$pos}
  if (chrs[k]=='C02'){Gaxis.samp<-data.samp$pos+230218}
  if (chrs[k]=='C03'){Gaxis.samp<-data.samp$pos+230218+813184}
  if (chrs[k]=='C04'){Gaxis.samp<-data.samp$pos+230218+813184+316620}
  if (chrs[k]=='C05'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933}
  if (chrs[k]=='C06'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874}
  if (chrs[k]=='C07'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161}
  if (chrs[k]=='C08'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940}
  if (chrs[k]=='C09'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643}
  if (chrs[k]=='C10'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888}
  if (chrs[k]=='C11'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751}
  if (chrs[k]=='C12'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816} 
  if (chrs[k]=='C13'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816+1078177}
  if (chrs[k]=='C14'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816+1078177+924431}
  if (chrs[k]=='C15'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816+1078177+924431+784333}
  if (chrs[k]=='C16'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816+1078177+924431+784333+1091291}
  Gaxis<-c(Gaxis,Gaxis.samp) 
}
MB=Gaxis/1e6


plot(MB,B8_K_final$p,
     xlab="",
     ylab="",
     main="",
     ylim=c(0,30),
     cex.main=3, #size of text for main title
     cex.lab=2, #size of text for axis labels
     type="n",	#type = none: don't make plot yet, need to make rectangles first
     axes=FALSE  #we will do the axes by hand next
)

#Here draw are gray rectangles that delineate chromosomes
rect(0.23,-10,1.04,80,col="grey80",lty=0)  
rect(1.36,-10,2.89,80,col="grey80",lty=0)  
rect(3.46,-10,3.73,80,col="grey80",lty=0)
rect(4.82,-10,5.39,80,col="grey80",lty=0)  
rect(5.83,-10,6.57,80,col="grey80",lty=0) 
rect(7.24,-10,8.32,80,col="grey80",lty=0)
rect(9.24,-10,10.03,80,col="grey80",lty=0)
rect(11.12,-10,12.07,80,col="grey80",lty=0) 

#Here draw points or lines.  Often it is useful to play around with point type, font size, color...
points(MB,B8_K_final$p, col="black",pch = 20, cex = 0.70,type="l")

####Candidate regions from Cubillos et al. 
rect(candidate$start,25,candidate$end,27,col="red",lty=0)

####Candidate regions from Burke 2014 
rect(candidate_2$start,25,candidate_2$end,27,col="#009E73",lty=0)

#99.9 quantile
#abline(h=quantile(B8_K_final$p,probs=.999, na.rm=TRUE),v=NA, col="red")


box() 
#Now draw axes back in (you have more flexibility this way)
axis(2,cex.axis=1.5, las = 1)

legend("right",inset=-.05, "K8",bty="n", cex=2)

legend("left",inset=-.08, "(D)",bty="n", cex=2)



#12 founders

B12 <- read.table("B12_Pearson.txt", header = TRUE)

B12$chr <- as.character(B12$chr)
replacement <- c('C01','C02','C03','C04','C05','C06','C07','C08','C09','C10','C11','C12','C13','C14','C15','C16')
old <- as.vector(unique(B12$chr))
for (i in 1:16){
  
  B12[B12==old[i]] <- replacement[i]
}


#S12
B12_S <- na.omit(B12[,1:3])
chr_index <- unique(B12_S$chr)

B12_S_final <- data.frame()
for (i in 1:16){
  
  temp <- subset(B12_S, chr == chr_index[i])
  
  temp2 <- splineAnalyze(-log10(temp$S12_p),temp$pos,smoothness=2000,
                         plotRaw=FALSE,plotWindows=FALSE,method=3 )
  
  p <- temp2[["windowData"]][["MeanY"]]
  pos <- temp2[["windowData"]][["WindowStop"]]
  chr <- rep( chr_index[i], length(p))
  temp3 <- data.frame(chr,pos,p )
  B12_S_final <- rbind(B12_S_final, temp3)
}

#Make a genome-wide scale so the entire genome can go on a single plot
#Add the max pos of C01 to C02, C01+C02 to C03, etc
#Asssuming you are starting with a data frame called xx with one colname = pos, and another colname = var1 (some variable you want to plot)
Gaxis <- numeric(length=0)
chrs<-c('C01','C02','C03','C04','C05','C06','C07','C08','C09','C10','C11','C12','C13','C14','C15','C16')

B12_S_final <- arrange(transform(B12_S_final,
                        chr=factor(chr,levels=chrs)),chr)

for (k in 1:length(chrs)){
  data.samp<-subset(B12_S_final,chr==chrs[k])
  if (chrs[k]=='C01'){Gaxis.samp<-data.samp$pos}
  if (chrs[k]=='C02'){Gaxis.samp<-data.samp$pos+230218}
  if (chrs[k]=='C03'){Gaxis.samp<-data.samp$pos+230218+813184}
  if (chrs[k]=='C04'){Gaxis.samp<-data.samp$pos+230218+813184+316620}
  if (chrs[k]=='C05'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933}
  if (chrs[k]=='C06'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874}
  if (chrs[k]=='C07'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161}
  if (chrs[k]=='C08'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940}
  if (chrs[k]=='C09'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643}
  if (chrs[k]=='C10'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888}
  if (chrs[k]=='C11'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751}
  if (chrs[k]=='C12'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816} 
  if (chrs[k]=='C13'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816+1078177}
  if (chrs[k]=='C14'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816+1078177+924431}
  if (chrs[k]=='C15'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816+1078177+924431+784333}
  if (chrs[k]=='C16'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816+1078177+924431+784333+1091291}
  Gaxis<-c(Gaxis,Gaxis.samp) 
}
MB=Gaxis/1e6

#S12
plot(MB,B12_S_final$p,
     xlab="",
     ylab="",
     main="",
     ylim=c(0,30),
     cex.main=3, #size of text for main title
     cex.lab=2, #size of text for axis labels
     type="n",	#type = none: don't make plot yet, need to make rectangles first
     axes=FALSE  #we will do the axes by hand next
)

#Here draw are gray rectangles that delineate chromosomes
rect(0.23,-10,1.04,80,col="grey80",lty=0)  
rect(1.36,-10,2.89,80,col="grey80",lty=0)  
rect(3.46,-10,3.73,80,col="grey80",lty=0)
rect(4.82,-10,5.39,80,col="grey80",lty=0)  
rect(5.83,-10,6.57,80,col="grey80",lty=0) 
rect(7.24,-10,8.32,80,col="grey80",lty=0)
rect(9.24,-10,10.03,80,col="grey80",lty=0)
rect(11.12,-10,12.07,80,col="grey80",lty=0) 

#Here draw points or lines.  Often it is useful to play around with point type, font size, color...
points(MB,B12_S_final$p, col="black", pch = 20, cex = 0.60, type="l")
box() 

####Candidate regions from Cubillos et al. 
rect(candidate$start,25,candidate$end,27,col="red",lty=0)

####Candidate regions from Burke 2014 
rect(candidate_2$start,25,candidate_2$end,27,col="#009E73",lty=0)

#99.9 quantile
#abline(h=quantile(B12_S_final$p,probs=.999, na.rm=TRUE),v=NA, col="red")

#Now draw axes back in (you have more flexibility this way)
axis(2,cex.axis=1.5, las = 1)
legend("right",inset=-.070, "S12",bty="n", cex=2)

legend("left",inset=-.080, "(E)",bty="n", cex=2)

#K12
B12_K <- na.omit(B12[,c(1,2,4)])
chr_index <- unique(B12_K$chr)

B12_K_final <- data.frame()
for (i in 1:16){
  
  temp <- subset(B12_K, chr == chr_index[i])
  
  temp2 <- splineAnalyze(-log10(temp$K12_p),temp$pos,smoothness=2000,
                         plotRaw=FALSE,plotWindows=FALSE,method=3 )
  
  p <- temp2[["windowData"]][["MeanY"]]
  pos <- temp2[["windowData"]][["WindowStop"]]
  chr <- rep( chr_index[i], length(p))
  temp3 <- data.frame(chr,pos,p )
  B12_K_final <- rbind(B12_K_final, temp3)
}

#Make a genome-wide scale so the entire genome can go on a single plot
#Add the max pos of C01 to C02, C01+C02 to C03, etc
#Asssuming you are starting with a data frame called xx with one colname = pos, and another colname = var1 (some variable you want to plot)
Gaxis <- numeric(length=0)
chrs<-c('C01','C02','C03','C04','C05','C06','C07','C08','C09','C10','C11','C12','C13','C14','C15','C16')

B12_K_final <- arrange(transform(B12_K_final,
                                 chr=factor(chr,levels=chrs)),chr)

for (k in 1:length(chrs)){
  data.samp<-subset(B12_K_final,chr==chrs[k])
  if (chrs[k]=='C01'){Gaxis.samp<-data.samp$pos}
  if (chrs[k]=='C02'){Gaxis.samp<-data.samp$pos+230218}
  if (chrs[k]=='C03'){Gaxis.samp<-data.samp$pos+230218+813184}
  if (chrs[k]=='C04'){Gaxis.samp<-data.samp$pos+230218+813184+316620}
  if (chrs[k]=='C05'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933}
  if (chrs[k]=='C06'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874}
  if (chrs[k]=='C07'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161}
  if (chrs[k]=='C08'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940}
  if (chrs[k]=='C09'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643}
  if (chrs[k]=='C10'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888}
  if (chrs[k]=='C11'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751}
  if (chrs[k]=='C12'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816} 
  if (chrs[k]=='C13'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816+1078177}
  if (chrs[k]=='C14'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816+1078177+924431}
  if (chrs[k]=='C15'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816+1078177+924431+784333}
  if (chrs[k]=='C16'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816+1078177+924431+784333+1091291}
  Gaxis<-c(Gaxis,Gaxis.samp) 
}
MB=Gaxis/1e6

plot(MB,B12_K_final$p,
     xlab="",
     ylab="",
     main="",
     ylim=c(0,50),
     cex.main=3, #size of text for main title
     cex.lab=2, #size of text for axis labels
     type="n",	#type = none: don't make plot yet, need to make rectangles first
     axes=FALSE  #we will do the axes by hand next
)

#Here draw are gray rectangles that delineate chromosomes
rect(0.23,-10,1.04,80,col="grey80",lty=0)  
rect(1.36,-10,2.89,80,col="grey80",lty=0)  
rect(3.46,-10,3.73,80,col="grey80",lty=0)
rect(4.82,-10,5.39,80,col="grey80",lty=0)  
rect(5.83,-10,6.57,80,col="grey80",lty=0) 
rect(7.24,-10,8.32,80,col="grey80",lty=0)
rect(9.24,-10,10.03,80,col="grey80",lty=0)
rect(11.12,-10,12.07,80,col="grey80",lty=0) 

#Here draw points or lines.  Often it is useful to play around with point type, font size, color...
points(MB,B12_K_final$p, col="black",pch = 20, cex = 0.70, type="l")

####Candidate regions from Cubillos et al. 
rect(candidate$start,45,candidate$end,48,col="red",lty=0)

####Candidate regions from Burke 2014 
rect(candidate_2$start,45,candidate_2$end,48,col="#009E73",lty=0)

#99.9 quantile
#abline(h=quantile(B12_K_final$p,probs=.999, na.rm=TRUE),v=NA, col="red")


box() 
#Now draw axes back in (you have more flexibility this way)
axis(2,cex.axis=1.5, las = 1)
legend("right",inset=-.070, "K12",bty="n", cex=2)

legend("left",inset=-.080, "(F)",bty="n", cex=2)

#x axis
axis(1, at = c(0,12.07), labels=c(0,12.07),tick=FALSE, line= 2, cex.axis=1.5)

#add y and x axis label
mtext(text="-log(p-value)",side=2,line=0,outer=TRUE,cex=1.5)
mtext(text="Position (mb)",side=1,line=1.5,outer=TRUE,cex=1.5)

#add chromosome labels at midpts
mtext("C2",line = .5,side=1, at =0.635, cex=1.2)
mtext("C4",line = .5,side=1, at =2.125, cex=1.2)
mtext("C6",line = .5,side=1, at =3.595, cex=1.2)
mtext("C8",line = .5,side=1, at =5.105, cex=1.2)
mtext("C10",line = .5,side=1, at =6.2, cex=1.2)
mtext("C12",line = .5,side=1, at =7.78, cex=1.2)
mtext("C14",line = .5,side=1, at =9.635, cex=1.2)
mtext("C16",line = .5,side=1, at =11.595, cex=1.2)


