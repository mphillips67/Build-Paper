setwd("~/Dropbox/popbuild/Mark_Scripts/Haplotype_analysis")
library(reshape2)
library(stringr)
library(tidyverse)

names <- c('chr','pos','DBVPG6765','DBVPG6044','BC187','SK1','L_1374','YJM975','YPS128','Y12','273614N','L_1528','UWOPS052173','YJM981')
#K12
K12 <- read.table("12k_haps.txt", header= TRUE)
founders <- as.vector(str_split_fixed(K12[1,5], ";",12))
K12_freqs <- cbind(K12[,2:3],colsplit(K12$adjfounderfreqs, ";", founders))
for (i in 3:ncol(K12_freqs)){
  K12_freqs <- subset(K12_freqs, K12_freqs[i] >= 0 & K12_freqs[i] <= 1 )
}

K12_freqs$chr <- as.character(K12_freqs$chr)

replacement <- c('C10','C11','C12','C13','C14','C15','C16','C01','C02','C03','C04','C05','C06','C07','C08','C09') 
old <- as.vector(unique(K12_freqs$chr))
for (i in 1:16){
  
  K12_freqs[K12_freqs==old[i]] <- replacement[i]
}

names(K12_freqs) <- names

K12_freqs <- K12_freqs[,c(1,2,3,14,5,12,6,11,10,4,7,9,8,13)]

#Make a genome-wide scale so the entire genome can go on a single plot
#Add the max pos of C01 to C02, C01+C02 to C03, etc
#Asssuming you are starting with a data frame called xx with one colname = pos, and another colname = var1 (some variable you want to plot)
Gaxis <- numeric(length=0)
chrs<-c('C01','C02','C03','C04','C05','C06','C07','C08','C09','C10','C11','C12','C13','C14','C15','C16')

K12_freqs <- arrange(transform(K12_freqs,
                               chr=factor(chr,levels=chrs)),chr)

for (k in 1:length(chrs)){
  data.samp<-subset(K12_freqs,chr==chrs[k])
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


#mfrow allows you to make multipaneled figures, mar changes the size of the margins, mpg changes the placement of the axis labels
par(mfrow=c(12,1),mar=c(1,3,1,0),mgp=c(3,1,0), oma=c(1,0,1,0))
#plots

for (i in 3:ncol(K12_freqs)){
  
  plot(MB,K12_freqs[,i],
       xlab="",
       ylab="",
       main="",
       ylim=c(0,1),
       cex.main=3, #size of text for main title
       cex.lab=2, #size of text for axis labels
       type="n",	#type = none: don't make plot yet, need to make rectangles first
       axes=FALSE  #we will do the axes by hand next
  )
  
  title(ylab = colnames(K12_freqs[i]), cex.lab = 1.4,
        line = 1)
  
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
  points(MB,K12_freqs[,i],  type = "l", col = "blue" )
  #legend(11,3.4,colnames(K12_freqs[i]),bty="n", cex=3)
  box() 
}

#add chromosome labels at midpts
mtext("C2",line = 1,side=1, at =0.635, cex=2)
mtext("C4",line = 1,side=1, at =2.125, cex=2)
mtext("C6",line = 1,side=1, at =3.595, cex=2)
mtext("C8",line = 1,side=1, at =5.105, cex=2)
mtext("C10",line = 1,side=1, at =6.2, cex=2)
mtext("C12",line = 1,side=1, at =7.78, cex=2)
mtext("C14",line = 1,side=1, at =9.635, cex=2)
mtext("C16",line = 1,side=1, at =11.595, cex=2)

#add y axis label
#mtext(text="Haplotype Frequency",side=2,line=0,outer=TRUE,cex=1.5)




#S12
S12 <- read.table("12S_haps.txt", header= TRUE)
founders <- as.vector(str_split_fixed(S12[1,5], ";",12))
S12_freqs <- cbind(S12[,2:3],colsplit(S12$adjfounderfreqs, ";", founders))
for (i in 3:ncol(S12_freqs)){
  S12_freqs <- subset(S12_freqs, S12_freqs[i] >= 0 & S12_freqs[i] <= 1 )
}

S12_freqs$chr <- as.character(S12_freqs$chr)

replacement <- c('C10','C11','C12','C13','C14','C15','C16','C01','C02','C03','C04','C05','C06','C07','C08','C09') 
old <- as.vector(unique(S12_freqs$chr))
for (i in 1:16){
  
  S12_freqs[S12_freqs==old[i]] <- replacement[i]
}

names(S12_freqs) <- names

S12_freqs <- S12_freqs[,c(1,2,3,14,5,12,6,11,10,4,7,9,8,13)]

#Make a genome-wide scale so the entire genome can go on a single plot
#Add the max pos of C01 to C02, C01+C02 to C03, etc
#Asssuming you are starting with a data frame called xx with one colname = pos, and another colname = var1 (some variable you want to plot)
Gaxis <- numeric(length=0)
chrs<-c('C01','C02','C03','C04','C05','C06','C07','C08','C09','C10','C11','C12','C13','C14','C15','C16')

S12_freqs <- arrange(transform(S12_freqs,
                               chr=factor(chr,levels=chrs)),chr)

for (k in 1:length(chrs)){
  data.samp<-subset(S12_freqs,chr==chrs[k])
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


#mfrow allows you to make multipaneled figures, mar changes the size of the margins, mpg changes the placement of the axis labels
par(mfrow=c(12,1),mar=c(1,3,1,0),mgp=c(3,1,0), oma=c(1,0,1,0))
#plots

for (i in 3:ncol(S12_freqs)){
  
  plot(MB,S12_freqs[,i],
       xlab="",
       ylab="",
       main="",
       ylim=c(0,1),
       cex.main=3, #size of text for main title
       cex.lab=2, #size of text for axis labels
       type="n",	#type = none: don't make plot yet, need to make rectangles first
       axes=FALSE  #we will do the axes by hand next
  )
  
  title(ylab = colnames(K12_freqs[i]), cex.lab = 1.4,
        line = 1)
  
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
  points(MB,S12_freqs[,i],  type = "l", col = "blue" )
  #legend(11,3.4,colnames(S12_freqs[i]),bty="n", cex=3)
  box() 
}


#add chromosome labels at midpts
mtext("C2",line = 1,side=1, at =0.635, cex=2)
mtext("C4",line = 1,side=1, at =2.125, cex=2)
mtext("C6",line = 1,side=1, at =3.595, cex=2)
mtext("C8",line = 1,side=1, at =5.105, cex=2)
mtext("C10",line = 1,side=1, at =6.2, cex=2)
mtext("C12",line = 1,side=1, at =7.78, cex=2)
mtext("C14",line = 1,side=1, at =9.635, cex=2)
mtext("C16",line = 1,side=1, at =11.595, cex=2)


#add y axis label
#mtext(text="Haplotype Frequency",side=2,line=0,outer=TRUE,cex=1)







