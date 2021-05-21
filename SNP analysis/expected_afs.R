setwd("~/Dropbox/popbuild/Mark_Scripts/AFS_Het_Ne/")
library(tidyverse)
library(viridis)
library(forcats)
library(ggplot2)
library(plyr)


#B4
all_snps <- read.table("PossibleSNPs_4F.txt", header = TRUE)
B4 <-read.table("~/Dropbox/popbuild/Mark_Scripts/Haplotype_analysis/4B/Inputfile_4.txt", header=TRUE)
B4 <- merge(B4, all_snps)

#curate SNPS based on founders
for (i in seq(5,11,by=2)){
  B4 <- subset(B4, B4[i] == 0 |B4[i] == 1 )
} ###only fixed sites in founders

#pull relevant columns and rename
B4 <- B4[,c(1:4,13:ncol(B4))]
type <- rep(c("4K", "4S"),each =3)
tp <- c("0","6","12")
maf <- paste("maf_", type, "_", tp, sep="")
cov <- paste("cov_", type, "_", tp, sep="")
names(B4)[c(seq(5,ncol(B4),by=2))] <- maf
names(B4)[c(seq(6,ncol(B4),by=2))] <- cov
alt_nuc <- B4[,maf]
new_name <- paste(type, "_", tp, sep="")
names(alt_nuc) <- new_name

#expected afs

possible <- read.table("PossibleSNPs_4F_Full.txt", header = TRUE)
expected_4 <- (apply(possible[,c(seq(5,ncol(possible),by=2))],1,sum))/4

#B8
all_snps <- read.table("PossibleSNPs_8F.txt", header = TRUE)
B8 <-read.table("~/Dropbox/popbuild/Mark_Scripts/Haplotype_analysis/8B/Inputfile_8.txt", header=TRUE)
B8 <- merge(B8, all_snps)

#curate SNPS based on founders
for (i in seq(5,19,by=2)){
  B8 <- subset(B8, B8[i] == 0 |B8[i] == 1 )
} ###only fixed sites in founders

#pull relevant columns and rename
B8 <- B8[,c(1:4,21:ncol(B8))]
type <- rep(c("8K", "8S"),each =3)
tp <- c("0","6","12")
maf <- paste("maf_", type, "_", tp, sep="")
cov <- paste("cov_", type, "_", tp, sep="")
names(B8)[c(seq(5,ncol(B8),by=2))] <- maf
names(B8)[c(seq(6,ncol(B8),by=2))] <- cov

#expected afs
possible <- read.table("PossibleSNPs_8F_Full.txt", header = TRUE)
expected_8 <- (apply(possible[,c(seq(5,ncol(possible),by=2))],1,sum))/8



#B12
all_snps <- read.table("PossibleSNPs_12F.txt", header = TRUE)
B12 <-read.table("~/Dropbox/popbuild/Mark_Scripts/Haplotype_analysis/12B/Inputfile_12.txt", header=TRUE)
B12 <- merge(B12, all_snps)

#curate SNPS based on founders
for (i in seq(5,27,by=2)){
  B12 <- subset(B12, B12[i] == 0 |B12[i] == 1 )
} ###only fixed sites in founders

#pull relevant columns and rename
B12 <- B12[,c(1:4,29:ncol(B12))]
type <- rep(c("12K", "12S"),each =3)
tp <- c("0","6","12")
maf <- paste("maf_", type, "_", tp, sep="")
cov <- paste("cov_", type, "_", tp, sep="")
names(B12)[c(seq(5,ncol(B12),by=2))] <- maf
names(B12)[c(seq(6,ncol(B12),by=2))] <- cov

#expected afs
possible <- read.table("PossibleSNPs_12F_Full.txt", header = TRUE)
expected_12 <- (apply(possible[,c(seq(5,ncol(possible),by=2))],1,sum))/12

#plot
names <- c("text", "value","color")
temp1 <- rep("A. Four Founders",length(expected_4))
temp2 <- rep("blue", length(expected_4))
e4 <- data.frame(cbind(temp1,expected_4,temp2))
names(e4) <- names


temp1 <- rep("B. Eight Founders",length(expected_8))
temp2 <- rep("green", length(expected_8))
e8 <- data.frame(cbind(temp1,expected_8,temp2))
names(e8) <- names

temp1 <- rep("C. Twelve Founders",length(expected_12))
temp2 <- rep("red", length(expected_12))
e12 <- data.frame(cbind(temp1,expected_12,temp2))
names(e12) <- names

plot_data <- rbind(e4,e8,e12)
plot_data$value <- as.numeric(plot_data$value)

new_order <- c("A. Four Founders", "B. Eight Founders", "C. Twelve Founders")
plot_data2 <- arrange(transform(plot_data,
                               text=factor(text,levels=new_order)),text)
  
p <- plot_data2 %>%
  ggplot( aes(x=value, color=color, fill=color)) +
  geom_histogram() +
  theme_bw()+
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 10,face="bold"),
    #panel.background = element_rect(fill='grey'),
    strip.background = element_rect(fill="grey")
  ) +
  xlab("SNP Frequency") +
  ylab("Counts") +
  facet_wrap(~text, ncol = 1)
p

              
              
              
              
