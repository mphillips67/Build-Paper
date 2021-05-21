setwd("~/Dropbox/popbuild/Mark_Scripts/AFS_Het_Ne/")
library(tidyverse)
library(viridis)
library(forcats)
library(ggplot2)
library(plyr)

#B4

all_snps4 <- read.table("PossibleSNPs_4F.txt", header = TRUE)

B4 <-read.table("~/Dropbox/popbuild/Mark_Scripts/Haplotype_analysis/4B/Inputfile_4.txt", header=TRUE)
B4 <- merge(B4, all_snps4)

#curate SNPS based on founders
for (i in seq(5,11,by=2)){
  B4 <- subset(B4, B4[i] == 0 |B4[i] == 1 )
} ###only fixed sites in founders

#pull relevant columns and rename
B4 <- B4[,c(1:4,13:ncol(B4))]
type <- rep(c("K4", "S4"),each =3)
tp <- c("0","6","12")
maf <- paste("maf_", type, "_", tp, sep="")
cov <- paste("cov_", type, "_", tp, sep="")
names(B4)[c(seq(5,ncol(B4),by=2))] <- maf
names(B4)[c(seq(6,ncol(B4),by=2))] <- cov
alt_nuc <- B4[,maf]
new_name <- paste(type, "_", tp, sep="")
names(alt_nuc) <- new_name

data_4 <- alt_nuc %>%
  gather(key="text", value="value") %>%
  mutate(text = gsub("\\.", " ",text)) %>%
  mutate(value = round(as.numeric(value),5))

color <- rep("blue", nrow(data_4))
data_4 <- cbind(data_4,color)

#B8
all_snps8 <- read.table("PossibleSNPs_12F.txt", header = TRUE)

B8 <-read.table("~/Dropbox/popbuild/Mark_Scripts/Haplotype_analysis/8B/Inputfile_8.txt", header=TRUE)
B8 <- merge(B8, all_snps8)

#curate SNPS based on founders
for (i in seq(5,19,by=2)){
  B8 <- subset(B8, B8[i] == 0 |B8[i] == 1 )
} ###only fixed sites in founders

#pull relevant columns and rename
B8 <- B8[,c(1:4,21:ncol(B8))]
type <- rep(c("K8", "S8"),each =3)
tp <- c("0","6","12")
maf <- paste("maf_", type, "_", tp, sep="")
cov <- paste("cov_", type, "_", tp, sep="")
names(B8)[c(seq(5,ncol(B8),by=2))] <- maf
names(B8)[c(seq(6,ncol(B8),by=2))] <- cov


alt_nuc <- B8[,maf]
new_name <- paste(type, "_", tp, sep="")
names(alt_nuc) <- new_name

data_8 <- alt_nuc %>%
  gather(key="text", value="value") %>%
  mutate(text = gsub("\\.", " ",text)) %>%
  mutate(value = round(as.numeric(value),5))

color <- rep("green", nrow(data_8))
data_8 <- cbind(data_8,color)

#B12

all_snps12 <- read.table("PossibleSNPs_12F.txt", header = TRUE)

B12 <-read.table("~/Dropbox/popbuild/Mark_Scripts/Haplotype_analysis/12B/Inputfile_12.txt", header=TRUE)
B12 <- merge(B12, all_snps12)

#curate SNPS based on founders
for (i in seq(5,27,by=2)){
  B12 <- subset(B12, B12[i] == 0 |B12[i] == 1 )
} ###only fixed sites in founders

#pull relevant columns and rename
B12 <- B12[,c(1:4,29:ncol(B12))]
type <- rep(c("K12", "S12"),each =3)
tp <- c("0","6","12")
maf <- paste("maf_", type, "_", tp, sep="")
cov <- paste("cov_", type, "_", tp, sep="")
names(B12)[c(seq(5,ncol(B12),by=2))] <- maf
names(B12)[c(seq(6,ncol(B12),by=2))] <- cov


alt_nuc <- B12[,maf]
new_name <- paste(type, "_", tp, sep="")
names(alt_nuc) <- new_name

data_12 <- alt_nuc %>%
  gather(key="text", value="value") %>%
  mutate(text = gsub("\\.", " ",text)) %>%
  mutate(value = round(as.numeric(value),5))

color <- rep("red", nrow(data_12))
data_12 <- cbind(data_12,color)

#plot
data_all <- rbind(data_4, data_8,data_12)
type <- rep(c("K4", "S4","K8","S8","K12","S12"),each =3)
tp <- c("0","6","12")
new_order <- paste(type, "_", tp, sep="")
data_all2 <- arrange(transform(data_all,
                           text=factor(text,levels=new_order)),text)



p <- data_all2 %>%
  ggplot( aes(x=value, color=color, fill=color)) +
  geom_histogram(bins = 20) +
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
  facet_wrap(~text, ncol = 3)
p


