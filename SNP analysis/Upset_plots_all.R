setwd("~/Dropbox/popbuild/Mark_Scripts/AFS_Het_Ne/")
library(UpSetR)
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
alt_nuc <- cbind(alt_nuc,B4[,1:2])

#SNP counts
k <- alt_nuc[,c(1:3,7:8)]
k <- subset(k, k[,1] > 0.02 & k[,1] < .98)
SNPS_4k_0 <- k[,4:5]
SNPS_4k_6 <- subset(k, k[,2] > 0.02  & k[,2] < .98)[,4:5]
SNPS_4k_12 <- subset(k, k[,3] > 0.02  & k[,3] < .98)[,4:5]


s <- alt_nuc[,4:8]
s <- subset(s, s[,1] > 0.02  & s[,1] < .98)
SNPS_4s_0 <- s[,4:5]
SNPS_4s_6 <- subset(s, s[,2] > 0.02  & s[,2] < .98)[,4:5]
SNPS_4s_12 <- subset(s, s[,3] > 0.02  & s[,3] < .98)[,4:5]


listInput2 <- list(Possible = paste(all_snps$CHROM, all_snps$POS), S4_0 =  paste(SNPS_4s_0$CHROM, SNPS_4s_0$POS), S4_6 = paste(SNPS_4s_6$CHROM, SNPS_4s_6$POS), S4_12 = paste(SNPS_4s_12$CHROM, SNPS_4s_12$POS), K4_0= paste(SNPS_4k_0$CHROM, SNPS_4k_0$POS), K4_6 = paste(SNPS_4k_6$CHROM, SNPS_4k_6$POS), K4_12 = paste(SNPS_4k_12$CHROM, SNPS_4k_12$POS))

upset(fromList(listInput2), sets.x.label = "SNP Sets", mainbar.y.label = "", keep.order = TRUE, mainbar.y.max= 55000, intersection = list(list("Possible"),list("Possible","S4_0","S4_6", "S4_12" ,"K4_0", "K4_6", "K4_12"),list("Possible","S4_0","S4_6", "S4_12" ), list("Possible","K4_0", "K4_6", "K4_12")), order.by = "degree",  text.scale = c(1,3,2,2,2,3), set_size.scale_max  = 156000)



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


#AFS plots
alt_nuc <- B8[,maf]
new_name <- paste(type, "_", tp, sep="")
names(alt_nuc) <- new_name
alt_nuc <- cbind(alt_nuc,B8[,1:2])

#SNP counts
k <- alt_nuc[,c(1:3,7:8)]
k <- subset(k, k[,1] > 0.02 & k[,1] < .98)
SNPS_8k_0 <- k[,4:5]
SNPS_8k_6 <- subset(k, k[,2] > 0.02  & k[,2] < .98)[,4:5]
SNPS_8k_12 <- subset(k, k[,3] > 0.02  & k[,3] < .98)[,4:5]


s <- alt_nuc[,4:8]
s <- subset(s, s[,1] > 0.02  & s[,1] < .98)
SNPS_8s_0 <- s[,4:5]
SNPS_8s_6 <- subset(s, s[,2] > 0.02  & s[,2] < .98)[,4:5]
SNPS_8s_12 <- subset(s, s[,3] > 0.02  & s[,3] < .98)[,4:5]

listInput2 <- list(Possible = paste(all_snps$CHROM, all_snps$POS), S8_0 =  paste(SNPS_8s_0$CHROM, SNPS_8s_0$POS), S8_6 = paste(SNPS_8s_6$CHROM, SNPS_8s_6$POS), S8_12 = paste(SNPS_8s_12$CHROM, SNPS_8s_12$POS), K8_0= paste(SNPS_8k_0$CHROM, SNPS_8k_0$POS), K8_6 = paste(SNPS_8k_6$CHROM, SNPS_8k_6$POS), K8_12 = paste(SNPS_8k_12$CHROM, SNPS_8k_12$POS))

 upset(fromList(listInput2), sets.x.label = "SNP Sets", mainbar.y.label = "", keep.order = TRUE, mainbar.y.max= 55000, intersection = list(list("Possible"),list("Possible","S8_0","S8_6", "S8_12" ,"K8_0", "K8_6", "K8_12"),list("Possible","S8_0","S8_6", "S8_12" ), list("Possible","K8_0", "K8_6", "K8_12")), order.by = "degree", text.scale = c(1,3,2,2,2,3), set_size.scale_max  = 156000)


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


#AFS plots
alt_nuc <- B12[,maf]
new_name <- paste(type, "_", tp, sep="")
names(alt_nuc) <- new_name
alt_nuc <- cbind(alt_nuc,B12[,1:2])

#SNP counts
k <- alt_nuc[,c(1:3,7:8)]
k <- subset(k, k[,1] > 0.02 & k[,1] < .98)
SNPS_12k_0 <- k[,4:5]
SNPS_12k_6 <- subset(k, k[,2] > 0.02  & k[,2] < .98)[,4:5]
SNPS_12k_12 <- subset(k, k[,3] > 0.02  & k[,3] < .98)[,4:5]


s <- alt_nuc[,4:8]
s <- subset(s, s[,1] > 0.02  & s[,1] < .98)
SNPS_12s_0 <- s[,4:5]
SNPS_12s_6 <- subset(s, s[,2] > 0.02  & s[,2] < .98)[,4:5]
SNPS_12s_12 <- subset(s, s[,3] > 0.02  & s[,3] < .98)[,4:5]


listInput2 <- list(Possible = paste(all_snps$CHROM, all_snps$POS), S12_0 =  paste(SNPS_12s_0$CHROM, SNPS_12s_0$POS), S12_6 = paste(SNPS_12s_6$CHROM, SNPS_12s_6$POS), S12_12 = paste(SNPS_12s_12$CHROM, SNPS_12s_12$POS), K12_0= paste(SNPS_12k_0$CHROM, SNPS_12k_0$POS), K12_6 = paste(SNPS_12k_6$CHROM, SNPS_12k_6$POS), K12_12 = paste(SNPS_12k_12$CHROM, SNPS_12k_12$POS))

upset(fromList(listInput2), sets.x.label = "SNP Sets", mainbar.y.label = "", keep.order = TRUE, mainbar.y.max= 55000, intersection = list(list("Possible"),list("Possible","S12_0","S12_6", "S12_12" ,"K12_0", "K12_6", "K12_12"),list("Possible","S12_0","S12_6", "S12_12" ), list("Possible","K12_0", "K12_6", "K12_12")), order.by = "degree",text.scale = c(1,3,2,2,2,3), set_size.scale_max  = 156000)

