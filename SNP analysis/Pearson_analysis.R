setwd("~/Dropbox/popbuild/Mark_Scripts/AFS_Het_Ne")
library(poolSeq)
#####4 founders 
B4 <- read.table("~/Dropbox/popbuild/Mark_Scripts/Haplotype_analysis/4B/Inputfile_4.txt", header=TRUE)

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
ref <- paste("ref_", type, "_", tp, sep="")
names(B4)[c(seq(5,ncol(B4),by=2))] <- maf
names(B4)[c(seq(6,ncol(B4),by=2))] <- cov
total_nuc <- B4[,cov]
alt_nuc <- B4[,maf] * total_nuc
ref_nuc <- (1-B4[,maf]) * total_nuc
names(ref_nuc) <- ref
chr <- B4$CHROM
pos <- B4$POS
B4_v2 <- data.frame(chr, pos, alt_nuc, ref_nuc)

#S4

S4_p <- chi.sq.test(A0=B4_v2$ref_S4_0, a0=B4_v2$maf_S4_0, At=B4_v2$ref_S4_12, at=B4_v2$maf_S4_12)

#K4
K4_p <- chi.sq.test(A0=B4_v2$ref_K4_0, a0=B4_v2$maf_K4_0, At=B4_v2$ref_K4_12, at=B4_v2$maf_K4_12)

#combine and print
B4_final <- data.frame(chr,pos, S4_p,K4_p)
write.table(B4_final,"B4_Pearson.txt", quote= FALSE, row.names= FALSE,col.names= TRUE, sep="\t")


####8 Founders
B8 <-read.table("~/Dropbox/popbuild/Mark_Scripts/Haplotype_analysis/8B/Inputfile_8.txt", header=TRUE)

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
ref <- paste("ref_", type, "_", tp, sep="")
names(B8)[c(seq(5,ncol(B8),by=2))] <- maf
names(B8)[c(seq(6,ncol(B8),by=2))] <- cov

total_nuc <- B8[,cov]
alt_nuc <- B8[,maf] * total_nuc
ref_nuc <- (1-B8[,maf]) * total_nuc
names(ref_nuc) <- ref
chr <- B8$CHROM
pos <- B8$POS
B8_v2 <- data.frame(chr, pos, alt_nuc, ref_nuc)

#S8

S8_p <- chi.sq.test(A0=B8_v2$ref_S8_0, a0=B8_v2$maf_S8_0, At=B8_v2$ref_S8_12, at=B8_v2$maf_S8_12)

#K8
K8_p <- chi.sq.test(A0=B8_v2$ref_K8_0, a0=B8_v2$maf_K8_0, At=B8_v2$ref_K8_12, at=B8_v2$maf_K8_12)

#combine and print
B8_final <- data.frame(chr,pos, S8_p,K8_p)
write.table(B8_final,"B8_Pearson.txt", quote= FALSE, row.names= FALSE,col.names= TRUE, sep="\t")


####12 Founders
B12 <-read.table("~/Dropbox/popbuild/Mark_Scripts/Haplotype_analysis/12B/Inputfile_12.txt", header=TRUE)

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
ref <- paste("ref_", type, "_", tp, sep="")
names(B12)[c(seq(5,ncol(B12),by=2))] <- maf
names(B12)[c(seq(6,ncol(B12),by=2))] <- cov
total_nuc <- B12[,cov]
alt_nuc <- B12[,maf] * total_nuc
ref_nuc <- (1-B12[,maf]) * total_nuc
names(ref_nuc) <- ref
chr <- B12$CHROM
pos <- B12$POS
B12_v2 <- data.frame(chr, pos, alt_nuc, ref_nuc)

#S12
S12_p <- chi.sq.test(A0=B12_v2$ref_S12_0, a0=B12_v2$maf_S12_0, At=B12_v2$ref_S12_12, at=B12_v2$maf_S12_12)

#K12
K12_p <- chi.sq.test(A0=B12_v2$ref_K12_0, a0=B12_v2$maf_K12_0, At=B12_v2$ref_K12_12, at=B12_v2$maf_K12_12)

#combine and print
B12_final <- data.frame(chr,pos, S12_p,K12_p)
write.table(B12_final,"B12_Pearson.txt", quote= FALSE, row.names= FALSE,col.names= TRUE, sep="\t")


