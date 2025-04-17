

library(dada2); packageVersion("dada2")
library(ShortRead); packageVersion("ShortRead")
library(dplyr); packageVersion("dplyr")
library(tidyr); packageVersion("tidyr")
library(Hmisc); packageVersion("Hmisc")
library(ggplot2); packageVersion("ggplot2")
library(plotly); packageVersion("plotly")

cutadapt <- "/usr/local/bin/cutadapt"
system2(cutadapt, args = "--version")
table.fp <- file.path("/home/sskoldas/TEZ/Disease_Stage/Dataset/v4", "tabletax")
dir.create(table.fp)

FWD <- "GTGCCAGCMGCCGCGGTAA"
REV <- "GGACTACHVGGGTWTCTAAT"


#PRJEB28175-------------------------------------------------------------------------------
path1 <- "/home/sskoldas/TEZ/Disease_Stage/Dataset/v4/PRJEB28175/raw/"
head(list.files(path1))
fn1 <- list.files(path1, pattern = "fastq.gz", full.names = TRUE)
filter.fp1 <- file.path("/home/sskoldas/TEZ/Disease_Stage/Dataset/v4/PRJEB28175", "filter")
trimmed.fp1 <- file.path("/home/sskoldas/TEZ/Disease_Stage/Dataset/v4/PRJEB28175", "trimmed")
fn.filtN1 <- file.path(path1, "filtN", basename(fn1))

# filter the data
filtfn1 <-file.path(path1, "filtered", basename(fn1)) 
filt_out1 <- filterAndTrim(fwd = fn1, filt = filtfn1, maxEE=1.5, truncQ = 8, truncLen = 150, rm.phix = TRUE, maxN = 0,
                          compress = TRUE, verbose = TRUE, multithread = TRUE)
# Learn the Error Rates
filtfn1 <- list.files(path = "/home/sskoldas/TEZ/Disease_Stage/Dataset/v4/PRJEB28175/raw/filtered",
                      pattern = "fastq.gz", full.names = TRUE)

sample.names1 <-sapply(strsplit(basename(filtfn1), ".fastq.gz", fixed = TRUE), `[`,1)
head(sample.names1)
names(filtfn1) <- sample.names1

set.seed(123)
err1 <- learnErrors(filtfn1, multithread = TRUE, verbose = 1)
# Dereplication, sequence inference
dadafn1 <- vector("list", length(sample.names1))
names(dadafn1) <- sample.names1
for (i in sample.names1){
  derep1 <- derepFastq(filtfn1[[i]])
  dadafn1[[i]] <- dada(derep1, err = err1, multithread = TRUE)
}
seqtab1 <- makeSequenceTable(dadafn1)

# PRJNA339813-------------------------------------------------------------------------------------------------
data.fp2 <- "/home/sskoldas/TEZ/Disease_Stage/Dataset/v4/PRJNA339813/raw"
list.files(data.fp2)
filter.fp2 <- file.path("/home/sskoldas/TEZ/Disease_Stage/Dataset/v4/PRJNA339813", "filter")
trimmed.fp2 <- file.path("/home/sskoldas/TEZ/Disease_Stage/Dataset/v4/PRJNA339813", "trimmed")
fnFs2 <- sort(list.files(data.fp2, pattern = "_R1.fastq.gz", full.names = TRUE))
fnRs2 <- sort(list.files(data.fp2, pattern = "_R2.fastq.gz", full.names = TRUE))
fnFs.filtN2 <- file.path(data.fp2, "filtN", basename(fnFs2))
fnRs.filtN2 <- file.path(data.fp2, "filtN", basename(fnRs2))
filterAndTrim(fnFs2, fnFs.filtN2, fnRs2, fnRs.filtN2, maxN = 0, multithread = TRUE)
if(!dir.exists(trimmed.fp2)) dir.create(trimmed.fp2)
fnFs.cut2 <- file.path(trimmed.fp2, basename(fnFs2))
fnRs.cut2 <- file.path(trimmed.fp2, basename(fnRs2))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# create the cutadapt flags
R1.flags <- paste("-g", FWD, "-a", REV.RC, "--minimum-length 100")
R2.flags <- paste("-G", REV, "-A", FWD.RC, "--minimum-length 100")

# run cutadapt
for(i in seq_along(fnFs2)){
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2,
                             "-o", fnFs.cut2[i], "-p", fnRs.cut2[i],
                             fnFs.filtN2[i], fnRs.filtN2[i]))
}

dir.create(filter.fp2)
subF.fp2 <- file.path(filter.fp2, "preprocessed_F")
subR.fp2 <- file.path(filter.fp2, "preprocessed_R")
dir.create(subF.fp2)
dir.create(subR.fp2)
fnFs.Q2 <- file.path(subF.fp2, basename(fnFs2))
fnRs.Q2 <- file.path(subR.fp2, basename(fnRs2))
file.rename(from = fnFs.cut2, to = fnFs.Q2)
file.rename(from = fnRs.cut2, to = fnRs.Q2)
filtpathF2 <- file.path(subF.fp2, "filtered")
filtpathR2 <- file.path(subR.fp2, "filtered")
fastqFs2 <- sort(list.files(subF.fp2, pattern = "fastq.gz"))
fastqRs2 <- sort(list.files(subR.fp2, pattern = "fastq.gz"))
if(length(fastqFs2) != length(fastqRs2)) stop("Forward and reverse files do not match!")

# filter the data
filt_out2 <- filterAndTrim(fwd = file.path(subF.fp2, fastqFs2), filt = file.path(filtpathF2, fastqFs2),
                          rev = file.path(subR.fp2, fastqRs2), filt.rev = file.path(filtpathR2, fastqRs2),
                          truncLen = c(245,230), maxEE = 1.5, truncQ = 10, maxN = 0, rm.phix = TRUE,
                          compress = TRUE, verbose = TRUE, multithread = TRUE)
# infer sequence variants
filtFs2 <- list.files(filtpathF2, pattern = "fastq.gz", full.names = TRUE)
filtRs2 <- list.files(filtpathR2, pattern = "fastq.gz", full.names = TRUE)
sample.names2 <- sapply(strsplit(basename(filtFs2), "_"), `[`,1)
sample.namesR2 <- sapply(strsplit(basename(filtRs2), "_"), `[`,1)

if(!identical(sample.names2, sample.namesR2)) stop("forward and reverse files do not match!")
names(filtFs2) <- sample.names2
names(filtRs2) <- sample.names2

set.seed(123)
errF2 <- learnErrors(filtFs2, multithread = TRUE)
errR2 <- learnErrors(filtRs2, multithread = TRUE)

mergers2 <- vector("list", length(sample.names2))
names(mergers2) <- sample.names2
ddF2 <- vector("list", length(sample.names2))
names(ddF2) <- sample.names2
ddR2 <- vector("list", length(sample.names2))
names(ddR2) <- sample.names2

for(i in sample.names2){
  derepF2 <- derepFastq(filtFs2[[i]])
  dadaF2 <- dada(derepF2, err = errF2, multithread = TRUE)
  ddF2[[i]] <- dadaF2
  
  derepR2 <- derepFastq(filtRs2[[i]])
  dadaR2 <- dada(derepR2, err = errR2, multithread = TRUE)
  ddR2[[i]] <- dadaR2
  
  merger2 <- mergePairs(ddF2[[i]], derepF2, ddR2[[i]], derepR2)
  mergers2[[i]] <- merger2
}

rm(derepF2); rm(derepR2)
seqtab2 <- makeSequenceTable(mergers2)

# PRJNA374847------------------------------------------------------------------------------------------
data.fp3 <- "/home/sskoldas/TEZ/Disease_Stage/Dataset/v4/PRJNA374847/raw/"
list.files(data.fp3)
filter.fp3 <- file.path("/home/sskoldas/TEZ/Disease_Stage/Dataset/v4/PRJNA374847", "filter")
trimmed.fp3 <- file.path("/home/sskoldas/TEZ/Disease_Stage/Dataset/v4/PRJNA374847", "trimmed")
fnFs3 <- sort(list.files(data.fp3, pattern = "_1.fastq.gz", full.names = TRUE))
fnRs3 <- sort(list.files(data.fp3, pattern = "_2.fastq.gz", full.names = TRUE))
fnFs.filtN3 <- file.path(data.fp3, "filtN", basename(fnFs3))
fnRs.filtN3 <- file.path(data.fp3, "filtN", basename(fnRs3))
filterAndTrim(fnFs3, fnFs.filtN3, fnRs3, fnRs.filtN3, maxN = 0, multithread = TRUE)

if(!dir.exists(trimmed.fp3)) dir.create(trimmed.fp3)
fnFs.cut3 <- file.path(trimmed.fp3, basename(fnFs3))
fnRs.cut3 <- file.path(trimmed.fp3, basename(fnRs3))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# create the cutadapt flags
R1.flags <- paste("-g", FWD, "-a", REV.RC, "--minimum-length 100")
R2.flags <- paste("-G", REV, "-A", FWD.RC, "--minimum-length 100")
# run cutadapt
for(i in seq_along(fnFs3)){
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2,
                             "-o", fnFs.cut3[i], "-p", fnRs.cut3[i],
                             fnFs.filtN3[i], fnRs.filtN3[i]))
}

dir.create(filter.fp3)
subF.fp3 <- file.path(filter.fp3, "preprocessed_F")
subR.fp3 <- file.path(filter.fp3, "preprocessed_R")
dir.create(subF.fp3)
dir.create(subR.fp3)

fnFs.Q3 <- file.path(subF.fp3, basename(fnFs3))
fnRs.Q3 <- file.path(subR.fp3, basename(fnRs3))
file.rename(from = fnFs.cut3, to = fnFs.Q3)
file.rename(from = fnRs.cut3, to = fnRs.Q3)

filtpathF3 <- file.path(subF.fp3, "filtered")
filtpathR3 <- file.path(subR.fp3, "filtered")
fastqFs3 <- sort(list.files(subF.fp3, pattern = "fastq.gz"))
fastqRs3 <- sort(list.files(subR.fp3, pattern = "fastq.gz"))

if(length(fastqFs3) != length(fastqRs3)) stop("Forward and reverse files do not match!")

# filter the data
filt_out3 <- filterAndTrim(fwd = file.path(subF.fp3, fastqFs3), filt = file.path(filtpathF3, fastqFs3),
                          rev = file.path(subR.fp3, fastqRs3), filt.rev = file.path(filtpathR3, fastqRs3),
                          truncLen = c(200,150), maxEE = c(3,3.5), truncQ = 2, maxN = 0, rm.phix = TRUE,
                          compress = TRUE, verbose = TRUE, multithread = TRUE)

# infer sequence variants
filtFs3 <- list.files(filtpathF3, pattern = "fastq.gz", full.names = TRUE)
filtRs3 <- list.files(filtpathR3, pattern = "fastq.gz", full.names = TRUE)

sample.names3 <- sapply(strsplit(basename(filtFs3), "_"), `[`,1)
sample.namesR3 <- sapply(strsplit(basename(filtRs3), "_"), `[`,1)

if(!identical(sample.names3, sample.namesR3)) stop("forward and reverse files do not match!")
names(filtFs3) <- sample.names3
names(filtRs3) <- sample.names3

set.seed(123)
errF3 <- learnErrors(filtFs3, multithread = TRUE)
errR3 <- learnErrors(filtRs3, multithread = TRUE)

mergers3 <- vector("list", length(sample.names3))
names(mergers3) <- sample.names3
ddF3 <- vector("list", length(sample.names3))
names(ddF3) <- sample.names3
ddR3 <- vector("list", length(sample.names3))
names(ddR3) <- sample.names3

for(i in sample.names3){
  derepF3 <- derepFastq(filtFs3[[i]])
  dadaF3 <- dada(derepF3, err = errF3, multithread = TRUE)
  ddF3[[i]] <- dadaF3
  
  derepR3 <- derepFastq(filtRs3[[i]])
  dadaR3 <- dada(derepR3, err = errR3, multithread = TRUE)
  ddR3[[i]] <- dadaR3
  
  merger3 <- mergePairs(ddF3[[i]], derepF3, ddR3[[i]], derepR3)
  mergers3[[i]] <- merger3
}

rm(derepF3); rm(derepR3)
seqtab3 <- makeSequenceTable(mergers3)

# PRJNA491861 ---------------------------------------------------------------------------------------------
data.fp4 <- "/home/sskoldas/TEZ/Disease_Stage/Dataset/v4/PRJNA491861/raw/"
list.files(data.fp4)
filter.fp4 <- file.path("/home/sskoldas/TEZ/Disease_Stage/Dataset/v4/PRJNA491861", "filter")
trimmed.fp4 <- file.path("/home/sskoldas/TEZ/Disease_Stage/Dataset/v4/PRJNA491861", "trimmed")
fnFs4 <- sort(list.files(data.fp4, pattern = "_R1.fastq.gz", full.names = TRUE))
fnRs4 <- sort(list.files(data.fp4, pattern = "_R2.fastq.gz", full.names = TRUE))
fnFs.filtN4 <- file.path(data.fp4, "filtN", basename(fnFs4))
fnRs.filtN4 <- file.path(data.fp4, "filtN", basename(fnRs4))
filterAndTrim(fnFs4, fnFs.filtN4, fnRs4, fnRs.filtN4, maxN = 0, multithread = TRUE)

if(!dir.exists(trimmed.fp4)) dir.create(trimmed.fp4)
fnFs.cut4 <- file.path(trimmed.fp4, basename(fnFs4))
fnRs.cut4 <- file.path(trimmed.fp4, basename(fnRs4))
FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# create the cutadapt flags
R1.flags <- paste("-g", FWD, "-a", REV.RC, "--minimum-length 100")
R2.flags <- paste("-G", REV, "-A", FWD.RC, "--minimum-length 100")

# run cutadapt
for(i in seq_along(fnFs4)){
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2,
                             "-o", fnFs.cut4[i], "-p", fnRs.cut4[i],
                             fnFs.filtN4[i], fnRs.filtN4[i]))
}

dir.create(filter.fp4)
subF.fp4 <- file.path(filter.fp4, "preprocessed_F")
subR.fp4 <- file.path(filter.fp4, "preprocessed_R")
dir.create(subF.fp4)
dir.create(subR.fp4)

fnFs.Q4 <- file.path(subF.fp4, basename(fnFs4))
fnRs.Q4 <- file.path(subR.fp4, basename(fnRs4))
file.rename(from = fnFs.cut4, to = fnFs.Q4)
file.rename(from = fnRs.cut4, to = fnRs.Q4)

filtpathF4 <- file.path(subF.fp4, "filtered")
filtpathR4 <- file.path(subR.fp4, "filtered")
fastqFs4 <- sort(list.files(subF.fp4, pattern = "fastq.gz"))
fastqRs4 <- sort(list.files(subR.fp4, pattern = "fastq.gz"))

if(length(fastqFs4) != length(fastqRs4)) stop("Forward and reverse files do not match!")

# filter the data
filt_out4 <- filterAndTrim(fwd = file.path(subF.fp4, fastqFs4), filt = file.path(filtpathF4, fastqFs4),
                          rev = file.path(subR.fp4, fastqRs4), filt.rev = file.path(filtpathR4, fastqRs4),
                          truncLen = c(200,190), maxEE = c(1.5,2), truncQ = 5, maxN = 0, rm.phix = TRUE,
                          compress = TRUE, verbose = TRUE, multithread = TRUE)
# infer sequence variants
filtFs4 <- list.files(filtpathF4, pattern = "fastq.gz", full.names = TRUE)
filtRs4 <- list.files(filtpathR4, pattern = "fastq.gz", full.names = TRUE)

sample.names4 <- sapply(strsplit(basename(filtFs4), "_"), `[`,1)
sample.namesR4 <- sapply(strsplit(basename(filtRs4), "_"), `[`,1)

if(!identical(sample.names4, sample.namesR4)) stop("forward and reverse files do not match!")
names(filtFs4) <- sample.names4
names(filtRs4) <- sample.names4

set.seed(123)
errF4 <- learnErrors(filtFs4, multithread = TRUE)
errR4 <- learnErrors(filtRs4, multithread = TRUE)

mergers4 <- vector("list", length(sample.names4))
names(mergers4) <- sample.names4
ddF4 <- vector("list", length(sample.names4))
names(ddF4) <- sample.names4
ddR4 <- vector("list", length(sample.names4))
names(ddR4) <- sample.names4

for(i in sample.names4){
  derepF4 <- derepFastq(filtFs4[[i]])
  dadaF4 <- dada(derepF4, err = errF4, multithread = TRUE)
  ddF4[[i]] <- dadaF4
  
  derepR4 <- derepFastq(filtRs4[[i]])
  dadaR4 <- dada(derepR4, err = errR4, multithread = TRUE)
  ddR4[[i]] <- dadaR4
  
  merger4 <- mergePairs(ddF4[[i]], derepF4, ddR4[[i]], derepR4)
  mergers4[[i]] <- merger4
}
rm(derepF4); rm(derepR4)
seqtab4 <- makeSequenceTable(mergers4)


# PRJNA520924-----------------------------------------------------------------------------------------------------------
data.fp5 <- "/home/sskoldas/TEZ/Disease_Stage/Dataset/v4/PRJNA520924/raw/"
list.files(data.fp5)
filter.fp5 <- file.path("/home/sskoldas/TEZ/Disease_Stage/Dataset/v4/PRJNA520924", "filter")
trimmed.fp5 <- file.path("/home/sskoldas/TEZ/Disease_Stage/Dataset/v4/PRJNA520924", "trimmed")
fnFs5 <- sort(list.files(data.fp5, pattern = "_1.fastq.gz", full.names = TRUE))
fnRs5 <- sort(list.files(data.fp5, pattern = "_2.fastq.gz", full.names = TRUE))
fnFs.filtN5 <- file.path(data.fp5, "filtN", basename(fnFs5))
fnRs.filtN5 <- file.path(data.fp5, "filtN", basename(fnRs5))
filterAndTrim(fnFs5, fnFs.filtN5, fnRs5, fnRs.filtN5, maxN = 0, multithread = TRUE)
if(!dir.exists(trimmed.fp5)) dir.create(trimmed.fp5)
fnFs.cut5 <- file.path(trimmed.fp5, basename(fnFs5))
fnRs.cut5 <- file.path(trimmed.fp5, basename(fnRs5))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

# create the cutadapt flags
R1.flags <- paste("-g", FWD, "-a", REV.RC, "--minimum-length 100")
R2.flags <- paste("-G", REV, "-A", FWD.RC, "--minimum-length 100")

# run cutadapt
for(i in seq_along(fnFs5)){
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2,
                             "-o", fnFs.cut5[i], "-p", fnRs.cut5[i],
                             fnFs.filtN5[i], fnRs.filtN5[i]))
}

dir.create(filter.fp5)
subF.fp5 <- file.path(filter.fp5, "preprocessed_F")
subR.fp5 <- file.path(filter.fp5, "preprocessed_R")
dir.create(subF.fp5)
dir.create(subR.fp5)

fnFs.Q5 <- file.path(subF.fp5, basename(fnFs5))
fnRs.Q5 <- file.path(subR.fp5, basename(fnRs5))
file.rename(from = fnFs.cut5, to = fnFs.Q5)
file.rename(from = fnRs.cut5, to = fnRs.Q5)

filtpathF5 <- file.path(subF.fp5, "filtered")
filtpathR5 <- file.path(subR.fp5, "filtered")
fastqFs5 <- sort(list.files(subF.fp5, pattern = "fastq.gz"))
fastqRs5 <- sort(list.files(subR.fp5, pattern = "fastq.gz"))

if(length(fastqFs5) != length(fastqRs5)) stop("Forward and reverse files do not match!")

# filter the data
filt_out5 <- filterAndTrim(fwd = file.path(subF.fp5, fastqFs5), filt = file.path(filtpathF5, fastqFs5),
                          rev = file.path(subR.fp5, fastqRs5), filt.rev = file.path(filtpathR5, fastqRs5),
                          truncLen = c(220,160), maxEE = c(1.6,2), truncQ = 8, maxN = 0, rm.phix = TRUE,
                          compress = TRUE, verbose = TRUE, multithread = TRUE)

# infer sequence variants
filtFs5 <- list.files(filtpathF5, pattern = "fastq.gz", full.names = TRUE)
filtRs5 <- list.files(filtpathR5, pattern = "fastq.gz", full.names = TRUE)

sample.names5 <- sapply(strsplit(basename(filtFs5), "_"), `[`,1)
sample.namesR5 <- sapply(strsplit(basename(filtRs5), "_"), `[`,1)

if(!identical(sample.names5, sample.namesR5)) stop("forward and reverse files do not match!")
names(filtFs5) <- sample.names5
names(filtRs5) <- sample.names5

set.seed(123)
errF5 <- learnErrors(filtFs5, multithread = TRUE)
errR5 <- learnErrors(filtRs5, multithread = TRUE)

mergers5 <- vector("list", length(sample.names5))
names(mergers5) <- sample.names5
ddF5 <- vector("list", length(sample.names5))
names(ddF5) <- sample.names5
ddR5 <- vector("list", length(sample.names5))
names(ddR5) <- sample.names5

for(i in sample.names5){
  derepF5 <- derepFastq(filtFs5[[i]])
  dadaF5 <- dada(derepF5, err = errF5, multithread = TRUE)
  ddF5[[i]] <- dadaF5
  
  derepR5 <- derepFastq(filtRs5[[i]])
  dadaR5 <- dada(derepR5, err = errR5, multithread = TRUE)
  ddR5[[i]] <- dadaR5
  
  merger5 <- mergePairs(ddF5[[i]], derepF5, ddR5[[i]], derepR5)
  mergers5[[i]] <- merger5
}

rm(derepF5); rm(derepR5)
seqtab5 <- makeSequenceTable(mergers5)

# PRJNA611611-------------------------------------------------------------------------------------------------
data.fp6 <- "/home/sskoldas/TEZ/Disease_Stage/Dataset/v4/PRJNA611611/raw/"
list.files(data.fp6)
filter.fp6 <- file.path("/home/sskoldas/TEZ/Disease_Stage/Dataset/v4/PRJNA611611", "filter")
trimmed.fp6 <- file.path("/home/sskoldas/TEZ/Disease_Stage/Dataset/v4/PRJNA611611", "trimmed")
fnFs6 <- sort(list.files(data.fp6, pattern = "_1.fastq.gz", full.names = TRUE))
fnRs6 <- sort(list.files(data.fp6, pattern = "_2.fastq.gz", full.names = TRUE))
fnFs.filtN6 <- file.path(data.fp6, "filtN", basename(fnFs6))
fnRs.filtN6 <- file.path(data.fp6, "filtN", basename(fnRs6))
filterAndTrim(fnFs6, fnFs.filtN6, fnRs6, fnRs.filtN6, maxN = 0, multithread = TRUE)
if(!dir.exists(trimmed.fp6)) dir.create(trimmed.fp6)
fnFs.cut6 <- file.path(trimmed.fp6, basename(fnFs6))
fnRs.cut6 <- file.path(trimmed.fp6, basename(fnRs6))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

# create the cutadapt flags
R1.flags <- paste("-g", FWD, "-a", REV.RC, "--minimum-length 100")
R2.flags <- paste("-G", REV, "-A", FWD.RC, "--minimum-length 100")

# run cutadapt
for(i in seq_along(fnFs6)){
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2,
                             "-o", fnFs.cut6[i], "-p", fnRs.cut6[i],
                             fnFs.filtN6[i], fnRs.filtN6[i]))
}

dir.create(filter.fp6)
subF.fp6 <- file.path(filter.fp6, "preprocessed_F")
subR.fp6 <- file.path(filter.fp6, "preprocessed_R")
dir.create(subF.fp6)
dir.create(subR.fp6)

fnFs.Q6 <- file.path(subF.fp6, basename(fnFs6))
fnRs.Q6 <- file.path(subR.fp6, basename(fnRs6))
file.rename(from = fnFs.cut6, to = fnFs.Q6)
file.rename(from = fnRs.cut6, to = fnRs.Q6)

filtpathF6 <- file.path(subF.fp6, "filtered")
filtpathR6 <- file.path(subR.fp6, "filtered")
fastqFs6 <- sort(list.files(subF.fp6, pattern = "fastq.gz"))
fastqRs6 <- sort(list.files(subR.fp6, pattern = "fastq.gz"))

if(length(fastqFs6) != length(fastqRs6)) stop("Forward and reverse files do not match!")

# filter the data
filt_out6 <- filterAndTrim(fwd = file.path(subF.fp6, fastqFs6), filt = file.path(filtpathF6, fastqFs6),
                          rev = file.path(subR.fp6, fastqRs6), filt.rev = file.path(filtpathR6, fastqRs6),
                          truncLen = c(220,160), maxEE = c(1.5,2), truncQ = 6, maxN = 0, rm.phix = TRUE,
                          compress = TRUE, verbose = TRUE, multithread = TRUE)

# infer sequence variants
filtFs6 <- list.files(filtpathF6, pattern = "fastq.gz", full.names = TRUE)
filtRs6 <- list.files(filtpathR6, pattern = "fastq.gz", full.names = TRUE)

sample.names6 <- sapply(strsplit(basename(filtFs6), "_"), `[`,1)
sample.namesR6 <- sapply(strsplit(basename(filtRs6), "_"), `[`,1)

if(!identical(sample.names6, sample.namesR6)) stop("forward and reverse files do not match!")
names(filtFs6) <- sample.names6
names(filtRs6) <- sample.names6

set.seed(123)
errF6 <- learnErrors(filtFs6, multithread = TRUE)
errR6 <- learnErrors(filtRs6, multithread = TRUE)

mergers6 <- vector("list", length(sample.names6))
names(mergers6) <- sample.names6
ddF6 <- vector("list", length(sample.names6))
names(ddF6) <- sample.names6
ddR6 <- vector("list", length(sample.names6))
names(ddR6) <- sample.names6

for(i in sample.names6){
  derepF6 <- derepFastq(filtFs6[[i]])
  dadaF6 <- dada(derepF6, err = errF6, multithread = TRUE)
  ddF6[[i]] <- dadaF6
  
  derepR6 <- derepFastq(filtRs6[[i]])
  dadaR6 <- dada(derepR6, err = errR6, multithread = TRUE)
  ddR6[[i]] <- dadaR6
  
  merger6 <- mergePairs(ddF6[[i]], derepF6, ddR6[[i]], derepR6)
  mergers6[[i]] <- merger6
}

rm(derepF6); rm(derepR6)
seqtab6 <- makeSequenceTable(mergers6)

# PRJNA733203-----------------------------------------------------------------------------------------------
data.fp7 <- "/home/sskoldas/TEZ/Disease_Stage/Dataset/v4/PRJNA733203/raw/"
list.files(data.fp7)
filter.fp7 <- file.path("/home/sskoldas/TEZ/Disease_Stage/Dataset/v4/PRJNA733203", "filter")
trimmed.fp7 <- file.path("/home/sskoldas/TEZ/Disease_Stage/Dataset/v4/PRJNA733203", "trimmed")
fnFs7 <- sort(list.files(data.fp7, pattern = "_1.fastq.gz", full.names = TRUE))
fnRs7 <- sort(list.files(data.fp7, pattern = "_2.fastq.gz", full.names = TRUE))
fnFs.filtN7 <- file.path(data.fp7, "filtN", basename(fnFs7))
fnRs.filtN7 <- file.path(data.fp7, "filtN", basename(fnRs7))
filterAndTrim(fnFs7, fnFs.filtN7, fnRs7, fnRs.filtN7, maxN = 0, multithread = TRUE)
if(!dir.exists(trimmed.fp7)) dir.create(trimmed.fp7)
fnFs.cut7 <- file.path(trimmed.fp7, basename(fnFs7))
fnRs.cut7 <- file.path(trimmed.fp7, basename(fnRs7))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

# create the cutadapt flags
R1.flags <- paste("-g", FWD, "-a", REV.RC, "--minimum-length 100")
R2.flags <- paste("-G", REV, "-A", FWD.RC, "--minimum-length 100")

# run cutadapt
for(i in seq_along(fnFs7)){
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2,
                             "-o", fnFs.cut7[i], "-p", fnRs.cut7[i],
                             fnFs.filtN7[i], fnRs.filtN7[i]))
}

dir.create(filter.fp7)
subF.fp7 <- file.path(filter.fp7, "preprocessed_F")
subR.fp7 <- file.path(filter.fp7, "preprocessed_R")
dir.create(subF.fp7)
dir.create(subR.fp7)

fnFs.Q7 <- file.path(subF.fp7, basename(fnFs7))
fnRs.Q7 <- file.path(subR.fp7, basename(fnRs7))
file.rename(from = fnFs.cut7, to = fnFs.Q7)
file.rename(from = fnRs.cut7, to = fnRs.Q7)

filtpathF7 <- file.path(subF.fp7, "filtered")
filtpathR7 <- file.path(subR.fp7, "filtered")
fastqFs7 <- sort(list.files(subF.fp7, pattern = "fastq.gz"))
fastqRs7 <- sort(list.files(subR.fp7, pattern = "fastq.gz"))

if(length(fastqFs7) != length(fastqRs7)) stop("Forward and reverse files do not match!")

# filter the data
filt_out7 <- filterAndTrim(fwd = file.path(subF.fp7, fastqFs7), filt = file.path(filtpathF7, fastqFs7),
                          rev = file.path(subR.fp7, fastqRs7), filt.rev = file.path(filtpathR7, fastqRs7),
                          truncLen = c(220,170), maxEE = c(1.6,2), truncQ = 5, maxN = 0, rm.phix = TRUE,
                          compress = TRUE, verbose = TRUE, multithread = TRUE)

# infer sequence variants
filtFs7 <- list.files(filtpathF7, pattern = "fastq.gz", full.names = TRUE)
filtRs7 <- list.files(filtpathR7, pattern = "fastq.gz", full.names = TRUE)

sample.names7 <- sapply(strsplit(basename(filtFs7), "_"), `[`,1)
sample.namesR7 <- sapply(strsplit(basename(filtRs7), "_"), `[`,1)

if(!identical(sample.names7, sample.namesR7)) stop("forward and reverse files do not match!")
names(filtFs7) <- sample.names7
names(filtRs7) <- sample.names7

set.seed(123)
errF7 <- learnErrors(filtFs7, multithread = TRUE)
errR7 <- learnErrors(filtRs7, multithread = TRUE)

mergers7 <- vector("list", length(sample.names7))
names(mergers7) <- sample.names7
ddF7 <- vector("list", length(sample.names7))
names(ddF7) <- sample.names7
ddR7 <- vector("list", length(sample.names7))
names(ddR7) <- sample.names7

for(i in sample.names7){
  derepF7 <- derepFastq(filtFs7[[i]])
  dadaF7 <- dada(derepF7, err = errF7, multithread = TRUE)
  ddF7[[i]] <- dadaF7
  
  derepR7 <- derepFastq(filtRs7[[i]])
  dadaR7 <- dada(derepR7, err = errR7, multithread = TRUE)
  ddR7[[i]] <- dadaR7
  
  merger7 <- mergePairs(ddF7[[i]], derepF7, ddR7[[i]], derepR7)
  mergers7[[i]] <- merger7
}

rm(derepF7); rm(derepR7)
seqtab7 <- makeSequenceTable(mergers7)



# PRJNA529405-----------------------------------------------------------------------------------------------
data.fp8 <- "/home/sskoldas/TEZ/Disease_Stage/Dataset/v4/PRJNA529405/raw/"
list.files(data.fp8)
filter.fp8 <- file.path("/home/sskoldas/TEZ/Disease_Stage/Dataset/v4/PRJNA529405", "filter")
trimmed.fp8 <- file.path("/home/sskoldas/TEZ/Disease_Stage/Dataset/v4/PRJNA529405", "trimmed")
fnFs8 <- sort(list.files(data.fp8, pattern = "_1.fastq.gz", full.names = TRUE))
fnRs8 <- sort(list.files(data.fp8, pattern = "_2.fastq.gz", full.names = TRUE))
fnFs.filtN8 <- file.path(data.fp8, "filtN", basename(fnFs8))
fnRs.filtN8 <- file.path(data.fp8, "filtN", basename(fnRs8))
filterAndTrim(fnFs8, fnFs.filtN8, fnRs8, fnRs.filtN8, maxN = 0, multithread = TRUE)
if(!dir.exists(trimmed.fp8)) dir.create(trimmed.fp8)
fnFs.cut8 <- file.path(trimmed.fp8, basename(fnFs8))
fnRs.cut8 <- file.path(trimmed.fp8, basename(fnRs8))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

# create the cutadapt flags
R1.flags <- paste("-g", FWD, "-a", REV.RC, "--minimum-length 100")
R2.flags <- paste("-G", REV, "-A", FWD.RC, "--minimum-length 100")

# run cutadapt
for(i in seq_along(fnFs8)){
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2,
                             "-o", fnFs.cut8[i], "-p", fnRs.cut8[i],
                             fnFs.filtN8[i], fnRs.filtN8[i]))
}

dir.create(filter.fp8)
subF.fp8 <- file.path(filter.fp8, "preprocessed_F")
subR.fp8 <- file.path(filter.fp8, "preprocessed_R")
dir.create(subF.fp8)
dir.create(subR.fp8)

fnFs.Q8 <- file.path(subF.fp8, basename(fnFs8))
fnRs.Q8 <- file.path(subR.fp8, basename(fnRs8))
file.rename(from = fnFs.cut8, to = fnFs.Q8)
file.rename(from = fnRs.cut8, to = fnRs.Q8)

filtpathF8 <- file.path(subF.fp8, "filtered")
filtpathR8 <- file.path(subR.fp8, "filtered")
fastqFs8 <- sort(list.files(subF.fp8, pattern = "fastq.gz"))
fastqRs8 <- sort(list.files(subR.fp8, pattern = "fastq.gz"))

if(length(fastqFs8) != length(fastqRs8)) stop("Forward and reverse files do not match!")

# filter the data
filt_out8 <- filterAndTrim(fwd = file.path(subF.fp8, fastqFs8), filt = file.path(filtpathF8, fastqFs8),
                           rev = file.path(subR.fp8, fastqRs8), filt.rev = file.path(filtpathR8, fastqRs8),
                           truncLen = c(240,170), maxEE = c(1.8,2), truncQ = 5, maxN = 0, rm.phix = TRUE,
                           compress = TRUE, verbose = TRUE, multithread = TRUE)


filt_out8 %>% data.frame() %>%
  mutate(Samples = rownames(.), percent_kept = 100*(reads.out/reads.in)) %>%
  select(Samples, everything()) %>%
  summarise(min_remaining = paste0(round(min(percent_kept),2),"%"),
            median_remaining = paste0(round(median(percent_kept),2),"%"),
            mean_remaining = paste0(round(mean(percent_kept),2),"%"),
            max_remaining = paste0(round(max(percent_kept),2),"%"))


# infer sequence variants
filtFs8 <- list.files(filtpathF8, pattern = "fastq.gz", full.names = TRUE)
filtRs8 <- list.files(filtpathR8, pattern = "fastq.gz", full.names = TRUE)

sample.names8 <- sapply(strsplit(basename(filtFs8), "_"), `[`,1)
sample.namesR8 <- sapply(strsplit(basename(filtRs8), "_"), `[`,1)

if(!identical(sample.names8, sample.namesR8)) stop("forward and reverse files do not match!")
names(filtFs8) <- sample.names8
names(filtRs8) <- sample.names8

set.seed(123)
errF8 <- learnErrors(filtFs8, multithread = TRUE)
errR8 <- learnErrors(filtRs8, multithread = TRUE)

mergers8 <- vector("list", length(sample.names8))
names(mergers8) <- sample.names8
ddF8 <- vector("list", length(sample.names8))
names(ddF8) <- sample.names8
ddR8 <- vector("list", length(sample.names8))
names(ddR8) <- sample.names8

for(i in sample.names8){
  derepF8 <- derepFastq(filtFs8[[i]])
  dadaF8 <- dada(derepF8, err = errF8, multithread = TRUE)
  ddF8[[i]] <- dadaF8
  
  derepR8 <- derepFastq(filtRs8[[i]])
  dadaR8 <- dada(derepR8, err = errR8, multithread = TRUE)
  ddR8[[i]] <- dadaR8
  
  merger8 <- mergePairs(ddF8[[i]], derepF8, ddR8[[i]], derepR8)
  mergers8[[i]] <- merger8
}

rm(derepF8); rm(derepR8)
seqtab8 <- makeSequenceTable(mergers8)


#--------------------------------------------------------------------------------
mergetab <- mergeSequenceTables(seqtab1,
                                seqtab2,
                                seqtab3,
                                seqtab4,
                                seqtab5,
                                seqtab6,
                                seqtab7,
                                seqtab8)

saveRDS(mergetab, paste0(table.fp, "/mergetab.rds"))
# remove chimeras and assign taxonomy
st.all <- readRDS(paste0(table.fp, "/mergetab.rds"))

seqtab.nochim <- removeBimeraDenovo(st.all, method = "consensus", multithread=TRUE)
100*sum(seqtab.nochim)/sum(mergetab)

tax <- assignTaxonomy(seqtab.nochim, "~/db/rdp_train_set_18.fa", tryRC = TRUE, multithread = TRUE, minBoot = 50)
saveRDS(seqtab.nochim, paste0(table.fp, "/mergetab_nochim_final.rds"))
saveRDS(tax, paste0(table.fp, "/tax_final.rds"))

seqtab.t <- as.data.frame(t(seqtab.nochim))
rep_set_ASVs <- as.data.frame(rownames(seqtab.t))
rep_set_ASVs <- mutate(rep_set_ASVs, ASV_ID = 1:n()) 
rep_set_ASVs$ASV_ID <- sub("^", "ASV_", rep_set_ASVs$ASV_ID)
rep_set_ASVs$ASV <- rep_set_ASVs$`rownames(seqtab.t)`
rep_set_ASVs$`rownames(seqtab.t)` <- NULL

rownames(seqtab.t) <- rep_set_ASVs$ASV_ID
taxonomy <- as.data.frame(tax)
taxonomy$ASV <- as.factor(rownames(taxonomy))
taxonomy <- merge(rep_set_ASVs, taxonomy, by="ASV")
rownames(taxonomy) <- taxonomy$ASV_ID
taxonomy_for_mctoolsr <- unite_(taxonomy, "taxonomy",
                                c("Kingdom","Phylum","Class","Order","Family","Genus","ASV_ID"),
                                sep = ";")

writeRepSetFasta <- function(data, filename){
  fastaLines = c()
  for(rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"name"], sep = "")))
    fastaLines = c(fastaLines, as.character(data[rowNum, "seq"]))
  }
  fileConn <- file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

taxonomy_for_fasta <- taxonomy %>%
  unite("TaxString", c("Kingdom","Phylum","Class", "Order","Family","Genus","ASV_ID"),
        sep = ";", remove = FALSE) %>%
  unite("name", c("ASV_ID", "TaxString"),
        sep = " ", remove = TRUE) %>% 
  select(ASV, name) %>%
  rename(seq = ASV)

writeRepSetFasta(taxonomy_for_fasta, paste0(table.fp, "/repset.fasta"))

seqtab_wTax <- merge(seqtab.t, taxonomy_for_mctoolsr, by=0)
seqtab_wTax$ASV <- NULL
out_fp <- paste0(table.fp,"/seqtab_wTax_mctoolsr.txt")
names(seqtab_wTax)[1] = "#ASV_ID"
write("#Exported for mctoolsr", out_fp)
suppressWarnings(write.table(seqtab_wTax, out_fp, sep = "\t", row.names = FALSE, append = TRUE))

write.table(seqtab.t, file = paste0(table.fp, "/seqtab_final.txt"),
            sep = "\t", row.names = TRUE, col.names = NA)
write.table(tax, file = paste0(table.fp, "/tax_final.txt"), 
            sep = "\t", row.names = TRUE, col.names = NA)



