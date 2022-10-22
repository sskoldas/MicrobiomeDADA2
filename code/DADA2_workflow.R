# This pipeline is prepared for MiSeq single-end fastq files that have been demultiplexed
# by sample and from which the adapters and barcodes have already been removed.
# The end products is an amplicon sequence variant (ASV) table and also assigned to the taxonomy, 
# and demonstrated into `phyloseq` R package for analysis of 16S -V4 region- microbiome data.
# for more information , please check this link out : https://benjjneb.github.io/dada2/tutorial.html


library(dada2) ; packageVersion("dada2")
library(ggplot2) ; packageVersion("ggplot2")
library(phyloseq) ; packageVersion("phyloseq")
library(phangorn) ; packageVersion("phangorn")
library(DECIPHER) ; packageVersion("DECIPHER")
library(tidyverse) ; packageVersion("tidyverse")

# Set the path to the fastq files:
path <- "./raw_data"
head(list.files(path))

# File name parsing: our filename's format is SAMPLENAME.fastq
raw_read <- sort(list.files(path, pattern = ".fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(raw_read), ".fastq"),
                       `[`, 1) #extract the first element of a subset
head(sample.names)

# Check the quality of first six reads
plotQualityProfile(raw_read[1:6])

# Place filtered files in filtered/ subdirectory
filtered_path <- file.path(path, "filtered")
filtered_read <- file.path(filtered_path, paste0(sample.names, "_filtered.fastq.gz"))

# Filter and Trimming:
out <- filterAndTrim(raw_read, filtered_read, truncLen = 240, maxN = 0,
              maxEE = 3,  truncQ = 2, rm.phix = TRUE, compress = TRUE,
              multithread = TRUE)
head(out)
plotQualityProfile(filtered_read[1:6])

# Learn the error rates: The DADA2 algorithm depends on a parametric error model
# and every amplicon dataset has a slightly different error rate. 
err <- learnErrors(filtered_read, multithread = TRUE)
plotErrors(err, nominalQ = TRUE) + theme_minimal()

# Dereplication combines all identical sequencing reads into "uniques sequences"
# with a corresponding "abundance"
derep <- derepFastq(filtered_read, verbose = TRUE)
names(derep) <- sample.names

# Apply the core sequence-variant inference algorithm to the dereplicated data
dada_reads <-dada(derep, err = err, multithread = TRUE)
dada_reads[[1]]

# Construct sequence table
seq_table <- makeSequenceTable(dada_reads)
dim(seq_table)

# Remove chimeras
seq_table.nochim <- removeBimeraDenovo(seq_table, method = "consensus",
                                       multithread=TRUE, verbose = TRUE)
dim(seq_table.nochim)

# Final check the progress so far
get_n <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dada_reads, get_n), rowSums(seq_table),
               rowSums(seq_table.nochim))
colnames(track) <- c("input", "filtered", "denoised", "tabled", "nochim")
rownames(track) <- sample.names
head(track)

# Assign Taxonomy :
taxa <- assignTaxonomy(seq_table.nochim,
                       "./database/rdp_train_set_18.dada2.fa.gz",
                       multithread = TRUE)
taxa <- addSpecies(taxa, "./database/rdp_species_assignment_18.dada2.fa.gz")
taxa_print <- taxa
rownames(taxa_print) <- NULL
head(taxa_print)

# Phylogenetic Tree
sequences <- getSequences(seq_table)
names(sequences) <- sequences
alignment <- AlignSeqs(DNAStringSet(sequences), anchor=NA)

phang_align <- phyDat(as(alignment, "matrix"), type = "DNA")
dm <- dist.ml(phang_align)
treeNJ <- NJ(dm)
fit = pml(treeNJ, data = phang_align)

fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model = "GTR", optInv = TRUE, optGamma = TRUE,
                    rearrangement = "stochastic",
                    control = pml.control(trace = 0))
detach("package:phangorn", unload = TRUE)

# Phyloseq
sample_data <-read.csv("./raw_data/sample_data.tsv", header = TRUE, sep = "\t", row.names = "Run")
physeq <- phyloseq(otu_table(seq_table.nochim, taxa_are_rows = FALSE),
                   sample_data(sample_data),
                   tax_table(taxa),
                   phy_tree(fitGTR$tree))
physeq <- prune_samples(sample_names(physeq) != "Mock", physeq)
plot_richness(physeq, x="Treatment", measures = c("Shannon", "Fisher")) +
  theme_minimal()

top20 <- names(sort(taxa_sums(physeq), decreasing = TRUE))[1:20]
physeq_top20 <- transform_sample_counts(physeq, function(OTU) OTU/sum(OTU))
physeq_top20 <- prune_taxa(top20, physeq_top20)

