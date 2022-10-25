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
library(Biostrings) ; packageVersion("Biostrings")
library(kableExtra) ; packageVersion("kableExtra")
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
set.seed(1234)
plotQualityProfile(raw_read)

# Place filtered files in filtered/ subdirectory
filtered_path <- file.path(path, "filtered")
filtered_read <- file.path(filtered_path, paste0(sample.names, "_filtered.fastq.gz"))

# Filter and Trimming:
out <- filterAndTrim(raw_read, filtered_read, truncLen = 250, maxN = 0,
              maxEE = 2,  truncQ = 2, rm.phix = TRUE, compress = TRUE,
              multithread = TRUE)
head(out)
plotQualityProfile(filtered_read[1:6])

# Learn the error rates: The DADA2 algorithm depends on a parametric error model
# and every amplicon dataset has a slightly different error rate. 
err <- learnErrors(filtered_read, nbases = 1e7, multithread = TRUE)
plotErrors(err, nominalQ = TRUE) + theme_minimal()

# Dereplication combines all identical sequencing reads into "uniques sequences"
# with a corresponding "abundance"
derep <- derepFastq(filtered_read, verbose = TRUE)
names(derep) <- sample.names

# Apply the core sequence-variant inference algorithm to the dereplicated data
dada_reads <-dada(derep, err = err, multithread = TRUE)
dada_reads[[1]]

# Construct sequence table: This table is a matrix with each row representing the samples, columns are the various ASVs, and each cell shows the number of that specific ASV within each sample.
seq_table <- makeSequenceTable(dada_reads)
dim(seq_table)

# Remove chimeras: dada2 will align each ASV to the other ASVs, and if an ASV’s left and right side align to two separate more abundant ASVs, the it will be flagged as a chimera and removed.
seq_table.nochim <- removeBimeraDenovo(seq_table, method = "consensus",
                                       multithread=TRUE, verbose = TRUE)
dim(seq_table.nochim)
# Write sequence table to disk
saveRDS(seq_table.nochim, "./output/seq_table.nochim")

# Final check the progress so far
get_n <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dada_reads, get_n), rowSums(seq_table),
               rowSums(seq_table.nochim))
colnames(track) <- c("input", "filtered", "denoised", "tabled", "nochim")
rownames(track) <- sample.names
kableExtra::kable(track)

# Assign Taxonomy :
## rdp training set
taxa <- assignTaxonomy(seq_table.nochim,
                       "./database/rdp_train_set_18.dada2.fa.gz",
                       multithread = TRUE)
taxa_rdp <- addSpecies(taxa, "./database/rdp_species_assignment_18.dada2.fa.gz")
saveRDS(taxa_rdp, "./output/taxa.rdp") #Write taxa to disk

## silva database 
taxa.silva <- assignTaxonomy(seq_table.nochim,
                             "./database/silva_nr99_v138.1_train_set.dada2.fa.gz",
                             multithread = TRUE)
taxa.silva <- addSpecies(taxa.silva, "./database/silva_species_assignment_v138.1.dada2.fa.gz")
saveRDS(taxa.silva, "./output/taxa.silva") #Write taxa to disk

taxa_print <- taxa_rdp
rownames(taxa_print) <- NULL
saveRDS(taxa_print, "./output/taxa_print.rdp") #Write taxa to disk

taxa_print.silva <- taxa.silva
rownames(taxa_print.silva) <- NULL
saveRDS(taxa_print.silva, "./output/taxa_print.silva") #Write taxa to disk

### ALTERNATIVES: `IdTaxa` taxonomic classification via DECIPHER
#dna <- DNAStringSet(getSequences(seq_table.nochim)) # create a DNA string set from the ASVs
#load("./database/RDP_v18-mod_July2020.RData")
#ids <- IdTaxa(dna, trainingSet, strand = "top", processors = NULL, verbose = FALSE)
#ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
#taxid <- t(sapply(ids, function(x) {
#  m <- match(ranks, x$rank)
#  taxa <- x$taxon[m]
#  taxa[startsWith(taxa, "unclassified_")] <- NA
#  taxa
#}))
#colnames(taxid) <- ranks ; rownames(taxid) <- getSequences(seq_table.nochim)

# Phylogenetic Tree
sequences <- getSequences(seq_table.nochim)
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
theme_set(theme_bw())
sample_data <-read.csv("./raw_data/sample_data.tsv", header = TRUE, sep = "\t", row.names = "Run")
physeq <- phyloseq(otu_table(seq_table.nochim, taxa_are_rows = FALSE),
                   sample_data(sample_data),
                   tax_table(taxa_rdp),
                   phy_tree(fitGTR$tree))
physeq <- prune_samples(sample_names(physeq) != "Mock", physeq)

# Rename ASVs to "ASV1, ASV2..."
# Store the DNA sequences of our ASVs in the refseq slot of the phyloseq object
dna <- Biostrings::DNAStringSet(taxa_names(physeq))
names(dna) <- taxa_names(physeq)
physeq <- merge_phyloseq(physeq, dna)
taxa_names(physeq) <- paste0("ASV", seq(ntaxa(physeq)))
physeq

rank_names(physeq)
table(tax_table(physeq)[, "Phylum"])

plot_richness(physeq, x="Treatment", measures = c("Shannon", "Fisher", "Simpson"))
physeq.prop <- transform_sample_counts(physeq, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(physeq.prop, method = "NMDS", distance = "bray")
plot_ordination(physeq.prop, ord.nmds.bray, title = "Bray NMDS")

top20 <- names(sort(taxa_sums(physeq), decreasing = TRUE))[1:20]
physeq_top20 <- transform_sample_counts(physeq, function(OTU) OTU/sum(OTU))
physeq_top20 <- prune_taxa(top20, physeq_top20)
plot_bar(physeq_top20, x="Treatment", fill = "Genus") + 
  facet_wrap(~When, scales = "free_x")




###EXAMPLE of creating a metadata for phyloseq
#samples.out <- rownames(seq_table.nochim)
#subject <- sapply(strsplit(samples.out, "D"), `[`,1)
#gender <- substr(subject,1,1)
#subject <- substr(subject,2,999)
#day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
#samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
#samdf$When <- "Early"
#samdf$When[samdf$Day>100] <- "Late"
#rownames(samdf) <- samples.out
#head(seq_table.nochim)
