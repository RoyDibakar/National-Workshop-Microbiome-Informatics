library(dada2)
#Getting ready
# First set directory, where you have downloaded the raw FASTQ files
path <- "C:/Users/User/Desktop/Conference_Microbiome_Informatics/Raw_FASTQ_Files" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnFs
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
fnRs

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names

#Inspect read quality profiles
#We start by visualizing the quality profiles of the forward reads
plotQualityProfile(fnFs[1:2])
#Now we visualize the quality profile of the reverse reads
plotQualityProfile(fnRs[1:2])

#Filter and trim
#Assign the filenames for the filtered fastq.gz files.

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
filtFs
filtRs

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160), maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out)

#Learn the Error Rates
errF <- learnErrors(filtFs, multithread=FALSE)
errR <- learnErrors(filtRs, multithread=FALSE)
#plotErrors(errF, nominalQ=TRUE)

#Sample Inference
dadaFs <- dada(filtFs, err=errF, multithread=FALSE)
dadaRs <- dada(filtRs, err=errR, multithread=FALSE)

dadaFs[[1]]

#Merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

#View(mergers$F3D0)
#Construct sequence table
seqtab <- makeSequenceTable(mergers)
# View(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE, verbose=TRUE)
#dim(seqtab.nochim)
#sum(seqtab.nochim)/sum(seqtab)

#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#Assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138.1_train_set.fa.gz", multithread=FALSE)

#Let’s inspect the taxonomic assignments:
View(taxa)  
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

#Downstream analysis
library(phyloseq)

samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
gender <- substr(subject,1,1)
subject <- substr(subject,2,999)
day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"
rownames(samdf) <- samples.out

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), sample_data(samdf), tax_table(taxa))

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

View(ps@otu_table)
View(ps@tax_table)
View(ps@sam_data)

#Off to microeco
library(file2meco)
library(microeco)
library(ggplot2)
library(ggdendro)
library(ggpubr)
library(magrittr)

meco_object <- phyloseq2meco(ps)
meco_object$tidy_dataset()
#Filter data
meco_object$filter_pollution(taxa = c("mitochondria", "chloroplast"))

#Alpha Diversity
t1 <- trans_alpha$new(dataset = meco_object, group = "When")
t1$cal_diff(method = "wilcox")
# y_increase: increased height for each label
t1$plot_alpha(measure = "Shannon", shape = "When", y_start = 0.1, y_increase = 0.1, xtext_size = 15) + theme(axis.text.y = element_text(size = 12, vjust = 0.5, face = "bold"), axis.text.x = element_text(size = 16, face = "bold"),panel.background = element_blank(), axis.title = element_text(size = 12, face = "bold"), axis.line = element_line(colour = "black"), legend.text = element_text(size = 12))
ggsave("alpha.tiff", height = 6, width = 6, dpi = 700)

#Beta Diversity.
meco_object$cal_betadiv()
t2 <- trans_beta$new(dataset = meco_object, group = "When", measure = "bray")
# PCoA, PCA, DCA and NMDS are available
t2$cal_ordination(method = "PCoA")
# plot the PCoA result with confidence ellipse
t2$plot_ordination(plot_color = "When", plot_shape = "When", plot_type = c("point", "ellipse")) + theme(axis.text.y = element_text(size = 12, vjust = 0.5, face = "bold"), axis.text.x = element_text(size = 16, face = "bold"),panel.background = element_blank(), axis.title = element_text(size = 12, face = "bold"), axis.line = element_line(colour = "black"), legend.text = element_text(size = 12))
ggsave("beta.tiff", height = 6, width = 6, dpi = 700)

# #Clustering
# meco_object$sample_table %<>% subset(When %in% c("Early", "Late"))
# meco_object$tidy_dataset()
# t3 <- trans_beta$new(dataset = meco_object, group = "When")
# # use replace_name to set the label name, group parameter used to set the color
# t3$plot_clustering(group = "When", replace_name = c("When"))
# ggsave("cluster.tiff", height = 6, width = 6, dpi = 700)

#Microbial composition
t4 <- trans_abund$new(dataset = meco_object, taxrank = "Genus", ntaxa = 10)
t4$plot_bar(others_color = "grey70", facet = "When", xtext_keep = FALSE, legend_text_italic = TRUE)
ggsave("abundance_1.tiff", height = 6, width = 10, dpi = 700)

# The groupmean parameter can be used to obtain the group-mean barplot.
t5 <- trans_abund$new(dataset = meco_object, taxrank = "Genus", ntaxa = 10, groupmean = "When")
g1 <- t5$plot_bar(bar_full = FALSE, legend_text_italic = TRUE)+theme(axis.text.y = element_text(size = 12, vjust = 0.5, face = "bold"), axis.text.x = element_text(size = 16, face = "bold"),panel.background = element_blank(), axis.title = element_text(size = 12, face = "bold"), axis.line = element_line(colour = "black"), legend.text = element_text(size = 12))
g1
ggsave("Abundance.tiff", height = 6, width = 6, dpi = 700)

# #Heatmap
# t6 <- trans_abund$new(dataset = meco_object, taxrank = "Genus", ntaxa = 15)
# g2 <- t6$plot_heatmap(facet = "When", xtext_keep = FALSE, withmargin = FALSE, plot_breaks = c(0.01, 0.1, 1, 10))
# g2

# #Boxplot
# t7 <- trans_abund$new(dataset = meco_object, taxrank = "Genus", ntaxa = 15)
# t7$plot_box(group = "When", xtext_angle = 90) +theme(axis.text.y = element_text(size = 12, vjust = 0.5, face = "bold"), axis.text.x = element_text(size = 10, face = "bold"),panel.background = element_blank(), axis.title = element_text(size = 12, face = "bold"), axis.line = element_line(colour = "black"), legend.text = element_text(size = 12))

#Differential analysis
#LEfSE
t8 <- trans_diff$new(dataset = meco_object, method = "lefse", group = "When", alpha = 0.05, lefse_subgroup = NULL, taxa_level = "Genus", p_adjust_method = "fdr")
# see t1$res_diff for the result
# From v0.8.0, threshold is used for the LDA score selection.
g3 <- t8$plot_diff_bar(threshold = 2)
g3
ggsave("LDA.tiff", height = 6, width = 10, dpi = 700)

#ALDEx2_kw
t9 <- trans_diff$new(dataset = meco_object, method = "ALDEx2_kw", group = "When", filter_thres = 0.01, taxa_level = "Genus")
t9$plot_diff_abund(add_sig = TRUE, use_number = 1:12, simplify_names = FALSE)
ggsave("ALDEx2_kw.tiff", height = 8, width = 12, dpi = 700)


