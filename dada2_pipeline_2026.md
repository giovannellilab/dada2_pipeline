################################################################################
## Pipeline for the analysis of sequences using Dada2. Designed for medium size
## dadasets, can also be used with large dataset where multiple runs are merged
## together. The specific steps used for merging of multiple runs are highlighted
## below. For information contact Donato Giovannelli at the University of
## Naples Federico II - donato.giovannelli@unina.it - 15 Sept 2026
################################################################################

# The analysis starts from within a working forlder where one or more folders
# containing the raw fastq files are stored. The different raw reads folder in
# this pipeline are called raw_1 raw_2 etc.., and cotain the fastq.gz files.

## Loading the library used for the pipeline
library(dada2); packageVersion("dada2") # Amplicon Sequencing Variants inference
library(phyloseq) # Data analysis and visualisation
library(DECIPHER) # Multiple sequence alignment and phylogenetic analysis
library(ape) # importing and handling phylogenetic trees

## STATING THE DADA2 PIPELINE
## This pipeline is based on the pipeline published in:
## Callahan et al. 2016. Bioconductor Workflow for Microbiome Data Analysis: from raw reads to community analyses. F1000Res 5, 1492.
## Several steps have been modified and adapted to the specific analyses carried
## out in the Giovannelli Lab. Please be sure that waht we do makes sense and
## applies to your specific project. This pipeline is shared as is.





#########################################################################################
#### ANALYSIS BLOCK 1 ###################################################################
########################################################################################

path1 <- "raw_1" # CHANGE it to the directory containing the fastq files after unzipping
list.files(path1) # Verify the file list

# separare forward e reverse reads. Change the pattern if necessary
fnFs1 <- sort(list.files(path1, pattern="_R1_001.fastq", full.names = TRUE))
fnRs1 <- sort(list.files(path1, pattern="_R2_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
# The Name is selected as ending at the FIRST underscore. This can be changed
# according to the sample name selected
sample.names1 <- sapply(strsplit(basename(fnFs1), "_"), `[`, 1)

# Place filtered files in path1/filtered/ subdirectory
filtFs1 <- file.path(path1, "filtered", paste0(sample.names1, "_F_filt.fastq.gz"))
filtRs1 <- file.path(path1, "filtered", paste0(sample.names1, "_R_filt.fastq.gz"))

## TRIMMING
# Checking the sequence quality. Different samples can be inspected by changing the
# numbers within the brackets. Refer to the original publication by Callahan et al. 2016
# for interpretation

plotQualityProfile(fnFs1[1:2]) # Forward sequences
plotQualityProfile(fnRs1[1:2]) # reverse sequences

# Select where the sequences should be truncated. usually where the median quality
# (the solid organge line in the plots) starts to drop

out1 <- filterAndTrim(fnFs1, filtFs1, fnRs1, filtRs1, truncLen=c(245,165),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              trimLeft=17, compress=TRUE, multithread=TRUE)

# Verify how many reads are lost. It is important to try different parameters
# (truncLen and maxEE, truncQ) to preserve a large number of sequences, without compormising
# on the quality of the reads. In case of low quality very few sequences will
# pair in the later steps and been retained. This step is the most important
# in the entire Dada2 pipeline. feel free to try several parameter settings and
# choose the best settings.Refer to Callahan et al. 2016 for details.
head(out1)

# Once the filter and trimming has been completed to satisfaction, proceede to
# the error lerning. This is computationally intensive
# To solve mitigate the effects fo the binned quality score now implemented from the sequencing facilities, new error models have been suggested to work around it: https://github.com/benjjneb/dada2/issues/1307

loessErrfun_mod <- function (trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow = 0, ncol = length(qq))
  for (nti in c("A", "C", "G", "T")) {
    for (ntj in c("A", "C", "G", "T")) {
      if (nti != ntj) {
        errs <- trans[paste0(nti, "2", ntj), ]
        tot <- colSums(trans[paste0(nti, "2", c("A",
                                                "C", "G", "T")), ])
        rlogp <- log10((errs + 1)/tot)
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q = qq, errs = errs, tot = tot,
                         rlogp = rlogp)
        mod.lo <- loess(rlogp ~ q, df, weights = log10(tot),span = 2)
        pred <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred) > maxrli] <- pred[[maxrli]]
        pred[seq_along(pred) < minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      }
    }
  }
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-07
  est[est > MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est < MIN_ERROR_RATE] <- MIN_ERROR_RATE
  err <- rbind(1 - colSums(est[1:3, ]), est[1:3, ], est[4,
                                                        ], 1 - colSums(est[4:6, ]), est[5:6, ], est[7:8, ], 1 -
                 colSums(est[7:9, ]), est[9, ], est[10:12, ], 1 - colSums(est[10:12,
                                                                              ]))
  rownames(err) <- paste0(rep(c("A", "C", "G", "T"), each = 4),
                          "2", c("A", "C", "G", "T"))
  colnames(err) <- colnames(trans)
  return(err)
}


errF1 <-learnErrors(filtFs1, errorEstimationFunction=loessErrfun_mod, multithread=TRUE)
errR1 <- learnErrors(filtRs1, errorEstimationFunction=loessErrfun_mod,multithread=TRUE)


plotErrors(errF1, nominalQ=TRUE)

# Dereplication step
derepFs1 <- derepFastq(filtFs1, verbose=TRUE)
derepRs1 <- derepFastq(filtRs1, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs1) <- sample.names1
names(derepRs1) <- sample.names1

# Infer the ASVs using the error prifiles ;earnined in the previous steps
dadaFs1 <- dada(derepFs1, err=errF1, pool=TRUE, errorEstimationFunction = loessErrfun_mod, multithread=TRUE) 
dadaRs1 <- dada(derepRs1, err=errR1, pool=TRUE, errorEstimationFunction = loessErrfun_mod, multithread=TRUE)

# Mate Pairing the reads
mergers1 <- mergePairs(dadaFs1, derepFs1, dadaRs1, derepRs1, verbose=TRUE)

# Create the ASV table from run_1
seqtab1 <- makeSequenceTable(mergers1)
dim(seqtab1) # Gives you info on the number of Samples and ASVs identified in run_1

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab1)))

# Drop the tail sequences that are either too short or too long compared to the majority
# Chage the seq numbers within brackets to match the min and max lenght desired
seqtab1 <- seqtab1[,nchar(colnames(seqtab1)) %in% seq(371,383)]

# Inspect distribution of sequence lengths after dropping the outliers
table(nchar(getSequences(seqtab1)))
dim(seqtab1) # Gives you info on the number of Samples and ASVs after tail sequence dropping

# Sanity check for run_1
getN <- function(x) sum(getUniques(x))
track1 <- cbind(out1, sapply(dadaFs1, getN), sapply(dadaRs1, getN), sapply(mergers1, getN))
colnames(track1) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
rownames(track1) <- sample.names1
write.csv("run1_asv_stats.csv", track1) # Export a csv file as archive
head(track1) # check the number of reads retained at each step

#### END OF ANALYSIS BLOCK 1 ##############################################################
## From here on you can either move to ANALYSIS BLOCK 3 if you have a single run to work with,
## or continue to ANALYSIS BLOCK 2 below if you need to analyze and merge multiple
## runs together
###########################################################################################





###########################################################################################
#### ANALYSIS BLOCK 2 #####################################################################
###########################################################################################

## Repeat analysis block 1 for every run you have, changing all the names from 1 to two
## for example run_1 to run_2 for the path, sample.names1 to sample.names2 or seqtab1 to
## seqtab2. ## For simplicity I have already transformed and copied the commands for run_2
## removing the comments below this line. The REPEAT block is highligthed wiht indents.
## Two comments are present only near the command where it is fundamental to interact with the script.

#### REPEAT FROM HERE <--------------------------------------------------------!
    path2 <- "raw_2"
    list.files(path2)
    fnFs2 <- sort(list.files(path2, pattern="_R1_001.fastq", full.names = TRUE))
    fnRs2 <- sort(list.files(path2, pattern="_R2_001.fastq", full.names = TRUE))
    sample.names2 <- sapply(strsplit(basename(fnFs2), "_"), `[`, 1)
    filtFs2 <- file.path(path2, "filtered", paste0(sample.names2, "_F_filt.fastq.gz"))
    filtRs2 <- file.path(path2, "filtered", paste0(sample.names2, "_R_filt.fastq.gz"))
    plotQualityProfile(fnFs2[1:2])
    plotQualityProfile(fnRs2[1:2])
    # Select where the sequences should be truncated iteratively if necessary
    out2 <- filterAndTrim(fnFs2, filtFs2, fnRs2, filtRs2, truncLen=c(245,165),
            maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
            trimLeft=17, compress=TRUE, multithread=TRUE)
    head(out2)
    errF2 <- learnErrors(filtFs2, errorEstimationFunction=loessErrfun_mod, multithread=TRUE, randomize=TRUE)
    errR2 <- learnErrors(filtRs2, errorEstimationFunction=loessErrfun_mod, multithread=TRUE, randomize=TRUE)
    plotErrors(errF2, nominalQ=TRUE)
    derepFs2 <- derepFastq(filtFs2, verbose=TRUE)
    derepRs2 <- derepFastq(filtRs2, verbose=TRUE)
    names(derepFs2) <- sample.names2
    names(derepRs2) <- sample.names2
    dadaFs2 <- dada(derepFs2, err=errF2, pool="pseudo", errorEstimationFunction = loessErrfun_mod, multithread=TRUE)
    dadaRs2 <- dada(derepRs2, err=errR2, pool="pseudo", errorEstimationFunction = loessErrfun_mod, multithread=TRUE)
    mergers2 <- mergePairs(dadaFs2, derepFs2, dadaRs2, derepRs2, verbose=TRUE)
    seqtab2 <- makeSequenceTable(mergers2)
    dim(seqtab2)
    table(nchar(getSequences(seqtab2)))
    seqtab2 <- seqtab2[,nchar(colnames(seqtab2)) %in% seq(371,383)] # Drop the tail sequences
    table(nchar(getSequences(seqtab2)))
    dim(seqtab2)
    track2 <- cbind(out2, sapply(dadaFs2, getN), sapply(dadaRs2, getN), sapply(mergers2, getN))
    colnames(track2) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
    rownames(track2) <- sample.names2
    write.csv("run2_asv_stats.csv", track2) # Export a csv file as archive
    head(track2) # check the number of reads retained at each step for run_2
#### TO HERE <------------------------------------------------------------------!
## Repeat As many time as necessary changing the number in the names. The follow the lines below

## Merging of multiple runs. Add as many seqtab.nochim with the approapriate number
## separated by a comma as needed.
seqtab <- mergeSequenceTables(seqtab1, seqtab2) # Merging of different runs

#### END OF ANALYSIS BLOCK 2 ##############################################################
## Continue the analysis on ANALYSIS BLOCK 3 below
###########################################################################################





###########################################################################################
#### ANALYSIS BLOCK 3 #####################################################################
###########################################################################################

## From here on we refer to the seqtab.nochim object otained from the merger done in
## ANALYSIS BLOCK 2. Change the name to seqtab.nochim1 if you have a single run or
## feel free to renave you objest to seqtab.nochim to avoid changing the code below

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab) # Gives you the percentage of sequences recovered after chimera removal
## This is a ggod check points. Even if a lot of ASVs have been removed, the majority of reads sould
## be retained. Usually >0.80 (aka 80%) are retained

## Assign Taxonomy. Point to where the silva database actually is
taxa <- assignTaxonomy(seqtab.nochim, "~/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "~/silva_species_assignment_v132.fa.gz")
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
write.csv(taxa.print, "taxa_print.csv") # For inspection in bash or excel

## Phylogenetic tree building. MSA and Phangorn used in the original publication
## scale quadratically with the number of ASV, and became quickly unisable. I have moved to
## AlignSeqs in the Dechipher package and FastTree in bash for making the phylogenetic tree
## Before continuing be sure that the fasttree execuatble is installed in your system and
## visible in your PATH

seqs <- getSequences(seqtab.nochim) ## Get the sequences
names(seqs) <- seqs
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA) # Generating the multiple sequence alignment
writeXStringSet(alignment, "alignment.fasta", format="fasta") # Exporting the alignment to an external file

## Build the tree using FastTree in bash using the GTR model. See more detail about
## FastTree at http://www.microbesonline.org/fasttree/
system('fasttree -gtr -nt alignment.fasta > alignment.tree', intern=TRUE)

tree <- read.tree(file = "alignment.tree") # Reading back the tree into the R environment

## Sanity check. The two numbers should match!
dim(seqtab.nochim)
tree$Nnode+2

## Import the environmental data for the project. Be sure that the first column
## of the csv file contains the sample names as the EXACTLY appear on the seqtab.nochim
## object

env_data <- read.csv("environmental_data.csv", header=T, sep=",", row.names=1)

summary(env_data) # check if some variables need to be converted in to factors using as.factor()

################################################################################
## Create the final phyloseq object
################################################################################

prok_data <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = F),
                      phy_tree(tree),
                      sample_data(env_data),
                      tax_table(taxa))

prok_data # Info on the phyloseq obect

################################################################################################
################## END OF THE DADA2 PIPELINE FOR THE GIOVANNELLI LAB############################
################################################################################################
## From here you can just keep analyzing your dataset. Take a look at our basic phyloseq
## pipeline or get the analysis file from the repository of one of our paper on GitHub
################################################################################################
