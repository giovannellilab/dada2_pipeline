#loading libraries
library(phyloseq)
library(tidyverse)
library(ape)
library(microbiome)
library(ggthemes) # additional themes fro ggplot2
library(ggpubr)
library(vegan)
library(repr)
library(ggpmisc) #to use stat_poly_eq
library(RColorBrewer) # nice color options
library(gridExtra) # gridding plots
library(viridis)
library(ggrepel)
library(reshape2)
library(knitr)
set.seed(10000)
options(repr.plot.width=16, repr.plot.height=12)

theme_glab <- function(base_size = 25,
                    base_family = "",
                    base_line_size = base_size / 180,
                    base_rect_size = base_size / 180) {
   
    font <- "Helvetica" #assign font family up front
   
    theme_bw(base_size = base_size,
                base_family = base_family,
                base_line_size = base_line_size) %+replace%
    theme(
        legend.background =  element_blank(),
        legend.title =       element_text(color = rgb(100, 100, 100, maxColorValue = 255),
                                          size = rel(0.65),
                                         hjust = 0),
        legend.text =        element_text(color = rgb(100, 100, 100, maxColorValue = 255),
                                          size = rel(0.65)),
        legend.key.size =    unit(0.8, "lines"),
     
      plot.title = element_text(
        color = rgb(100, 100, 100, maxColorValue = 255),
        hjust = 0),
       
      axis.title = element_text(
        color = rgb(100, 100, 100, maxColorValue = 255),
        size = rel(0.65)),
      axis.text = element_text(
        color = rgb(100, 100, 100, maxColorValue = 255),
        size = rel(0.65)),
       
      plot.caption = element_text(
        color = rgb(100, 100, 100, maxColorValue = 255),
        size = rel(0.7),
        hjust = 1),
       
      panel.grid.major = element_blank(),  
      panel.grid.minor = element_blank(),  
      panel.border = element_rect(fill = NA, colour = rgb(100, 100, 100, maxColorValue = 255)),

     
      complete = TRUE
    )
}

# loading dataset
env_data<-read.csv("dataset/env_data.csv",header=T,sep=",",row.names=2) #if you used 'tab' as separator, use "\t" instead of ","
head(env_data)

#Since all columns are read as character, change in numeric the one you need
env_data<-cbind(env_data[,1:9],env_data[,10:32] %>% mutate_if(is.character,as.numeric))
env_data

#loading .rds files
seqtab_nochim<-readRDS("rds/seqtab_nochim.rds")
taxa<-readRDS("rds/taxa.rds")
tree<-readRDS("rds/tree.rds")

prok_data_raw<-phyloseq(sample_data(env_data),
                        otu_table(seqtab_nochim, taxa_are_rows = F),
                        tax_table(taxa),
                        phy_tree(tree))
prok_data_raw
message("Total number of reads:")
sum(readcount(prok_data_raw))

# Define a function to remove negative controls

prune_negatives = function(physeq, negs, samps) {
  negs.n1 = prune_taxa(taxa_sums(negs)>=1, negs)
  samps.n1 = prune_taxa(taxa_sums(samps)>=1, samps)
  allTaxa <- names(sort(taxa_sums(physeq),TRUE))
  negtaxa <- names(sort(taxa_sums(negs.n1),TRUE))
  taxa.noneg <- allTaxa[!(allTaxa %in% negtaxa)]
  return(prune_taxa(taxa.noneg,samps.n1))
}

# Removing the ASV found in the negative controls from all the samples¶
blanks = subset_samples(prok_data_raw, sample_names(prok_data_raw) == "blk")
samples = subset_samples(prok_data_raw, sample_names(prok_data_raw) != "blk")
prok_data = prune_negatives(prok_data_raw,blanks,samples)

# Check the stats after removing the blanks
prok_data
sum(readcount(prok_data))
(sum(readcount(prok_data))/sum(readcount(prok_data_raw)))*100# check percentage of reads after blank removal

# Clean up unwanted sequences from Eukarya, mitochrondria and chloroplast
prok_data <- subset_taxa(prok_data,  (Kingdom != "Eukaryota") | is.na(Kingdom))
prok_data <- subset_taxa(prok_data, (Order!="Chloroplast") | is.na(Order))
prok_data <- subset_taxa(prok_data, (Family!="Mitochondria") | is.na(Family))
prok_data
sum(readcount(prok_data))
(sum(readcount(prok_data))/sum(readcount(prok_data_raw)))*100# check percentage of reads after cleaning step II

## Removing the known DNA Extraction contaminants from Sheik et al., 2018
# Assuming you are starting from a phyloseq object called prok_data_raw
prok_data_contam <- subset_taxa(prok_data,  (Genus != "Afipia") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Aquabacterium") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Asticcacaulis") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Aurantimonas") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Beijerinckia") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Bosea") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Bradyrhizobium") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Brevundimonas") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Caulobacter") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Craurococcus") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Devosia") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Hoefleae") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Mesorhizobium") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Methylobacterium") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Novosphingobium") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Ochrobactrum") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Paracoccus") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Pedomicrobium") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Phyllobacterium") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Rhizobium") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Sphingobium") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Sphingomonas") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Sphingopyxis") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Acidovorax") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Azoarcus") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Azospira") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Burkholderia") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Comamonas") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Cupriavidus") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Curvibacter") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Delftiae") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Duganella") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Herbaspirillum") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Janthinobacterium") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Kingella") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Leptothrix") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Limnobacter") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Massilia") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Methylophilus") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Methyloversatilis") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Oxalobacter") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Pelomonas") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Polaromonas") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Neisseria") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Ralstonia") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Schlegelella") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Sulfuritalea") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Undibacterium") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Variovorax") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Acinetobactera") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Enhydrobacter") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Enterobacter") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Escherichia") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Nevskia") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Pasteurella") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Pseudoxanthomonas") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Psychrobacter") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Stenotrophomonas") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Xanthomonas") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Aeromicrobium") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Actinomyces") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Arthrobacter") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Beutenbergia") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Brevibacterium") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Corynebacterium") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Curtobacterium") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Dietzia") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Janibacter") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Kocuria") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Microbacterium") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Micrococcus") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Microlunatus") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Patulibacter") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Propionibacterium") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Rhodococcus") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Tsukamurella") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Chryseobacterium") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Dyadobacter") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Flavobacterium") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Hydrotalea") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Niastella") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Parabacteroides") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Pedobacter") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Prevotella") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Wautersiella") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Deinococcus") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Abiotrophia") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Bacillus") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Brevibacillus") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Brochothrix") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Facklamia") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Olivibacter") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Lactobacillus") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Paenibacillus") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Ruminococcus") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Staphylococcus") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Streptococcus") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Veillonella") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Fusobacterium") | is.na(Genus))
prok_data_contam
message("reads saved from prok_data_raw after contaminants removal")
sum(readcount(prok_data_contam))/sum(readcount(prok_data_raw)) # reads saved from prok_data_raw after contaminants removal 
message("reads saved from prok_data_prune_mit after Eukaryota, Chloroplast, and Mitochondria removal")
sum(readcount(prok_data_contam))/sum(readcount(prok_data)) # reads saved from prok_data_prune_mit after Eukaryota, Chloroplast, and Mitochondria removal

# Removing the potential human pathogens and contaminants
# This step needs to be evaluated with attention since many of these genera
# might be relevant in many environmental settings. Usually it is better to compare
# before-after removal to see what and how much you are removing. Feel free to
# experiment with the different groups and evaluate the results.

prok_data_humpath <- subset_taxa(prok_data_contam, (Genus != "Abiotrophia") |  is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Achromobacter") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Acinetobacter") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Actinobacillus") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Arcanobacterium") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Babesia") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Bifidobacterium") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Bartonella") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Bordetella") |  is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Borrelia") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Brodetella") |  is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Brucella") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Capnocytophaga") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Chlamydia") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Citrobacter") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Comamonas") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Corynebacterium_1") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Corynebacterium") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Coxiella") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Cronobacter") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Cutibacterium") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Dermatophilus") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Ehrlichia") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Order != "Enterobacteriales") | is.na(Order))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Enterococcus") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Erysipelothrix") |  is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Escherichia") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Escherichia/Shigella") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Francisella") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Gardnerella") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Granulicatella") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Haemophilus") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Hafnia") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Helicobacter") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Klebsiella") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Kocuria") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Lactococcus") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Lactobacillus") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Lawsonia") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Legionella") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Leptospira") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Listeria") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Merkel_cell") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Micrococcus") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Morganella") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Mycoviridis") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Neisseria") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Nocardia") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Pasteurella") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Plesiomonas") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Propionibacterium") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Proteus") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Providencia") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Pseudomonas") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Rhodococcus") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Rickettsiae") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Roseomonas") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Rothia") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Salmonella") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Serratia") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Shewanella") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Shigella") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Sphaerophorus") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Staphylococcus") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Stenotrophomonas") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Streptococcus") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Treponema") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Vibrio") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Yersinia") | is.na(Genus))
prok_data_humpath
sum(readcount(prok_data_humpath))
message("reads saved from prok_data_raw after contaminants removal")
sum(readcount(prok_data_humpath))/sum(readcount(prok_data_raw)) # reads saved from prok_data_raw after contaminants removal 
message("reads saved from prok_data_contam after DNA KIT EXTRAC. contaminants removal")
sum(readcount(prok_data_humpath))/sum(readcount(prok_data_contam)) # reads saved from prok_data_contam after DNA KIT EXTRAC. contaminants removal

#Sanity check
prok_data2 = filter_taxa(prok_data_humpath, function(x) sum(x) > 0, TRUE)
prok_data2
sum(readcount(prok_data2))
(sum(readcount(prok_data2))/sum(readcount(prok_data_raw)))*100

# Normalize the counts across the different samples by converting the abundance to relative abundance 
# and multiply by the median library size
prok_ndata <- transform_sample_counts(prok_data2, function(x) ((x / sum(x))*median(readcount(prok_data2))))
prok_ndata
sum(readcount(prok_ndata))

# Transform normalized abundance to relative abundances for plotting and some stats
prok_ra = transform_sample_counts(prok_ndata, function(x){x / sum(x)})
prok_ra

#estimating rarefaction curves
#svg("../plots/svg/rarefaction_curves.svg", width=12,height=12)
rarecurve(as.matrix(data.frame(otu_table(prok_data_raw))), col='black', step=100 , lwd=3, ylab="ASVs", label=T)
#dev.off()

#Checking column names
print(colnames(env_data))

plot_richness(prok_data_raw, measures=c("Observed", "Shannon", "Simpson"), x="type") + 
geom_boxplot(aes(fill=type), size=0.25) +  
scale_fill_viridis("Sample Type",option = "viridis", discrete=T, labels=c('Biofilm', 'Blk', 'Fluid', 'Sediment')) +
xlab("") + 
scale_x_discrete(labels=c('Biofilm', 'Blk', 'Fluid', 'Sediment')) + #Change or modify x axis tick labels
theme_glab() + theme(aspect.ratio = 0.75)

alpha_div_summary<-cbind(
    data.frame(sample_sums(prok_data_raw)),
   (estimate_richness(prok_data_raw, split = TRUE, measures = "Shannon")),
    (estimate_richness(prok_data_raw, split = TRUE, measures = "Simpson")),
    (estimate_richness(prok_data_raw, split = TRUE, measures = "Chao1")),
    (estimate_richness(prok_data_raw, split = TRUE, measures = "Observed"))
                        )
alpha_div_summary <- alpha_div_summary[,-4]
alpha_div_summary<-cbind(alpha_div_summary,alpha_div_summary[,4]/alpha_div_summary[,3])
colnames(alpha_div_summary)<-c("Reads","Shannon","Simpson","Chao1","Observed","Coverage")
alpha_div_summary[order(alpha_div_summary$Reads), ]

library(gt)
library(webshot2)

#Check https://r-graph-gallery.com/package/gt.html
alpha_div_table<-alpha_div_summary[order(alpha_div_summary$Reads), ] %>%
                  gt() %>%
                    tab_header(title = md("Supplementary Table xxx"),
                    subtitle = md("Alpha Diversity"))


gtsave(alpha_div_table, "tables/alpha_div_table.pdf")

## Agglomerate at a specific taxonomic level at the Genus level
prok_ra_genus = tax_glom(prok_ra, "Genus", NArm = FALSE)
prok_ra_family = tax_glom(prok_ra, "Family", NArm = FALSE)
prok_ra_order = tax_glom(prok_ra, "Order", NArm = FALSE)
prok_ra_class = tax_glom(prok_ra, "Class", NArm = FALSE)
prok_ra_phyla = tax_glom(prok_ra, "Phylum", NArm = FALSE)

#Saving a .csv with combined taxonomy and ASVs abundances for the selected taxonomic level
summary_object <-prok_ra_phyla 
summary_otu <- as.matrix(t(otu_table(summary_object)))
summary_tax <- as.matrix(tax_table(summary_object))
summary_comb <- cbind(summary_otu, summary_tax)
write.csv(summary_comb, "results/summary/summary_prok_ra_phylum.csv")

library(MicEco)

#Subset the relative abundance at 1%
pseq.phylum_perc1 <- ps_prune(prok_ra_phyla, min.abundance = 0.01)
pseq.class_perc1 <- ps_prune(prok_ra_class, min.abundance = 0.01)
pseq.order_perc1 <- ps_prune(prok_ra_order, min.abundance = 0.01)
pseq.family_perc1 <- ps_prune(prok_ra_family, min.abundance = 0.01)
pseq.gen_perc1 <- ps_prune(prok_ra_genus, min.abundance = 0.01)

#Its core purpose is to clean a taxonomy table (tax_table) by replacing missing (NA) values 
#with the name of the most specific known taxonomic classification available for that entry, 
#automatically adding the corresponding rank as a prefix.
pseq.phylum_perc1<-ps_tax_clean(pseq.phylum_perc1)
pseq.class_perc1<-ps_tax_clean(pseq.class_perc1)
pseq.order_perc1<-ps_tax_clean(pseq.order_perc1)
pseq.family_perc1<-ps_tax_clean(pseq.family_perc1)
pseq.gen_perc1<-ps_tax_clean(pseq.gen_perc1)

#Correct all the NA_NA with NA, avoiding issues making the barplots
tax_phylum <- tax_table(pseq.phylum_perc1) %>% as.data.frame()
tax_phylum[tax_phylum == "NA_NA"] <- NA
tax_table(pseq.phylum_perc1) <- tax_table(as.matrix(tax_phylum))

tax_class <- tax_table(pseq.class_perc1) %>% as.data.frame()
tax_class[tax_class == "NA_NA"] <- NA
tax_table(pseq.class_perc1) <- tax_table(as.matrix(tax_class))

tax_order <- tax_table(pseq.order_perc1) %>% as.data.frame()
tax_order[tax_order == "NA_NA"] <- NA
tax_table(pseq.order_perc1) <- tax_table(as.matrix(tax_order))

tax_family <- tax_table(pseq.family_perc1) %>% as.data.frame()
tax_family[tax_family == "NA_NA"] <- NA
tax_table(pseq.family_perc1) <- tax_table(as.matrix(tax_family))

tax_genus <- tax_table(pseq.gen_perc1) %>% as.data.frame()
tax_genus[tax_genus == "NA_NA"] <- NA
tax_table(pseq.gen_perc1) <- tax_table(as.matrix(tax_genus))

#Barplot
plot_bar(pseq.phylum_perc1, fill ="Phylum", x='code') +
  scale_fill_manual(na.value="black","Phylum",values = colorRampPalette(brewer.pal(9, "PuBuGn"))(13)) + #Adapt the palette PuBuGn to the number of taxa you have in the phyloseq object you are plotting (i.e, 13 in my case)
  ylab("Relative Abundance (%)") + 
  scale_y_continuous(labels=c(0,25,50,75,100))  +
  theme_glab() + 
  theme(strip.text.x = element_text(size=14),
        legend.position = "bottom", aspect.ratio = 0.5)
    
plot_bar(pseq.class_perc1, fill ="Class", x='code') +
  scale_y_continuous(labels=c(0,25,50,75,100)) + ylab("Relative Abundance (%)") +
  scale_fill_manual(na.value="black","Class",values = colorRampPalette(brewer.pal(11, "Spectral"))(15)) + #Adapt the palette Spectral to the number of taxa you have in the phyloseq object you are plotting (i.e, 15 in my case) 
  theme_glab() + 
  theme(strip.text.x = element_text(size=14),
        legend.position = "bottom", aspect.ratio = 0.5)
    
plot_bar(pseq.order_perc1, fill ="Order", x='code') +
  scale_y_continuous(labels=c(0,25,50,75,100)) + ylab("Relative Abundance (%)") +
  scale_fill_manual(na.value="black","Order",values = colorRampPalette(brewer.pal(11, "PuOr"))(20)) + #Adapt the palette PuOr to the number of taxa you have in the phyloseq object you are plotting (i.e, 20 in my case)
  theme_glab() + 
  theme(strip.text.x = element_text(size=14),
        legend.position = "bottom", aspect.ratio = 0.5)

plot_bar(pseq.family_perc1, fill ="Family", x='code') +
  scale_y_continuous(labels=c(0,25,50,75,100)) + ylab("Relative Abundance (%)") +
  scale_fill_manual(na.value="black","Family",values = colorRampPalette(brewer.pal(11, "BrBG"))(22)) + #Adapt the palette BrBG to the number of taxa you have in the phyloseq object you are plotting (i.e, 22 in my case)
  theme_glab() + 
  theme(strip.text.x = element_text(size=14),
        legend.position = "bottom", aspect.ratio = 0.5)

plot_bar(pseq.gen_perc1, fill ="Genus", x='code') +
  scale_y_continuous(labels=c(0,25,50,75,100)) + ylab("Relative Abundance (%)") +
  scale_fill_manual(na.value="black","Genus",values = colorRampPalette(brewer.pal(11, "RdBu"))(23)) + #Adapt the palette RdBu to the number of taxa you have in the phyloseq object you are plotting (i.e, 23 in my case)
  theme_glab() + 
  theme(strip.text.x = element_text(size=14),
        legend.position = "bottom", aspect.ratio = 0.5)

### NMDS Jaccard similarity index: Weighted and Unweighted

prok_dist_wjac <- phyloseq::distance(prok_ndata, method = "jaccard")
prok_dist_unjac <- phyloseq::distance(prok_ndata, method = "jaccard", binary = TRUE)

# https://rdrr.io/bioc/phyloseq/man/distance.html
prok_nmds_jw <- ordinate(prok_ndata,prok_dist_wjac, method = "NMDS",trymax=10000)

# https://rdrr.io/bioc/phyloseq/man/distance.html
prok_nmds_juw <- ordinate(prok_ndata,prok_dist_unjac, method = "NMDS",trymax=10000)

plot_ordination(prok_ndata, prok_nmds_jw, type="samples",title="nMDS weighted Jaccard similarity") +
geom_text(aes(label= code), size=5, hjust=-.1,vjust=2) +
geom_point(aes(fill=type,shape=island),size=8,color="black",stroke=0.3) +
scale_fill_viridis("Sample Type",option = "viridis", discrete=T, labels=c('Biofilm', 'Fluid', 'Sediment')) +
scale_shape_manual("Island",values=c(21,25,22,24,23)) +
theme_glab() + theme(legend.position = "right") +
guides(fill = guide_legend(override.aes = list(shape = 21) ),
            shape = guide_legend(override.aes = list(fill = "black")))

#ggsave("../../xxx.svg", width=14, height=12)

plot_ordination(prok_ndata, prok_nmds_juw, type="samples",title="nMDS unweighted Jaccard similarity") +
geom_text(aes(label= code), size=5, hjust=-.1,vjust=2) +
geom_point(aes(fill=type,shape=island),size=8,color="black",stroke=0.3) +
scale_fill_viridis("Sample Type",option = "viridis", discrete=T, labels=c('Biofilm', 'Fluid', 'Sediment')) +
scale_shape_manual("Island",values=c(21,25,22,24,23)) +
theme_glab() + theme(legend.position = "right") +
guides(fill = guide_legend(override.aes = list(shape = 21) ),
            shape = guide_legend(override.aes = list(fill = "black")))

#ggsave("../../xxx.svg", width=14, height=12)

### NMDS Unifrac similarity index: Weighted and Unweighted
### To use this index it is mandatory having a tree within the phyloseq object

prok_dist_wunif <- distance(prok_ndata, method = "wunifrac")
prok_dist_ununif <- distance(prok_ndata, method = "uunifrac")

# https://rdrr.io/bioc/phyloseq/man/distance.html
prok_nmds_wuni <- ordinate(prok_ndata,prok_dist_wunif, method = "NMDS",trymax=9999)

# https://rdrr.io/bioc/phyloseq/man/distance.html
prok_nmds_ununif <- ordinate(prok_ndata,prok_dist_ununif, method = "NMDS",trymax=9999)

nmds_wUnif_p<-plot_ordination(prok_ndata, prok_nmds_wuni, type="samples",title="NMDS - Weighted UniFrac") +
geom_text(aes(label= code), size=5, hjust=-.1,vjust=2) +
geom_point(aes(fill=type,shape=island),size=8,color="black",stroke=0.3) +
scale_fill_viridis("Sample Type",option = "viridis", discrete=T, labels=c('Biofilm', 'Fluid', 'Sediment')) +
scale_shape_manual("Island",values=c(21,25,22,24,23)) +
theme_glab() + theme(legend.position = "right",theme(strip.text.x = element_text(size=14),
                aspect.ratio = 0.75)) +
guides(fill = guide_legend(override.aes = list(shape = 21) ),
            shape = guide_legend(override.aes = list(fill = "black")))

nmds_uUnif_p<-plot_ordination(prok_ndata, prok_nmds_ununif, type="samples",title="NMDS - Unweighted UniFrac") +
geom_text(aes(label= code), size=5, hjust=-.1,vjust=2) +
geom_point(aes(fill=type,shape=island),size=8,color="black",stroke=0.3) +
scale_fill_viridis("Sample Type",option = "viridis", discrete=T, labels=c('Biofilm', 'Fluid', 'Sediment')) +
scale_shape_manual("Island",values=c(21,25,22,24,23)) +
theme_glab() + theme(legend.position = "right",theme(strip.text.x = element_text(size=14),
                aspect.ratio = 0.75)) +
guides(fill = guide_legend(override.aes = list(shape = 21) ),
            shape = guide_legend(override.aes = list(fill = "black")))


ggarrange(nmds_wUnif_p,NA, nmds_uUnif_p,NA,
          ncol=2, nrow=2,
          align = c("hv"),
          labels= c('A)',NA, 'B)'),
          font.label = list(size = 20, color = "black", face = "plain", family = 'Helvetica'),
          legend = "none",
  common.legend = TRUE       
         )
#ggsave("../plots/svg/nmds_w_uw_jacard.svg", width=16, height=16)





save.image()
