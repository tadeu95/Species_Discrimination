library("Biostrings")
library(stringr)
library(dplyr)
library(tidyverse)
library(seqRFLP)
library(plyr)
library(ape)
library(MonoPhy)
library(readr)
library(tidyr)

#Set working directory
setwd("C:/User/folder_with_MSA")

#Load the  multiple sequence alignment (MSA) fasta obtained using MAFFT
fastaFile = readDNAStringSet("MSA_from_MAFFT.fasta")
seq_name = names(fastaFile)
sequence = paste(fastaFile)
df <- data.frame(seq_name, sequence)
df$number<-str_count(df$sequence, pattern="[A-Z]")
plot(df$number)

#Remove duplicate sequences (i.e., the queries)
df2<-df %>% 
  group_by(seq_name) %>% 
  slice_max(number)
df3<-ddply(df2,.(seq_name),function(x) x[sample(nrow(x),1),])
df4<-separate(df3,seq_name,into = c("ID","species"),sep = "\\|",remove = FALSE,extra = "merge")
mean(df4$number)
plot(df4$number)
length(unique(df4$species))

#Create new fasta file and realign the sequences if needed
dataframe2fas(df4[,c("seq_name","sequence")],"MSA.fasta")

fastaFile = readDNAStringSet("MSA.fasta")
seq_name = names(fastaFile)
sequence = paste(fastaFile)
df <- data.frame(seq_name, sequence)
df$number<-str_count(df$sequence, pattern="[A-Z]")
mean(df$number)
df2<-separate(df,seq_name,into = c("ID","species"),sep = "\\|",remove = FALSE,extra = "merge")
length(unique(df2$species))

#If the datasets are too large, remove duplicate haplotypes for each species from the MSA
#Keep  two repeated haplotypes per species if you do not want to keep singletons from being intruders in another species' clade
df_filtered2 <- df2 %>%
  group_by(species, sequence) %>%
  slice_head(n=2) %>%
  ungroup()
length(unique(df_filtered2$species))
for_fasta<-data.frame(df_filtered2[,c("seq_name","sequence")])
dataframe2fas(for_fasta,"MSA_with_2_repeated_haplotypes_per_species.fasta")

###############################################################################
##############################
#Construct a phylogenetic tree for each of your multiple sequence alignments with the inclusion of one outgroup sequence for re-rooting
#You can use MEGA for example, but you can also use a tree constructed within R or any other software, preferably saved in newick format
##############################
###############################################################################

#Create taxonomy data set for monophyly test using the fasta containing the MSA that was used for the construction of the tree
fastaFile = readDNAStringSet("MSA.fasta")
seq_name = names(fastaFile)
sequence = paste(fastaFile)
df <- data.frame(seq_name, sequence)
df$number<-str_count(df$sequence, pattern="[A-Z]")
mean(df$number)
df2<-separate(df,seq_name,into = c("ID","species"),sep = "\\|",remove = FALSE,extra = "merge")
####If the tree was constructed in MEGA (10 or 11), replace spaces with underscores in the taxonomy data set, because MEGA automatically replaces spaces in tip names with underscores
####create the taxonomy data sets for the levels you want (i.e., species, genus, etc...)
df2$seq_name<-gsub(" ","_",df2$seq_name)
taxonomy_species<-df2[,c("seq_name","species")]

#Load the tree
#The tree must include an outgroup, so it is rerooted and the monophyly test works.
Tree <- ape::read.tree("Tree.nwk")
#Reroot the tree using an outgroup sequence
Tree2<-root(Tree,outgroup = "GB000001|Outgroup",resolve.root=T)

#Assess the monophyly at the species level
mono_species<-AssessMonophyly(Tree2,taxonomy = taxonomy_species,verbosity=100)
summary_species<-mono_species$species$summary
result_species<-mono_species$species$result
summary_species<- tibble::rownames_to_column(summary_species, "VALUE")
result_species<- tibble::rownames_to_column(result_species, "VALUE")
names(result_species)[6]<-"Number_Intruders"
names(result_species)[8]<-"Number_Outliers"

#Check intruder singletons
intruder_species <- result_species$Intruders
species_list <- strsplit(intruder_species, ",")
species_vector <- unlist(species_list)
species_vector <- trimws(species_vector)
unique_species <- data.frame(unique(species_vector))
names(unique_species)<-"names"
result_species2<-result_species %>%
  mutate(Monophyly = ifelse((Monophyly=="Monotypic") & (VALUE %in% unique_species$names),"Monotypic_Intruder", Monophyly))
new_row <- data.frame(VALUE = "Monotypic_intruder", Taxa=nrow(result_species2[result_species2$Monophyly=="Monotypic_Intruder",]),Tips=nrow(result_species2[result_species2$Monophyly=="Monotypic_Intruder",]))
summary_species <- rbind(summary_species, new_row)
summary_species[4,2]<-as.numeric(summary_species[4,2])-as.numeric(summary_species[7,2])
summary_species[3,2]<-as.numeric(summary_species[3,2])+as.numeric(summary_species[7,2])
summary_species <- summary_species[-nrow(summary_species), ]

#Save the results as tsv files by copying to clipboard and pasting for example in Excel
write.table(result_species2, file="clipboard-16384", sep="\t", row.names=FALSE)
write.table(summary_species, "clipboard", sep="\t", row.names=FALSE)

#Or save them by creating tsv files in the directory
write_tsv(result_species2,"Monophyly_result.tsv")
write_tsv(summary_species,"Monophyly_summary.tsv")

#Check the reasons behind non-monophyly

#Load the "result" table from the monophyly analysis if needed
result_species2 <- read.delim("Monophyly_result.tsv")

nomonophyly<-result_species2[result_species2$Monophyly=="No" | result_species2$Monophyly=="Monotypic_Intruder",]
nrow(nomonophyly)

#Load the fasta file that was used for tree construction if needed 
fastaFile = readDNAStringSet("MSA.fasta")
seq_name = names(fastaFile)
sequence = paste(fastaFile)
df <- data.frame(seq_name, sequence)
df$number<-str_count(df$sequence, pattern="[A-Z]")
mean(df$number)
df2<-separate(df,seq_name,into = c("ID","species"),sep = "\\|",remove = FALSE,extra = "merge")
df2<-df2[df2$species %in% nomonophyly$VALUE,]

# Find species that share a haplotype with at least one other species
df2$shared_seqs <- NA 
for (i in 1:nrow(df2)) {
  current_seq <- df2$sequence[i]
  current_species <- df2$species[i]
  shared_seqs <- df2$seq_name[df2$sequence %in% current_seq & df2$species != current_species]
  
  if (length(shared_seqs) > 0) {
    df2$shared_seqs[i] <- paste(shared_seqs, collapse = ", ")
  } else {
    df2$shared_seqs[i] <- "None"
  }
}
df2$number_shared <- sapply(df2$shared_seqs, function(x) length(strsplit(x, ",\\s*")[[1]]))
df3<-df2 %>%
  mutate(number_shared = ifelse(shared_seqs=="None",0,number_shared))
species_with_shared<-df3[df3$number_shared>0,]
species_with_shared<-data.frame(unique(species_with_shared$species))  

#Percentage of species with at least one shared haplotype from another species 
nrow(species_with_shared)/nrow(nomonophyly)*100
species_with_shared$genus<- str_extract(species_with_shared$unique.species_with_shared.species., '[A-Za-z]+')
species_per_genus<-data.frame(table(species_with_shared$genus))
names(species_per_genus)[1]<-"genus"
all<-left_join(species_with_shared,species_per_genus,by="genus")

#Calculate percentage of non-monophyletic species with intruders and/or outliers
nomonophyly_real<-result_species2[result_species2$Monophyly=="No",]
nrow(nomonophyly_real)

#percentage with both intruders and outliers
percentage_intruders_outliers<-nrow(nomonophyly_real[nomonophyly_real$Number_Intruders>0 & nomonophyly_real$Number_Outliers>0,])/nrow(nomonophyly_real)*100
percentage_intruders_outliers

#percentage with intruders
percentage_intruders<-nrow(nomonophyly_real[nomonophyly_real$Number_Intruders>0 & nomonophyly_real$Number_Outliers==0,])/nrow(nomonophyly_real)*100
percentage_intruders

#percentage with outliers
percentage_outliers<-nrow(nomonophyly_real[nomonophyly_real$Number_Intruders==0 & nomonophyly_real$Number_Outliers>0,])/nrow(nomonophyly_real)*100
percentage_outliers


