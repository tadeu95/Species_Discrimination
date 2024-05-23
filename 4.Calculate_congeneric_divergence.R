library("Biostrings")
library(stringr)
library(tidyr)
library(dplyr)
library(seqRFLP)
library(msa)
library(ape)

#Set working directory
setwd("C:/User/folder_with_MSA")

#Load MSA (load the version before removal of duplicate haplotypes)
fastaFile = readDNAStringSet("MSA.fasta")
seq_name = names(fastaFile)
sequence = paste(fastaFile)
df <- data.frame(seq_name, sequence)
df$number<-str_count(df$sequence, pattern="[A-Z]")
mean(df$number)
df2<-separate(df,seq_name,into = c("ID","species"),sep = "\\|",remove = FALSE,extra = "merge")
df2$genus<- str_extract(df2$species, '[A-Za-z]+')

#Keep only genera with at least two different species
df3<-df2[c("species","genus")]
df4<-unique(df3)
speciespergenus<-data.frame(table(df4$genus))
at_least_2_species<-speciespergenus[speciespergenus$Freq>1,]
df5<-df2[df2$genus %in% at_least_2_species$Var1,]
length(unique(df5$species))
length(unique(df5$genus))
nrow(df5)

#Remove duplicate haplotypes per species
#This step is optional because if you're doing the analysis for records that include all your target genomic regions, you should keep duplicate haplotypes per species
df5 <- df5 %>%
  group_by(species, sequence) %>%
  slice(1) %>%
  ungroup()
nrow(df5)
df5$sequence<-gsub("-","",df5$sequence)

#Create new fasta where every genus has at least two species and, if chosen, the duplicate haplotypes per species were removed
dataframe2fas(data.frame(df5)[,c("seq_name","sequence")],"Atleast_two_species_per_genus.fasta")

#Read back the fasta file after creating it
myseq <- readDNAStringSet("Atleast_two_species_per_genus.fasta", format = "fasta")
species_names <- sapply(strsplit(names(myseq), "\\|"), function(x) x[2])
genus_names<- str_extract(species_names, '[A-Za-z]+')
genus_names<-unique(genus_names)

#Calculate congeneric divergence
result <- data.frame(genus = character(),
                     avg_distance = numeric(),
                     row_count = integer(),
                     stringsAsFactors = FALSE)

#Optionally save the pairwise distances if you're doing the analysis for records that include all your target genomic regions
#pairwise_distances <- data.frame(species1 = character(),
#                                species2 = character(),
#                               distance = numeric(),
#                              stringsAsFactors = FALSE)

for (genus in genus_names) {
  print(paste0("Genus:",genus))
  subset_indices <- grep(paste0(genus, " "), names(myseq))
  subset_myseq <- myseq[subset_indices, ]
  aligned_seqs <- msa(subset_myseq,method = "ClustalW",type="DNA")
  aligned_stringset <- as(aligned_seqs, "DNAStringSet")
  aligned_matrix <- as.matrix(aligned_stringset)
  mydnabin <- as.DNAbin(aligned_matrix)
  mydist <- dist.dna(mydnabin,model="raw",pairwise.deletion = TRUE)
  df_dist <- setNames(reshape2::melt(as.matrix(mydist)), c('rows', 'vars', 'values'))
  names(df_dist) <- c("seqname1", "seqname2", "dist")
  print(paste0("Number of seqs:",length(unique(df_dist$seqname1))))
  df_dist$seqname1 <- as.character(df_dist$seqname1)
  df_dist$seqname2 <- as.character(df_dist$seqname2)
  df_dist<-separate(df_dist,seqname1,into = c("ID1","species1"),sep = "\\|",remove = FALSE,extra = "merge")
  df_dist<-separate(df_dist,seqname2,into = c("ID2","species2"),sep = "\\|",remove = FALSE,extra = "merge")
  congeneric <- df_dist[df_dist$species1 != df_dist$species2, ]
  congeneric <- congeneric[congeneric$species1 < congeneric$species2, ]
  genus_avg_distance <- mean(congeneric$dist) * 100
  print(paste0("Average distance:",genus_avg_distance))
  genus_row_count <- nrow(congeneric)
  result_row <- data.frame(genus = genus,
                           avg_distance = genus_avg_distance,
                           row_count = genus_row_count)
  #pairwise_distances <- rbind(pairwise_distances, congeneric)
  result <- rbind(result, result_row)
}

#Check the number of pairwise comparisons, mean, min and max average congeneric distances, and standard error of the mean
sum(result$row_count)
mean(result$avg_distance)
min(result$avg_distance)
max(result$avg_distance)
sd <- sd(result$avg_distance)
n <- length(result$avg_distance)
se <- sd/sqrt(n)
se
nrow(result)
length(unique(species_names))
length(myseq)

#Copy the results to clipboard to paste in excel
write.table(result, "clipboard", sep="\t", row.names=FALSE)

#If you used the option of saving the pairwise distances
#write.table(pairwise_distances, "clipboard", sep="\t", row.names=FALSE)
#write_tsv(pairwise_distances,"Scombri_intersection_fullCOI_pairwise.txt")
