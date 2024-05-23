library("Biostrings")
library(tidyr)
library(stringr)
library(worms)
library(dplyr)
library(readr)
library(seqRFLP)

#NOTE: Depending on the target taxonomic groups, several steps of this script should be skipped or changed

#Set working directory with GenBank records
setwd("C:/User/folder_with_GenBank_records")

#Read fasta files with GenBank records. Combine two fastas if two data mining approaches were used
fastaFile = readDNAStringSet("GenBank_records_datamining1.fasta")
seq_name = names(fastaFile)
sequence = paste(fastaFile)
df <- data.frame(seq_name, sequence)
fastaFile = readDNAStringSet("GenBank_records_datamining2.fasta")
seq_name = names(fastaFile)
sequence = paste(fastaFile)
df1.5 <- data.frame(seq_name, sequence)
df1.6<-rbind(df,df1.5)

#Alternatively, read unique file with GenBank records if only one data mining approach was used:
#fastaFile = readDNAStringSet("GenBank_records.fasta")
#seq_name = names(fastaFile)
#sequence = paste(fastaFile)
#df1.7 <- data.frame(seq_name, sequence)

#Data preprocessing
df2<-separate(df1.6,seq_name,into = c("ID","species"),sep = "\\s",remove = FALSE,extra = "merge")
df3<-df2 %>% 
  extract(species, into = c("species_name", "other_info"), "^(\\S+\\s+\\S+)\\s+(.*)")
df3<-df3[,-1]
df3$number<-str_count(df3$sequence, pattern="[A-Z]")
plot(df3$number)
names(df3)[2]<-"species"
names(df3)[4]<-"sequence"
df4<-df3 %>% 
  group_by(ID) %>% 
  slice_max(number)
df5<-ddply(df4,.(ID),function(x) x[sample(nrow(x),1),])
#View(data.frame(table(df5$ID)))

#Remove COI-like sequences
coi_info<-df5[,c("ID","species","other_info","sequence","number")]
df5_info<-coi_info[coi_info$ID %in% df5$ID ,]
#Confirm if the records truly comprise "COI-like" sequences by manually inspecting
coi_like<-df5_info[grep("like",df5_info$other_info),]
df5.5<-df5[!(df5$ID %in% coi_like$ID),]
#Otherwise use this approach
#df5.5<-df5[!grepl("cytochrome oxidase subunit 1-like",df5$other_info),]
#df5.5<-df5.5[!grepl("CO1-like",df5.5$other_info),]
#df5.5<-df5.5[!grepl("COI-like",df5.5$other_info),]
#df5.5<-df5.5[!grepl("COI like",df5.5$other_info),]

#Remove "predicted" and "unverified" records, ambiguous characters in the species names and records not identified to the species level
df6<-df5.5[!grepl("PREDICTED", df5.5$species),]
df6$species<-gsub("\\(","",df6$species)
df6$species<-gsub("\\'","",df6$species)
ambiguous<-df6[grepl("\\.", df6$species),]
ambiguous2<-ambiguous %>% 
  mutate(counter = str_count(ambiguous$species, pattern = "\\w+"))
df7<-df6[!(df6$ID %in% ambiguous2$ID),]
df7<-df7[!grepl("UNVERIFIED", df7$species),]

#Keep only marine species and standardize taxonomy, "new_name" is the correct species name according to WoRMS
worms<-wormsbynames(as.character(unique(df7$species)), marine_only=FALSE, ids=TRUE)
worms2<-worms%>%
  mutate(new_name = ifelse(status!="accepted" & !is.na(status),valid_name,
                           ifelse(status=="accepted" & !is.na(status),name,name)))
colnames(worms2)[colnames(worms2) == "name"] <- "species"
worms2<-worms2[worms2$isMarine==1,]
#If you want fresh water and brackish species, for example, use instead:
#worms2<-worms2[worms2$isFreshWater==1|worms2$isBrackish==1,]
worms3<-worms2[,c("new_name","species","class","order","family")]
worms3<-na.omit(worms3)
df8<-df7[df7$species %in% worms3$species,]
df9<-left_join(df8,worms3,by="species")
df9<-df9[,c("ID","species","other_info","sequence","number","new_name","class","order","family")]
df9$genus<- str_extract(df9$new_name, '[A-Za-z]+')

#Save data in tabular format
write_tsv(df9,"DF.tsv")
df9<-read.delim("DF.tsv")

#Create fasta for all sequences
df9$name<-paste(df9$ID, df9$new_name, sep=" | ")
df9<-data.frame(df9$name,df9$sequence)
names(df9)<-c("name","sequence")
dataframe2fas(df9,"sequences.fasta")

#Subset into target taxonomic groups
orders<-data.frame(table(df9$order))
df10<-df9[df9$order=="Pleuronectiformes" | df9$order=="Carangiformes" | df9$order=="Carangaria incertae sedis",]

#Create subset fasta from which the target amplicons/regions will be extracted in the MAFFT alignment
df10$name<-paste(df10$ID,df10$new_name,df10$family,sep="|")
dataframe2fas(df10[,c("name","sequence")],file="full-length_sequences_subset.fasta")

#Search the mitogenome data set to select queries for local alignment
mitogenomes<-(df10[df10$number>15000,c("ID","species","other_info")])

#Check the most well-represented families in the data set to choose mitochondrial genome query sequences for local alignment
families<-df10[,c("species","family")]
families2<-unique(families)
View(table(families2$family))
df11<-df10[(df10$family=="Carangidae" & df10$number>15000) | (df10$family=="Pleuronectidae" & df10$number>15000)|
              (df10$family=="Paralichthyidae" & df10$number>15000)|(df10$family=="Soleidae" & df10$number>15000)|
              (df10$family=="Cynoglossidae" & df10$number>15000)|(df10$family=="Bothidae" & df10$number>15000),]

#Select, for example, two mitochondrial genomes from each family, if available
df12<-ddply(df11,.(family),function(x) x[sample(nrow(x),2),])
df13<-df12[,c("ID","new_name","family","sequence")]
df13$name<-paste(df13$ID,df13$new_name,df13$family,sep="|")
#This following fasta contains the sequences to be trimmed by the primer annealing regions in python,
#to obtain the metabarcoding amplicons for local and multiple sequence alignment in MAFFT 
dataframe2fas(df113[,c("name","sequence")],file="sequences_to_cut_in_python.fasta")

###############################################################################
##############################
#Use the python script "Cut_Sequences_By_Primer_Annealing_Region.py" to cut sequences in "sequences_to_cut_in_python.fasta"
#Afterwards, use the query sequences obtained to extract the target regions from "full-length_sequences_subset.fasta"
#Do this resorting to the following MAFFT calculation: https://mafft.cbrc.jp/alignment/server/specificregion-last.html
#NOTE: For each genomic region you analyze, you should add an outgroup sequence to the alignment in order to determine the monophyhletic species
##############################
###############################################################################