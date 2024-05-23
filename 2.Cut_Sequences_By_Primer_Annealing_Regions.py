# -*- coding: utf-8 -*-
from Bio import SeqIO
from Bio.Data.IUPACData import ambiguous_dna_values
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from difflib import SequenceMatcher
from itertools import product

#Load "sequences_to_cut_in_python.fasta" generated in the R script "Preprocessing_GenBank_data.R"
#This script will trimm the sequences in order to obtain the query regions to be extracted from "full-length_sequences_subset.fasta"
#NOTE: For each genomic region you analyze, you should add an outgroup sequence to the alignment in order to determine the monophyhletic species
records = list(SeqIO.parse("sequences_to_cut_in_python.fasta", "fasta"))

list_seqs = []
list_names = []

for record in records:
    list_seqs.append(str(record.seq))
    list_names.append(str(record.description))
names_seqs = dict(zip(list_names, list_seqs))
total_seqs = len(records)

#The primer pairs to be searched for. Please type them in 5' to 3'
forward_primers = {"MiFish-U_F":"GTCGGTAAAACTCGTGCCAGC","Teleo_F": "ACACCGCCCGTCACTCT","Berry-fish_F":"GACCCTATGGAGCTTTAGAC","Ac16S_F":"CCTTTTGCATCATGATTTAGC"}
reverse_primers = {"MiFish-U_R":"CAAACTGGGATTAGATACCCCACTATG","Teleo_R":"CATGGTAAGTGTACCGGAAG","Berry-fish_R":"AGTTACYHTAGGGATAACAGCG","Ac16S_R":"GCCTAAAAGCAGCCACCTG"}

#You can adjust the degree of similarity between the primers and the annealing regions (from 0.3 to 1) in the "Cut_Regions" function
#NOTE: adjustments should be made taking into account the target regions you are working with
res_list = []

def similar(a, b):
    return SequenceMatcher(None, a, b).ratio()

def find_similar(needle, haystack, backwards = False, min_diff_ratio = 0.30):
    n_len = len(needle)
    if backwards:
        r = range(len(haystack)-n_len, -1, -1)
    else:
        r = range(len(haystack) - n_len)
    for i in r:
        substr = haystack[i:i + n_len]
        if similar(needle, substr) >= min_diff_ratio:
            if not backwards:
                result = haystack[i + n_len:]
            else:
                result = haystack[:i]
            return result
    return []

def extend_ambiguous_dna(seq):
   d = ambiguous_dna_values  #From IUPACData
   return  list(map("".join, product(*map(d.get, seq))))

results_fasta = []

def Cut_Regions(list_seqs, forward_primers, reverse_primers):
    forwardP = []
    reverseP = []
    seq_length = []
    seq_species = []
    for counter, s in enumerate(list_seqs):
        upatedString = s
        for x, (index, i) in enumerate(forward_primers.items()):
            i = extend_ambiguous_dna(i)
            forward_found = False
            for f in i:
                if forward_found:
                    break
                temp = find_similar(f, s, False, 0.95) #Adjust the degree of similarity here
                reverse_found = False
                if len(temp) != 0:
                    upatedString = temp
                    forward_str = temp
                    temp = []
                    print(index, "FORWARD", list_names[counter])
                    forward_found = True
                    for y, (indexr, t) in enumerate(reverse_primers.items()):
                        t = extend_ambiguous_dna(t)
                        if reverse_found:
                            break
                        if x == y:
                            for r in t:
                                temp = find_similar(r, forward_str, True,0.95) #Adjust the degree of similarity here
                                if len(temp) != 0:
                                    upatedString = temp
                                    temp = []
                                    print("%s\n%s\n%i" % (list_names[counter], upatedString, len(upatedString)))
                                    #Collect data
                                    forwardP.append(index)
                                    reverseP.append(indexr)
                                    seq_length.append(len(upatedString))
                                    seq_species.append(" ".join(list_names[counter].split(" ")[1:3]))
                                    ##Save fasta
                                    rec = SeqRecord(Seq(upatedString), id=list_names[counter], description=" ".join(list_names[counter].split(" ")[1:3]) + " " + index + " " + indexr)
                                    res_list.append(rec)
                                    reverse_found = True
                                    break
                            else:
                                print("Reverse primer not found")  
                        else:
                            continue
                else:
                    print("Forward primer not found") 
    return ##df

Cut_Regions(list_seqs, forward_primers, reverse_primers)

print("Total sequences in record:", len(records))

#Save resulting fasta file with trimmed sequences and extract these regions in MAFFT: https://mafft.cbrc.jp/alignment/server/specificregion-last.html
SeqIO.write(res_list, "Output.fasta", "fasta")


