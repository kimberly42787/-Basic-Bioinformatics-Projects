# Main Testing File

# Import everything from the Toolkit_Biopython File(Function File)
from Toolkit_Biopython import *

# Import gene sequence from a Fasta File
    # Using Biopython package to parse through the Fasta File and extract the ID, length, and sequence
from Bio import SeqIO
from collections import Counter
import Bio
from Bio.Seq import Seq



# # Test out program with the TP53 sequence using a FASTA file

gene_seq = []

for seq_record in SeqIO.parse("/Users/kim/Desktop/repos/-Basic-Bioinformatics-Projects/DNA to Protein Data Visualization/TP53.fna", "fasta"):
    gene_seq.append(seq_record)
    print("Id: " + seq_record.id  + " \t " + "Length: " + str("{:,d}".format(len(seq_record))) )
    print(repr(seq_record.seq) + "\n")

type(gene_seq)



# # Assign the sequence to make it easier to call later. Since there's two sequences, we asssign each to its own variable

seq1 = gene_seq[0].seq # First ID
seq2 = gene_seq[1].seq # Second ID


# # Reverse Complement of DNA strand:
# print(reverseComp(seq1))
# print(reverseComp(seq1))


# # Print out GC content:

# print("GC content of Seq1: " + str(gc(seq1)) + "\n")
# print("GC content of Seq1: " + str(gc(seq2)) + "\n")



# # Assign nuc composition to a variable. Print the nucleotide composition of both sequences

seq1_nuc = (nucleotides_composition(seq1))
seq2_nuc = (nucleotides_composition(seq2))

print('\n' + 'Seq 1 Composition: ' + str(nucleotides_composition(seq1)))
print('\n' + 'Seq 2 Composition: ' + str(nucleotides_composition(seq2)))


# # Turn each nuc composition dictionary to a dataframe. This will be used later for data visualization. Print out to check if it's correct
seq1_df = nuc_df(seq1_nuc)
seq2_df = nuc_df(seq2_nuc)

# print(seq1_df)
# print(seq2_df)


# # Create a column that identifies which sequence the composition is from. Will be useful for our graph 
seq1_df['Sequence'] = 'Sequence 1' 
seq2_df['Sequence'] = 'Sequence 2'
merge_df = pd.concat([seq1_df,seq2_df]) # Merge two df together
# print(merge_df)


# # Print out kmers (trimers - 3)
kmer_seq1 = kmers_list(seq1, 3)
kmer_seq2 = kmers_list(seq2, 3)

# print(kmer_seq1)
# print(kmer_seq2)



# # Get the # of occurrence for each kmers 
seq1_kmer_count = count_dict(kmer_seq1)
seq1_kmer_count = count_dict(kmer_seq2)

# # Change it into a dataframe
seq1_kmer_df = kmer_df(seq1_kmer_count)
seq2_kmer_df = kmer_df(seq1_kmer_count)

seq1_kmer_df['Sequence'] = 'Sequence 1'
seq2_kmer_df['Sequence'] = 'Sequence 2'
merge_kmer = pd.concat([seq1_kmer_df,seq2_kmer_df]) # Merge the two df together so we can compare the two sequences
# print(merge_kmer)



# # Transcribe the DNA seq into mRNA

TP53_mRNA1 = (transcription(seq1))
TP53_mRNA2 = (transcription(seq2))

# print(TP53_mRNA1)

# # Translate mRNA into protein
protein1 = (translation1(TP53_mRNA1))
protein2 = (translation1(TP53_mRNA2))
# print(protein1)

# Get the number of occurence for each codons 
protein_count1 = Counter(protein1)
protein_count1_df = codon_df(protein_count1)

protein_count2 = Counter(protein2)
protein_count2_df = codon_df(protein_count2)


aa_chains1 = protein1.split('*') # List
aa_chains2 = protein1.split('*')

# print((aa_chains1))
# print((aa_chains2))


aa_length1 = aa_chain_length(aa_chains1) #List
aa_length2 = aa_chain_length(aa_chains2) #List

# Create a dictionary using the two lists (aa_chain and aa_length)
aa_dict1 = aa_chain_dict(aa_chains1,aa_length1)
aa_dict2 = aa_chain_dict(aa_chains1,aa_length2)

# # Create a dataframe
aa_df1 = aa_df(aa_dict1)
aa_df2 = aa_df(aa_dict2)

# Assign a new column to assign each dataframe to their specific sequence
aa_df1['Sequence'] = 'Sequence 1' 
aa_df2['Sequence'] = 'Sequence 2'
merge_aa= pd.concat([aa_df1,aa_df2]) # Merge the two df together so we can compare the two sequences
print(merge_aa)













