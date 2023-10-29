# All of the functions will be stored here
    # Import all of the tables and variables we need
from Data_Table import *

# Import all of the packages
import pandas as pd # data processing, CSV file I/O (e.g. pd.read_csv)
from collections import Counter



### Nucleotide Composition 
# Get the composition of each nucleotide in a sequence
def nucleotides_composition(seq):
    nucleotides = {'A': 0,
                    'G': 0,
                    'C': 0,
                    'T': 0}
    for bases in seq: 
        nucleotides[bases] += 1
    return nucleotides

# Turn the nucleotide composition into a dataframe to compare the two sequences
    # They are currently a dictionary
def nuc_df(seq):
    seq = pd.DataFrame.from_dict(seq, orient = 'index')
    seq = seq.reset_index()
    seq = seq.rename(columns={"index": "Nucleotide", 0: "Composition"})
    return seq



### Template Strand - Reverse Complement of DNA strand 
def reverseComp(seq):
    seq = seq[::-1]
    return''.join([Comp[nuc] for nuc in seq])



### GC Content
# Calculate GC content
def gc(seq):
    gc_count = 0
    seq = str(seq)
    for bases in seq: 
        gc_count = ((seq.count('G') + seq.count('C'))/ len(seq) * 100)
        return gc_count



### Kmers
# Get all of the kmers 
def kmers_list(seq, ksize):
    seq = str(seq)
    kmers = []
    n_kmers = len(seq) - ksize + 1

    for bases in range(n_kmers):
        kmer = seq[bases:bases + ksize]
        kmers.append(kmer)
    return kmers

# Change data type into dictionary to be used for data visualization
def count_dict(kmer):   
    kmers_count = Counter(kmer)
    kmers_count = dict(kmers_count)
    return kmers_count

# Change kmer data type into a dataframe
def kmer_df(seq):
    seq = pd.DataFrame.from_dict(seq, orient = 'index')
    seq = seq.reset_index()
    seq = seq.rename(columns={"index": "Kmer", 0: "# of Occurrence"})
    #seq.sort_values(by=seq[seq.columns[1]], ascending=False),
    seq.style.bar(subset=seq[seq.columns[1]], color='#').background_gradient(cmap='Reds')
    return seq



### Transcription - DNA --> RNA
def transcription(seq):
    return seq.replace('T','U')



### Translation - mRNA --> protein
def translation1(seq, init_pos=0):
    return ''.join([RNA_PROTEIN_code[seq[pos:pos+3]] for pos in range(init_pos, len(seq) - 2, 3)])


# Turn into a dataframe
def codon_df(seq):
    seq = pd.DataFrame.from_dict(seq, orient = 'index')
    seq = seq.reset_index()
    seq = seq.rename(columns={"index": "Amino Acids", 0: "Composition"})
    return seq


### Protein - Amino Acid Chains
def aa_chain_length(seq): 
    aa_length = []
    for i in seq:
        aa_length.append(len(i))
    return(aa_length)

# Turn the two lists (aa_chain and length of aa_chain into a dictionary)
def aa_chain_dict(dict1, dict2):
    aa_dict = {}
    for i in range(len(dict1)):
        aa_dict[dict1[i]] = dict2[i]
    return (aa_dict)
    
# Turn into a dataframe
def aa_df(seq):
    seq = pd.DataFrame.from_dict(seq, orient = 'index')
    seq = seq.reset_index()
    seq = seq.rename(columns={"index": "Amino Acid Chain", 0: "Length"})
    return seq



