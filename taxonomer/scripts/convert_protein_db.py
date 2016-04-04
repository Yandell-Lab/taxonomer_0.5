from Bio import SeqIO
import argparse
import sys

parser = argparse.ArgumentParser(description="Change amino acids to predetermined codons")
parser.add_argument("protein_fasta", type=str, help="fasta protein file")
parser.add_argument("output", type=str, help="output file name")
parser.add_argument("--afterburner", help="flag for afterburner", action="store_true",default=False)
args = parser.parse_args()


bases = ['T', 'C', 'A', 'G']
codons = [a+b+c for a in bases for b in bases for c in bases]
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
codon_table = dict(zip(codons, amino_acids))
seen_aa = {}
aa_to_codon = {}
for cdn,aa in codon_table.items():
    if aa in seen_aa:
        continue
    aa_to_codon[aa] = cdn
    seen_aa[aa] = cdn
print "\n"
print aa_to_codon

if args.afterburner:
    #collapsed_alphabet = {'S':'S','G':'S','V':'V','M':'V','L':'V','I':'V','C':'V','F':'F','Y':'F','N':'N','E':'N','Z':'N','B':'N','D':'N','K':'N','Q':'N','R':'N','T':'T','*':'*','A':'A','W':'W','X':'X','P':'P','H':'H'}
    collapsed_alphabet = {'A': 'A', 'C': 'C', 'E': 'E', 'D': 'D', 'G': 'G', 'F': 'F', 'I': 'I', 'H': 'H', 'K': 'K', 'M': 'M', 'L': 'I', 'N': 'N', 'Q': 'E', 'P': 'P', 'S': 'S', 'R': 'K', 'T': 'S', 'W': 'W', 'V': 'I', 'Y': 'F', '*':'*'} #kmeans 14
    for aa in aa_to_codon:
        aa_to_codon[aa] = aa_to_codon[collapsed_alphabet[aa]]
    print "\nafterburner codons"
    print aa_to_codon
    print "\n"

out = open(args.output,"w")

for seq_record in SeqIO.parse(args.protein_fasta, "fasta"):
    seq = ""
    for aa in seq_record.seq:
        if aa in aa_to_codon:
            seq += aa_to_codon[aa]
        else:
            seq += "NNN"
            #seq += aa_to_codon["*"]
    out.write(">%s\n%s\n"%(seq_record.id,seq))

out.close()
