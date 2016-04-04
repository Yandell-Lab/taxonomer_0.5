import argparse
from Bio import SeqIO
import math

parser = argparse.ArgumentParser(description="Calculate effective lengths of sequences")
parser.add_argument("fasta_file", type=str, help="fasta file of reference sequences")
parser.add_argument("cleaned_output", type=str, help="cleaned fasta output")
args = parser.parse_args()

def clean_letter(l):
    if l == "A" or l == "a":
        return l
    elif l == "T" or l == "t":
        return l
    elif l == "G" or l == "g":
        return l
    elif l == "C" or l == "c":
        return l
    else:
        return "N"


out = open(args.cleaned_output,'w')
handle = open(args.fasta_file,'r')
for record in SeqIO.parse(handle, "fasta"):
    out.write(">%s\n"%(record.id))
    out.write("%s\n"%( "".join([clean_letter(l) for l in str(record.seq)]) ))
handle.close()
out.close()
