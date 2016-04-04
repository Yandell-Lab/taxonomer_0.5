import sys
import argparse
from Bio import SeqIO
import collections
import math

parser = argparse.ArgumentParser(description="Get transcript covariates")
parser.add_argument("reference_fasta", type=str, help="reference fasta file of transcripts")
parser.add_argument("output_file", type=str, help="output file -- tab delimited with transcript covariates")
args = parser.parse_args()

bases = ['T', 'C', 'A', 'G']
dinucleotides = [a+b for a in bases for b in bases]

def gc_content(seq):
    counts = collections.Counter(list(seq))
    gc_count = 0
    if 'G' in counts:
        gc_count += counts['G']
    if 'C' in counts:
        gc_count += counts['C']
    return float(gc_count)/float(len(seq))

out = open(args.output_file,'w')
handle = open(args.reference_fasta,'r')
for record in SeqIO.parse(handle, "fasta"):
    s = str(record.seq)
    length = len(s)
    gcc = gc_content(s)
    out.write("%s\t%f\t%f\t"%(record.id,length,gcc))
    dicounts = collections.Counter([s[i:i+2] for i in xrange(len(s)-1)])
    freqs = []
    for di in dinucleotides:
        if di in dicounts:
            freqs.append( str( float(dicounts[di])/ float(len(s)-2) ) )
        else:
            freqs.append("0.0")
    out.write("%s\n"%("\t".join(freqs)))
out.close()
handle.close()
