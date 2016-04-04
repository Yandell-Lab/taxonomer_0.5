import os
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description="Count kmer by reference")
parser.add_argument("reference", type=str, help="reference fasta sequences")
parser.add_argument("kl",type=int, help="kmer length")
parser.add_argument("output", type=str, help="output file of kmer counts")
parser.add_argument("-canonical",help="flag for canonical kmers",action="store_true",default=False)
args = parser.parse_args()

ref_count = 0
handle = open(args.reference,'r')
for record in SeqIO.parse(handle, "fasta"):
    ref_count += 1
    print ref_count
    ref_out = open("tmp_ref.fa",'w')
    ref_out.write(">%s\n%s\n"%(str(record.id),str(record.seq)))
    ref_out.close()
    if args.canonical == False:
        os.system("java -jar ~/kanalyze-0.9.7/kanalyze.jar count -k %d -m fasta -o tmp_counts.fa -f fasta tmp_ref.fa"%(args.kl))
    else:
        os.system("java -jar ~/kanalyze-0.9.7/kanalyze.jar count -k %d -rcanonical -m fasta -o tmp_counts.fa -f fasta tmp_ref.fa"%(args.kl))
    os.system("cat tmp_counts.fa >> kmers.fa")
handle.close()

if args.canonical == False:
    os.system("java -jar ~/kanalyze-0.9.7/kanalyze.jar count -k %d -m hex -o %s -f fasta kmers.fa"%(args.kl,args.output))
else:
    os.system("java -jar ~/kanalyze-0.9.7/kanalyze.jar count -k %d -rcanonical -m hex -o %s -f fasta kmers.fa"%(args.kl,args.output))
