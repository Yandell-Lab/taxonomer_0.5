import classify_READpy
import argparse
from Bio import SeqIO
import math

parser = argparse.ArgumentParser(description="Calculate effective lengths of sequences")
parser.add_argument("fasta_file", type=str, help="fasta file of reference sequences")
parser.add_argument("db_prefix", type=str, help="taxonomer database file prefix")
parser.add_argument("tri_file", type=str, help=".tri file that gives taxid relationships")
parser.add_argument("sti_file", type=str, help=".sti file for sequence id / taxid relationships")
parser.add_argument("output_file", type=str, help="output file -- tab delimited with seqid and effective length")
args = parser.parse_args()

cr = classify_READpy.classify_READpy()
cr.load_classification_dbs(args.db_prefix+".mi", args.db_prefix+".tbi", args.db_prefix+".rsi",args.tri_file,args.sti_file,0,0,0)
kl = cr.get_kmer_length()
total_kc = cr.get_total_kmer_count()

seqid_taxid = {}
cur_seqid = None
f = open(args.sti_file,'r')
for line in f:
    if line[0] == ">":
        cur_seqid = line[1:].strip()
    else:
        seqid_taxid[cur_seqid] = line.strip()
f.close()

def calculate_effective_length(seq):
    hobs = 0
    hexp = 0
    expected_freq = 1.0 / total_kc

    for kmer in [seq[i:i+kl] for i in range(len(seq) + 1 - kl)]:
        count, num_taxids = cr.get_kmer_info(kmer)
        if count > 0:
            rel_freq = float(count) / float(total_kc)
            hobs += rel_freq * math.log(rel_freq)
            hexp += expected_freq * math.log(expected_freq)

    el = 0.0
    try:
        el = (hexp / hobs)*float(len(seq))
    except:
        pass

    return el


out = open(args.output_file,'w')
handle = open(args.fasta_file,'r')
for record in SeqIO.parse(handle, "fasta"):
    el = calculate_effective_length(str(record.seq))
    out.write("%s\t%f\n"%(seqid_taxid[record.id],el))
handle.close()
out.close()
