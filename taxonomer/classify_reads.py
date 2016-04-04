import os
path = os.path.dirname(os.path.abspath(__file__))
import sys
sys.path.append("%s/scripts/cython/"%(path))
import argparse
import time
import classify_READpy

parser = argparse.ArgumentParser(description="Classify fasta / fastq file of sequences")
parser.add_argument("db_prefix", type=str, help="file name for database prefix")
parser.add_argument("tri_file", type=str, help=".tri file that gives taxid relationships")
parser.add_argument("sti_file", type=str, help=".sti file that gives taxid / seqid relationships")
#parser.add_argument("taxid_confidence", type=str, help="taxid confidence scores")
parser.add_argument("reads", type=str, help="fasta/fastq file of sequences to classify")
parser.add_argument("output", type=str, help="output file name", default = None)
parser.add_argument("-np", type=int, help="number of processors", default = 1)
parser.add_argument("-load_db", help="set to false to turn off loading db into memory before queries", action="store_false", default=True)
parser.add_argument("--ties", type=int, help="1 = output ties, 0 = do not output ties (default=0)",default=0)
parser.add_argument("--protein","-prot", type=int, help="1 = classify in protein space, DB *MUST* be build for this", default=0)
parser.add_argument("--afterburner", "-aft", type=int, help = "1=use afterburner search, DB *MUST* be built for this", default=0)
parser.add_argument("--kmer_cutoff", type=int, help="set kmer cutoff to change default database build cutoff", default=None)
args = parser.parse_args()

if args.afterburner == 1:
    args.protein = 1

cr = classify_READpy.classify_READpy()
cr.load_classification_dbs(args.db_prefix+".mi", args.db_prefix+".tbi", args.db_prefix+".rsi",args.tri_file,args.sti_file,int(args.load_db),args.ties,args.protein,args.afterburner)
if args.kmer_cutoff != None:
    cr.set_kmer_cutoff(args.kmer_cutoff)

start = time.time()
cr.classify_reads(args.reads, args.np, args.output)
print "time to classify reads: %f"%(time.time() - start)
