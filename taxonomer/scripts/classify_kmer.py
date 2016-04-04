import os
path = os.path.dirname(os.path.abspath(__file__))
import sys
sys.path.append("%s/cython/"%(path))
import classify_READpy
import argparse
import time

parser = argparse.ArgumentParser(description="Build Taxonomer database for sequence classification")
parser.add_argument("db_prefix", type=str, help="file name for database prefix")
parser.add_argument("tri_file", type=str, help=".tri file that gives taxid relationships")
parser.add_argument("sti_file", type=str, help=".sti file that gives taxid / seqid relationships")
parser.add_argument("reads", type=str, help="comma separated lsit of read sequences to be classified")
parser.add_argument("-load_db", help="set to false to turn off loading db into memory before queries", action="store_false", default=True)
parser.add_argument("--ties", type=int, help="1 = output ties, 0 = do not output ties (default=0)",default=0)
parser.add_argument("--protein","-prot", type=int, help="1 = classify in protein space, DB *MUST* be built for this", default=0)
parser.add_argument("--afterburner", "-aft", type=int, help = "1=use afterburner search, DB *MUST* be built for this", default=0)
args = parser.parse_args()

if args.afterburner == 1:
    args.protein = 1

cr = classify_READpy.classify_READpy()
cr.load_classification_dbs(args.db_prefix+".mi", args.db_prefix+".tbi", args.db_prefix+".rsi",args.tri_file,args.sti_file,int(args.load_db),args.ties,args.protein,args.afterburner)
kl = cr.get_kmer_length()

seqs = args.reads.strip().split(",")
for seq in seqs:
    j = 0
    for kmer in [seq[i:i+kl] for i in xrange(len(seq) - (kl-1))]:
        start = time.time()
        count, num_taxids, rev_count, rev_ntaxids = cr.get_kmer_info(kmer)
        if count > 0 or num_taxids > 0 or rev_count > 0 or rev_ntaxids > 0:
            print "%s\t%d\t%d\t%d\t%d\t%f\t%d"%(kmer, count, num_taxids, rev_count, rev_ntaxids,time.time() - start,j)
        #if j % 3 == 0:
        #    print kmer,j
        j += 1
