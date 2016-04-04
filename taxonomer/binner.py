import os
path = os.path.dirname(os.path.abspath(__file__))
import sys
sys.path.append("%s/scripts/cython/"%(path))
import binner_PY
import argparse
import time

parser = argparse.ArgumentParser(description="Classify fasta / fastq file of sequences")
parser.add_argument("db_prefix", type=str, help="file name for database prefix")
parser.add_argument("reads", type=str, help="fasta/fastq file of sequences to classify")
parser.add_argument("output", type=str, help="output file name", default = None)
parser.add_argument("-np", type=int, help="number of processors", default = 1)
parser.add_argument("-load_db", help="set to turn off loading db into memory before queries", action="store_false", default=True)
args = parser.parse_args()

bp = binner_PY.binner_PY()
start = time.time()
bp.bin_reads(args.reads,args.np,int(args.load_db),args.output,args.db_prefix+".bmi", args.db_prefix+".btbi")
print "time to load db and classify reads: %f"%(time.time() - start)
