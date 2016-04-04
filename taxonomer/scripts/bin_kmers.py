import binner_PY
import argparse
import time
import sys

parser = argparse.ArgumentParser(description="Classify fasta / fastq file of sequences")
parser.add_argument("db_prefix", type=str, help="file name for database prefix")
parser.add_argument("read", type=str, help="read sequence")
parser.add_argument("output", type=str, help="output file name", default = None)
parser.add_argument("-load_db", help="set to turn off loading db into memory before queries", action="store_false", default=True)
args = parser.parse_args()

bp = binner_PY.binner_PY()
bp.bin_kmers(args.read, len(args.read), int(args.load_db), args.db_prefix+".bmi", args.db_prefix+".btbi")
