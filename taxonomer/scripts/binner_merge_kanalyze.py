import os
path = os.path.dirname(os.path.abspath(__file__))
import sys
sys.path.append("%s/cython/"%(path))
import argparse
import binner_mergepy

parser = argparse.ArgumentParser(description="Merge kmer count files for subsequent binning")
parser.add_argument("kmer_counts_1", type=str, help="name of kmer count file in kanalyze hex output format")
parser.add_argument("kmer_counts_2", type=str, help="name of kmer count file to merge, must also be in kanalyze hex output format")
parser.add_argument("-id_1", type=int, help="integer value used to identify first database, must be a multiple of 2", default=0)
parser.add_argument("-id_2", type=int, help="integer value used to identify second databse, must be a multiple of 2", default=0)
parser.add_argument("output", type=str, help="merged output file")
args = parser.parse_args()

binner_mergepy.merge_counts(args.kmer_counts_1,args.id_1,args.kmer_counts_2,args.id_2,args.output)
