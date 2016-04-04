import os
path = os.path.dirname(os.path.abspath(__file__))
import sys
sys.path.append("%s/cython/"%(path))
import classify_READpy
import argparse
import time

parser = argparse.ArgumentParser(description="Get confidence scores for db sequences")
parser.add_argument("db_prefix", type=str, help="file name for database prefix")
parser.add_argument("tri_file", type=str, help=".tri file that gives taxid relationships")
parser.add_argument("sti_file", type=str, help=".sti file for sequence id / taxid relationships")
parser.add_argument("db_file", type=str, help="fasta file of database reference sequences")
parser.add_argument("read_length", type=int, help="approximate read length for confidence calculation")
parser.add_argument("-output", type=str, help="output file name", default = None)
parser.add_argument("-np", type=int, help="number of processors", default = 1)
args = parser.parse_args()

cr = classify_READpy.classify_READpy()
cr.load_classification_dbs(args.db_prefix+".mi", args.db_prefix+".tbi", args.db_prefix+".rsi",args.tri_file)
cr.classify_db(args.db_file, args.sti_file, args.read_length, args.output, args.np)
