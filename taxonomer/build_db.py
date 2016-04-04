import os
path = os.path.dirname(os.path.abspath(__file__))
import sys
sys.path.append("%s/scripts/cython/"%(path))
import build_DBpy
import argparse

parser = argparse.ArgumentParser(description="Build Taxonomer database for sequence classification")
parser.add_argument("db_prefix", type=str, help="file name for database prefix")
parser.add_argument("kanalyze_input", type=str, help="kanalyze file of fasta reference")
parser.add_argument("sti_file", type=str, help=".sti file for sequence id / taxid relationships")
parser.add_argument("fasta_file", type=str, help="fasta reference file to be built")
parser.add_argument("-kl",type=int,help="kmer length (default=31)",default=31)
parser.add_argument("-binner",help="flag to build binner database",action="store_true",default=False)
parser.add_argument("-tc","--tax_cutoff",type=int,help="taxid cutoff (default=500)",default=500)
parser.add_argument("-kc","--kmer_cutoff",type=int,help="kmer cutoff (defulat=10000)",default=10000)
parser.add_argument("--protein","-prot", type=int, help="1 = build db for protein searches", default=0)
parser.add_argument("--afterburner", "-aft", type=int, help = "1=build db for afterburner searches", default=0)
args = parser.parse_args()

if args.afterburner == 1:
    args.protein = 1

if args.protein == 1:
    if args.kl % 3 != 0 or args.kl < 12 or args.kl > 30:
        print "when building a protein database, the kmer length must be a multiple of 3 that is >= 12 and <= 30"
        sys.exit(1)
else:
    if args.kl < 8 or args.kl > 31:
        print "when build a database, the kmer length must be >= 8 and <= 31"
        sys.exit(1)

bdb = build_DBpy.build_DBpy()
bdb.set_db_prefix(args.db_prefix)
bdb.set_kanalyze_input(args.kanalyze_input)
bdb.set_sti_file(args.sti_file)
bdb.set_fasta_file(args.fasta_file)
if args.binner:
    bdb.build_binner_dbs(args.kl)
else:
    bdb.build_dbs(args.kl, args.tax_cutoff, args.kmer_cutoff, args.protein, args.afterburner)
