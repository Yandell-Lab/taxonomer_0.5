import os
path = os.path.dirname(os.path.abspath(__file__))
import sys
sys.path.append("%s/cython/"%(path))
import argparse
import itertools
import multiprocessing as mp
import binner_splitpy

parser = argparse.ArgumentParser(description="Split binner output into specified output and give summary")
parser.add_argument("-b", "--binner_results", help="<file>\tthe output from the binner.", action="store",
                    required=True)

parser.add_argument("-k", "--database_key", help="<file>  a tab seperated list of the database (int) and description \
                    (string). Each line should cooresponds to a new database. if not provide db numbers will be used.",
                    action="store", required=False, default=None)

parser.add_argument("--use_default_key", help="<file>\tIf specified will use the default database key for \
                    binner_db_2.0.jf. [Default=True]", action="store_false", required=False, default=True)

parser.add_argument("-d", "--bin_db_output", help="<int>\tthe id of the binner database required for output",
                    action="store", required=True, type=int, default=None)

parser.add_argument("-c", "--kmer_cutoff", help="<int>\tnumber of kmers required to indicated a sequence belongs to \
                    selected bin. [DEFAULT=11]", action="store", type=int, required=False, default=11)

parser.add_argument("-o", "--base_output", help="<string>\tthe base file name to output binned reads. \
                    <.fa OR _1.fa, _2.fa> will be automatically appended.", action="store", default=None,
                    required=True)

parser.add_argument("-r", "--range_between_dbs", help="<int>\tnumber of kmers required to seperate read as belong to \
                    one bin over the other. dbs indicated with the --intersections_ok flag are igonored. reads within \
                    the range will be summerized as belonging to both bins. [Default=0]", action="store", default=0,
                    required=False)

parser.add_argument("-i", "--intersections_ok", help="<list,>\tcomma seperated list of databases are considered OK to \
                    intersect with --bin_db_output. This is required for database that are subset of other databases. \
                    (e.g. human transcripts and the human genome; Molusk LSU/SSU and other euk LSU/SSU)",
                    required=False, type=str, action="store", default=None)

parser.add_argument("-np", help="number of cups", required=False, type=int, action="store", default = 1)

args = parser.parse_args()


collisions = args.bin_db_output
if args.intersections_ok:
    for c in args.intersections_ok.split(","):
        collisions = collisions | long(c)

def process_line(line):
    if len(line.strip()) == 0:
        sys.stderr.write("empty line found")
        return ""
    else:
        return binner_splitpy.parse_line(line, args.kmer_cutoff, args.bin_db_output, collisions)

p = mp.Pool(processes = args.np)

out = open(args.base_output+".fa",'w')
f = open(args.binner_results,'r')
#results = itertools.imap(process_line, f)
results = p.imap(process_line, f, 500)
for r in results:
    if r:
        out.write(r)
out.close()
f.close()
