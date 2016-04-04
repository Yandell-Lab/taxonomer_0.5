import numpy as np
import argparse
import itertools
import random
import multiprocessing as mp
from Bio import SeqIO

parser = argparse.ArgumentParser(description="Break ties in taxonomer output, meant for transcript abundance estimation")
parser.add_argument("taxonomer_tied_output", type=str, help="taxonomer output with ties")
parser.add_argument("reference_fasta", type=str, help="reference fasta file")
parser.add_argument("sti", type=str, help=".sti file")
parser.add_argument("effective_lengths", type=str, help="effective lengths file")
parser.add_argument("read_length", type=float, help="read length")
parser.add_argument("output_base", type=str, help="output base -- taxonomer output format with tie breaks")
parser.add_argument("-ecc", type=float, help="effective coverage cutoff", default = 2.0)
parser.add_argument("-np", type=int, help="number of processors to use", default = 1)
args = parser.parse_args()

class tiesChunk:

    def __init__(self,file):
        self.f = open(file,'r')
        self.mapped_reads = 0
        self.line = self.f.readline()

    def __iter__(self):
        return self

    def next(self):
        lines = []
        while self.line:
            if self.line[0] == "U":
                self.line = self.f.readline()
                continue
            self.mapped_reads += 1
            lines = [self.line]
            data = self.line.strip().split("\t")
            nties = int(data[5]) - 1
            for i in range(nties):
                lines.append(self.f.readline())
            self.line = self.f.readline()
            return lines

        raise StopIteration

    def number_mapped_reads(self):
        return self.mapped_reads

effective_coverage = {} #{taxid}->effective coverage
unique_taxid_counts = {}

print "parsing sti"
cur_transcript = None
transcript_taxid = {}
taxid_transcript = {}
f = open(args.sti,'r')
for line in f:
    if line[0] == ">":
        cur_transcript = line[1:].strip()
    else:
        transcript_taxid[cur_transcript] = int(line.strip())
        taxid_transcript[int(line.strip())] = cur_transcript
f.close()

print "reference lengths"
transcript_length = {}
handle = open(args.reference_fasta,'r')
for record in SeqIO.parse(handle, "fasta"):
    tid = transcript_taxid[record.id]
    transcript_length[tid] = len(str(record.seq))
handle.close()

def breakTies(lines):
    if len(lines) == 1:
        return lines[0]
    else:
        total_effective_coverage = 0.0
        read_probabilities = []
        use_lines = []
        for line in lines:
            data = line.strip().split("\t")
            if int(data[2]) in effective_coverage:
                use_lines.append(line)
                read_probabilities.append(effective_coverage[int(data[2])])
                total_effective_coverage += read_probabilities[-1]
        if total_effective_coverage < .00001: #if total effective coverage is very small, return a random line / unclassified
            ri = random.randint(0,len(lines)-1)
            data = lines[ri].strip().split("\t")
            return lines[ri]
            #data = lines[ri].strip().split("\t")
            #new_line = "U\t%s\n"%(data[1])
            #return new_line
        if len(use_lines) == 1:
            return use_lines[0]
        read_probabilities = [x/total_effective_coverage for x in read_probabilities]
        random_num = random.random()
        total_rp = 0.0
        for i,rp in enumerate(read_probabilities):
            total_rp += rp
            if total_rp >= random_num:
                data = use_lines[i].strip().split("\t")
                return use_lines[i]
        ri = random.randint(0,len(lines)-1) #if assignment is not made, return a random line
        data = lines[ri].strip().split("\t")
        return lines[ri]
        #new_line = "U\t%s\n"%(data[1])
        #return new_line

effective_lengths = {} #{taxid}->effective lengths
f = open(args.effective_lengths,'r')
for line in f:
    data = line.strip().split("\t")
    effective_lengths[int(data[0])] = float(data[1])
f.close()

print "starting first pass of taxonomer output"
#unique_taxid_counts = {} #taxid->counts
f = open(args.taxonomer_tied_output,'r') #first pass get unique taxid counts
for line in f:
    if line[0] == "U":
        continue
    data = line.strip().split("\t")
    if int(data[5]) == 1:
        tid = int(data[2])
        if tid not in unique_taxid_counts:
            unique_taxid_counts[tid] = 0
        unique_taxid_counts[tid] += 1
f.close()
print "done with first pass"

#get effective coverage estimates
for taxid in unique_taxid_counts:
    if abs(effective_lengths[taxid]) > .00001:
        effective_coverage[taxid] = (float(unique_taxid_counts[taxid])/effective_lengths[taxid])/args.read_length
    else:
        effective_coverage[taxid] = 0.0


#resolve ties with effective coverage
#p = mp.Pool(processes = args.np)
#results = p.imap(breakTies, tiesChunk(args.taxonomer_tied_output),50)
print "breaking ties"
tc = tiesChunk(args.taxonomer_tied_output)
results = itertools.imap(breakTies, tc)

out = open(args.output_base,'w')
for r in results:
    out.write(r)
out.close()
print "done breaking ties"
