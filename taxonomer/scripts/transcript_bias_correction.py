from sklearn.linear_model import LinearRegression
from sklearn import preprocessing
from Bio import SeqIO
import numpy as np
import sys
import argparse
import time


parser = argparse.ArgumentParser(description="Classify input reads")
parser.add_argument("taxonomer_input", type=str, help="tiebroken taxonomer input file name")
parser.add_argument("transcript_covariates", type=str, help="transcript specific covariates")
parser.add_argument("reference_fasta", type=str, help="reference fasta sequences")
parser.add_argument("sti", type=str, help=".sti file")
parser.add_argument("output", type=str, help="output file")
parser.add_argument("-rpkm", action="store_true", help="adjust rpkm values instead of counts")
args = parser.parse_args()

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

transcript_counts = {} #taxid -> counts
f = open(args.taxonomer_input,'r')
mapped_reads = 0
for line in f:
    if line[0] == "U":
        continue
    mapped_reads += 1
    data = line.strip().split("\t")
    #tid = int(data[2])
    tnid = data[3] #transcript name
    if tnid not in transcript_taxid:
    	print "problem, transcript %s not found  in sti"%(tnid)
    	continue
    tid = transcript_taxid[tnid]
    if tid not in transcript_counts:
        transcript_counts[tid] = 0
    transcript_counts[tid] += 1
f.close()

print "reference lengths"
transcript_length = {}
handle = open(args.reference_fasta,'r')
for record in SeqIO.parse(handle, "fasta"):
    tid = transcript_taxid[record.id]
    transcript_length[tid] = len(str(record.seq))
handle.close()

rpkms = {} #transcript -> rpkm
for tid,tcount in transcript_counts.items():
    rpkm = float(tcount)/((mapped_reads / 1000000.0)*(transcript_length[tid] / 1000.0))
    rpkms[tid] = rpkm

all_x = [] #covariates
all_y = [] #counts
trans_ids = [] #transcript_ids

f = open(args.transcript_covariates,'r')
for line in f:
    data = line.strip().split("\t")
    tid = transcript_taxid[data[0]]
    if args.rpkm:
        if tid in rpkms:
            all_y.append(rpkms[tid])
            all_x.append([float(x) for x in data[1:]])
            trans_ids.append(tid)
    else:
        if tid in transcript_counts:
            all_y.append(transcript_counts[tid])
            all_x.append([float(x) for x in data[1:]])
            trans_ids.append(tid)
f.close()

all_x = np.array(all_x)
all_y = np.log2(all_y)

print "normalizing data"
#my_normalization = preprocessing.StandardScaler().fit(all_x)
#all_x = my_normalization.transform(all_x)
print "fitting linear regression"
clf = LinearRegression(fit_intercept=True)
clf = clf.fit(all_x, all_y)
print "making predictions"
predictions = clf.predict(all_x)

count_mean = np.mean(all_y)
residuals = all_y - predictions
corrected_counts = residuals + count_mean

out = open(args.output,'w')
for i in xrange(len(all_y)):
    out.write("%s\t%d\t%d\t%f\t%f\n"%(taxid_transcript[trans_ids[i]],trans_ids[i],transcript_counts[trans_ids[i]],rpkms[trans_ids[i]],2**corrected_counts[i]))
out.close()
