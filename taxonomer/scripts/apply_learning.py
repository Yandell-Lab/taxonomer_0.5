from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn import preprocessing
import sqlite3 as sql
import numpy
import sys
import argparse
import time


parser = argparse.ArgumentParser(description="Classify input reads")
parser.add_argument("taxonomer_input", type=str, help="taxonomer input file name")
parser.add_argument("learn_input", type=str, help="file to learn from")
parser.add_argument("tri_file", type=str, help="tri file")
parser.add_argument("output", type=str, help="output file")
args = parser.parse_args()

f = open(args.tri_file,'r')
taxid_rel = {} #child->parent
cur_child = None
for line in f:
	if line[0] == ">":
		cur_child = line[1:].strip()
	else:
		taxid_rel[cur_child] = line.strip()
f.close()

print "building random forest..."
all_x = []
all_y = []
f = open(args.learn_input,'r')
for line in f:
	data = line.strip().split("\t")
	all_y.append(int(data[0]))
	all_x.append([float(x) for x in data[1:]])
f.close()

all_x = numpy.array(all_x)
all_y = numpy.array(all_y)

#separate into test and train sets and fit the normalization on the training data
#normalize the training data
#normalize the test data with same parameters as training data

indices = numpy.random.permutation(len(all_y))
train_cut = int(len(all_y)*.8)
train_x = all_x[indices[indices < train_cut]]
train_y = all_y[indices[indices < train_cut]]

test_x = all_x[indices[indices >= train_cut]]
test_y = all_y[indices[indices >= train_cut]]

my_normalization = preprocessing.StandardScaler().fit(train_x)
train_x = my_normalization.transform(train_x)
test_x = my_normalization.transform(test_x)

#train random forest classification
#print "fitting forest"
clf = RandomForestClassifier(n_estimators=15)
clf = clf.fit(train_x, train_y)
print "mean accuracy on test set: %f"%(clf.score(test_x,test_y))

print "processing reads"
f = open(args.taxonomer_input,'r')
out = open(args.output,'w')
for line in f:
	if line[0] == "U":
		out.write(line)
		continue
	data = line.strip().split("\t")
	read_data = [float(x) for x in data[4:]]
	tran_vals = my_normalization.transform(read_data)
	probs = clf.predict_proba(tran_vals)
	new_taxid = data[2]
	if probs[0][1] < .7:
		new_taxid = taxid_rel[new_taxid]
	out.write("%s\t%s\t%s\t%f\n"%("\t".join(data[:2]),new_taxid,"\t".join(data[3:]),probs[0][1]))

f.close()
