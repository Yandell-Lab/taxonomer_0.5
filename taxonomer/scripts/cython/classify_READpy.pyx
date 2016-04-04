cimport classifyREAD
from cython.view cimport array as cvarray
import numpy as np
import re
from libc.stdlib cimport malloc, free
#import matplotlib
#matplotlib.use("Agg")
#import matplotlib.pyplot as plt

cdef struct read_tax_scores:
	unsigned long long * taxids
	float * scores
	int * counts
	char * seqid
	unsigned long long * lca
	int * ntaxids

cdef class classify_READpy:

	def __cinit__(self):
		pass

	cpdef load_classification_dbs(self, char * mi_file, char * db_file, char * rsi_file, char * tri_file, char * sti_file, int load_mem, int ties, int protein, int afterburner):
		classifyREAD.load_classification_dbs(mi_file, db_file, rsi_file, tri_file, sti_file, load_mem, ties, protein, afterburner)

	cpdef classify_read(self, char * read, long int read_len):
		cdef read_tax_scores * results = classifyREAD.classify_read(read, read_len)
		cdef int tmp_ntaxids = results.ntaxids[0]
		if tmp_ntaxids < 1:
			return 0, 0, None, None, None
		taxids_np = np.asarray(<unsigned long long[:tmp_ntaxids]> results.taxids)
		scores_np = np.asarray(<float[:tmp_ntaxids]> results.scores)
		counts_np = np.asarray(<int[:tmp_ntaxids]> results.counts)
		#n, bins, patches = plt.hist(scores_np, 50, normed=1, facecolor='green', alpha=0.5)
		#plt.savefig("myfig")
		return results.ntaxids[0], results.lca[0], taxids_np, scores_np, counts_np

	cpdef get_kmer_info(self, char * kmer):
		cdef unsigned long long * info = classifyREAD.get_kmer_info(kmer)
		return info[0],info[1],info[2],info[3]

	def classify_reads(self, char * read_file, int np, char * output_file):
		classifyREAD.classify_reads(read_file, np, output_file)

	def classify_reads_stream_parallel(self, socket_io):
		#only classify fasta for now
		mysplit = re.compile("\s+")
		sequence = ""
		seqid = ""
		cdef int num_seqs = 0
		for line in socket_io:
			socket_io.write(line)
			#print line.strip()
			if line[0] == ">":
				if len(sequence) > 1:
					classifyREAD.add_reads(sequence,seqid,len(sequence),len(seqid),0)
					num_seqs += 1
				seqid = mysplit.split(line)[0][1:]
				sequence = ""
			else:
				sequence += line.strip()
		if len(sequence) > 1:
			classifyREAD.add_reads(sequence,seqid,len(sequence),len(seqid),0)
		classifyREAD.add_reads(NULL,NULL,0,0,1)

	def classify_read_from_stream(self, char * read_name, char * read_seq, long int read_len):
		cdef char * read_output = classifyREAD.classify_read_stream(read_seq,read_name,read_len)
		pyread_out = read_output + ""
		free(read_output)
		return pyread_out #note:  read_output was allocated but never freed, need to manage memory after return

	def classify_db(self, char * db_file, char * sti_file, int read_len, char * output_file, int np):
		classifyREAD.classify_db(db_file, sti_file, read_len, output_file, np)

	def get_kmer_length(self):
		return classifyREAD.get_kmer_length()

	def get_total_kmer_count(self):
		return classifyREAD.get_total_kmer_count()

	def get_unique_kmer_count(self):
		return classifyREAD.get_unique_kmer_count()

	def set_kmer_cutoff(self, int kc):
		classifyREAD.set_kmer_cutoff(kc)

	def set_number_processors(self, int np):
		classifyREAD.set_number_processors(np)
