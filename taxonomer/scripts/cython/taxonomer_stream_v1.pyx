import classify_READpy
import binner_PY
import re
from libc.stdlib cimport malloc, free
import itertools
#import matplotlib
#matplotlib.use("Agg")
#import matplotlib.pyplot as plt

cdef class taxonomer_stream_v1:

	def __cinit__(self):
		pass

	cpdef taxonomer_v1(self, binner_classification_objs, stream_io):
		#binner_classification_objs -- binner, bacteria, fungi, protonomer, afterburner, human transcripts
		seqid = None
		mysplit = re.compile("\s+")
		sequence = ""
		for line in stream_io:
			#just fasta for now
			if line[0] == ">":
				if len(sequence) > 0:
					read_output = binner_classification_objs[2].classify_read_from_stream(seqid, sequence, len(sequence))
					stream_io.write(read_output)

				seqid = mysplit.split(line)[0][1:]
				sequence = ""
			else:
				sequence += line.strip()

