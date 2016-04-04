cimport binnerPY

cdef class binner_PY:

	def __cinit__(self):
		pass

	cpdef bin_reads(self, char * read_file, int np, int load_mem, char * output_file, char * mi_file, char * db_file):
		binnerPY.load_binner_dbs(mi_file, db_file, load_mem)
		binnerPY.bin_reads(read_file, np, output_file)

	cpdef bin_kmers(self, char * read, int read_len, int load_mem, char * mi_file, char * db_file):
		binnerPY.load_binner_dbs(mi_file, db_file, load_mem)
		binnerPY.bin_kmers(read, read_len)

	cpdef load_binner_db(self, int load_mem, char * mi_file, char * db_file):
		binnerPY.load_binner_dbs(mi_file,db_file,load_mem)

	cpdef bin_read_stream(self, char * read_seq, int read_len, int bin_cutoff):
		return binnerPY.bin_read_stream(read_seq, read_len, bin_cutoff)
