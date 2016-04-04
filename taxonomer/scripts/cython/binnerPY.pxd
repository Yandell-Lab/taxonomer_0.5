cdef extern from "../../src/binner.h":
  void bin_reads(char * read_file, int np, char * output_file)
  void load_binner_dbs(char * mi_file, char * db_file, int load_mem)
  void bin_kmers(char * read_seq, int read_len)
  unsigned long long bin_read_stream(char * read_seq, int read_len, int bin_cutoff)
