cdef extern from "../../src/classify.h":

	struct read_tax_scores:
		unsigned long long * taxids
		float * scores
		int * counts
		char * seqid
		unsigned long long * lca
		int * ntaxids

	void classify_reads(char * read_file, int np, char * output_file)
	char * classify_read_stream(char * read, char * read_name, long int read_length)
	void add_reads(char * read, char * read_name, long int read_length, int name_length, int classify)
	read_tax_scores * classify_read(char * read, long int read_len)
	void load_classification_dbs(char * mi_file, char * db_file, char * rsi_file, char * tri_file, char * sti_file, int load_mem, int ties, int protein, int afterburner)
	unsigned long long * get_kmer_info(char * kmer)
	void classify_db(char * db_file, char * sti_file, int read_len, char * output_file, int np)
	int get_kmer_length()
	unsigned long long get_total_kmer_count()
	unsigned long long get_unique_kmer_count()
	void set_kmer_cutoff(int kc)
	void set_number_processors(int np)
