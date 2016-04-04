cdef extern from "../../src/build_db.h":

	void build_dbs(int kl, int tc, int kc, int protein, int afterburner)
	void build_binner_dbs(int kl)
	void set_db_prefix(char * tdb)
	void set_kanalyze_input(char * tka)
	void set_sti_file(char * tsti)
	void set_fasta_file(char * tff)
