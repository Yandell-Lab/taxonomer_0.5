cimport buildDB

cdef class build_DBpy:

	def __cinit__(self):
		pass

	cpdef build_dbs(self, int kl, int tc, int kc, int protein, int afterburner):
		buildDB.build_dbs(kl, tc, kc, protein, afterburner)

	cpdef build_binner_dbs(self, int kl):
		buildDB.build_binner_dbs(kl)

	cpdef set_db_prefix(self, char * tdb):
		buildDB.set_db_prefix(tdb)

	cpdef set_kanalyze_input(self, char * tka):
		buildDB.set_kanalyze_input(tka)

	cpdef set_sti_file(self, char * tsti):
		buildDB.set_sti_file(tsti)

	cpdef set_fasta_file(self, char * tff):
		buildDB.set_fasta_file(tff)
