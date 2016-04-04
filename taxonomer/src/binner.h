//
//  binner.h
//  taxonomer
//
//  Created by steven on 11/21/14.
//  Copyright (c) 2014 yandell lab. All rights reserved.
//

#ifndef __taxonomer__binner__
#define __taxonomer__binner__

#include "taxonomer_headers.h"
#include "load_dbs.h"
#include "kmer_utils.h"
#include "query_db.h"

extern struct tax_info tax_globals;
extern struct seqid_taxid_single * seqid_taxid_rel;
extern struct taxid_seqid * tax_to_seqs;
extern struct taxid_confidence * tax_c_scores;
//static const int BIN_MAX = 65535;

struct bin_info {
    uint64_t * bins;
    int * counts;
    int nbins;
};

void bin_reads(char * read_file, int np, char * output_file);
uint64_t bin_read_stream(char * read_seq, int read_len, int bin_cutoff);
void bin_kmers(char * read_seq, int read_len);
void load_binner_dbs(char * mi_file, char * db_file, int load_mem);

#endif /* defined(__taxonomer__binner__) */
