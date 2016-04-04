//
//  classify.h
//  taxonomer
//
//  Created by Steven Flygare on 10/6/14.
//  Copyright (c) 2014 yandell lab. All rights reserved.
//

#ifndef __taxonomer__classify__
#define __taxonomer__classify__

#include "taxonomer_headers.h"
#include "load_dbs.h"
#include "kmer_utils.h"
#include "query_db.h"
#include "process_results.h"

extern struct tax_info tax_globals;
extern struct seqid_taxid_single * seqid_taxid_rel;
extern struct taxid_seqid * tax_to_seqs;
extern struct taxid_confidence * tax_c_scores;

struct read_tax_scores {
    uint64_t * taxids;
    float * scores;
    int * counts;
    char * seqid;
    uint64_t * lca;
    int * ntaxids;
};

struct read_info {
    char * read;
    uint64_t read_len;
    char * name;
};

void classify_reads(char * read_file, int np, char * output_file);
struct read_tax_scores * classify_read(char * read, uint64_t read_len);
void classify_db(char * db_file, char * sti_file, int read_len, char * output_file, int np);
char * classify_read_stream(char * read, char * read_name, uint64_t read_length);
void add_reads(char * read, char * read_name, uint64_t read_length, int name_length, int classify);
void load_classification_dbs(char * mi_file, char * db_file, char * rsi_file, char * tri_file, char * sti_file, int load_mem, int ties, int protein, int afterburner);
int64_t score_sort(struct read_tax_info * a, struct read_tax_info * b);
float score_sort_kw(struct read_tax_info * a, struct read_tax_info * b);
uint64_t * get_kmer_info(char * kmer);
int get_kmer_length();
uint64_t get_total_kmer_count();
uint64_t get_unique_kmer_count();
void set_kmer_cutoff(int kc);
void set_number_processors(int np);

#endif /* defined(__taxonomer__classify__) */
