//
//  kmer_counts.h
//  taxonomer
//
//  Created by Steven Flygare on 10/3/14.
//  Copyright (c) 2014 yandell lab. All rights reserved.
//

#ifndef __taxonomer__kmer_counts__
#define __taxonomer__kmer_counts__

#include "taxonomer_headers.h"
#include "kmer_utils.h"
#include "load_dbs.h"
#include "query_db.h"

KHASH_MAP_INIT_INT64(m64, int);

struct min_kmer_db{
    uint64_t offset;
    uint64_t kmer_count;
};

//global variables
struct min_kmer_db ** mi;
extern struct tax_info tax_globals;
extern struct seqid_taxid_single * seqid_taxid_rel;
unsigned char * mmkmer;
unsigned char * mmrsi;

void build_index(char * db_prefix, char * kmer_db);
void build_binner_index(char * db_prefix, char * kmer_db);
void populate_rsi(char * fasta_file, char * db_prefix);

static void process_kanalyze_line(char * line, uint64_t * kmer, uint64_t * kmer_count);
static void get_kco(uint64_t kmer, uint64_t * info);
static void populate_kmer_counts(char * kmer_db, int * minimizer_kcount, char * db_prefix);
static void populate_binner_kmer_counts(char * kmer_db, int * minimizer_kcount, char * db_prefix);
static void get_kmer_counts(khash_t(m64) * kc, char * ref_seq, uint64_t ref_length, uint64_t * ref_kmer_count);
static void add_refs_rsi(khash_t(m64) * kc, uint64_t taxid, uint64_t ref_kmer_count);

#endif /* defined(__taxonomer__kmer_counts__) */
