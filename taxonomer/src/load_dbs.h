//
//  load_dbs.h
//  taxonomer
//
//  Created by Steven Flygare on 10/3/14.
//  Copyright (c) 2014 yandell lab. All rights reserved.
//

#ifndef __taxonomer__load_dbs__
#define __taxonomer__load_dbs__

#include "taxonomer_headers.h"
#include "kmer_counts.h"

#define SEQID_SIZE 100

struct seqid_taxid_single{
    char seqid[SEQID_SIZE];
    uint64_t taxid;
    UT_hash_handle hh;
};

struct taxid_seqid{
    uint64_t taxid;
    char ** seqids;
    int num_seqids;
    int max_seqs;
    UT_hash_handle hh;
};

struct taxid_confidence{
    uint64_t taxid;
    float confidence;
    UT_hash_handle hh;
};

extern struct tax_info tax_globals;
extern struct min_kmer_db ** mi;
extern unsigned char * mmkmer;
extern unsigned char * mmrsi;
struct seqid_taxid_single * seqid_taxid_rel;
struct taxid_seqid * tax_to_seqs;
struct taxid_confidence * tax_c_scores;

void load_minimizer_index(char * minimizer_file);
void load_kmer_db(char * kmer_db);
void load_rsi_db(char * rsi_file);
void load_seqid_taxid_rel(char * seqid_taxid_file);
void load_taxid_confidence(char * taxid_confidence_file);

#endif /* defined(__taxonomer__load_dbs__) */
