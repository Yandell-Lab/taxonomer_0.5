//
//  query_db.h
//  taxonomer
//
//  Created by Steven Flygare on 10/3/14.
//  Copyright (c) 2014 yandell lab. All rights reserved.
//

#ifndef __taxonomer__query_db__
#define __taxonomer__query_db__

#include "taxonomer_headers.h"
#include "kmer_counts.h"
#include "kmer_utils.h"

struct read_tax_info{
    uint64_t taxid;
    float score;
    int counts;
    UT_hash_handle hh;
};

extern struct min_kmer_db ** mi;
extern struct tax_info tax_globals;
extern unsigned char * mmkmer;
extern unsigned char * mmrsi;
static int KMER_PER_MINIMIZER = 500000;

uint64_t query_kc(uint64_t kmer);
uint64_t query_kc2(uint64_t kmer, struct min_kmer_db * tmi);

void query_kco(uint64_t kmer, uint64_t * info);
void query_kcrn(uint64_t kmer, uint64_t * info);

uint64_t kmer_binary_search(uint64_t * offsets, uint64_t offset_size, uint64_t kmer);

int process_read2(char * read_seq, uint64_t read_len, struct read_tax_info ** tax_hits);

void process_read(char * read_seq, uint64_t read_len, struct read_tax_info ** tax_hits);

void query_kmer(uint64_t kmer, struct min_kmer_db * tmi, uint64_t ** offsets, uint64_t * offset_len, struct read_tax_info ** tax_hits);

int query_kmer2(uint64_t kmer, struct min_kmer_db * tmi, uint64_t ** offsets, uint64_t * offset_len, struct read_tax_info ** tax_hits);

uint64_t query_minimizer(uint64_t kmer, uint64_t ckmer, struct min_kmer_db * tmi, uint64_t ** offsets, uint64_t * offset_len);

#endif /* defined(__taxonomer__query_db__) */
