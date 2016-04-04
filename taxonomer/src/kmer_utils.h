//
//  kmer_utils.h
//  taxonomer
//
//  Created by steven on 10/3/14.
//  Copyright (c) 2014 yandell lab. All rights reserved.
//

#ifndef __taxonomer__kmer_utils__
#define __taxonomer__kmer_utils__

#include "taxonomer_headers.h"

struct tax_info {
    int KMER_FLAG;
    uint64_t INDEX2_XOR_MASK;
    int KMER_NOT_FOUND;
    uint8_t kl;
    uint8_t ml;
    uint64_t total_kmer_count;
    uint64_t unique_kmer_count;
    int tax_cutoff;
    int kmer_cutoff;
    int output_ties;
    int protein;
    int afterburner;
    int load_mem;
    int np;
};

//global variable
struct tax_info tax_globals;

uint64_t kmer_minimizer(uint64_t);
uint64_t canonical_representation(uint64_t, uint8_t);
uint64_t canonical_representation2(uint64_t);
uint64_t reverse_complement(uint64_t, uint8_t);
int next_kmer(char *, uint64_t *);
int next_protein_kmer(char * read_seq, uint64_t * kmer);
void set_kl(uint8_t);
void set_ml(uint8_t);
void set_tax_globals();

#endif /* defined(__taxonomer__kmer_utils__) */
