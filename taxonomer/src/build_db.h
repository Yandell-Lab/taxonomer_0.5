//
//  build_db.h
//  taxonomer
//
//  Created by Steven Flygare on 10/6/14.
//  Copyright (c) 2014 yandell lab. All rights reserved.
//

#ifndef __taxonomer__build_db__
#define __taxonomer__build_db__

#include "taxonomer_headers.h"
#include "kmer_utils.h"
#include "kmer_counts.h"
#include "load_dbs.h"

static char db_prefix[100];
static char kanalyze_input[100];
static char sti_file[100];
static char fasta_file[100];

void build_dbs(int kl, int tc, int kc, int protein, int afterburner);
void build_binner_dbs(int kl);
void set_db_prefix(char * tdb);
void set_kanalyze_input(char * tka);
void set_sti_file(char * tsti);
void set_fasta_file(char * tff);

#endif /* defined(__taxonomer__build_db__) */
