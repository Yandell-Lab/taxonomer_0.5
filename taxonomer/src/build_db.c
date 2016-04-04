//
//  build_db.c
//  taxonomer
//
//  Created by Steven Flygare on 10/7/14.
//  Copyright (c) 2014 yandell lab. All rights reserved.
//

#include "build_db.h"

void build_dbs(int kl, int tc, int kc, int protein, int afterburner){
    set_tax_globals();
    int ml = (int)kl / 2;
    tax_globals.kl = kl;
    tax_globals.ml = ml;
    tax_globals.tax_cutoff = tc;
    tax_globals.kmer_cutoff = kc;
    tax_globals.protein = protein;
    tax_globals.afterburner = afterburner;
    build_index(db_prefix, kanalyze_input);
    printf("done building index\n");
    load_seqid_taxid_rel(sti_file);
    populate_rsi(fasta_file, db_prefix);
}

void build_binner_dbs(int kl){
    set_tax_globals();
    int ml = (int)kl / 2;
    tax_globals.kl = kl;
    tax_globals.ml = ml;
    build_binner_index(db_prefix, kanalyze_input);
}

void set_db_prefix(char * tdb){
    memset(db_prefix, '\0', sizeof(char)*100);
    strcpy(db_prefix, tdb);
}

void set_kanalyze_input(char * tka){
    memset(kanalyze_input, '\0', sizeof(char)*100);
    strcpy(kanalyze_input, tka);
}

void set_sti_file(char * tsti){
    memset(sti_file, '\0', sizeof(char)*100);
    strcpy(sti_file, tsti);
}

void set_fasta_file(char * tff){
    memset(fasta_file, '\0', sizeof(char)*100);
    strcpy(fasta_file, tff);
}
