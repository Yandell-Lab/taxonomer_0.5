//
//  main.c
//  taxonomer
//
//  Created by steven on 10/3/14.
//  Copyright (c) 2014 yandell lab. All rights reserved.
//

#include "taxonomer_headers.h"
//#include "../include/build_db.h"
#include "classify.h"
//#include "../include/binner.h"

extern struct tax_info tax_globals;
KSEQ_INIT(gzFile, gzread)

int main(int argc, const char * argv[]) {
    
    /*
    tax_globals.kmer_cutoff = 10000;
    tax_globals.tax_cutoff = 500;
    tax_globals.kl = 31;
    tax_globals.ml = 15;
    tax_globals.KMER_FLAG = -1;
    tax_globals.KMER_NOT_FOUND = 7;
    tax_globals.INDEX2_XOR_MASK = 0xe37e28c4271b5a2dULL;
    */
     
    /*
    char * mi_file = "/Users/steven/Desktop/gg_dbs/tmp.mi";
    char * tbi_file = "/Users/steven/Desktop/gg_dbs/tmp.tbi";
    char * rsi_file = "/Users/steven/Desktop/gg_dbs/tmp.rsi";
    char * tri_file = "/Users/steven/Lab/github/TaxaMer/test_data/gg_subset/taxid_rel_fa.tri";
    load_classification_dbs(mi_file, tbi_file, rsi_file, tri_file);
    uint64_t * class_res =  classify_read("CGTAAAGCGCACGCAGGCGGTTTGTTAAGTCAGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCTGATACTGGCAAGCTTGAGTCTCGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTA", 134);
    printf("%llu\t%llu\n",class_res[0],class_res[1]);
    */
    
    //<mi file> <tbi file> <rsi file> <tri file> <read file>
    
    char * mi_file = (char *)argv[1];
    char * tbi_file = (char *)argv[2];
    char * rsi_file = (char *)argv[3];
    char * tri_file = (char *)argv[4];
    char * sti_file = (char *)argv[5];
    //char * taxid_scores = (char *)argv[6];
    char * read_file = (char *)argv[6];
    load_classification_dbs(mi_file, tbi_file, rsi_file, tri_file, sti_file, 1, 0, 0, 0);
    //tax_globals.output_ties = 1;
    classify_reads(read_file, 1, "/Users/steven/Desktop/debug_small/tmp.out");
    //struct read_tax_scores * rts = classify_read("AACAGGATTAGATACCCTGGTAGTCCACGCC", 31);
    //int i=0;
    
    //load_binner_dbs("/Users/steven/Desktop/debug_small/tmp_binner.bmi", "/Users/steven/Desktop/debug_small/tmp_binner.btbi", 0);
    //bin_reads("/Users/steven/Desktop/debug_small/tmp.fa", 4, "/Users/steven/Desktop/debug_small/tmp.out");
    //load_binner_dbs(argv[1], argv[2], 0);
    //bin_reads(argv[3], 1, argv[4]);
     
    /*
    gzFile seq_tax = gzopen(argv[5], "r");
    kseq_t * seq_t = kseq_init(seq_tax);
    int l = 0;
    
    while ((l = kseq_read(seq_t)) >= 0) {
        uint64_t * class_res = classify_read(seq_t->seq.s, seq_t->seq.l);
    }
    gzclose(seq_tax);
    */
    
    return 0;
}
