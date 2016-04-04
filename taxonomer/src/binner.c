//
//  binner.c
//  taxonomer
//
//  Created by steven on 11/21/14.
//  Copyright (c) 2014 yandell lab. All rights reserved.
//

#include "binner.h"

#define READS_PER_BLOCK 1000
#define READ_OUTPUT_SIZE 500

KSEQ_INIT(gzFile, gzread)

struct read_info {
    char * read;
    uint64_t read_len;
    char * name;
};

struct bin_counts {
    uint64_t bin;
    int count;
    UT_hash_handle hh;
};

void bin_kmers(char * read_seq, int read_len){
    int i = 0;
    struct min_kmer_db tmi;
    tmi.kmer_count = 0;
    tmi.offset = 0;
    uint64_t kmer_mask = ~0;
    kmer_mask >>= sizeof(kmer_mask) * 8 - (tax_globals.kl * 2);
    while (i + tax_globals.kl <= read_len) {
        uint64_t kmer = 0;
        int flag = next_kmer(read_seq + i, &kmer);
        if (flag == tax_globals.KMER_FLAG) {
            kmer &= kmer_mask;
            //uint64_t bin = query_kc(kmer);
            uint64_t bin = query_kc2(kmer, &tmi);
            char tmpstr[tax_globals.kl + 1];
            memset(tmpstr, '\0', sizeof(char)*(tax_globals.kl + 1));
            strncpy(tmpstr, read_seq + i, tax_globals.kl);
            printf("%s\t%llu\t%llu\t%llu\t%d\n",tmpstr,kmer,canonical_representation2(kmer),bin,i);
            i++;
        }
        else {
            i += (flag+1);
        }
    }
}

/*
void bin_read2(char * read_seq, int read_len, struct bin_info * bi){
    //bi->bins = (uint64_t *)malloc(sizeof(uint64_t)*(BIN_MAX+1));
    bi->counts = (int *)calloc(BIN_MAX+1, sizeof(int));
    bi->nbins = 0;
    int i = 0;
    uint64_t * offsets = (uint64_t *)malloc(sizeof(uint64_t)*KMER_PER_MINIMIZER);
    uint64_t offset_len = KMER_PER_MINIMIZER;
    struct min_kmer_db tmi;
    tmi.kmer_count = 0;
    tmi.offset = 0;
    uint64_t kmer_mask = ~0;
    kmer_mask >>= sizeof(kmer_mask) * 8 - (tax_globals.kl * 2);
    while (i + tax_globals.kl <= read_len) {
        uint64_t kmer = 0;
        int flag = next_kmer(read_seq + i, &kmer);
        if (flag == tax_globals.KMER_FLAG) {
            kmer &= kmer_mask;
            uint64_t bin = query_kc2(kmer, &tmi);
            if (bi->counts[bin] == 0) {
                bi->nbins += 1;
            }
            bi->counts[bin] += 1;
            i++;
        }
        else {
            i += (flag+1);
        }
    }
    free(offsets);
    
}
*/

uint64_t bin_read_stream(char * read_seq, int read_len, int bin_cutoff){
    struct bin_counts * bcs = NULL;
    struct bin_counts * tbcs = NULL;
    int i = 0;
    struct min_kmer_db tmi;
    tmi.kmer_count = 0;
    tmi.offset = 0;
    uint64_t kmer_mask = ~0;
    kmer_mask >>= sizeof(kmer_mask) * 8 - (tax_globals.kl * 2);
    while (i + tax_globals.kl <= read_len) {
        uint64_t kmer = 0;
        int flag = next_kmer(read_seq + i, &kmer);
        if (flag == tax_globals.KMER_FLAG) {
            kmer &= kmer_mask;
            //uint64_t bin = query_kc(kmer);
            uint64_t bin = query_kc2(kmer, &tmi);
            //printf("%llu\t%llu\t%d\n",kmer,bin,i);
            HASH_FIND_INT(bcs, &bin, tbcs);
            if (tbcs == NULL) {
                tbcs = (struct bin_counts *)malloc(sizeof(struct bin_counts));
                tbcs->bin = bin;
                tbcs->count = 1;
                HASH_ADD_INT(bcs, bin, tbcs);
            }
            else{
                tbcs->count++;
            }
            i++;
        }
        else {
            i += (flag+1);
        }
    }
    //free(offsets);
    struct bin_counts * tmp = NULL;
    int index = 0;
    int zero_counts = 0;
    uint64_t max_bin = 0;
    int max_counts = 0;
    HASH_ITER(hh, bcs, tbcs, tmp) {
        if (tbcs->bin == 0) {
            zero_counts = tbcs->count;
        }
        else if (tbcs->count > max_counts){
            max_counts = tbcs->count;
            max_bin = tbcs->bin;
        }
        else if (tbcs->count == max_counts){
            max_bin = max_bin | tbcs->bin;
        }
        index++;
        free(tbcs);
    }
    
    return max_counts >= bin_cutoff ? max_bin : 0;
}


void bin_read(char * read_seq, int read_len, struct bin_info * bi){
    struct bin_counts * bcs = NULL;
    struct bin_counts * tbcs = NULL;
    int i = 0;
    struct min_kmer_db tmi;
    tmi.kmer_count = 0;
    tmi.offset = 0;
    uint64_t kmer_mask = ~0;
    kmer_mask >>= sizeof(kmer_mask) * 8 - (tax_globals.kl * 2);
    while (i + tax_globals.kl <= read_len) {
        uint64_t kmer = 0;
        int flag = next_kmer(read_seq + i, &kmer);
        if (flag == tax_globals.KMER_FLAG) {
            kmer &= kmer_mask;
            //uint64_t bin = query_kc(kmer);
            uint64_t bin = query_kc2(kmer, &tmi);
            //printf("%llu\t%llu\t%d\n",kmer,bin,i);
            HASH_FIND_INT(bcs, &bin, tbcs);
            if (tbcs == NULL) {
                tbcs = (struct bin_counts *)malloc(sizeof(struct bin_counts));
                tbcs->bin = bin;
                tbcs->count = 1;
                HASH_ADD_INT(bcs, bin, tbcs);
            }
            else{
                tbcs->count++;
            }
            i++;
        }
        else {
            i += (flag+1);
        }
    }
    //free(offsets);
    int nbins = HASH_COUNT(bcs);
    bi->nbins = nbins;
    if (nbins < 1) {
        bi->nbins = 0;
        return;
    }
    bi->bins = (uint64_t *)malloc(sizeof(uint64_t)*nbins);
    bi->counts = (int *)malloc(sizeof(int)*nbins);
    struct bin_counts * tmp = NULL;
    int index = 0;
    HASH_ITER(hh, bcs, tbcs, tmp) {
        bi->bins[index] = tbcs->bin;
        bi->counts[index] = tbcs->count;
        index++;
        free(tbcs);
    }
}

void bin_read_chunk(struct read_info ** reads, int num_reads, char output[READS_PER_BLOCK][READ_OUTPUT_SIZE]){
    int i= 0;
#pragma omp parallel for
    for (i=0; i < num_reads; i++) {
        struct bin_info bi;
        char tmp_output[400];
        memset(tmp_output, '\0', sizeof(char)*400);
        //printf("%s\t%llu\n",reads[i]->name,reads[i]->read_len);
        bin_read(reads[i]->read, reads[i]->read_len, &bi);
        /*if (bi.nbins < 1) {
            sprintf(output[i], "0-0\t");
            continue;
        }
        int j = 0;
        for (j=0; j < BIN_MAX+1; j++) {
            if (bi.counts[j] > 0) {
                char small_tmp[50];
                memset(small_tmp, '\0', 50*sizeof(char));
                sprintf(small_tmp, "%d-%d;", j,bi.counts[j]);
                strcat(tmp_output, small_tmp);
            }
        }
        sprintf(output[i], "%s\t%s\t", reads[i]->name, tmp_output);
        free(bi.counts);
        */
        if (bi.nbins > 0) {
            char small_tmp2[50];
            memset(small_tmp2, '\0', 50*sizeof(char));
            sprintf(small_tmp2, "%llu-%d", bi.bins[0],bi.counts[0]); //print without semi-colon in front
            strcpy(tmp_output, small_tmp2);
            int j = 1;
            for (j = 1; j < bi.nbins; j++) {
                char small_tmp[50];
                memset(small_tmp, '\0', 50*sizeof(char));
                sprintf(small_tmp, ";%llu-%d", bi.bins[j],bi.counts[j]);
                strcat(tmp_output, small_tmp);
            }
            sprintf(output[i], "%s\t%s", reads[i]->name, tmp_output);
            free(bi.bins);
            free(bi.counts);
        }
        else {
            sprintf(output[i], "%s\t0-0", reads[i]->name);
        }
    }
}

void bin_reads(char * read_file, int np, char * output_file){
#ifdef _OPENMP
    omp_set_num_threads(np);
#endif
    
    FILE * out = fopen(output_file, "w");
    
    gzFile seq_index_file = gzopen(read_file, "r");
    kseq_t * seq_t = kseq_init(seq_index_file);
    int l = 0;
    struct read_info ** reads = (struct read_info **)malloc(sizeof(struct read_info *)*READS_PER_BLOCK);
    int num_reads = 0;
    uint64_t reads_read = 0;
    uint64_t total_processed = 0;
    while ((l = kseq_read(seq_t)) >= 0) {
        reads_read++;
        if (num_reads < READS_PER_BLOCK) {
            struct read_info * myseq = (struct read_info *)malloc(sizeof(struct read_info));
            myseq->read = (char *)calloc(seq_t->seq.l, sizeof(char));
            memcpy(myseq->read, seq_t->seq.s, sizeof(char)*seq_t->seq.l);
            myseq->name = (char *)calloc(seq_t->name.l + 1, sizeof(char));
            memcpy(myseq->name, seq_t->name.s, sizeof(char)*seq_t->name.l);
            myseq->read_len = seq_t->seq.l;
            reads[num_reads] = myseq;
            num_reads++;
        }
        else {
            
            char output[READS_PER_BLOCK][READ_OUTPUT_SIZE];
            bin_read_chunk(reads, num_reads, output);
            int i = 0;
            for (i=0; i < num_reads; i++) {
                fprintf(out, "%s\t%s\n", output[i], reads[i]->read);
                free(reads[i]->read);
                free(reads[i]->name);
                free(reads[i]);
            }
            
            total_processed += num_reads;
            num_reads = 0;
            struct read_info * myseq = (struct read_info *)malloc(sizeof(struct read_info));
            myseq->read = (char *)calloc(seq_t->seq.l, sizeof(char));
            memcpy(myseq->read, seq_t->seq.s, sizeof(char)*seq_t->seq.l);
            myseq->name = (char *)calloc(seq_t->name.l + 1, sizeof(char));
            memcpy(myseq->name, seq_t->name.s, sizeof(char)*seq_t->name.l);
            myseq->read_len = seq_t->seq.l;
            reads[num_reads] = myseq;
            num_reads++;
        }
    }
    
    char output[READS_PER_BLOCK][READ_OUTPUT_SIZE];
    bin_read_chunk(reads, num_reads, output);
    int i = 0;
    for (i=0; i < num_reads; i++) {
        fprintf(out, "%s\t%s\n", output[i], reads[i]->read);
        free(reads[i]->read);
        free(reads[i]->name);
        free(reads[i]);
    }
    total_processed += num_reads;
    printf("total reads processed: %llu\n",total_processed);
    
    gzclose(seq_index_file);
    kseq_destroy(seq_t);
    free(reads);
    fclose(out);
}

void load_binner_dbs(char * mi_file, char * db_file, int load_mem){
    set_tax_globals();
    tax_globals.load_mem = load_mem;
    printf("loading kmer index...\n");
    load_kmer_db(db_file);
    printf("done loading kmer index\n");
    printf("loading minimizer index...");
    load_minimizer_index(mi_file);
    printf("done\n");
    //printf("loading taxid scores...");
    //load_taxid_confidence(taxid_confidence);
    //printf("done\n");
}

