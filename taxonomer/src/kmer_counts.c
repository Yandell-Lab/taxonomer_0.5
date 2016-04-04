//
//  kmer_counts.c
//  taxonomer
//
//  Created by Steven Flygare on 10/3/14.
//  Copyright (c) 2014 yandell lab. All rights reserved.
//

#include "kmer_counts.h"

KSEQ_INIT(gzFile, gzread)

static void process_kanalyze_line(char * line, uint64_t * kmer, uint64_t * kmer_count){
    char * pch;
    char * saveptr;
    pch = strtok_r (line,"\t",&saveptr);
    int tok_count = 0;
    while (pch != NULL)
    {
        //std::cout << pch << std::endl;
        if (tok_count == 0) {
            *kmer = strtoull(pch, NULL, 16);
        }
        else if (tok_count == 1) {
            *kmer_count = strtoull(pch, NULL, 10);
        }
        else{
            printf("problem with line %s , shouldn't have more than 2 tokens\n",line);
        }
        pch = strtok_r (NULL, "\t", &saveptr);
        tok_count++;
    }

}


void build_index(char * db_prefix, char * kmer_db){
    FILE * infile = fopen(kmer_db, "r");
    if (! infile) {
        printf("file %s cannot be read\n", kmer_db);
        return;
    }
    
    char line[100];
    int * minimizer_kcount = calloc((uint64_t)pow(4,tax_globals.ml), sizeof(int));
    uint64_t kmer, kmer_count;
    uint64_t unique_kmer_counts = 0;
    uint64_t total_kmer_counts = 0;
    while(fgets(line, 100, infile)){
        unique_kmer_counts++;
        process_kanalyze_line(line, &kmer, &kmer_count);
        total_kmer_counts += kmer_count;
        uint64_t minimizer = kmer_minimizer(kmer);
        minimizer_kcount[minimizer]++;
    }
    fclose(infile);
    
    tax_globals.total_kmer_count = total_kmer_counts;
    tax_globals.unique_kmer_count = unique_kmer_counts;
    
    char tmini[100];
    memset(tmini, '\0', 100);
    strcpy(tmini, db_prefix);
    strcat(tmini, ".mi");
    
    char tbkmer[100];
    memset(tbkmer, '\0', 100);
    strcpy(tbkmer, db_prefix);
    strcat(tbkmer, ".tbi");
    
    FILE * mini = fopen( tmini ,"wb");
    FILE * bkmer = fopen( tbkmer , "wb");
    
    //write head to bkmer file
    fwrite(&tax_globals.kl, sizeof(uint8_t), 1, bkmer);
    fwrite(&tax_globals.ml, sizeof(uint8_t), 1, bkmer);
    fwrite(&unique_kmer_counts, sizeof(uint64_t), 1, bkmer);
    fwrite(&total_kmer_counts, sizeof(uint64_t), 1, bkmer);
    fwrite(&tax_globals.kmer_cutoff, sizeof(int), 1, bkmer);
    fwrite(&tax_globals.tax_cutoff, sizeof(int), 1, bkmer);
    
    uint64_t offset = 2*sizeof(uint8_t) + 2*sizeof(uint64_t) + 2*sizeof(int);
    uint64_t ref_score_offset = 0;
    uint64_t total_minimizers = (uint64_t)pow(4, tax_globals.ml);
    mi = (struct min_kmer_db **)calloc((uint64_t)pow(4, tax_globals.ml), sizeof(struct min_kmer_db *));
    uint64_t tkmer = 0;
    
    //write basic organization to both minimizer and kmer binary files
    uint64_t i = 0;
    for (i=0; i < total_minimizers; i++) {
        if (minimizer_kcount[i] == 0) {
            continue;
        }
        uint64_t * offsets = (uint64_t *)malloc(sizeof(uint64_t)*minimizer_kcount[i]);
        uint64_t offset_index = 0;
        uLong c = 0;
        
        int j=0;
        for (j = 0; j < minimizer_kcount[i]; j++) {
            
            offsets[offset_index] = offset;
            offset_index += 1;
            
            c = fwrite(&(tkmer), sizeof(uint64_t), 1, bkmer); //writing kmer
            if (c != 1) {
                printf("incorrect number of elements written1\n");
            }
            offset += sizeof(uint64_t);
            
            c = fwrite(&(tkmer), sizeof(uint64_t), 1, bkmer); //writing placeholder for count
            if (c != 1) {
                printf("incorrect number of elements written2\n");
            }
            offset += sizeof(uint64_t);
            
            c = fwrite(&(ref_score_offset), sizeof(uint64_t), 1, bkmer); //writing offset to refs / scores
            if (c != 1) {
                printf("incorrect number of elements written3\n");
            }
            offset += sizeof(uint64_t);
        }
        
        mi[i] = (struct min_kmer_db *)malloc(sizeof(struct min_kmer_db));
        mi[i]->kmer_count = offset_index;
        mi[i]->offset = offset;
        
        c = fwrite(&i, sizeof(uint64_t), 1, mini); //writing minimizer
        if (c != 1) {
            printf("incorrect number of elements written5\n");
        }
        
        c = fwrite(&offset_index, sizeof(uint64_t), 1, mini); //writing number of kmers
        if (c != 1) {
            printf("incorrect number of elements written6\n");
        }
        
        c = fwrite(&offset, sizeof(uint64_t), 1, mini); //writing offset to start of kmer offsets
        if (c != 1) {
            printf("incorrect number of elements written7\n");
        }
        
        c = fwrite(offsets, sizeof(uint64_t), offset_index, bkmer); //writing kmer offsets
        if (c != offset_index) {
            printf("incorrect number of elements written8\n");
        }
        offset += sizeof(uint64_t)*offset_index;
        
        free(offsets);
        offsets = NULL;
    }
    
    fclose(mini);
    fclose(bkmer);
    memset(minimizer_kcount, 0, sizeof(int)*(uint64_t)pow(4, tax_globals.ml));
    load_kmer_db(tbkmer);
    populate_kmer_counts(kmer_db, minimizer_kcount, db_prefix);
    
    //free memory held by kmers and minimizers
    free(minimizer_kcount);
    minimizer_kcount = NULL;
    
}

void build_binner_index(char * db_prefix, char * kmer_db){
    FILE * infile = fopen(kmer_db, "r");
    if (! infile) {
        printf("file %s cannot be read\n", kmer_db);
        return;
    }
    
    char line[100];
    int * minimizer_kcount = calloc((uint64_t)pow(4,tax_globals.ml), sizeof(int));
    uint64_t kmer, kmer_count;
    uint64_t unique_kmer_counts = 0;
    uint64_t total_kmer_counts = 0;
    while(fgets(line, 100, infile)){
        unique_kmer_counts++;
        process_kanalyze_line(line, &kmer, &kmer_count);
        total_kmer_counts += kmer_count;
        uint64_t minimizer = kmer_minimizer(kmer);
        minimizer_kcount[minimizer]++;
    }
    fclose(infile);
    
    tax_globals.total_kmer_count = total_kmer_counts;
    tax_globals.unique_kmer_count = unique_kmer_counts;
    
    char tmini[100];
    memset(tmini, '\0', 100);
    strcpy(tmini, db_prefix);
    strcat(tmini, ".bmi");
    
    char tbkmer[100];
    memset(tbkmer, '\0', 100);
    strcpy(tbkmer, db_prefix);
    strcat(tbkmer, ".btbi");
    
    FILE * mini = fopen( tmini ,"wb");
    FILE * bkmer = fopen( tbkmer , "wb");
    
    //write head to bkmer file
    fwrite(&tax_globals.kl, sizeof(uint8_t), 1, bkmer);
    fwrite(&tax_globals.ml, sizeof(uint8_t), 1, bkmer);
    fwrite(&unique_kmer_counts, sizeof(uint64_t), 1, bkmer);
    fwrite(&total_kmer_counts, sizeof(uint64_t), 1, bkmer);
    fwrite(&tax_globals.kmer_cutoff, sizeof(int), 1, bkmer);
    fwrite(&tax_globals.tax_cutoff, sizeof(int), 1, bkmer);
    
    uint64_t offset = 2*sizeof(uint8_t) + 2*sizeof(uint64_t) + 2*sizeof(int);
    uint64_t ref_score_offset = 0;
    uint64_t total_minimizers = (uint64_t)pow(4, tax_globals.ml);
    mi = (struct min_kmer_db **)calloc((uint64_t)pow(4, tax_globals.ml), sizeof(struct min_kmer_db *));
    uint64_t tkmer = 0;
    
    //write basic organization to both minimizer and kmer binary files
    uint64_t i = 0;
    for (i=0; i < total_minimizers; i++) {
        if (minimizer_kcount[i] == 0) {
            continue;
        }
        uint64_t * offsets = (uint64_t *)malloc(sizeof(uint64_t)*minimizer_kcount[i]);
        uint64_t offset_index = 0;
        uLong c = 0;
        
        int j=0;
        for (j = 0; j < minimizer_kcount[i]; j++) {
            
            offsets[offset_index] = offset;
            offset_index += 1;
            
            c = fwrite(&(tkmer), sizeof(uint64_t), 1, bkmer); //writing kmer
            if (c != 1) {
                printf("incorrect number of elements written1\n");
            }
            offset += sizeof(uint64_t);
            
            c = fwrite(&(tkmer), sizeof(uint64_t), 1, bkmer); //writing placeholder for count
            if (c != 1) {
                printf("incorrect number of elements written2\n");
            }
            offset += sizeof(uint64_t);
            
            /*
            c = fwrite(&(ref_score_offset), sizeof(uint64_t), 1, bkmer); //writing offset to refs / scores
            if (c != 1) {
                printf("incorrect number of elements written3\n");
            }
            offset += sizeof(uint64_t);
            */
        }
        
        mi[i] = (struct min_kmer_db *)malloc(sizeof(struct min_kmer_db));
        mi[i]->kmer_count = offset_index;
        mi[i]->offset = offset;
        
        c = fwrite(&i, sizeof(uint64_t), 1, mini); //writing minimizer
        if (c != 1) {
            printf("incorrect number of elements written5\n");
        }
        
        c = fwrite(&offset_index, sizeof(uint64_t), 1, mini); //writing number of kmers
        if (c != 1) {
            printf("incorrect number of elements written6\n");
        }
        
        c = fwrite(&offset, sizeof(uint64_t), 1, mini); //writing offset to start of kmer offsets
        if (c != 1) {
            printf("incorrect number of elements written7\n");
        }
        
        c = fwrite(offsets, sizeof(uint64_t), offset_index, bkmer); //writing kmer offsets
        if (c != offset_index) {
            printf("incorrect number of elements written8\n");
        }
        offset += sizeof(uint64_t)*offset_index;
        
        free(offsets);
        offsets = NULL;
    }
    
    fclose(mini);
    fclose(bkmer);
    memset(minimizer_kcount, 0, sizeof(int)*(uint64_t)pow(4, tax_globals.ml));
    load_kmer_db(tbkmer);
    populate_binner_kmer_counts(kmer_db, minimizer_kcount, db_prefix);
    
    //free memory held by kmers and minimizers
    free(minimizer_kcount);
    minimizer_kcount = NULL;
    
}

static void populate_binner_kmer_counts(char * kmer_db, int * minimizer_kcount, char * db_prefix){
    FILE * infile = fopen(kmer_db, "r");
    if (! infile) {
        printf("file %s cannot be read\n", kmer_db);
        return;
    }
    
    struct min_kmer_db * tmk = NULL;
    char line[100];
    uint64_t kmer, kmer_count;
    while(fgets(line, 100, infile)){
        process_kanalyze_line(line, &kmer, &kmer_count);
        uint64_t minimizer = kmer_minimizer(kmer);
        tmk = mi[minimizer];
        if (tmk != NULL) {
            uint64_t koffset = 0;
            int kmer_num = minimizer_kcount[minimizer];
            memcpy(&koffset, mmkmer+tmk->offset+sizeof(uint64_t)*kmer_num, sizeof(uint64_t));
            memcpy(mmkmer+koffset, &kmer, sizeof(uint64_t));
            memcpy(mmkmer+koffset+sizeof(uint64_t), &kmer_count, sizeof(uint64_t));
            minimizer_kcount[minimizer]++;
        }
        else {
            printf("not going to use kmer line: %s\n",line);
        }
    }
    //printf("%llu,%llu\n", kmer,kmer_count); //last line in kmer file
    //printf("bytes written to rsi file: %llu\n", rsi_offset);
    fclose(infile);
}


static void populate_kmer_counts(char * kmer_db, int * minimizer_kcount, char * db_prefix){
    FILE * infile = fopen(kmer_db, "r");
    if (! infile) {
        printf("file %s cannot be read\n", kmer_db);
        return;
    }
    
    //start rsi file
    //creating scaffold here for rsi file to avoid reading through kmers again
    char trs[100];
    memset(trs, '\0', 100);
    strcpy(trs, db_prefix);
    strcat(trs, ".rsi");
    FILE * brefscore = fopen(trs, "wb");
    
    uint64_t rsi_offset = 0;
    uint64_t tzero = 0;
    uint64_t trefs[tax_globals.tax_cutoff];
    memset(trefs, 1, sizeof(uint64_t)*tax_globals.tax_cutoff);
    float tscores[tax_globals.tax_cutoff];
    memset(tscores, 1, sizeof(float)*tax_globals.tax_cutoff);
    
    struct min_kmer_db * tmk = NULL;
    char line[100];
    uint64_t kmer, kmer_count;
    while(fgets(line, 100, infile)){
        process_kanalyze_line(line, &kmer, &kmer_count);
        uint64_t minimizer = kmer_minimizer(kmer);
        tmk = mi[minimizer];
        if (tmk != NULL) {
            uint64_t koffset = 0;
            int kmer_num = minimizer_kcount[minimizer];
            memcpy(&koffset, mmkmer+tmk->offset+sizeof(uint64_t)*kmer_num, sizeof(uint64_t));
            memcpy(mmkmer+koffset, &kmer, sizeof(uint64_t));
            memcpy(mmkmer+koffset+sizeof(uint64_t), &kmer_count, sizeof(uint64_t));
            memcpy(mmkmer+koffset+sizeof(uint64_t)*2, &rsi_offset, sizeof(uint64_t));
            minimizer_kcount[minimizer]++;
            
            uLong c = fwrite(&tzero, sizeof(uint64_t), 1, brefscore); //writing placeholder for ref count
            if (c != 1) {
                printf("incorrect number of elements written2\n");
            }
            rsi_offset += sizeof(uint64_t);
            
            uint64_t num_to_write = kmer_count < tax_globals.tax_cutoff ? kmer_count : tax_globals.tax_cutoff;
            
            c = fwrite(trefs, sizeof(uint64_t), num_to_write, brefscore); //writing placeholders for refs
            if (c != num_to_write) {
                printf("incorrect number of elements written2\n");
            }
            rsi_offset += sizeof(uint64_t)*num_to_write;
            
            c = fwrite(tscores, sizeof(float), num_to_write, brefscore); //writing placeholders for scores
            if (c != num_to_write) {
                printf("incorrect number of elements written2\n");
            }
            rsi_offset += sizeof(float)*num_to_write;
        }
        else {
            printf("not going to use kmer line: %s\n",line);
        }
    }
    //printf("%llu,%llu\n", kmer,kmer_count); //last line in kmer file
    //printf("bytes written to rsi file: %llu\n", rsi_offset);
    fclose(infile);
    fclose(brefscore);
}

void populate_rsi(char * fasta_file, char * db_prefix){
    char trs[100];
    memset(trs, '\0', 100);
    strcpy(trs, db_prefix);
    strcat(trs, ".rsi");
    
    load_rsi_db(trs);
    
    gzFile seq_index_file = gzopen(fasta_file, "r");
    kseq_t * seq_t = kseq_init(seq_index_file);
    int l = 0;
    struct seqid_taxid_single * stsr;
    int ref_count = 0;
    
    while ((l = kseq_read(seq_t)) >= 0) {
        //find reference, if not, print error and continue
        ref_count++;
        if (ref_count % 1000 == 0) {
            printf("references processed: %d\n",ref_count);
        }
        
        HASH_FIND_STR(seqid_taxid_rel, seq_t->name.s, stsr);
        if (stsr == NULL) {
            printf("seqid %s has no associated taxid, will skip reference\n", seq_t->name.s);
            continue;
        }
        uint64_t taxid = stsr->taxid;
        khash_t(m64) * kc = kh_init(m64);
        uint64_t ref_kmer_count = 0;
        get_kmer_counts(kc, seq_t->seq.s, seq_t->seq.l, &ref_kmer_count);
        add_refs_rsi(kc, taxid, ref_kmer_count);
        kh_destroy_m64(kc);
    }
    gzclose(seq_index_file);
    kseq_destroy(seq_t);
}

static void get_kmer_counts(khash_t(m64) * kc, char * ref_seq, uint64_t ref_length, uint64_t * ref_kmer_count){
    uint64_t kmer_mask = ~0;
    kmer_mask >>= sizeof(kmer_mask) * 8 - (tax_globals.kl * 2);
    int i =0;
    khint_t k;
    int is_missing, ret;
    while (i + tax_globals.kl <= ref_length) {
        uint64_t kmer = 0;
        int flag = tax_globals.protein == 0 ? next_kmer(ref_seq + i, &kmer) : next_protein_kmer(ref_seq + i, &kmer);
        if (flag == tax_globals.KMER_FLAG) {
            (*ref_kmer_count)++;
            kmer &= kmer_mask;
            uint64_t canonical_kmer = tax_globals.protein == 0 ? canonical_representation2(kmer) : kmer;
            //uint64_t canonical_kmer = canonical_representation2(kmer);
            k = kh_get(m64, kc, canonical_kmer);
            is_missing = (k == kh_end(kc));
            if (is_missing) {
                k = kh_put(m64, kc, canonical_kmer, &ret);
                kh_value(kc, k) = 1;
            }
            else{
                kc->vals[k]++;
            }
            
            i++;
        }
        else {
            i += (flag+1);
        }
    }
    
}

static void get_kco(uint64_t kmer, uint64_t * info){
    uint64_t minimizer = kmer_minimizer(kmer);
    struct min_kmer_db * tmk = mi[minimizer];
    if (tmk == NULL) {
        return;
    }
    uint64_t * offsets = (uint64_t *)malloc(sizeof(uint64_t)*tmk->kmer_count);
    memcpy(offsets, mmkmer+tmk->offset, tmk->kmer_count*sizeof(uint64_t));
    uint64_t pos = kmer_binary_search(offsets, tmk->kmer_count, kmer);
    free(offsets);
    if (pos == tax_globals.KMER_NOT_FOUND) {
        return;
    }
    
    memcpy(info, mmkmer+pos+sizeof(uint64_t), sizeof(uint64_t));
    memcpy(info+1, mmkmer+pos+sizeof(uint64_t)+sizeof(uint64_t), sizeof(uint64_t));
}

static void add_refs_rsi(khash_t(m64) * kc, uint64_t taxid, uint64_t ref_kmer_count){
    uint64_t info[2];
    uint64_t total_kmer_count = tax_globals.total_kmer_count;
    khint_t k;
    for (k = kh_begin(kc); k != kh_end(kc); k++) {
        if (kh_exist(kc, k)) {
            memset(info, 0, sizeof(uint64_t)*2);
            get_kco(kc->keys[k], info);
            if (info[0] <= 0) {
                continue;
            }
            float llp = log((float)kc->vals[k] / (float)info[0]) - log((float)info[0] / (float)(total_kmer_count));
            uint64_t num_refs = 0;
            uint64_t this_cutoff = info[0] < tax_globals.tax_cutoff ? info[0] : tax_globals.tax_cutoff;
            //printf("%llu\t%llu\t%lu\n",info[0],info[1],kc->keys[k]);
            memcpy(&num_refs, mmrsi + info[1], sizeof(uint64_t));
            if (num_refs < this_cutoff) {
                memcpy(mmrsi+info[1]+sizeof(uint64_t)+sizeof(uint64_t)*num_refs, &taxid, sizeof(uint64_t));
                memcpy(mmrsi+info[1]+sizeof(uint64_t)+sizeof(uint64_t)*this_cutoff+sizeof(float)*num_refs, &llp, sizeof(float));
                num_refs++;
                memcpy(mmrsi + info[1], &num_refs, sizeof(uint64_t));
            }

        }
    }
    
}


