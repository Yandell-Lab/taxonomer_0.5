//
//  load_dbs.c
//  taxonomer
//
//  Created by Steven Flygare on 10/3/14.
//  Copyright (c) 2014 yandell lab. All rights reserved.
//

#include "load_dbs.h"

KSEQ_INIT(gzFile, gzread);

void load_minimizer_index(char * minimizer_file){
    
    mi = (struct min_kmer_db **)calloc((uint64_t)pow(4, tax_globals.ml), sizeof(struct min_kmer_db *));
    uint64_t offset = 0;
    int line_size = 24; //3 * 8bytes
    int fdSrc = open(minimizer_file,O_RDONLY,0);
    struct stat st;
    fstat(fdSrc, &st);
    unsigned char * mmdata = (unsigned char *)mmap(NULL, st.st_size, PROT_READ, MAP_SHARED, fdSrc, 0);
    if(mmdata == MAP_FAILED){
        printf("index map failed\n");
        exit(1);
    }
    uint64_t ldata[3];
    struct min_kmer_db * tmi;
    while (offset <= st.st_size - line_size) {
        memcpy(ldata, mmdata+offset, line_size);
        offset += line_size;
        tmi = (struct min_kmer_db *)malloc(sizeof(struct min_kmer_db));
        tmi->kmer_count = ldata[1];
        tmi->offset = ldata[2];
        mi[ldata[0]] = tmi;
        
    }
    munmap((void *)mmdata, st.st_size);
    mmdata = NULL;
    close(fdSrc);
    
}

void load_kmer_db(char * kmer_db){
    
    int fdSrc = open(kmer_db,O_RDWR,0);
    struct stat st;
    fstat(fdSrc, &st);
    mmkmer = (unsigned char *)mmap(NULL, st.st_size, PROT_READ | PROT_WRITE, MAP_SHARED, fdSrc, 0);
    if(mmkmer == MAP_FAILED){
        printf("byte map failed\n");
        exit(1);
    }
    if (tax_globals.load_mem == 1) {
        size_t filesize = st.st_size;
        size_t page_size = getpagesize();
        unsigned char * buf[page_size];
        size_t pos=0;
        for (pos=0; pos < filesize; pos += page_size) {
            size_t this_page_size = filesize - pos;
            if (this_page_size > page_size){
                this_page_size = page_size;
            }
            memcpy(buf, mmkmer + pos, this_page_size);
        }
    }
    memcpy(&tax_globals.kl, mmkmer, sizeof(uint8_t));
    printf("kmer length: %d\n",tax_globals.kl);
    memcpy(&tax_globals.ml, mmkmer+sizeof(uint8_t), sizeof(uint8_t));
    printf("minimizer length: %d\n",tax_globals.ml);
    memcpy(&tax_globals.unique_kmer_count, mmkmer+sizeof(uint8_t)*2, sizeof(uint64_t));
    printf("unique kmer count: %llu\n",tax_globals.unique_kmer_count);
    memcpy(&tax_globals.total_kmer_count, mmkmer+sizeof(uint8_t)*2+sizeof(uint64_t), sizeof(uint64_t));
    printf("total kmer count: %llu\n",tax_globals.total_kmer_count);
    memcpy(&tax_globals.kmer_cutoff, mmkmer+sizeof(uint8_t)*2+sizeof(uint64_t)*2, sizeof(int));
    printf("kmer cutoff when constructing db: %d\n",tax_globals.kmer_cutoff);
    memcpy(&tax_globals.tax_cutoff, mmkmer+sizeof(uint8_t)*2+sizeof(uint64_t)*2+sizeof(int), sizeof(int));
    printf("tax cutoff when constructing db: %d\n",tax_globals.tax_cutoff);
    
    
}

void load_rsi_db(char * rsi_file){
 
    int fdSrc = open(rsi_file,O_RDWR,0);
    struct stat st;
    fstat(fdSrc, &st);
    mmrsi = (unsigned char *)mmap(NULL, st.st_size, PROT_READ | PROT_WRITE, MAP_SHARED, fdSrc, 0);
    if(mmrsi == MAP_FAILED){
        printf("byte map failed\n");
        exit(1);
    }
    if (tax_globals.load_mem == 1) {
        size_t filesize = st.st_size;
        size_t page_size = getpagesize();
        unsigned char * buf[page_size];
        size_t pos=0;
        for (pos=0; pos < filesize; pos += page_size) {
            size_t this_page_size = filesize - pos;
            if (this_page_size > page_size){
                this_page_size = page_size;
            }
            memcpy(buf, mmrsi + pos, this_page_size);
        }
    }

}

void load_seqid_taxid_rel(char * seqid_taxid_file){
    seqid_taxid_rel = NULL;
    tax_to_seqs = NULL;
    gzFile seq_tax = gzopen(seqid_taxid_file, "r");
    kseq_t * seq_t = kseq_init(seq_tax);
    struct seqid_taxid_single * tstr;
    struct taxid_seqid * tts;
    int l = 0;
    
    while ((l = kseq_read(seq_t)) >= 0) {
        //add taxid related to seqid
        uint64_t tmp_taxid = 0;
        HASH_FIND_STR(seqid_taxid_rel, seq_t->name.s, tstr);
        if (tstr == NULL) {
            if (seq_t->name.l > SEQID_SIZE) {
                printf("seqid %s is too long, must be less than 100 characters.  Exiting...", seq_t->name.s);
                exit(1);
            }
            tstr = (struct seqid_taxid_single *)malloc(sizeof(struct seqid_taxid_single));
            memset(tstr->seqid, '\0', SEQID_SIZE*sizeof(char));
            strncpy(tstr->seqid, seq_t->name.s, seq_t->name.l);
            tmp_taxid = strtoull(seq_t->seq.s, NULL, 10);
            tstr->taxid = tmp_taxid;
            HASH_ADD_STR(seqid_taxid_rel, seqid, tstr);
            
        }
        else {
            printf("%s already seen in hash",seq_t->name.s);
        }
        
        //add seqid(s) related to taxid
        HASH_FIND(hh, tax_to_seqs, &tmp_taxid, sizeof(uint64_t), tts);
        if (tts == NULL) {
            tts = (struct taxid_seqid *)malloc(sizeof(struct taxid_seqid));
            tts->taxid = tmp_taxid;
            tts->num_seqids = 0;
            tts->max_seqs = 10;
            tts->seqids = (char **)malloc(sizeof(char *)*10);
            tts->seqids[tts->num_seqids] = (char *)calloc(SEQID_SIZE,sizeof(char));
            strncpy(tts->seqids[tts->num_seqids], seq_t->name.s, seq_t->name.l);
            tts->num_seqids++;
            HASH_ADD(hh, tax_to_seqs, taxid, sizeof(uint64_t), tts);
        }
        else {
            if (tts->num_seqids >= tts->max_seqs) {
                tts->max_seqs = 2 * tts->max_seqs;
                tts->seqids = (char **)realloc(tts->seqids, sizeof(char *)*tts->max_seqs);
            }
            tts->seqids[tts->num_seqids] = (char *)calloc(SEQID_SIZE,sizeof(char));
            strncpy(tts->seqids[tts->num_seqids], seq_t->name.s, seq_t->name.l);
            tts->num_seqids++;
        }
    }
    gzclose(seq_tax);
}

void load_taxid_confidence(char * taxid_confidence_file){
    tax_c_scores = NULL;
    struct taxid_confidence * ttc = NULL;
    gzFile tcf = gzopen(taxid_confidence_file, "r");
    kseq_t * seq_t = kseq_init(tcf);
    int l = 0;
    while ((l = kseq_read(seq_t)) >= 0) {
        uint64_t taxid = strtoull(seq_t->name.s, NULL, 10);
        uint64_t score = strtof(seq_t->seq.s, NULL);
        HASH_FIND(hh, tax_c_scores, &taxid, sizeof(uint64_t), ttc);
        if (ttc != NULL) {
            printf("problem, taxid %llu shouldn't already have a confidence\n",taxid);
        }
        else {
            ttc = (struct taxid_confidence *)malloc(sizeof(struct taxid_confidence));
            ttc->taxid = taxid;
            ttc->confidence = score;
            HASH_ADD(hh, tax_c_scores, taxid, sizeof(uint64_t), ttc);
        }
    }
    gzclose(tcf);
}




