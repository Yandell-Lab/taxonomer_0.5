//
//  classify.c
//  taxonomer
//
//  Created by Steven Flygare on 10/6/14.
//  Copyright (c) 2014 yandell lab. All rights reserved.
//

#include "classify.h"

#define READS_PER_BLOCK 10000
#define READ_OUTPUT_SIZE 200

KSEQ_INIT(gzFile, gzread)

char rev_base[] = {'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'T', 'N', 'G', 'N', 'N', 'N', 'C', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'A', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 't', 'N', 'g', 'N', 'N', 'N', 'c', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'a', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N'};


void classify_read_output_ties(struct read_info ** reads, int num_reads, char ** tied_output){
    int i= 0;
#pragma omp parallel for
    for (i=0; i < num_reads; i++) {
        struct read_tax_info * rti = NULL;
        process_read(reads[i]->read, reads[i]->read_len, &rti);
        
        if (tax_globals.protein == 1) {
            //get reverse complement
            uint64_t rl = reads[i]->read_len;
            char rev_seq[rl+1];
            memset(rev_seq, '\0', rl+1);
            int rev_pos = 0;
            int j = 0;
            for (j = rl - 1; j >= 0; j--) {
                rev_seq[rev_pos] = rev_base[(int)reads[i]->read[j]];
                rev_pos++;
            }
            process_read(rev_seq, rl, &rti);
        }
        
        struct read_output_info roi = {0,0, 0.0, 0.0, 0.0};
        uint64_t * max_taxids = NULL;
        ties_from_taxids(rti, &roi, &max_taxids);
        int output_read_buffer_size = 500;
        int cur_output_size = 0;
        tied_output[i] = (char *)calloc(output_read_buffer_size, sizeof(char));
        int j = 0;
        struct taxid_seqid * tts = NULL;
        if (roi.num_taxids > 0) {
            for (j = 0; j < roi.num_taxids; j++) {
                //printf("C\t%s\t%llu\t%llu\n",reads[i]->name,max_taxids[j],roi.num_taxids);
                char tmp_output[READ_OUTPUT_SIZE];
                HASH_FIND(hh, tax_to_seqs, &max_taxids[j], sizeof(uint64_t), tts);
                if (tts == NULL) {
                    printf("problem, taxid %llu not found in structure that contains corresponding seqids, will skip\n",max_taxids[j]);
                    continue;
                }
                
                int tmp_output_size = sprintf(tmp_output, "C\t%s\t%llu\t%s\t%d\t%llu\t%f\t%f\t%llu\n",reads[i]->name,max_taxids[j],tts->seqids[rand() % tts->num_seqids], roi.tax_depth, roi.num_taxids, roi.max_score, roi.mean_score, reads[i]->read_len);
                if (j == 0) {
                    strcpy(tied_output[i], tmp_output);
                }
                else {
                    strcat(tied_output[i], tmp_output);
                }
                cur_output_size += tmp_output_size;
                cur_output_size++;
                if (output_read_buffer_size - cur_output_size <= READ_OUTPUT_SIZE) {
                    output_read_buffer_size = 2 * output_read_buffer_size;
                    tied_output[i] = (char *)realloc(tied_output[i], sizeof(char)*output_read_buffer_size);
                }
            }
        }
        else {
            sprintf(tied_output[i], "U\t%s\n",reads[i]->name);
        }
        
        free(max_taxids);
        struct read_tax_info *s, *tmp;
        HASH_ITER(hh, rti, s, tmp){ //need to free allocated structs
            HASH_DEL(rti, s);
            free(s);
        }
    }
}

char * classify_read_stream(char * read, char * read_name, uint64_t read_length){
    struct read_tax_info * rti = NULL;
    process_read(read, read_length, &rti);
    
    if (tax_globals.protein == 1) {
        //get reverse complement
        uint64_t rl = read_length;
        char rev_seq[rl+1];
        memset(rev_seq, '\0', rl+1);
        int rev_pos = 0;
        int j = 0;
        for (j = rl - 1; j >= 0; j--) {
            rev_seq[rev_pos] = rev_base[(int)read[j]];
            rev_pos++;
        }
        process_read(rev_seq, rl, &rti);
    }
    
    
    struct taxid_confidence * ttc = NULL;
    struct read_tax_info *s, *tmp;
    /*HASH_ITER(hh, rti, s, tmp){ //apply weight
     HASH_FIND(hh, tax_c_scores, &s->taxid, sizeof(uint64_t), ttc);
     if (ttc == NULL) {
     continue;
     }
     else {
     s->score = s->score * ttc->confidence;
     }
     }*/
    char * output = calloc(READ_OUTPUT_SIZE,sizeof(char));
    //char output[READ_OUTPUT_SIZE];
   //memset(output,'\0', sizeof(char)*READ_OUTPUT_SIZE);
    struct read_output_info roi = {0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0};
    uint64_t * max_taxids = NULL;
    lca_from_taxids(rti, &roi, &max_taxids);
    //lca_from_taxids_kw(rti, &roi, &max_taxids);
    struct taxid_seqid * tts = NULL;
    if (roi.lca > 0) {
        HASH_FIND(hh, tax_to_seqs, &max_taxids[0], sizeof(uint64_t), tts);
        if (tts == NULL) {
            printf("problem, taxid %llu not found in structure that contains corresponding seqids, will skip\n",max_taxids[0]);
            sprintf(output, "U\t%s\n",read_name);
        }
        else {
            sprintf(output, "C\t%s\t%llu\t%s\t%d\t%llu\t%f\t%f\t%llu\n",read_name,roi.lca,tts->seqids[rand() % tts->num_seqids], roi.tax_depth, roi.num_taxids, roi.max_score, roi.mean_score,read_length);
        }
    }
    else {
        sprintf(output, "U\t%s\n",read_name);
    }
    
    //fprintf(stdout, "%s", output);
    
    free(max_taxids);
    HASH_ITER(hh, rti, s, tmp){ //need to free allocated structs
        HASH_DEL(rti, s);
        free(s);
    }
    
    return output;

}

void add_reads(char * read, char * read_name, uint64_t read_length, int name_length, int classify){
    static struct read_info * reads[READS_PER_BLOCK];
    static int num_reads = 0;
    
    //printf("%s,%llu,%d,%d\n",read_name,read_length,name_length,num_reads);
    
    if (num_reads >= READS_PER_BLOCK || (classify == 1 && num_reads > 0)) { //enable classification if no more reads and not larger than READS_PER_BLOCK
        int i = 0;
#ifdef _OPENMP
        omp_set_num_threads(tax_globals.np);
#endif

#pragma omp parallel for
        for (i=0; i < num_reads; i++) {
            classify_read_stream(reads[i]->read, reads[i]->name, reads[i]->read_len);
            free(reads[i]->read);
            free(reads[i]->name);
            free(reads[i]);
        }
        
        num_reads = 0;
    }
    
    if (classify == 1) { //return without adding reads
        return;
    }
    
    reads[num_reads] = (struct read_info *)malloc(sizeof(struct read_info));
    reads[num_reads]->read = (char *)calloc(read_length+1,sizeof(char));
    memcpy(reads[num_reads]->read, read, read_length*sizeof(char));
    reads[num_reads]->name = (char *)calloc(name_length+1,sizeof(char));
    memcpy(reads[num_reads]->name, read_name, name_length*sizeof(char));
    reads[num_reads]->read_len = read_length;
    num_reads++;
    
    
}

void classify_read_chunk(struct read_info ** reads, int num_reads, char output[READS_PER_BLOCK][READ_OUTPUT_SIZE]){
    int i= 0;
    #pragma omp parallel for
    for (i=0; i < num_reads; i++) {
        struct read_tax_info * rti = NULL;
        //printf("%s\t%llu\n",reads[i]->name,reads[i]->read_len);
        process_read(reads[i]->read, reads[i]->read_len, &rti);
        
        if (tax_globals.protein == 1) {
            //get reverse complement
            uint64_t rl = reads[i]->read_len;
            char rev_seq[rl+1];
            memset(rev_seq, '\0', rl+1);
            int rev_pos = 0;
            int j = 0;
            for (j = rl - 1; j >= 0; j--) {
                rev_seq[rev_pos] = rev_base[(int)reads[i]->read[j]];
                rev_pos++;
            }
            process_read(rev_seq, rl, &rti);
        }
        
        
        struct taxid_confidence * ttc = NULL;
        struct read_tax_info *s, *tmp;
        /*HASH_ITER(hh, rti, s, tmp){ //apply weight
            HASH_FIND(hh, tax_c_scores, &s->taxid, sizeof(uint64_t), ttc);
            if (ttc == NULL) {
                continue;
            }
            else {
                s->score = s->score * ttc->confidence;
            }
        }*/
        struct read_output_info roi = {0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0};
        uint64_t * max_taxids = NULL;
        lca_from_taxids(rti, &roi, &max_taxids);
        //lca_from_taxids_kw(rti, &roi, &max_taxids);
        memset(output[i],'\0', sizeof(char)*READ_OUTPUT_SIZE);
        struct taxid_seqid * tts = NULL;
        if (roi.lca > 0) {
            HASH_FIND(hh, tax_to_seqs, &max_taxids[0], sizeof(uint64_t), tts);
            if (tts == NULL) {
                printf("problem, taxid %llu not found in structure that contains corresponding seqids, will skip\n",max_taxids[0]);
                continue;
            }
            sprintf(output[i], "C\t%s\t%llu\t%s\t%d\t%llu\t%f\t%f\t%llu\n",reads[i]->name,roi.lca,tts->seqids[rand() % tts->num_seqids], roi.tax_depth, roi.num_taxids, roi.max_score, roi.mean_score, reads[i]->read_len);
        }
        else {
            sprintf(output[i], "U\t%s\n",reads[i]->name);
        }
        
        free(max_taxids);
        HASH_ITER(hh, rti, s, tmp){ //need to free allocated structs
            HASH_DEL(rti, s);
            free(s);
        }
    }
}

void classify_reads(char * read_file, int np, char * output_file){
    srand (time(NULL));
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
            if (tax_globals.output_ties == 0) {
                char output[READS_PER_BLOCK][READ_OUTPUT_SIZE];
                classify_read_chunk(reads, num_reads, output);
                int i = 0;
                for (i=0; i < num_reads; i++) {
                    fprintf(out, "%s", output[i]);
                    free(reads[i]->read);
                    free(reads[i]->name);
                    free(reads[i]);
                }
            }
            else {
                char ** output = (char **)malloc(sizeof(char *)*READS_PER_BLOCK);
                classify_read_output_ties(reads, num_reads, output);
                int i = 0;
                for (i=0; i < num_reads; i++) {
                    fprintf(out, "%s", output[i]);
                    free(output[i]);
                    free(reads[i]->read);
                    free(reads[i]->name);
                    free(reads[i]);
                }
                free(output);
            }
            //free(output);
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
    
    if (tax_globals.output_ties == 0) {
        char output[READS_PER_BLOCK][READ_OUTPUT_SIZE];
        classify_read_chunk(reads, num_reads, output);
        int i = 0;
        for (i=0; i < num_reads; i++) {
            fprintf(out, "%s", output[i]);
            free(reads[i]->read);
            free(reads[i]->name);
            free(reads[i]);
        }
    }
    else {
        char ** output = (char **)malloc(sizeof(char *)*READS_PER_BLOCK);
        classify_read_output_ties(reads, num_reads, output);
        int i = 0;
        for (i=0; i < num_reads; i++) {
            fprintf(out, "%s", output[i]);
            free(output[i]);
            free(reads[i]->read);
            free(reads[i]->name);
            free(reads[i]);
        }
        free(output);
    }
    total_processed += num_reads;
    printf("total reads processed: %llu\n",total_processed);
    
    gzclose(seq_index_file);
    kseq_destroy(seq_t);
    free(reads);
    fclose(out);
}

int64_t score_sort(struct read_tax_info * a, struct read_tax_info * b){
    return a->score - b->score;
}

float score_sort_kw(struct read_tax_info * a, struct read_tax_info * b){
    return a->score/(float)a->counts - b->score/(float)b->counts;
}

struct read_tax_scores * classify_read(char * read, uint64_t read_len){
    //<classify> <db_prefix> <tri> <sti> <read>
    srand (time(NULL));
    struct read_tax_info * rti = NULL;
    process_read(read, read_len, &rti);
    
    if (tax_globals.protein == 1) {
        //get reverse complement
        uint64_t rl = read_len;
        char rev_seq[rl+1];
        memset(rev_seq, '\0', rl+1);
        int rev_pos = 0;
        int j = 0;
        for (j = rl - 1; j >= 0; j--) {
            rev_seq[rev_pos] = rev_base[(int)read[j]];
            rev_pos++;
        }
        process_read(rev_seq, rl, &rti);
    }

    struct read_output_info roi = {0,0, 0.0, 0.0, 0.0};
    uint64_t * max_taxids = NULL;
    lca_from_taxids(rti, &roi, &max_taxids);
    //free(max_taxids);
    struct read_tax_scores * rts = (struct read_tax_scores *)malloc(sizeof(struct read_tax_scores));
    rts->lca = (uint64_t *)calloc(1, sizeof(uint64_t));
    rts->ntaxids = (int *)calloc(1, sizeof(int));
    rts->seqid  = (char *)calloc(100, sizeof(char));
    memset(rts->seqid, '\0', sizeof(char)*100);
    *rts->lca = roi.lca;
    struct taxid_seqid * tts = NULL;
    if (roi.lca > 0) {
        HASH_FIND(hh, tax_to_seqs, &max_taxids[0], sizeof(uint64_t), tts);
        if (tts == NULL) {
            printf("problem, taxid %llu not found in structure that contains corresponding seqids, will skip\n",max_taxids[0]);
        }
        else{
            strcpy(rts->seqid, tts->seqids[rand() % tts->num_seqids]);
        }
    }
    free(max_taxids);
    *rts->ntaxids = HASH_COUNT(rti);
    rts->taxids = (uint64_t *)malloc(sizeof(uint64_t)*(*rts->ntaxids));
    rts->scores = (float *)malloc(sizeof(float)*(*rts->ntaxids));
    rts->counts = (int *)malloc(sizeof(int)*(*rts->ntaxids));
    struct read_tax_info *s, *tmp;
    HASH_SORT(rti, score_sort);
    //HASH_SORT(rti, score_sort_kw);
    int i = 0;
    HASH_ITER(hh, rti, s, tmp){ //need to free allocated structs
        //printf("%llu\t%f\t%d\n",s->taxid,s->score,s->counts);
        rts->taxids[i] = s->taxid;
        rts->scores[i] = s->score;
        rts->counts[i] = s->counts;
        
        i++;
        HASH_DEL(rti, s);
        free(s);
    }
    return rts;
}

void classify_db(char * db_file, char * sti_file, int read_len, char * output_file, int np){
    
#ifdef _OPENMP
    omp_set_num_threads(np);
#endif
    
    printf("loading .sti file...");
    load_seqid_taxid_rel(sti_file);
    printf("done\n");
    
    FILE * out = fopen(output_file, "w");
    
    gzFile db_seq = gzopen(db_file, "r");
    kseq_t * seq_t = kseq_init(db_seq);
    int l = 0;
    struct seqid_taxid_single * stsr;
    
    while ((l = kseq_read(seq_t)) >= 0) {
        HASH_FIND_STR(seqid_taxid_rel, seq_t->name.s, stsr);
        if (stsr == NULL) {
            printf("seqid %s has no associated taxid, will skip reference\n", seq_t->name.s);
            continue;
        }
        uint64_t taxid = stsr->taxid;
        uint64_t i = 0;
        int correct = 0;
        //#pragma omp parallel for
        for (i=0; i + read_len <= seq_t->seq.l; i++) {
            char read[read_len+1];
            read[read_len] = '\0';
            memcpy(read, (seq_t->seq.s)+i, read_len);
            //printf("read: %s\n",read);
            struct read_tax_info * rti = NULL;
            process_read(read, read_len, &rti);
            struct read_output_info roi = {0,0, 0.0, 0.0, 0.0};
            uint64_t * max_taxids = NULL;
            lca_from_taxids(rti, &roi, &max_taxids);
            free(max_taxids);
            int tmp_correct = roi.lca == taxid ? 1 : 0;
//#pragma omp atomic
            correct += tmp_correct;
            struct read_tax_info *s, *tmp;
            HASH_ITER(hh, rti, s, tmp){ //need to free allocated structs
                HASH_DEL(rti, s);
                free(s);
            }
        }
        
        fprintf(out, "%s\t%llu\t%f\n",seq_t->name.s,taxid,(float)correct/(float)i );
        
    }
    gzclose(db_seq);
    fclose(out);
}


uint64_t * get_kmer_info(char * kmer){
    uint64_t * res = (uint64_t *)calloc(4, sizeof(uint64_t));
    //uint64_t * info = (uint64_t *)calloc(2, sizeof(uint64_t));
    uint64_t info[2];
    memset(info, '\0', 2*sizeof(uint64_t));
    uint64_t k = 0;
    int ret = tax_globals.protein == 0 ? next_kmer(kmer, &k) : next_protein_kmer(kmer, &k);
    if (ret == tax_globals.KMER_FLAG) {
        query_kcrn(k, info);
        res[0] = info[0];
        res[1] = info[1];
    }
    if (tax_globals.protein == 1) {
        //get reverse complement
        uint64_t rl = tax_globals.kl;
        char rev_seq[rl+1];
        memset(rev_seq, '\0', rl+1);
        int rev_pos = 0;
        int j = 0;
        for (j = rl - 1; j >= 0; j--) {
            rev_seq[rev_pos] = rev_base[(int)kmer[j]];
            rev_pos++;
        }
        int ret = next_protein_kmer(rev_seq, &k);
        if (ret == tax_globals.KMER_FLAG) {
            memset(info, '\0', 2*sizeof(uint64_t));
            query_kcrn(k, info);
            res[2] = info[0];
            res[3] = info[1];
        }
        
    }
    
    return res;
}

void load_classification_dbs(char * mi_file, char * db_file, char * rsi_file, char * tri_file, char * sti_file,int load_mem, int ties, int protein, int afterburner){
    set_tax_globals();
    tax_globals.output_ties = ties;
    tax_globals.protein = protein;
    tax_globals.afterburner = afterburner;
    tax_globals.load_mem = load_mem;
    printf("loading kmer index...\n");
    load_kmer_db(db_file);
    printf("done loading kmer index\n");
    printf("loading minimizer index...");
    load_minimizer_index(mi_file);
    printf("done\n");
    printf("loading rsi db...");
    load_rsi_db(rsi_file);
    printf("done\n");
    printf("loading taxid relationships...");
    load_tri_file(tri_file);
    printf("done\n");
    printf("loading seqid taxid relationships...");
    load_seqid_taxid_rel(sti_file);
    printf("done\n");
    if (protein == 1 && tax_globals.kl%3 != 0) {
        printf("kmer length MUST be a multiple of 3 to classify in protein space\ndatabase must be built specifically for this, going to exit");
        exit(0);
    }
    
    //printf("loading taxid scores...");
    //load_taxid_confidence(taxid_confidence);
    //printf("done\n");
}

int get_kmer_length(){
    return tax_globals.kl;
}

uint64_t get_total_kmer_count(){
    return tax_globals.total_kmer_count;
}

uint64_t get_unique_kmer_count(){
    return tax_globals.unique_kmer_count;
}

void set_kmer_cutoff(int kc){
    tax_globals.kmer_cutoff = kc;
}

void set_number_processors(int np){
    tax_globals.np = np;
}



