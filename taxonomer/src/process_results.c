//
//  process_results.c
//  taxonomer
//
//  Created by steven on 10/7/14.
//  Copyright (c) 2014 yandell lab. All rights reserved.
//

#include "process_results.h"

KSEQ_INIT(gzFile, gzread)
KBTREE_INIT(path, uint64_t, kb_generic_cmp)
KHASH_MAP_INIT_INT64(m64_2, uint64_t)
khash_t(m64_2) * taxid_rel;

void load_tri_file(char * tri_file){
    
    taxid_rel = kh_init(m64_2);
    
    gzFile seq_index_file = gzopen(tri_file, "r");
    kseq_t * seq_t = kseq_init(seq_index_file);
    int l = 0;
    khint_t k;
    int is_missing,ret;
    while ((l = kseq_read(seq_t)) >= 0) {
        uint64_t child = strtoull(seq_t->name.s, NULL, 10);
        uint64_t parent = strtoull(seq_t->seq.s, NULL, 10);
        k = kh_get(m64_2, taxid_rel, child);
        is_missing = (k == kh_end(taxid_rel));
        if (is_missing) {
            k = kh_put(m64_2, taxid_rel, child, &ret);
            kh_value(taxid_rel, k) = parent;
        }
        else{
            printf("problem! taxid has more than one parent: %llu\n",child);
        }
    }
    gzclose(seq_index_file);
    kseq_destroy(seq_t);
}

static float score_sample_variance(struct read_tax_info * taxid_results, float mean){
    float variance = 0.0;
    if (HASH_COUNT(taxid_results) < 2){
        return variance;
    }
    struct read_tax_info *s, *tmp;
    float sum_squares = 0.0;
    HASH_ITER(hh, taxid_results, s, tmp) {
        sum_squares += ((s->score - mean) * (s->score - mean));
    }
    variance = (1.0 / (float)(HASH_COUNT(taxid_results)-1)) * sum_squares;
    if (fabs(variance) < .1) {
        return 0.0;
    }
    return sqrtf(variance);
}

static int taxid_depth(uint64_t taxid){
    int depth = 0;
    uint64_t tmp_taxid = taxid;
    khint_t k;
    while (tmp_taxid > 0) {
        depth++;
        k = kh_get(m64_2, taxid_rel, tmp_taxid);
        tmp_taxid = taxid_rel->vals[k];
    }
    return depth;
}

void lca_from_taxids(struct read_tax_info * taxid_results, struct read_output_info * roi, uint64_t ** max_taxids){
    struct read_tax_info *s, *tmp;
    int initial_ties = 500;
    *max_taxids = (uint64_t *)malloc(sizeof(uint64_t)*initial_ties);
    uint64_t size = 0;
    double sum_score = 0;
    double ln_sum_score = 0;
    float max_score = -99999999;
    HASH_ITER(hh, taxid_results, s, tmp) {
        //printf("%llu\t%f\t%d\n",s->taxid,s->score,s->counts);
        sum_score += s->score;
        ln_sum_score += log(s->score);
        if (fabs(s->score - max_score) < .1) {
            if (size >= initial_ties) {
                initial_ties = 2*initial_ties;
                uint64_t * tmp_ptr = (uint64_t *)realloc(*max_taxids, initial_ties*sizeof(uint64_t));
                *max_taxids = tmp_ptr;
                tmp_ptr = NULL;
            }
            (*max_taxids)[size] = s->taxid;
            size++;
        }
        else if (s->score > max_score) {
            max_score = s->score;
            size = 0;
            (*max_taxids)[size] = s->taxid;
            size++;
        }
    }
    roi->mean_score = sum_score / (float)HASH_COUNT(taxid_results);
    //roi->var_score = score_sample_variance(taxid_results, roi->mean_score);
    //roi->var_score = 0.0;
    roi->max_score = max_score;
    roi->num_taxids = size;
    roi->lca = get_lca(*max_taxids, size);
    roi->tax_depth = taxid_depth(roi->lca);
    //calculate gamma parameters
    //float tmps = log( (1.0 / (float)HASH_COUNT(taxid_results)) * sum_score ) - (1.0 / (float)HASH_COUNT(taxid_results)) * ln_sum_score;
    //roi->kappa = (3 - tmps + sqrt( (tmps-3)*(tmps-3) + 24*tmps) ) / (12*tmps);
    //roi->theta = sum_score / (roi->kappa * (float)HASH_COUNT(taxid_results));
    //free(max_taxids);
    
}

void lca_from_taxids_kw(struct read_tax_info * taxid_results, struct read_output_info * roi, uint64_t ** max_taxids){
    struct read_tax_info *s, *tmp;
    int initial_ties = 500;
    *max_taxids = (uint64_t *)malloc(sizeof(uint64_t)*initial_ties);
    uint64_t size = 0;
    double sum_score = 0;
    double ln_sum_score = 0;
    float max_score = -99999999;
    HASH_ITER(hh, taxid_results, s, tmp) {
        //printf("%llu\t%f\t%d\n",s->taxid,s->score,s->counts);
        double kw = s->score / s->counts;
        sum_score += kw;
        ln_sum_score += log(kw);
        if (fabs(kw - max_score) < .1) {
            if (size >= initial_ties) {
                initial_ties = 2*initial_ties;
                uint64_t * tmp_ptr = (uint64_t *)realloc(*max_taxids, initial_ties*sizeof(uint64_t));
                *max_taxids = tmp_ptr;
                tmp_ptr = NULL;
            }
            (*max_taxids)[size] = s->taxid;
            size++;
        }
        else if (kw > max_score) {
            max_score = kw;
            size = 0;
            (*max_taxids)[size] = s->taxid;
            size++;
        }
    }
    roi->mean_score = sum_score / (float)HASH_COUNT(taxid_results);
    roi->max_score = max_score;
    roi->num_taxids = size;
    roi->lca = get_lca(*max_taxids, size);
    roi->tax_depth = taxid_depth(roi->lca);
    //free(max_taxids);
    
}


void ties_from_taxids(struct read_tax_info * taxid_results, struct read_output_info * roi, uint64_t ** max_taxids){
    struct read_tax_info *s, *tmp;
    int initial_ties = 500;
    *max_taxids = (uint64_t *)malloc(sizeof(uint64_t)*initial_ties);
    uint64_t size = 0;
    double sum_score = 0;
    double ln_sum_score = 0;
    float max_score = -99999999;
    HASH_ITER(hh, taxid_results, s, tmp) {
        sum_score += s->score;
        ln_sum_score += log(s->score);
        if (fabs(s->score - max_score) < .1) {
            if (size >= initial_ties) {
                initial_ties = 2*initial_ties;
                uint64_t * tmp_ptr = (uint64_t *)realloc(*max_taxids, initial_ties*sizeof(uint64_t));
                *max_taxids = tmp_ptr;
                tmp_ptr = NULL;
            }
            (*max_taxids)[size] = s->taxid;
            size++;
        }
        else if (s->score > max_score) {
            max_score = s->score;
            size = 0;
            (*max_taxids)[size] = s->taxid;
            size++;
        }
    }
    roi->mean_score = sum_score / (float)HASH_COUNT(taxid_results);
    //roi->var_score = score_sample_variance(taxid_results, roi->mean_score);
    roi->var_score = 0.0;
    roi->max_score = max_score;
    roi->num_taxids = size;
    roi->lca = get_lca(*max_taxids, size);
    roi->tax_depth = taxid_depth(roi->lca);
}


uint64_t get_lca(uint64_t * max_taxids, uint64_t size){
    if (size == 0) {
        return 0;
    }
    else if (size == 1){
        return max_taxids[0];
    }
    uint64_t lca = max_taxids[0];
    int i = 1;
    for (i = 1; i < size; i++) {
        lca = get_lca2(lca, max_taxids[i]);
    }
    
    return lca;
}

static uint64_t get_lca2(uint64_t a, uint64_t b){
    if (a == 0 || b == 0){
        return a ? a : b;
    }
    
    khint_t k;
    kbtree_t(path) * a_path = kb_init(path, KB_DEFAULT_SIZE);
    while (a > 0) {
        kb_put(path, a_path, a);
        k = kh_get(m64_2, taxid_rel, a);
        a = taxid_rel->vals[k];
    }
    while (b > 0) {
        if (kb_get(path, a_path, b)) {
            kb_destroy(path, a_path);
            return b;
        }
        k = kh_get(m64_2, taxid_rel, b);
        b = taxid_rel->vals[k];
    }
    kb_destroy(path, a_path);
    return 1;
}


