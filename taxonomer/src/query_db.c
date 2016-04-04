//
//  query_db.c
//  taxonomer
//
//  Created by Steven Flygare on 10/3/14.
//  Copyright (c) 2014 yandell lab. All rights reserved.
//

#include "query_db.h"

uint64_t query_kc(uint64_t kmer){
    uint64_t canonical_kmer = canonical_representation2(kmer);
    //uint64_t canonical_kmer = kmer;
    uint64_t minimizer = kmer_minimizer(canonical_kmer);
    struct min_kmer_db * tmk = mi[minimizer];
    if (tmk == NULL) {
        return 0;
    }
    uint64_t * offsets = (uint64_t *)malloc(sizeof(uint64_t)*tmk->kmer_count);
    memcpy(offsets, mmkmer+tmk->offset, tmk->kmer_count*sizeof(uint64_t));
    uint64_t pos = kmer_binary_search(offsets, tmk->kmer_count, canonical_kmer);
    free(offsets);
    if (pos == tax_globals.KMER_NOT_FOUND) {
        return 0;
    }
    uint64_t kmer_count = 0;
    memcpy(&kmer_count, mmkmer+pos+sizeof(uint64_t), sizeof(uint64_t));
    return kmer_count;
}

uint64_t query_kc2(uint64_t kmer, struct min_kmer_db * tmi){
    uint64_t * offsets = NULL;
    uint64_t offset_len = 0;
    uint64_t canonical_kmer = tax_globals.protein == 0 ? canonical_representation2(kmer) : kmer;
    //uint64_t canonical_kmer = canonical_representation2(kmer);
    uint64_t pos = tax_globals.KMER_NOT_FOUND;
    if (tmi->kmer_count == 0) { //check if minimizer is not set, if not, set it and get kmer byte offsets
        pos = query_minimizer(kmer, canonical_kmer, tmi, &offsets, &offset_len);
    }
    else if (tmi->kmer_count > 0) { //if minimizer is set, do binary search to see if kmer is in range, if not then recompute minimizer and search
        /*if (canonical_kmer >= *(uint64_t *)(mmkmer+offsets[0]) && canonical_kmer <= *(uint64_t *)(mmkmer+offsets[tmi->kmer_count - 1])) {
         pos = kmer_binary_search(offsets,tmi->kmer_count,canonical_kmer);
         }*/
        pos = kmer_binary_search(offsets,tmi->kmer_count,canonical_kmer);
        if (pos == tax_globals.KMER_NOT_FOUND) {
            pos = query_minimizer(kmer, canonical_kmer, tmi, &offsets, &offset_len);
        }
    }
    
    if (pos != tax_globals.KMER_NOT_FOUND) {
        uint64_t * kmer_count = (uint64_t *)(mmkmer+pos+sizeof(uint64_t));
        return *kmer_count;
    }
    else {
        return 0;
    }
}


void query_kco(uint64_t kmer, uint64_t * info){
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
    memcpy(&info[0], mmkmer+pos+sizeof(uint64_t), sizeof(uint64_t));
    memcpy(&info[1], mmkmer+pos+sizeof(uint64_t)*2, sizeof(uint64_t));
}

void query_kcrn(uint64_t kmer, uint64_t * info){
    uint64_t minimizer = kmer_minimizer(kmer);
    struct min_kmer_db * tmk = mi[minimizer];
    if (tmk == NULL) {
        return;
    }
    uint64_t canonical_kmer = tax_globals.protein == 0 ? canonical_representation2(kmer) : kmer;
    uint64_t * offsets = (uint64_t *)malloc(sizeof(uint64_t)*tmk->kmer_count);
    memcpy(offsets, mmkmer+tmk->offset, tmk->kmer_count*sizeof(uint64_t));
    uint64_t pos = kmer_binary_search(offsets, tmk->kmer_count, canonical_kmer);
    free(offsets);
    if (pos == tax_globals.KMER_NOT_FOUND) {
        return;
    }
    uint64_t rsi_offset = 0;
    memcpy(&info[0], mmkmer+pos+sizeof(uint64_t), sizeof(uint64_t));
    memcpy(&rsi_offset, mmkmer+pos+sizeof(uint64_t)*2, sizeof(uint64_t));
    uint64_t * num_refs = (uint64_t *)(mmrsi+(rsi_offset));
    memcpy(&info[1], num_refs, sizeof(uint64_t));
}

uint64_t kmer_binary_search(uint64_t * offsets, uint64_t offset_size, uint64_t kmer){
    int64_t min = 0;
    int64_t max = offset_size - 1;
    while (max >= min) {
        int64_t mid = (min+max)/2;
        if (*(uint64_t*)(mmkmer+offsets[mid]) == kmer) {
            return offsets[mid];
        }
        else if (*(uint64_t*)(mmkmer+offsets[mid]) < kmer){
            min = mid+1;
        }
        else{
            max = mid-1;
        }
    }
    
    return tax_globals.KMER_NOT_FOUND;
}

int process_read2(char * read_seq, uint64_t read_len, struct read_tax_info ** tax_hits){
    int i = 0;
    uint64_t kmer_mask = ~0;
    kmer_mask >>= sizeof(kmer_mask) * 8 - (tax_globals.kl * 2);
    //uint64_t * offsets = (uint64_t *)malloc(sizeof(uint64_t)*KMER_PER_MINIMIZER);
    //uint64_t offset_len = KMER_PER_MINIMIZER;
    uint64_t * offsets = NULL;
    uint64_t offset_len = 0;
    struct min_kmer_db tmi;
    tmi.kmer_count = 0;
    tmi.offset = 0;
    int covered_count = 0;
    int prev_found = -1;
    while (i + tax_globals.kl <= read_len) {
        uint64_t kmer = 0;
        int flag = tax_globals.protein == 0 ? next_kmer(read_seq + i, &kmer) : next_protein_kmer(read_seq + i, &kmer);
        if (flag == tax_globals.KMER_FLAG) {
            kmer &= kmer_mask;
            int hit = query_kmer2(kmer, &tmi, &offsets, &offset_len, tax_hits);
            if (hit == 1) {
                if (prev_found < 0) {
                    covered_count += tax_globals.kl;
                    prev_found = i;
                }
                else {
                    if (i - prev_found > tax_globals.kl) {
                        covered_count += tax_globals.kl;
                    }
                    else {
                        covered_count += i - prev_found;
                    }
                    prev_found = i;
                }
            }
            i++;
        }
        else {
            i += (flag+1);
        }
    }
    return covered_count;
    //free(offsets);
}

void process_read(char * read_seq, uint64_t read_len, struct read_tax_info ** tax_hits){
    int i = 0;
    uint64_t kmer_mask = ~0;
    kmer_mask >>= sizeof(kmer_mask) * 8 - (tax_globals.kl * 2);
    //uint64_t * offsets = (uint64_t *)malloc(sizeof(uint64_t)*KMER_PER_MINIMIZER);
    //uint64_t offset_len = KMER_PER_MINIMIZER;
    uint64_t * offsets = NULL;
    uint64_t offset_len = 0;
    struct min_kmer_db tmi;
    tmi.kmer_count = 0;
    tmi.offset = 0;
    while (i + tax_globals.kl <= read_len) {
        uint64_t kmer = 0;
        int flag = tax_globals.protein == 0 ? next_kmer(read_seq + i, &kmer) : next_protein_kmer(read_seq + i, &kmer);
        if (flag == tax_globals.KMER_FLAG) {
            kmer &= kmer_mask;
            query_kmer(kmer, &tmi, &offsets, &offset_len, tax_hits);
            i++;
        }
        else {
            i += (flag+1);
        }
    }
    //free(offsets);
}


void query_kmer(uint64_t kmer, struct min_kmer_db * tmi, uint64_t ** offsets, uint64_t * offset_len, struct read_tax_info ** tax_hits){
    uint64_t canonical_kmer = tax_globals.protein == 0 ? canonical_representation2(kmer) : kmer;
    //uint64_t canonical_kmer = canonical_representation2(kmer);
    uint64_t pos = tax_globals.KMER_NOT_FOUND;
    if (tmi->kmer_count == 0) { //check if minimizer is not set, if not, set it and get kmer byte offsets
        pos = query_minimizer(kmer, canonical_kmer, tmi, offsets, offset_len);
    }
    else if (tmi->kmer_count > 0) { //if minimizer is set, do binary search to see if kmer is in range, if not then recompute minimizer and search
        /*if (canonical_kmer >= *(uint64_t *)(mmkmer+offsets[0]) && canonical_kmer <= *(uint64_t *)(mmkmer+offsets[tmi->kmer_count - 1])) {
            pos = kmer_binary_search(offsets,tmi->kmer_count,canonical_kmer);
        }*/
        pos = kmer_binary_search(*offsets,tmi->kmer_count,canonical_kmer);
        if (pos == tax_globals.KMER_NOT_FOUND) {
            pos = query_minimizer(kmer, canonical_kmer, tmi, offsets, offset_len);
        }
    }
    
    if (pos != tax_globals.KMER_NOT_FOUND) {
        uint64_t * kmer_count = (uint64_t *)(mmkmer+pos+sizeof(uint64_t));
        if (*kmer_count >= tax_globals.kmer_cutoff) { //enforce kmer cutoff
            return;
        }
        uint64_t num_ref_places = *kmer_count < tax_globals.tax_cutoff ? *kmer_count : tax_globals.tax_cutoff;
        uint64_t * rsi_offset = (uint64_t *)(mmkmer+pos+2*sizeof(uint64_t));
        uint64_t * num_refs = (uint64_t *)(mmrsi+(*rsi_offset));
        if (*num_refs >= tax_globals.tax_cutoff) { //enforce taxid number cutoff
            return;
        }
        //uint64_t * refs = (uint64_t *)malloc(sizeof(uint64_t) * (*num_refs));
        //memcpy(refs, mmrsi+(*rsi_offset)+sizeof(uint64_t), (*num_refs)*sizeof(uint64_t));
        uint64_t * refs = (uint64_t *) (mmrsi + (*rsi_offset) + sizeof(uint64_t));
        uint64_t total_offset = (*rsi_offset)+sizeof(uint64_t)+(num_ref_places)*sizeof(uint64_t);
        //float * scores = (float *)malloc(sizeof(float) * (*num_refs));
        //memcpy(scores, mmrsi+total_offset, (*num_refs)*sizeof(float));
        float * scores = (float *)(mmrsi+total_offset);
        int k = 0;
        for (k=0; k < *num_refs; k++) {
            //uint64_t tref = refs[k];
            //uint64_t ttax = db_info->usid_taxid[refs[k]];
            //std::cout << refs[k] << "\t" << scores[k] << std::endl;
            //tax_hits->operator[](refs[k]) += scores[k];
            //tax_counts->operator[](refs[k]) += 1.0;
            struct read_tax_info * trti = NULL;
            HASH_FIND(hh, *tax_hits, &refs[k], sizeof(uint64_t), trti);
            if (trti == NULL) {
                trti = (struct read_tax_info *)malloc(sizeof(struct read_tax_info));
                trti->taxid = refs[k];
                trti->score = 0.0;
                trti->counts = 0;
                HASH_ADD(hh, *tax_hits, taxid, sizeof(uint64_t), trti);
            }
            trti->score += scores[k];
            trti->counts += 1;
        }
        //free(refs);
        //free(scores);
        refs = NULL;
        scores = NULL;
    }
}

int query_kmer2(uint64_t kmer, struct min_kmer_db * tmi, uint64_t ** offsets, uint64_t * offset_len, struct read_tax_info ** tax_hits){
    uint64_t canonical_kmer = tax_globals.protein == 0 ? canonical_representation2(kmer) : kmer;
    //uint64_t canonical_kmer = canonical_representation2(kmer);
    uint64_t pos = tax_globals.KMER_NOT_FOUND;
    if (tmi->kmer_count == 0) { //check if minimizer is not set, if not, set it and get kmer byte offsets
        pos = query_minimizer(kmer, canonical_kmer, tmi, offsets, offset_len);
    }
    else if (tmi->kmer_count > 0) { //if minimizer is set, do binary search to see if kmer is in range, if not then recompute minimizer and search
        /*if (canonical_kmer >= *(uint64_t *)(mmkmer+offsets[0]) && canonical_kmer <= *(uint64_t *)(mmkmer+offsets[tmi->kmer_count - 1])) {
         pos = kmer_binary_search(offsets,tmi->kmer_count,canonical_kmer);
         }*/
        pos = kmer_binary_search(*offsets,tmi->kmer_count,canonical_kmer);
        if (pos == tax_globals.KMER_NOT_FOUND) {
            pos = query_minimizer(kmer, canonical_kmer, tmi, offsets, offset_len);
        }
    }
    
    if (pos != tax_globals.KMER_NOT_FOUND) {
        uint64_t * kmer_count = (uint64_t *)(mmkmer+pos+sizeof(uint64_t));
        if (*kmer_count >= tax_globals.kmer_cutoff) { //enforce kmer cutoff
            return 0;
        }
        uint64_t num_ref_places = *kmer_count < tax_globals.tax_cutoff ? *kmer_count : tax_globals.tax_cutoff;
        uint64_t * rsi_offset = (uint64_t *)(mmkmer+pos+2*sizeof(uint64_t));
        uint64_t * num_refs = (uint64_t *)(mmrsi+(*rsi_offset));
        if (*num_refs >= tax_globals.tax_cutoff) { //enforce taxid number cutoff
            return 0;
        }
        //uint64_t * refs = (uint64_t *)malloc(sizeof(uint64_t) * (*num_refs));
        //memcpy(refs, mmrsi+(*rsi_offset)+sizeof(uint64_t), (*num_refs)*sizeof(uint64_t));
        uint64_t * refs = (uint64_t *) (mmrsi + (*rsi_offset) + sizeof(uint64_t));
        uint64_t total_offset = (*rsi_offset)+sizeof(uint64_t)+(num_ref_places)*sizeof(uint64_t);
        //float * scores = (float *)malloc(sizeof(float) * (*num_refs));
        //memcpy(scores, mmrsi+total_offset, (*num_refs)*sizeof(float));
        float * scores = (float *)(mmrsi+total_offset);
        int k = 0;
        for (k=0; k < *num_refs; k++) {
            //uint64_t tref = refs[k];
            //uint64_t ttax = db_info->usid_taxid[refs[k]];
            //std::cout << refs[k] << "\t" << scores[k] << std::endl;
            //tax_hits->operator[](refs[k]) += scores[k];
            //tax_counts->operator[](refs[k]) += 1.0;
            struct read_tax_info * trti = NULL;
            HASH_FIND(hh, *tax_hits, &refs[k], sizeof(uint64_t), trti);
            if (trti == NULL) {
                trti = (struct read_tax_info *)malloc(sizeof(struct read_tax_info));
                trti->taxid = refs[k];
                trti->score = 0.0;
                trti->counts = 0;
                HASH_ADD(hh, *tax_hits, taxid, sizeof(uint64_t), trti);
            }
            trti->score += scores[k];
            trti->counts += 1;
        }
        //free(refs);
        //free(scores);
        refs = NULL;
        scores = NULL;
        return 1;
    }
    return 0;
}


uint64_t query_minimizer(uint64_t kmer, uint64_t ckmer, struct min_kmer_db * tmi, uint64_t ** offsets, uint64_t * offset_len){
    uint64_t pos = tax_globals.KMER_NOT_FOUND;
    uint64_t minimizer = kmer_minimizer(ckmer);
    //printf("query minimizer %llu\n",minimizer);
    if (mi[minimizer] == NULL) {
        tmi->kmer_count = 0;
        tmi->offset = 0;
    }
    else{
        tmi = mi[minimizer];
        //printf("number of kmers: %llu\t%llu\n",tmi->kmer_count,*offset_len);
        *offsets = (uint64_t *)(mmkmer + tmi->offset);
        *offset_len = tmi->kmer_count;
        /*if (tmi->kmer_count > *offset_len) {
            free(*offsets);
            *offsets = NULL;
            *offset_len = tmi->kmer_count;
            //printf("number of offsets (in if): %llu\n",*offset_len);
            *offsets = (uint64_t *)malloc(sizeof(uint64_t)*(*offset_len));
        }
        //printf("number of offsets: %llu\n",*offset_len);
        memcpy(*offsets, mmkmer+tmi->offset, tmi->kmer_count*sizeof(uint64_t));*/
        pos = kmer_binary_search(*offsets,tmi->kmer_count,ckmer);
    }
    //printf("exiting query minimizer\n");
    return pos;
}


