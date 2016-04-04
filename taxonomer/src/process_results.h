//
//  process_results.h
//  taxonomer
//
//  Created by steven on 10/7/14.
//  Copyright (c) 2014 yandell lab. All rights reserved.
//

#ifndef __taxonomer__process_results__
#define __taxonomer__process_results__

#include "taxonomer_headers.h"
#include "query_db.h"

struct read_output_info {
    uint64_t lca;
    uint64_t num_taxids;
    int tax_depth;
    float max_score;
    float mean_score;
    float var_score;
    float kappa; //parameters for gamma distribution
    float theta;
};

void load_tri_file(char * tri_file);
void lca_from_taxids(struct read_tax_info * taxid_results, struct read_output_info * roi, uint64_t ** max_taxids);
void lca_from_taxids_kw(struct read_tax_info * taxid_results, struct read_output_info * roi, uint64_t ** max_taxids);
void ties_from_taxids(struct read_tax_info * taxid_results, struct read_output_info * roi, uint64_t ** max_taxids);
uint64_t get_lca(uint64_t * max_taxids, uint64_t size);
static uint64_t get_lca2(uint64_t a, uint64_t b);
static float score_sample_variance(struct read_tax_info * taxid_results, float mean);
static int taxid_depth(uint64_t taxid);

#endif /* defined(__taxonomer__process_results__) */
