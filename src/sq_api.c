
#define _XOPEN_SOURCE 700
#include <string.h>
#include <slow5/slow5.h>

#include "squigulator.h"
#include "version.h"
#include "sq.h"
#include "ref.h"
#include "format.h"
#include "seq.h"
#include "rand.h"
#include "error.h"
#include "misc.h"





struct squig_s {
    core_t *core;
    squig_opt_t opt;
    // can later have the header attributes and enum labels here
};

struct squig_batch_s {
    db_t *db;
    squig_rec_t **rec;
    int n_rec;
};

core_t *init_core(opt_t opt, profile_t p, char *refname, char *output_file, char *fasta, char *paf, char *sam, char *trans_count);
void init_opt(opt_t *opt);
profile_t set_profile(char *prof_name, opt_t *opt);
void free_core(core_t *core);
void free_db(db_t* db);
db_t* init_db(core_t* core, int32_t n_rec);
void process_db(core_t* core,db_t* db);
void output_db(core_t* core, db_t* db);

void squig_init_opt(squig_opt_t *opts){
    opts->avg_rlen = 1000;
    opts->random_seed = 0;
    opts->num_threads = 8;
    opts->batch_size = 1000;
    opts->output_blow5 = "output.blow5";
    opts->profile = "dna-r9-prom";
}

squig_t *squig_init(char *ref_fasta, squig_opt_t *opts){
    squig_t *sq = (squig_t *)malloc(sizeof(squig_t));
    MALLOC_CHK(sq);

    opt_t opt;
    init_opt(&opt);

    opt.rlen = opts->avg_rlen;
    opt.num_thread = opts->num_threads;
    opt.batch_size = opts->batch_size;
    if (opts->random_seed == 0){
        opts->random_seed = realtime();
    }
    opt.seed = opts->random_seed;

    profile_t p = set_profile(opts->profile, &opt);

    char *output_file = opts->output_blow5;
    char *fasta = NULL;
    char *paf = NULL;
    char *sam = NULL;
    char *trans_count = NULL;


    sq->core = init_core(opt, p, ref_fasta, output_file, fasta, paf, sam, trans_count);

    sq->opt = *opts;

    return sq;
}

void squig_free(squig_t *sq){
    free_core(sq->core);
    free(sq);

}


void squig_free_batch(squig_batch_t *batch){
    for(int i = 0; i < batch->n_rec; i++){
        free(batch->rec[i]->read_id);
        free(batch->rec[i]->raw_signal);
        free(batch->rec[i]->channel_number);
        free(batch->rec[i]);
    }
    free(batch->rec);
    free_db(batch->db);
    free(batch);
}

squig_batch_t *squig_sim_batch(squig_t *sq){

    squig_batch_t *batch = (squig_batch_t *)malloc(sizeof(squig_batch_t));
    MALLOC_CHK(batch);
    batch->n_rec = sq->opt.batch_size;
    batch->db = init_db(sq->core, batch->n_rec);
    batch->rec = (squig_rec_t **)malloc(batch->n_rec * sizeof(squig_rec_t *));
    MALLOC_CHK(batch->rec);
    for(int i = 0; i < batch->n_rec; i++){
        batch->rec[i] = (squig_rec_t *)malloc(sizeof(squig_rec_t));
        MALLOC_CHK(batch->rec[i]);
    }
    process_db(sq->core, batch->db);
    output_db(sq->core, batch->db);
    return batch;
}

//this function is hacky, inefficient and repeats unnecessary work. Can be fixed later if performance becomes an issue,
// by integrating into core function. Not wanting to touch  those at the moment
squig_rec_t *squig_get_rec(squig_t *sq, squig_batch_t *batch, int i){

    if(i < 0 || i >= batch->n_rec){
        fprintf(stderr, "Index out of bounds. batch size: %d, requested index: %d\n", batch->n_rec, i);
        return NULL;
    }

    slow5_rec_t *rec = NULL;

    size_t bytes = batch->db->mem_bytes[i];
    char *mem = (char *)malloc(bytes);
    MALLOC_CHK(mem);

    if(sq->core->sp->format == SLOW5_FORMAT_ASCII){
        memcpy(mem, batch->db->mem_records[i], bytes);
        mem[bytes-1] = '\0';
    } else if (sq->core->sp->format == SLOW5_FORMAT_BINARY){
        memcpy(mem, batch->db->mem_records[i]+8, bytes-8);
    } else {
        fprintf(stderr, "Something is not right!\n");
        exit(EXIT_FAILURE);
    }

    if (slow5_decode(&mem, &batch->db->mem_bytes[i], &rec, sq->core->sp) < 0){
        fprintf(stderr,"Error in decoding record\n");
        exit(EXIT_FAILURE);
    }

    //from primary SLOW5 fields
    batch->rec[i]->read_id = strdup(rec->read_id);
    batch->rec[i]->read_group = rec->read_group;
    batch->rec[i]->digitisation = rec->digitisation;
    batch->rec[i]->offset = rec->offset;
    batch->rec[i]->range = rec->range;
    batch->rec[i]->sampling_rate = rec->sampling_rate;
    batch->rec[i]->len_raw_signal = rec->len_raw_signal;
    batch->rec[i]->raw_signal = rec->raw_signal; rec->raw_signal = NULL;

    //from auxilliary SLOW5 fields
    int ret=0;
    uint64_t len=0;
    batch->rec[i]->channel_number = strdup(slow5_aux_get_string(rec,"channel_number", &len, &ret)); NEG_CHK(ret);
    batch->rec[i]->median_before = slow5_aux_get_double(rec,"median_before", &ret); NEG_CHK(ret);
    batch->rec[i]->read_number = slow5_aux_get_int32(rec,"read_number", &ret); NEG_CHK(ret);
    batch->rec[i]->start_mux = slow5_aux_get_uint8(rec,"start_mux", &ret); NEG_CHK(ret);
    batch->rec[i]->start_time = slow5_aux_get_uint64(rec,"start_time", &ret); NEG_CHK(ret);
    batch->rec[i]->end_reason = slow5_aux_get_enum(rec,"end_reason", &ret); NEG_CHK(ret);

    free(mem);
    slow5_rec_free(rec);

    return batch->rec[i];
}
