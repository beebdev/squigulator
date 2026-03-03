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
    slow5_file_t *sp;
};

struct squig_batch_s {
    db_t *db;
};

// db_t* init_db(core_t* core, int32_t n_rec);
core_t *init_core(opt_t opt, profile_t p, char *refname, char *output_file, char *fasta, char *paf, char *sam, char *trans_count);
void init_opt(opt_t *opt);
profile_t set_profile(char *prof_name, opt_t *opt);
void free_core(core_t *core);
void free_db(db_t* db);
db_t* init_db(core_t* core, int32_t n_rec);
void process_db(core_t* core,db_t* db);
void output_db(core_t* core, db_t* db);

squig_t *squig_init(char *ref_fasta, char* output_blow5, char *profile){
    squig_t *sq = (squig_t *)malloc(sizeof(squig_t *));
    MALLOC_CHK(sq);

    opt_t opt;
    init_opt(&opt);

    profile_t p = set_profile(profile, &opt);

    char *output_file = output_blow5;
    char *fasta = NULL;
    char *paf = NULL;
    char *sam = NULL;
    char *trans_count = NULL;

    opt.seed = 1; //todo, fix this

    sq->core = init_core(opt, p, ref_fasta, output_file, fasta, paf, sam, trans_count);

    return sq;
}

void squig_free(squig_t *sq){
    free_core(sq->core);
    free(sq);

}


void squig_free_batch(squig_batch_t *batch){
    free_db(batch->db);
    free(batch);
}

squig_batch_t *squig_sim_batch(squig_t *sq, int num_reads){
    //todo check if the n<=c
    squig_batch_t *batch = (squig_batch_t *)malloc(sizeof(squig_batch_t *));
    batch->db = init_db(sq->core, num_reads);
    process_db(sq->core, batch->db);
    output_db(sq->core, batch->db);
    return batch;
}

slow5_rec_t *squig_decode(squig_t *sq, squig_batch_t *batch, int i){
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

    free(mem);

    return rec;
}
