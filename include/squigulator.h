#ifndef SQUIGULATOR_H
#define SQUIGULATOR_H

// squigulator C API which is still in a very early stages
/*
Currently, signal is in parallel dumped to a BLOW5, which can be disabled later if deemed a bottleneck
squig_get_rec can be optimised significantly, as it does some repeated work, that can have been better handled in the core
Need to hide symbols in the .a
Need to add some test cases
*/

#include <stdint.h>

typedef struct {
    //from primary SLOW5 fields
    char* read_id;
    uint32_t read_group;
    double digitisation;
    double offset;
    double range;
    double sampling_rate;
    uint64_t len_raw_signal;
    int16_t* raw_signal;

    //from auxilliary SLOW5 fields
    char *channel_number;
    double median_before;
    int32_t read_number;
    uint8_t start_mux;
    uint64_t start_time;
    uint8_t end_reason;

} squig_rec_t;

typedef struct {
    int32_t avg_rlen;
    int64_t random_seed;
    int32_t num_threads;
    int32_t batch_size;
    char *output_blow5;
    char *profile;

    //things to add later ...

} squig_opt_t;

typedef struct squig_s squig_t;
typedef struct squig_batch_s squig_batch_t;

void squig_init_opt(squig_opt_t *opts);
squig_t *squig_init(char *ref_fasta, squig_opt_t *opts);
void squig_free(squig_t *sq);
squig_batch_t* squig_sim_batch(squig_t *sq);
void squig_free_batch(squig_batch_t *batch);
squig_rec_t *squig_get_rec(squig_t *sq, squig_batch_t *batch, int i);


#endif