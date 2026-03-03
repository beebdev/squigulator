#ifndef SQUIGULATOR_H
#define SQUIGULATOR_H

#include <slow5/slow5.h>

typedef struct squig_s squig_t;
typedef struct squig_batch_s squig_batch_t;


squig_t *squig_init(char *ref_fasta, char* output_blow5, char *profile );

void squig_free(squig_t *sq);

squig_batch_t* squig_sim_batch(squig_t *sq, int num_reads);

void squig_free_batch(squig_batch_t *batch);

slow5_rec_t *squig_decode(squig_t *sq, squig_batch_t *batch, int i);

#endif