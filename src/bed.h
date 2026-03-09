/* @file bed.h
** BED file parsing for targeted read simulation.
*/

#ifndef SQ_BED_H
#define SQ_BED_H

#include <stdint.h>
#include "ref.h"

typedef struct {
    int32_t *seq_idx;    // index into ref_t arrays
    int32_t *starts;     // 0-based start
    int32_t *ends;       // 0-based end (exclusive)
    int n;               // number of intervals
    int64_t sum;         // total interval length (for uniform sampling)
} bed_t;

bed_t *load_bed(const char *filename, ref_t *ref);
void free_bed(bed_t *bed);

#endif
