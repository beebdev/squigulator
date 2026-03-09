/* @file  bed.c
** BED file parsing for targeted read simulation.
**
** Loads BED3 intervals (chrom, start, end) and resolves chromosome names
** to ref_t indices for use in targeted position selection.
*/

/*
MIT License

Copyright (c) 2023 Hasindu Gamaarachchi

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bed.h"
#include "error.h"

// find ref_t index for a chromosome name, -1 if not found
static int find_seq_idx(ref_t *ref, const char *chrom) {
    for (int i = 0; i < ref->num_ref; i++) {
        if (strcmp(ref->ref_names[i], chrom) == 0) {
            return i;
        }
    }
    return -1;
}

bed_t *load_bed(const char *filename, ref_t *ref) {
    FILE *fp = fopen(filename, "r");
    F_CHK(fp, filename);

    int capacity = 64;
    bed_t *bed = (bed_t *)malloc(sizeof(bed_t));
    MALLOC_CHK(bed);
    bed->seq_idx = (int32_t *)malloc(capacity * sizeof(int32_t));
    MALLOC_CHK(bed->seq_idx);
    bed->starts = (int32_t *)malloc(capacity * sizeof(int32_t));
    MALLOC_CHK(bed->starts);
    bed->ends = (int32_t *)malloc(capacity * sizeof(int32_t));
    MALLOC_CHK(bed->ends);
    bed->n = 0;
    bed->sum = 0;

    char line[4096];
    int line_num = 0;
    while (fgets(line, sizeof(line), fp)) {
        line_num++;
        // skip empty lines and comments
        if (line[0] == '#' || line[0] == '\n' || line[0] == '\r') {
            continue;
        }

        char chrom[256];
        int32_t start, end;
        int ret = sscanf(line, "%255s %d %d", chrom, &start, &end);
        if (ret < 3) {
            WARNING("BED line %d: expected at least 3 columns, skipping", line_num);
            continue;
        }
        if (start < 0 || end <= start) {
            WARNING("BED line %d: invalid interval %s:%d-%d, skipping", line_num, chrom, start, end);
            continue;
        }

        int idx = find_seq_idx(ref, chrom);
        if (idx < 0) {
            WARNING("BED line %d: chromosome '%s' not found in reference, skipping", line_num, chrom);
            continue;
        }

        // clamp end to reference length
        if (end > ref->ref_lengths[idx]) {
            WARNING("BED line %d: end %d exceeds reference length %d for %s, clamping",
                    line_num, end, ref->ref_lengths[idx], chrom);
            end = ref->ref_lengths[idx];
        }

        // grow arrays if needed
        if (bed->n >= capacity) {
            capacity *= 2;
            bed->seq_idx = (int32_t *)realloc(bed->seq_idx, capacity * sizeof(int32_t));
            MALLOC_CHK(bed->seq_idx);
            bed->starts = (int32_t *)realloc(bed->starts, capacity * sizeof(int32_t));
            MALLOC_CHK(bed->starts);
            bed->ends = (int32_t *)realloc(bed->ends, capacity * sizeof(int32_t));
            MALLOC_CHK(bed->ends);
        }

        bed->seq_idx[bed->n] = idx;
        bed->starts[bed->n] = start;
        bed->ends[bed->n] = end;
        bed->sum += (end - start);
        bed->n++;
    }

    fclose(fp);

    if (bed->n == 0) {
        ERROR("No valid BED intervals loaded from '%s'", filename);
        free(bed->seq_idx);
        free(bed->starts);
        free(bed->ends);
        free(bed);
        exit(EXIT_FAILURE);
    }

    INFO("load_bed: Loaded %d intervals with total length %ld bases from %s", bed->n, (long)bed->sum, filename);
    return bed;
}

void free_bed(bed_t *bed) {
    if (bed) {
        free(bed->seq_idx);
        free(bed->starts);
        free(bed->ends);
        free(bed);
    }
}
