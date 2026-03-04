/* @file  gensig.c
** squigulator, a nanopore signal simulator
**
** @@
******************************************************************************/

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

#include <assert.h>
#include "sq.h"
#include "seq.h"
#include "format.h"

char *attach_prefix(core_t *core, const char *read, int32_t *len, aln_t *aln);
void gen_prefix_dna(int16_t *raw_signal, int64_t* n, int64_t *c, profile_t *profile, double offset);
int16_t *gen_prefix_rna(core_t *core, int16_t *raw_signal, int64_t* n, int64_t *c, double offset, int tid, aln_t *aln);


static inline void slow5_hdr_add_check(const char *attr_name, slow5_hdr_t *header){
    if (slow5_hdr_add(attr_name, header) < 0){
        ERROR("Error adding %s attribute", attr_name);
        exit(EXIT_FAILURE);
    }
}

static inline void slow5_hdr_set_check(const char *attr_name, const char *value, slow5_hdr_t *header){
    if (slow5_hdr_set(attr_name, value, 0, header) < 0){
        ERROR("Error setting %s attribute in read group %d", attr_name, 0);
        exit(EXIT_FAILURE);
    }
}

void set_header_attributes(slow5_file_t *sp, int8_t rna, int8_t r10, int8_t prom, double sample_frequency){

    slow5_hdr_t *header=sp->header;

    //header group attributes

    slow5_hdr_add_check("asic_id", header);
    slow5_hdr_add_check("exp_start_time", header);
    slow5_hdr_add_check("experiment_name", header);
    slow5_hdr_add_check("experiment_type", header);
    slow5_hdr_add_check("flow_cell_id", header);
    slow5_hdr_add_check("flow_cell_product_code", header);

    slow5_hdr_add_check("protocol_run_id", header);
    slow5_hdr_add_check("protocol_start_time", header);
    slow5_hdr_add_check("run_id", header);
    slow5_hdr_add_check("sample_frequency", header);
    slow5_hdr_add_check("sample_id", header);
    slow5_hdr_add_check("sequencing_kit", header);
    slow5_hdr_add_check("sequencer_position", header);
    slow5_hdr_add_check("system_name", header);


    //dna/rna
    const char* experiment_type = rna ? "rna" : "genomic_dna" ;

    //flow cell
    const char* flow_cell = ".";
    if(r10){
        if(rna){
            flow_cell = prom ? "FLO-PRO004RA" : "FLO-MIN004RA";
        } else {
            flow_cell = prom ? "FLO-PRO114M" : "FLO-MIN114";
        }
    } else {
        flow_cell = prom ? "FLO-PRO002" : "FLO-MIN106";
    }

    //sequencing kit
    const char* kit = ".";
    if(rna){
        kit = r10 ? "sqk-rna004" : "sqk-rna002";
    } else{
        kit = r10 ? "sqk-lsk114" : "sqk-lsk109";
    }

    //sample_frequency
    if(sample_frequency<=0 || sample_frequency>1000000000){
        ERROR("%s","A weird sample frequency. It should be between 0 and 1000000000 Hz");
        exit(EXIT_FAILURE);
    }
    char sample_frequency_str[100];
    sprintf(sample_frequency_str, "%d", (int)sample_frequency);

    //set header attributes
    slow5_hdr_set_check("asic_id", "asic_id_0", header);
    slow5_hdr_set_check("exp_start_time", "2026-02-20T00:00:00Z", header);
    slow5_hdr_set_check("experiment_name", "experiment_0", header);
    slow5_hdr_set_check("experiment_type", experiment_type, header);
    slow5_hdr_set_check("flow_cell_id", "FAN00000", header);
    slow5_hdr_set_check("flow_cell_product_code", flow_cell, header);
    slow5_hdr_set_check("protocol_run_id", "protocol_run_0", header);
    slow5_hdr_set_check("protocol_start_time", "2026-02-20T00:00:00Z", header);
    slow5_hdr_set_check("run_id", "run_0", header);
    slow5_hdr_set_check("sample_frequency", sample_frequency_str, header);
    slow5_hdr_set_check("sample_id", "squigulator", header);
    slow5_hdr_set_check("sequencing_kit", kit, header);
    slow5_hdr_set_check("sequencer_position", "P0", header);
    slow5_hdr_set_check("system_name", "PCA000000", header);


}

void set_header_aux_fields(slow5_file_t *sp){

    //add auxilliary field: channel number
    if (slow5_aux_add("channel_number", SLOW5_STRING, sp->header) < 0){
        ERROR("%s","Error adding channel_number auxilliary field\n");
        exit(EXIT_FAILURE);
    }

    //add auxilliary field: median_before
    if (slow5_aux_add("median_before", SLOW5_DOUBLE, sp->header) < 0) {
        ERROR("%s","Error adding median_before auxilliary field\n");
        exit(EXIT_FAILURE);
    }

    //add auxilliary field: read_number
    if(slow5_aux_add("read_number", SLOW5_INT32_T, sp->header) < 0){
        ERROR("%s","Error adding read_number auxilliary field\n");
        exit(EXIT_FAILURE);
    }
    //add auxilliary field: start_mux
    if(slow5_aux_add("start_mux", SLOW5_UINT8_T, sp->header) < 0){
        ERROR("%s","Error adding start_mux auxilliary field\n");
        exit(EXIT_FAILURE);
    }
    //add auxilliary field: start_time
    if(slow5_aux_add("start_time", SLOW5_UINT64_T, sp->header) < 0){
        ERROR("%s","Error adding start_time auxilliary field\n");
        exit(EXIT_FAILURE);
    }
    //add auxilliary field: end_reason
    const char *enum_labels[]={"unknown","partial","mux_change","unblock_mux_change","data_service_unblock_mux_change","signal_positive","signal_negative"};
    uint8_t num_labels = 7;
    if (slow5_aux_add_enum("end_reason", enum_labels, num_labels, sp->header) < 0){
        ERROR("%s","Error adding end_reason auxilliary field\n");
        exit(EXIT_FAILURE);
    }

}

void set_record_primary_fields(profile_t *profile, slow5_rec_t *slow5_record, char *read_id, double offset, int64_t len_raw_signal, int16_t *raw_signal){

    slow5_record -> read_id = read_id;
    slow5_record-> read_id_len = strlen(slow5_record -> read_id);
    slow5_record -> read_group = 0;
    slow5_record -> digitisation = profile->digitisation;
    slow5_record -> offset = offset;
    slow5_record -> range = profile->range;
    slow5_record -> sampling_rate = profile->sample_rate;
    slow5_record -> len_raw_signal = len_raw_signal;
    slow5_record -> raw_signal = raw_signal;

}

void set_record_aux_fields(slow5_rec_t *slow5_record, slow5_file_t *sp, double median_before, int32_t read_number, uint64_t start_time){

    const char *channel_number = "0";
    uint8_t start_mux = 0;

    if(slow5_aux_set_string(slow5_record, "channel_number", channel_number, sp->header) < 0){
        ERROR("%s","Error setting channel_number auxilliary field\n");
        exit(EXIT_FAILURE);
    }

    if(slow5_aux_set(slow5_record, "median_before", &median_before, sp->header) < 0){
        ERROR("%s","Error setting median_before auxilliary field\n");
        exit(EXIT_FAILURE);
    }

    if(slow5_aux_set(slow5_record, "read_number", &read_number, sp->header) < 0){
        ERROR("%s","Error setting read_number auxilliary field\n");
        exit(EXIT_FAILURE);
    }

    if(slow5_aux_set(slow5_record, "start_mux", &start_mux, sp->header) < 0){
        ERROR("%s","Error setting start_mux auxilliary field\n");
        exit(EXIT_FAILURE);
    }

    if(slow5_aux_set(slow5_record, "start_time", &start_time, sp->header) < 0){
        ERROR("%s","Error setting start_time auxilliary field\n");
        exit(EXIT_FAILURE);
    }


    uint8_t end_reason = 0;
    if(slow5_aux_set(slow5_record, "end_reason", &end_reason, sp->header) < 0){
        ERROR("%s","Error setting end_reason auxilliary field\n");
        exit(EXIT_FAILURE);
    }


}


int16_t * gen_sig_core_seq(core_t *core, int16_t *raw_signal, int64_t* n, int64_t *c, double offset, const char *read, int32_t len, int tid, aln_t *aln){

    profile_t *profile = &core->profile;
    uint32_t kmer_size = core->kmer_size;
    model_t *pore_model = core->model;
    if(core->opt.meth_freq){
        pore_model = core->cpgmodel;
    }

    int8_t ideal = core->opt.flag & SQ_IDEAL;
    int8_t ideal_time = core->opt.flag & SQ_IDEAL_TIME;
    int8_t ideal_amp = core->opt.flag & SQ_IDEAL_AMP;
    int sps = (int)profile->dwell_mean;

    int64_t n_kmers = len-kmer_size+1;

    if(len<kmer_size){ //a hack
        n_kmers=5;
        read="ACGTACGTACGT";
    }

    if(aln) aln->sig_start = 0;

    for (int i=0; i< n_kmers; i++){
        uint32_t kmer_rank = get_kmer_rank(read+i, kmer_size);
        if(core->opt.meth_freq){
            kmer_rank = get_meth_kmer_rank(read+i, kmer_size);
        }
        if(!(ideal || ideal_time)){
            sps = round(nrng(core->rand_time[tid]));
            sps = sps < 1 ? -sps + 1 : sps;
            //fprintf(stderr,"%d %d %d %d\n",*n,n_kmers,*c,sps);
        }
        for(int j=0; j<sps; j++){
            if(*n==*c){
                *c *= 2;
                raw_signal = (int16_t *)realloc(raw_signal, *c*sizeof(int16_t));
            }
            float s = 0;
            if(ideal || ideal_amp){
                s=pore_model[kmer_rank].level_mean;
            } else {
                s=nrng(core->kmer_gen[tid][kmer_rank]);
            }
            raw_signal[*n] = s*(profile->digitisation)/(profile->range)-(offset);
            *n = *n+1;
        }
        if(aln) {
            if(aln->ss_n==aln->ss_c){
                aln->ss_c *= 2;
                aln->ss = (int32_t *)realloc(aln->ss, aln->ss_c*sizeof(int32_t));
                MALLOC_CHK(aln->ss);
            }
            aln->ss[aln->ss_n] = sps >=0 ? sps : 0;
            aln->ss_n++;
        }
    }

    if(aln) aln->sig_end = *n;

    return raw_signal;

}




int16_t *gen_sig_core(core_t *core, const char *read, int32_t len, double *offset, double *median_before, int64_t *len_raw_signal, int tid, aln_t *aln){

    profile_t *profile = &core->profile;
    uint32_t kmer_size = core->kmer_size;
    //model_t *pore_model = core->model;

    int8_t ideal = core->opt.flag & SQ_IDEAL;
    //int8_t ideal_time = core->opt.ideal_time;
    //int8_t ideal_amp = core->opt.ideal_amp;

    int64_t n_kmers = len < kmer_size ? 1 : len-kmer_size+1;

    int64_t n=0;
    int sps = (int)profile->dwell_mean;

    int64_t c = n_kmers * sps + 2000;
    int16_t *raw_signal = (int16_t *)malloc(c*sizeof(int16_t));
    MALLOC_CHK(raw_signal);

    if(ideal) {
        *offset = profile->offset_mean;
        *median_before = profile->median_before_mean;
    } else {
        *offset = nrng(core->rand_offset[tid]);
        *median_before = nrng(core->rand_median_before[tid]);
    }

    //todo adaptor sequence and prefix
    // if(core->opt.prefix && !core->opt.rna){
    //     gen_prefix_dna(raw_signal,&n,&c, profile, *offset);
    // }

    char *tmpread = NULL;
    if(core->opt.flag & SQ_PREFIX){
        read = tmpread = attach_prefix(core, read, &len, aln);
    }
    //fprintf(stderr, "read: %s\n", read);

    raw_signal = gen_sig_core_seq(core, raw_signal, &n, &c, *offset, read, len, tid, aln);

    if(core->opt.flag & SQ_PREFIX && core->opt.flag & SQ_RNA){
        raw_signal=gen_prefix_rna(core, raw_signal,&n,&c, *offset, tid, aln);
    }
    if(tmpread){
        free(tmpread);
    }
    assert(n<=c);

    *len_raw_signal = n;
    return raw_signal;
}


int16_t *gen_sig(core_t *core, const char *read, int32_t len, double *offset, double *median_before, int64_t *len_raw_signal, int8_t rna, int tid, aln_t *aln){
    int16_t *sig = gen_sig_core(core, read, len, offset, median_before, len_raw_signal, tid, aln);
    if(rna){
        for(int i=0; i<*len_raw_signal/2; i++){
            int16_t tmp = sig[i];
            sig[i] = sig[*len_raw_signal-1-i];
            sig[*len_raw_signal-1-i] = tmp;
        }
    }
    return sig;
}