// From the squigulator repository root directory:
// make build/libsquigulator.a
// gcc -Wall -O2 -I include/ examples/example.c build/libsquigulator.a -o example -lz -lm -lpthread
// ./example

#include <stdio.h>
#include <squigulator.h>

int main(){

    squig_opt_t opt;
    squig_init_opt(&opt);
    opt.output_blow5 = "test.blow5";
    opt.profile = "dna-r9-prom";
    opt.num_threads = 8;
    opt.batch_size = 1000;

    squig_t *sq = squig_init("test/nCoV-2019.reference.fasta", &opt);

    fprintf(stdout,"read_id\tread_group\tdigitisation\toffset\trange\tlen_raw_signal\tsum_raw_signal\n");

    for(int n=0; n<5; n++){ //5 batches of 1000  reads in each batch

        squig_batch_t *batch = squig_sim_batch(sq);

        for(int i=0; i < opt.batch_size; i++){
            squig_rec_t *rec = squig_get_rec(sq, batch, i); //thread safe, can be called in parallel
            fprintf(stdout, "%s\t%d\t%.0f\t%.0f\t%.0f\t%.0f\t%ld\t", rec->read_id, rec->read_group, rec->digitisation, rec->offset, rec->range, rec->sampling_rate, rec->len_raw_signal);
            int64_t raw_signal_sum = 0;
            for(int j=0; j < rec->len_raw_signal; j++){
                raw_signal_sum += rec->raw_signal[j];
            }
            fprintf(stdout, "%ld\n", raw_signal_sum);
        }

        squig_free_batch(batch);
    }

    squig_free(sq);

    return 0;
}