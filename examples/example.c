// gcc -Wall -I include/ examples/example.c build/libsquigulator.a -o example -lz -lm -lpthread

#include <squigulator.h>
#include <slow5/slow5.h>

int main(){
    squig_t *sq = squig_init("test/nCoV-2019.reference.fasta", "out.blow5","dna-r9-prom");
    squig_batch_t *batch = squig_sim_batch(sq, 10);

    for(int i=0; i<10; i++){
        slow5_rec_t *rec = squig_decode(sq, batch, i);
        fprintf(stdout, "%s\n", rec->read_id);
        slow5_rec_free(rec);
    }

    squig_free_batch(batch);
    squig_free(sq);
    return 0;
}