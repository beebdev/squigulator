#cython: language_level=3
#from libc.stdio cimport *
from libc.stdint cimport *
#from libc.stdlib cimport *

cdef extern from "squigulator.h":

	ctypedef struct squig_t:
		pass

	ctypedef struct squig_batch_t:
		pass

	ctypedef struct squig_rec_t:
		char* read_id;
		uint32_t read_group;
		double digitisation;
		double offset;
		double range;
		double sampling_rate;
		uint64_t len_raw_signal;
		int16_t* raw_signal;

		char *channel_number;
		double median_before;
		int32_t read_number;
		uint8_t start_mux;
		uint64_t start_time;
		uint8_t end_reason;

		pass

	ctypedef struct squig_opt_t:
		int32_t avg_rlen;
		int64_t random_seed;
		int32_t num_threads;
		int32_t batch_size;
		char *output_blow5;
		char *profile;
		pass


	# squigulator interface
	void squig_init_opt(squig_opt_t *opts);
	squig_t *squig_init(char *ref_fasta, squig_opt_t *opts);
	void squig_free(squig_t *sq);
	squig_batch_t* squig_sim_batch(squig_t *sq);
	void squig_free_batch(squig_batch_t *batch);
	squig_rec_t *squig_get_rec(squig_t *sq, squig_batch_t *batch, int i);

