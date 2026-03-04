# distutils: language = c
# cython: language_level=3
# cython: profile=True

import sys
# import time
# import logging
import copy
from libc.stdlib cimport malloc, free
from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free
from libc.string cimport strdup
cimport pysquigulator

# Import the Python-level symbols of numpy
import numpy as np
# Import the C-level symbols of numpy
cimport numpy as np
# Numpy must be initialized. When using numpy from C or Cython you must
# _always_ do that, or you will have segfaults
np.import_array()


cdef class sim:
    '''
    Creates a new pysquigulator object

    '''
    cdef squig_opt_t opt
    cdef squig_t *sq
    cdef char *reference
    cdef char *profile
    cdef char *output_blow5
    cdef squig_batch_t *batch

    cdef squig_rec_t *rec
    cdef np.npy_intp shape_seq[1]

    def __cinit__(self, reference, profile='dna-r9-prom', num_threads=8, batch_size=100, output_blow5='test.blow5'):
        '''
        C init
        '''
        self.sq = NULL
        self.batch = NULL

        self.reference = NULL
        self.profile = NULL
        self.output_blow5 = NULL

        self.reference = strdup(str.encode(reference))
        self.profile = strdup(str.encode(profile))
        self.output_blow5 = strdup(str.encode(output_blow5))

        squig_init_opt(&self.opt)

        self.opt.profile = self.profile
        self.opt.output_blow5 = self.output_blow5
        self.opt.num_threads = num_threads
        self.opt.batch_size = batch_size

        self.sq = squig_init(self.reference, &self.opt)

        if self.sq is NULL:
            self.logger.error("pysquigulator not initialised")


    def __init__(self, reference, profile='dna-r9-prom', num_threads=8, batch_size=100, output_blow5='test.blow5'):
        '''
        python init
        '''
        pass

    def __dealloc__(self):
        '''
        free memory
        '''

        if self.sq is not NULL:
            squig_free(self.sq)
        if self.reference is not NULL:
            free(self.reference)
        if self.profile is not NULL:
            free(self.profile)
        if self.output_blow5 is not NULL:
            free(self.output_blow5)



    def simulate_batch(self):
        '''
        process a batch of of signals
        '''
        row = {}

        self.batch = squig_sim_batch(self.sq)
        if self.batch is NULL:
            print("Failed to simulate batch")

        for i in range(self.opt.batch_size):
            row = {}
            self.rec = squig_get_rec(self.sq, self.batch, i)
            # todo check null

            if type(self.rec.read_id) is bytes:
                row['read_id'] = self.rec.read_id.decode()
            else:
                row['read_id'] = self.rec.read_id
            row['read_group'] = self.rec.read_group
            row['digitisation'] = self.rec.digitisation
            row['offset'] = self.rec.offset
            row['range'] = self.rec.range
            row['sampling_rate'] = self.rec.sampling_rate
            row['len_raw_signal'] = self.rec.len_raw_signal

            self.shape_seq[0] = <np.npy_intp> self.rec.len_raw_signal
            signal = copy.deepcopy(np.PyArray_SimpleNewFromData(1, self.shape_seq,
                        np.NPY_INT16, <void *> self.rec.raw_signal))
            np.PyArray_UpdateFlags(signal, signal.flags.num | np.NPY_ARRAY_OWNDATA)
            row['signal'] = signal

            if type(self.rec.channel_number) is bytes:
                row['channel_number'] = self.rec.channel_number.decode()
            else:
                row['channel_number'] = self.rec.channel_number

            row['median_before'] = self.rec.median_before
            row['read_number'] = self.rec.read_number
            row['start_mux'] = self.rec.start_mux
            row['start_time'] = self.rec.start_time

            yield row

        squig_free_batch(self.batch)
        self.batch = NULL
