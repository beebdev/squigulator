#!/usr/bin/env python3

# python3 setup.py install
# python3 python/example.py

import argparse
import sys
import numpy as np
import pysquigulator

def main():

    reference = 'test/nCoV-2019.reference.fasta'

    pysq = pysquigulator.sim(reference, profile='dna-r9-prom', num_threads=8, batch_size=100, output_blow5='test.blow5')

    recs = pysq.simulate_batch()

    for read in recs:
        print("read_id:", read['read_id'])
        print("read_group:", read['read_group'])
        print("digitisation:", read['digitisation'])
        print("offset:", read['offset'])
        print("range:", read['range'])
        print("sampling_rate:", read['sampling_rate'])
        print("len_raw_signal:", read['len_raw_signal'])
        print("signal:", read['signal'][:10])
        print("signal:", read['signal'][:10])
        print("================================")



if __name__ == '__main__':
    main()

