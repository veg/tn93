#! /usr/bin/env python3
'''
Calculate all pairwise distances between sequences from the same individual
'''

# imports
from gzip import open as gopen
from subprocess import check_output, DEVNULL, PIPE
from sys import argv
import argparse

# parse user args
def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', required=False, type=str, default='stdin', help="Input FASTA")
    args = parser.parse_args()
    return args

# load FASTA sequences from a file
def load_fasta(fn):
    # set up input file stream
    if fn.strip().lower() == 'stdin':        # standard input
        from sys import stdin as f
    elif fn.strip().lower().endswith('.gz'): # gzip-compressed file
        f = gopen(fn, 'rt')
    else:                                    # uncompressed file
        f = open(fn, 'r')

    # parse FASTA
    seqs = dict(); curr_ID = None; curr_seq = ''
    for line in f:
        l = line.strip()
        if len(l) == 0: # skip empty lines
            continue
        if l[0] == '>': # header line
            if len(curr_seq) != 0:
                seqs[curr_ID] = curr_seq
            curr_ID = l; curr_seq = ''
        else:           # sequence line
            curr_seq += l
    if len(curr_seq) != 0:
        seqs[curr_ID] = curr_seq
    return seqs

# return dictionary where seqs_per_person[person] = set of sequence IDs from `person`
# this is currently a dummy implementation that clones each sequence 10 times (AND MODIFIES SEQS ACCORDINGLY!!!)
# TODO in practice, replace this with a function that actually loads the true person-to-sequences mapping
def load_seqs_per_person_dummy(seqs):
    seqs_per_person = dict(); orig_keys = set(seqs.keys())
    for k in orig_keys:
        seqs_per_person[k] = {k} # first copy is the original
        for i in range(9):       # create 9 more copies
            k2 = '%s.CLONE.%d' % (k,i); seqs_per_person[k].add(k2); seqs[k2] = seqs[k]
    return seqs_per_person

# calculate all pairwise distances between sequences from the same individual (this is the main logic of this script)
# currently not parallelizing across people because tn93 parallelizes pairwise distance calculations,
# but we could easily parallelize across people as well
def calc_intrahost_dists(seqs, seqs_per_person):
    for person_num, person in enumerate(seqs_per_person):
        curr_fasta = '\n'.join('%s\n%s' % (k,seqs[k]) for k in seqs_per_person[person])
        o = check_output(['tn93', '-t', '0.03', '-l', '1'], stderr=DEVNULL, input=curr_fasta, encoding='ascii')
        if person_num == 0:
            print(o.strip())
        else:
            for l in o.splitlines()[1:]:
                print(l.strip())

# run tool
if __name__ == "__main__":
    args = parse_args()
    seqs = load_fasta(args.input)
    seqs_per_person = load_seqs_per_person_dummy(seqs) # TODO replace with true way to determine which seqs come from the same person
    calc_intrahost_dists(seqs, seqs_per_person)
