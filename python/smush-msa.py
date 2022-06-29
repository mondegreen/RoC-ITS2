from collections import Counter
import os
import multiprocessing as mp
from Bio import SeqIO, SeqRecord, Seq
from Bio.SeqIO import FastaIO
import argparse
import sys
import tempfile
import subprocess
import re

if sys.version_info[0] < 3:
    raise Exception("Python 3 or a more recent version is required.")

parser = argparse.ArgumentParser()
parser.add_argument('--msa', required=True, help='MSA containing processed subreads')
args = parser.parse_args()

def check_candidacy(sub, i):
    cur = []
    call = None
    
    for s in sub:
        if sub[s][i] == '-':
            continue
        elif sub[s][i] == call:
            cur.append(s)
        elif call == None:
            call = sub[s][i]
            cur.append(s)
        else:
            call = 'X'

    if call == None:
        raise Exception(f'did not call candidate for position {i}')

    return call, cur

def test_pure(sub, call, i):
    for s in sub:
        if sub[s][i] == '-':
            continue
        elif sub[s][i] != call:
            return False
    return True

def try_merge(sub, call, cur, i):
    # need to make a loop so we can check over more than one position
    for j in range(i+1, len(sub[cur[0]])):

        # if a base exists in any of the positions, we are done
        for s in cur:
            if sub[s][j] != '-':
                return

        # is j a candidate for merging?
        if not test_pure(sub, call, j):
            continue

        # i think we're good column j only has gaps and our base
        # and all the positions in i have a gap in j
        for s in cur:
            sub[s][j] = call
        for s in sub:
            del sub[s][i]
        return

if __name__ == '__main__':

    # a dictionary of our sub reads
    sub = {}

    # the length of our alignmet
    alnlen = 0

    # read in the fasta file
    for record in SeqIO.parse(args.msa, 'fasta'):
        seq = record.seq.upper()
        sub[record.id] = list(seq)
        # build a thing based on the length of sequences
        alnlen = len(seq)

    # template for good positions
    good = ''

    # start with the most trivial case, adjacent columns
    # with only one base that is the same base

    # work from the end so we can collapse columns
    for i in range(alnlen-1,-1,-1):

        # skip the first column
        if i+1 == alnlen:
            continue

        call, cur = check_candidacy(sub, i)

        if call == 'X':
            continue

        try_merge(sub, call, cur, i)
        
    records = []
    for s in sub:
        record = SeqRecord.SeqRecord(
            Seq.Seq(''.join(sub[s]).replace('n','-')),
            id=s,
            description=''
        )
        records.append(record)

    FastaIO.FastaWriter(sys.stdout, wrap=None).write_file(records)

        
