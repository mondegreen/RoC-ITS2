from collections import Counter
from Bio import SeqIO, SeqRecord, Seq
from Bio.SeqIO import FastaIO
import sys
import argparse
import os
import time

if sys.version_info[0] < 3:
    raise Exception("Python 3 or a more recent version is required.")

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--msa', help='MSA containing processed subreads')
parser.add_argument('--benchmark', action='store_true',
                    help='Print the time script starts and ends to the project log file')
parser.add_argument('--pipeline', action='store_true',
                    help='Pipeline base path [project/prefix/subread]')
parser.add_argument('--base', help='Output base. Required for pipeline or benchmark')
parser.add_argument('--read', required=True, help='Name of the read')
args = parser.parse_args()

# configure based on if part of a pipeline
if args.pipeline:
    if args.base is None:
        parser.error('You must provide --base with --pipeline')
    args.benchmark = True
    args.msa = f'{args.base}-smushed-improved.fa'
    if os.path.exists(f'{args.base}-make-consensus.benchmark'):
        quit()
    output = open(f'{args.base}-consensus.fa', 'w')

else:
    if args.msa is None:
        parser.error('You must provide either --msa or --pipeline')
    if args.benchmark and args.base is None:
        parser.error('You must provide --base with --benchmark')
    output = sys.stdout

def score_seq(seq, con):
    '''return portion of matches with supplied consensus'''
    match = 0
    for i,val in enumerate(seq):
        if val == con[i]:
            match += 1
        elif val == '-':
            match += .5
    return match/len(seq)

def make_consensus(seqs):
    '''return a consensus sequence for the array of sequences'''
    con = ''
    for i in range(len(seqs[0])):
        base = Counter()
        for s in seqs:
            base.update(s[i])
        con += base.most_common(1)[0][0]
    return con

def score_alignment(seqs, con):
    (good,all) = (0,0)
    for s in seqs:
        for i,val in enumerate(s):
            all += 1
            if val == con[i]:
                good += 1
            
    return good/all

def do_iteration(seqs):
    '''make an iteration of consensus seqs'''
    con = make_consensus(seqs)
    seqs2 = [s for s in seqs if score_seq(s,con) > .65]
    return (seqs2, con)

if __name__ == '__main__':

    if args.benchmark:
        start_time = time.time()

    seqs = []
    for record in SeqIO.parse(args.msa, 'fasta'):
        seqs.append(record.seq.upper())

    orig_len = len(seqs)

    (s2,consensus) = do_iteration(seqs)

    while len(s2) != len(seqs) and len(s2) > 3:
        seqs = s2
        (s2, consensus) = do_iteration(seqs)


    record = SeqRecord.SeqRecord(
        Seq.Seq(consensus.replace('-','')),
        id=args.read,
        description=f'/subreads={len(s2)} /consensus={score_alignment(s2,consensus):.4f} /removed-subreads={orig_len-len(s2)}'
)

    FastaIO.FastaWriter(output, wrap=None).write_file([record])

    # finish
    if args.benchmark:
        end_time = time.time()
        with open(f'{args.base}-make-consensus.benchmark', 'a') as bf:
            bf.write(f'{start_time}\n')
            bf.write(f'{end_time}\n')

