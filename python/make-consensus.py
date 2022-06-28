from collections import Counter
from Bio import SeqIO, SeqRecord, Seq
from Bio.SeqIO import FastaIO
import sys
import argparse

if sys.version_info[0] < 3:
    raise Exception("Python 3 or a more recent version is required.")

parser = argparse.ArgumentParser()
parser.add_argument('--msa', required=True, help='MSA containing processed subreads')
parser.add_argument('--read', required=True, help='Name of the read')
args = parser.parse_args()



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

FastaIO.FastaWriter(sys.stdout, wrap=None).write_file([record])
