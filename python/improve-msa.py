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
import time

if sys.version_info[0] < 3:
    raise Exception("Python 3 or a more recent version is required.")

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--msa', help='MSA containing processed subreads')
parser.add_argument('--max-cores', type=int, default=4, help='Maximum number of cores to use')
parser.add_argument('--scan', action='store_true', help='Print number of regions and exit')
parser.add_argument('--benchmark', action='store_true',
                    help='Print the time script starts and ends to the project log file')
parser.add_argument('--pipeline', action='store_true',
                    help='Pipeline base path [project/prefix/subread]')
parser.add_argument('--base', help='Output base. Required for pipeline or benchmark')
parser.add_argument('--word-size', type=int, default=12,
                    help='Word size to use')
args = parser.parse_args()

# configure based on if part of a pipeline
if args.pipeline:
    if args.base is None:
        parser.error('You must provide --base with --pipeline')
    args.max_cores = 1
    args.benchmark = True
    args.msa = f'{args.base}-corrected-smushed.best.fas'
    if os.path.exists(f'{args.base}-improve-msa.benchmark'):
        quit()
    output = open(f'{args.base}-smushed-improved.fa', 'w')

else:
    if args.msa is None:
        parser.error('You must provide either --msa or --pipeline')
    if args.benchmark and args.base is None:
        parser.error('You must provide --base with --benchmark')
    output = sys.stdout

class Region:
    '''A region across the alignent that will be evaluated'''

    def __init__(self, begin, end):
        self.begin = begin
        self.end = end
        self.dir = tempfile.TemporaryDirectory()
        self.input = self.dir.name + '/input'
        self.output = self.dir.name + '/output'

def prepare_alignment(region, sub):
    '''Process a region - create fasta file as input for probcons'''
    records = []

    for s in sub:
        seq = ''
           
        for i in range(region.begin, region.end):
            chr = sub[s][i]
            if i < region.end-args.word_size :
                sub[s][i] = ''
                seq += chr
            else :
                # subsequences that have gaps in the right anchor region are by definition
                # rare but need to be replaced by N's. This is important as later non-dashes
                # in the right anchor are removed to when the new msa is inserted.
                if chr == '-':
                    seq += 'n'
                else:
                    seq += chr

        record = SeqRecord.SeqRecord(
            Seq.Seq(seq.replace('-','')),
            id=' '.join([s, str(region.begin), str(region.end)]),
            description=''
        )
        records.append(record)
    FastaIO.FastaWriter(region.input, wrap=None).write_file(records)

def process_alignment(region, sub):
    '''After a region has been processed, incorporate the results'''
    for record in SeqIO.parse(region.output, 'fasta'):
        s = record.id
        res = list(record.seq)
        removed = 0

        for k in reversed(range(len(res))):
            if wp.search(res[k]):
                res[k] = ''
                removed += 1
                if removed == args.word_size:
                    break

        sub[s][region.begin] = ''.join(res).lower()

def run_job(fin, fout):
    '''run the job on one of our workers'''
    f = open(fout, 'w')
    subprocess.run(['probcons', '-c', '0', fin], stdout=f, check=True)

if __name__ == '__main__':

    if args.benchmark:
        start_time = time.time()

    # bracket our sequences with 22 ns on each side
    GGGS = 'n' * 22
    ANCHOR = 'G' * args.word_size

    # a dictionary of our sub reads
    sub = {}

    # the length of our alignmet
    alnlen = 0

    # read in the fasta file
    for record in SeqIO.parse(args.msa, 'fasta'):
        seq = GGGS + record.seq.upper() + GGGS
        sub[record.id] = list(seq)
        # build a thing based on the length of sequences
        alnlen = len(seq)

    # template for good positions
    good = ''

    for i in range(alnlen):
        base = Counter()
        for s in sub.values():
            base.update(s[i])
        best = base.most_common(1)[0]
        subcount = len(sub)

        # I changed the gap cutoff from >= 10% to >10% to simplify the logic
        if ( best[0] == '-' or best[1]/subcount < 0.9 ) :
            good += '-'
        else:
            good += 'G'

    # find regions of gaps, defined as gap to the first sequence of word_size G bases
    # we will keep a list of region coordinates
    regions = []

    # i will count our way across the alignment
    start = args.word_size

    # process all the gaps
    while (i := good.find('-', start)) > 0:

        # want an upstream anchor, default end is end of alignment
        begin = i - args.word_size
        end = alnlen

        # if there is an anchor, start after it
        if ( j := good.find(ANCHOR, i)) > 0:
            end = j + args.word_size

        # save the region
        regions.append(Region(begin,end))

        # update start
        start = end - args.word_size

    # regular expression for matching word
    wp = re.compile('\w')

    if args.scan:
        print(f'Found {len(regions)} regions with word-size {args.word_size}')
        sys.exit()

    sys.stderr.write(f'We have {len(regions)} regions to investigate')

    # farm out probcons to a worker pool, output is written to a file in a temp dir
    with mp.Pool(processes=args.max_cores) as pool:
        results = []
        for reg in regions:
            prepare_alignment(reg, sub)
            results.append(pool.apply_async(run_job, (reg.input, reg.output)))
        for res in results:
            res.get()

    for reg in regions:
        process_alignment(reg, sub)

    records = []
    for s in sub:
        record = SeqRecord.SeqRecord(
            Seq.Seq(''.join(sub[s]).replace('n','-')),
            id=s,
            description=''
        )
        records.append(record)

    FastaIO.FastaWriter(output, wrap=None).write_file(records)

    # finish
    if args.benchmark:
        end_time = time.time()
        with open(f'{args.base}-improve-msa.benchmark', 'a') as bf:
            bf.write(f'{start_time}\n')
            bf.write(f'{end_time}\n')
