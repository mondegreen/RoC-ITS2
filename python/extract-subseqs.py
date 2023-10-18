import sys
import os
import argparse
import time
import tempfile
import seaborn as sns
import matplotlib.pyplot as plt
from math import log
from statistics import mean, stdev

import re
from collections import defaultdict

from Bio import SearchIO, SeqIO, SeqRecord, Seq, AlignIO
from Bio.Align import AlignInfo
from Bio.Align.Applications import MafftCommandline
from io import StringIO

from dataclasses import dataclass

if sys.version_info[0] < 3:
    raise Exception("Python 3 or a more recent version is required.")

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--project', required=True, help='Basename for output files and directory name for subreads')
parser.add_argument('--pipeline', action='store_true', help='This script is being called as part of the pipeline')
parser.add_argument('--hmm', help='HMMSearch results')
parser.add_argument('--min-insert', default=1500, type=int, help='minimum insert size')
parser.add_argument('--max-insert', default=3500, type=int, help='maximum insert size')
parser.add_argument('--min-coverage', default=5, type=int,
                    help='minimum number of inserts required to proceed')
parser.add_argument('--max-coverage', default=100, type=int,
                    help='maximum number of inserts allowed to proceed')
parser.add_argument('--min-barcode-coverage', default=3, type=int,
                    help='minimun number of each barcode required to proceed')
parser.add_argument('--barcode-consensus-cutoff', default=.65, type=float,
                     help='minimum portion of read agreement to call base in barcode')
parser.add_argument('--reads', help='fastq file with reads')
parser.add_argument('--benchmark', action='store_true', help='Print the time script starts and ends to the project log file')

args = parser.parse_args()

# if pipeline isn't supplied, we need hmm and reads
if args.pipeline is None and args.hmm is None:
    parser.error('You must either supply an hmm file with --hmm or run as part of a pipeline --pipeline')
if args.pipeline is None and args.reads is None:
    parser.error('You must either supply a fasta file with --reads or run as part of a pipeline --pipeline')

# if a pipeline, we put it in a subdirectory
if args.pipeline:
    args.benchmark = True
    args.project = f'{args.project}/{args.project}'
    args.hmm = f'{args.project}-tbl.out'
    args.reads = f'{args.project}.fq'
    if os.path.exists(f'{args.project}-extract-subseqs.benchmark'):
        quit()

class Read:
    '''A collection of Subreads'''

    reads = []

    MIN_INSERT = args.min_insert
    MAX_INSERT = args.max_insert
    MIN_BARCODE_CON = args.barcode_consensus_cutoff
    MIN_COVERAGE = args.min_coverage
    MAX_COVERAGE = args.max_coverage
    MIN_BARCODE_COV = args.min_barcode_coverage

    def __init__(self, hit, is_reversed):
        self.hit = hit
        self.name = hit.id
        self.subreads = []
        self.subread_count = 0
        self.subreads_filtered = 0
        self.insert_lengths = []
        self.barcodes = { 5 : None, 3 : None }
        self.barcode_coverage = { 5 : 0, 3 : 0 }
        self.insert_length_mean = None
        self.insert_length_stdev = None
        self.passed = True
        self.fail_step = ''
        self.align = {}
        self.is_reversed = is_reversed

        # retrieve our sequence
        if self.name not in SEQS:
            raise Exception(f'{self.name} not in our fasta reads')
        self.seq = SEQS[self.name]
        self.length = len(self.seq)

        # are we reverse complemented?
        self.rev_ratio = self.__countReversed()
        if self.rev_ratio > .25:
            self.passed = False
            self.fail_step = 'Mixed Insert Orientation'
        self.__processHSPs()
        Read.reads.append(self)

    def __countReversed(self):
        '''Count the portion of HSPs that are reversed.'''
        hsps = 0
        rev_hsps = 0
        for hsp in self.hit.hsps:
            hsps += 1
            if hsp.is_reversed:
                rev_hsps += 1
        return rev_hsps / hsps

    def __processHSPs(self):
        '''process the read from left to right processing subreads'''

        start = 0
        end = self.length
        for hsp in self.hit.hsps:
            # do not split on reversed hsps, but tag current subread
            if hsp.is_reversed:
                continue
            else:
                sub = Subread(self, start, hsp)
                self.subreads.append(sub)
                self.subread_count += 1
                start = sub.end
        
    def __align_barcodes(self):
        '''align both the 5' and 3' barcode for the read'''
        for barcode in [5, 3]:
            # create fasta file for input
            seqs = []
            for s in self.subreads:
                if s.barcodes[barcode] is None:
                    continue
                record = SeqRecord.SeqRecord(
                    #Seq.Seq(s.barcodes[barcode]),
                    s.barcodes[barcode],
                    id=s.name,
                    description='')
                seqs.append(record)

            # if we don't have enough barcodes, skip this part
            self.barcode_coverage[barcode] = len(seqs)
            if self.barcode_coverage[barcode] < Read.MIN_BARCODE_COV:
                if self.passed:
                    self.passed = False
                    self.fail_step = 'Minimum Barcode Coverage'
                continue

            fp = tempfile.NamedTemporaryFile(mode = 'w')
            jnk = SeqIO.write(seqs, fp, 'fasta')
            fp.flush()

            # run MAFFT
            mafft_cline = MafftCommandline(input=fp.name, thread=4)
            stdout, stderr = mafft_cline()
            self.align[barcode] = AlignIO.read(StringIO(stdout), 'fasta')

            # build consensus
            ai = AlignInfo.SummaryInfo(self.align[barcode])
            self.barcodes[barcode] = str(ai.gap_consensus(threshold=Read.MIN_BARCODE_CON)).replace('-','')
            if len(self.barcodes[barcode]) > 5 and self.passed:
                self.passed = False
                self.fail_step = 'Barcode Mismatch'

    @classmethod
    def align_barcodes(cls):
        '''class method to align barcode for all the reads that haven't been excluded yet'''
        for r in [r for r in Read.reads if r.passed]:
            r.__align_barcodes()

    @classmethod
    def build_statistics(cls):
        for read in Read.reads:
            for s in read.subreads:
                # collect insert lengths and barcodes
                read.insert_lengths.append(s.insert_length)
            # insert_length distribution - exclude first insert
            if read.subread_count > 1:
                read.insert_length_mean = mean(read.insert_lengths[1:])
            if read.subread_count > 2:
                read.insert_length_stdev = stdev(read.insert_lengths[1:])

    @classmethod
    def plot_insert_distribution(cls,name):
        subreads = []
        for read in Read.reads:
            for sr in read.subreads:
                i_len = sr.insert_length
                if i_len < 1:
                    i_len = 1
                subreads.append(i_len)
        sns_hist = sns.histplot(subreads, log_scale = 10)
        fig = sns_hist.get_figure()
        fig.savefig(name)
        # matplotlib is garbage at the OOP, so cleanup
        fig.clf()

    @classmethod
    def plot_read_length_histogram(cls,name):
        reads = []
        for read in Read.reads:
            reads.append(read.length)
        sns_hist = sns.histplot(reads, log_scale = 10)
        fig = sns_hist.get_figure()
        fig.savefig(name)
        # matplotlib is garbage at the OOP, so cleanup
        fig.clf()

    @classmethod
    def plot_hsp_histogram(cls,name):
        reads = []
        for read in Read.reads:
            hsps = read.hit.hsps
            reads.append(len(hsps))
        sns_hist = sns.histplot(reads)
        fig = sns_hist.get_figure()
        fig.savefig(name)
        # matplotlib is garbage at the OOP, so cleanup
        fig.clf()

    @classmethod
    def plot_subread_histogram(cls,name):
        reads = []
        for read in Read.reads:
            reads.append(read.subread_count)
        sns_hist = sns.histplot(reads)
        fig = sns_hist.get_figure()
        fig.savefig(name)
        # matplotlib is garbage at the OOP, so cleanup
        fig.clf()

    @classmethod
    def plot_status_piechart(cls,name):
        status = defaultdict(int)
        seen = 0
        for read in Read.reads:
            seen += 1
            if read.passed == True:
                status['Pass QC'] += 1
            else:
                status[read.fail_step] += 1
        status['No HSPs'] = len(SEQS) - seen
        labels = []
        counts = []
        keys = list(status.keys())
        keys.sort()
        for x in keys:
            labels.append(x)
            counts.append(status[x])
        plt.pie(counts, labels = labels)
        plt.savefig(name)
        fig = plt.gcf()
        fig.clf()
            

    @classmethod
    def write_subread_report(cls, file):
        '''Write subread information to file'''
        f = open(file, 'w+')
        for r in Read.reads:
            for s in r.subreads:
                f.write(str(s)+'\n')
        f.close()

    @classmethod
    def write_read_report(cls, file):
        '''Write read information to file'''
        headers = ['Name', 'Read Length', 'Subread Count', '# of Subreads Filtered', 'Mean Insert Length',
                   'Insert Length SD', 'Passed QC?', 'Fail Step', '5\' Barcode Coverage',
                   '5\' Barcode', '3\' Barcode Coverage', '3\' Barcode', 'Reverse Oriented?', 'Portion of Wrongly Oriented HSPs']
        f = open(file, 'w+')
        f.write('\t'.join(headers) + '\n')
        for r in Read.reads:
            f.write(str(r)+'\n')
        f.close()

    def __str__(self):
        '''string override'''
        values = [self.name,
                  self.length,
                  self.subread_count,
                  self.subreads_filtered,
                  self.insert_length_mean,
                  self.insert_length_stdev,
                  self.passed,
                  self.fail_step,
                  self.barcode_coverage[5],
                  self.barcodes[5],
                  self.barcode_coverage[3],
                  self.barcodes[3],
                  self.is_reversed,
                  self.rev_ratio
                  ]

        return '\t'.join([f'{x}' for x in values])

    @classmethod
    def count_subreads(cls):
        counter = 0
        for r in Read.reads:
            for s in r.subreads:
                counter += 1
        return counter

    @classmethod
    def write_subreads(cls):

        for r in [r for r in Read.reads if r.passed]:
            sub_i = min(len(r.subreads), 5)
            index = r.name[0:2]
            path = args.project + '/' + index
            if not os.path.exists(path):
                os.makedirs(path)
            path += f'/{r.name}'
            records = []
            for s in r.subreads:
                record = SeqRecord.SeqRecord(
                    Seq.Seq(str(r.seq.seq[s.start : s.hsp.hit_start])),
                    id=s.name,
                    description=(f'/start={s.start} /end={s.end} /length={s.insert_length}'
                                 f' /its_barcode={s.barcodes[5]} /16s_barcode={s.barcodes[3]}'))
                record.letter_annotations['phred_quality'] = r.seq.letter_annotations['phred_quality'][s.start : s.hsp.hit_start]
                records.append(record)

            fq = open(path+'-subreads.fq', '+w')
            SeqIO.write(records, fq, 'fastq')
            fq.close()
            fa = open(path+'-subreads.fa', '+w')
            SeqIO.write(records, fa, 'fasta')
            fa.close()

    @classmethod
    def length_filter(cls):
        '''Remove subreads below the min or above the max across all reads'''
        for r in Read.reads:
            kept = [i for i in r.subreads if ( i.insert_length <= Read.MAX_INSERT
                                               and i.insert_length >= Read.MIN_INSERT )]

            kept_count = len(kept)
            r.subreads_filtered = r.subread_count - kept_count
            r.subread_count = kept_count
            r.subreads = kept

            # see if we failed
            if r.subread_count < Read.MIN_COVERAGE and ( r.passed or r.fail_step == 'Mixed Insert Orientation' ):
                r.passed = False
                r.fail_step = 'Minimum Coverage'

            if r.subread_count > Read.MAX_COVERAGE and ( r.passed or r.fail_step == 'Mixed Insert Orientation' ):
                r.passed = False
                r.fail_step = 'Maximum Coverage'

class Subread:
    '''A detected subread - separated by joint sequence'''

    def __init__(self, read, start, hsp):
        self.read = read
        self.start = start
        self.end = hsp.hit_end
        self.name = f'{read.name}-{self.start:09}-{self.end:09}'
        self.length = self.end - self.start
        self.insert_length = hsp.hit_start - self.start
        self.qc = None 
        self.hsp = hsp
        self.hit_length = hsp.query_end - hsp.query_start

        self.barcodes = { 5 : None, 3 : None }
        
        # check for backwards joint
        # shouldn't be any of these left
        if hsp.is_reversed:
            print('next hit is backwards')

        # check for 5' barcode
        if hsp.query_start == 0 and hsp.hit_start >= 5:
            self.barcodes[5] = read.seq.seq[hsp.hit_start - 5   : hsp.hit_start]

        # check for 3' barcode
        if hsp.query_end == 335 and hsp.hit_end + 5 <= read.length:
            self.barcodes[3] = read.seq.seq[hsp.hit_end : hsp.hit_end+5]

    def __str__(self):

        values = [self.read.name,
                  self.start,
                  self.end,
                  self.insert_length,
                  self.hit_length,
                  self.hsp.evalue,
                  self.barcodes[5],
                  self.barcodes[3]
                  ]

        return '\t'.join([f'{x}' for x in values])

@dataclass
class HSP:
    '''An HSP or envelope in the nhmmer'''
    query_start : int
    query_end : int
    hit_start : int
    hit_end : int
    evalue : str
    is_reversed : bool = False

    def __post_init__(self):
        if self.hit_end < self.hit_start:
            tmp = self.hit_end
            self.hit_end = self.hit_start
            self.hit_start = tmp
            self.is_reversed = True

class Hit:
    '''A Hit with nhmmer'''

    def __init__(self, id, hsps):
        self.id = id
        self.hsps = hsps

def parse_reads(file):
    return {v.id:v for v in SeqIO.parse(file, 'fastq')}

def parse_hmm(file):
    '''Parse the HMM tbl.out file and determine which reads are reverse complemented'''
    re_tab = re.compile('\s+')
    hsp_map = defaultdict(list)
    ori_map = defaultdict(lambda: {'+' : 0, '-' : 0 })

    # read in all the lines and save them and record orientation
    with open(file) as f:
        for line in f:
            if line[0] == '#':
                continue
            fields = re_tab.split(line)
            ori_map[fields[0]][fields[11]] += 1
            hsp_map[fields[0]].append([int(fields[4]), int(fields[5]), int(fields[6]),
                                       int(fields[7]), fields[11], fields[12]])

    # reverse reads where appropriate
    for read,cnt in ori_map.items():
        is_reversed = False

        if cnt['-'] > cnt['+']:
            is_reversed = True
            seq_len = len(SEQS[read])

            # rc the read
            SEQS[read] = SEQS[read].reverse_complement()
          
            # edit the hsp information for these
            for hsp in hsp_map[read]:
                hsp[2] = seq_len - hsp[2] + 1
                hsp[3] = seq_len - hsp[3] + 1
                hsp[4] = '+' if '-' else '-'

        hsps = []

        # now load the reads
        for x in hsp_map[read]:
            hsp = HSP(x[0]-1, x[1], x[2]-1, x[3], x[5])
            hsps.append(hsp)
        hsps.sort(key=lambda x: x.hit_start)
        hit = Hit(read, hsps)
        Read(hit, is_reversed)

def parse_hmm2(file):
    '''Instatiate Read and Subread objects'''
    re_tab = re.compile("\s+")
    hit_map = defaultdict(list)

    with open(file) as f:
        for line in f:
            if line[0] == '#':
                continue
            fields = re_tab.split(line)
            hsp = HSP(int(fields[4])-1, int(fields[5]), int(fields[6])-1, int(fields[7]), fields[12])            
            hit_map[fields[0]].append(hsp)

    for id, hsps in hit_map.items():
        #print(f'{id}')
        #print(hsps)
        hsps.sort(key=lambda x: x.hit_start)
        #print(hsps)
        hit = Hit(id, hsps)
        Read(hit)

            

#    [Read(hit) for hit in res]

SEQS = parse_reads(args.reads)

if __name__ == '__main__':

    if args.benchmark:
        start_time = time.time()

    parse_hmm(args.hmm)
    dir_name = args.project

    # create the dir if it doesnt exist
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)

    # write subread report
    Read.write_subread_report(dir_name+'-subread-list.tsv')

    # print raw subread length distribution
    Read.plot_insert_distribution(dir_name+'-subread-length-hist.png')

    # apply length and minimum coverage filter
    Read.length_filter()

    Read.plot_insert_distribution(dir_name+'-filtered-subread-length-hist.png')

    # align barcodes and filter on consensus quality and barcode count
    Read.align_barcodes()

    # filter on barcode consensus
    Read.build_statistics()

    # build read report
    Read.write_read_report(dir_name+'-read-list.tsv')

    # write fastx output
    Read.write_subreads()

    Read.plot_read_length_histogram(dir_name+'-read-length-hist.png')

    # plot status piechart
    Read.plot_status_piechart(dir_name+'-read-status-piechart.png')
    Read.plot_hsp_histogram(dir_name+'-hsp-count-hist.png')
    Read.plot_subread_histogram(dir_name+'-subread-count-hist.png')

    # build summary

    # print summary

    # finish
    if args.benchmark:
        end_time = time.time()
        with open(f'{args.project}-extract-subseqs.benchmark', 'a') as bf:
            bf.write(f'{start_time}\n')
            bf.write(f'{end_time}\n')

