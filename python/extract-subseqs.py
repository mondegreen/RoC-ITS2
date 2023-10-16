import sys
import os
import argparse
import tempfile
import seaborn as sns
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

parser = argparse.ArgumentParser()
parser.add_argument('--base', required=True, help='Basename for output files and directory name for subreads')
parser.add_argument('--hmm', required=True, help='HMMSearch results')
parser.add_argument('--min-insert', default=1500, type=int, help='minimum insert size')
parser.add_argument('--max-insert', default=3500, type=int, help='maximum insert size')
parser.add_argument('--min-coverage', default=5, type=int,
                     help='minimum number of inserts required to proceed (default = 5)')
parser.add_argument('--max-coverage', default=100, type=int,
                     help='maximum number of inserts required to proceed (default = 100)')
parser.add_argument('--min-barcode-coverage', default=3, type=int,
                    help='minimun number of each barcode required to proceed')
parser.add_argument('--barcode-consensus-cutoff', default=.65, type=float,
                     help='minimum portion of read agreement to call base in barcode')
parser.add_argument('--reads', required=True, help='fastq file with reads')

args = parser.parse_args()

class Read:
    '''A collection of Subreads'''

    reads = []

    MIN_INSERT = args.min_insert
    MAX_INSERT = args.max_insert
    MIN_BARCODE_CON = args.barcode_consensus_cutoff
    MIN_COVERAGE = args.min_coverage
    MAX_COVERAGE = args.max_coverage
    MIN_BARCODE_COV = args.min_barcode_coverage

    def __init__(self, hit):
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

        # retrieve our sequence
        if self.name not in SEQS:
            raise Exception(f'{self.name} not in our fasta reads')
        self.seq = SEQS[self.name]
        self.length = len(self.seq)

        self.__processHSPs()
        Read.reads.append(self)

    def __processHSPs(self):
        '''process the read from left to right processing subreads'''
        start = 0
        end = self.length
        for hsp in self.hit.hsps:
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
        f = open(file, 'w+')
        for r in Read.reads:
            f.write(str(r)+'\n')
        f.close()

    def __str__(self):
        '''string override'''
        values = [self.name,
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

        if not os.path.exists(args.base):
            os.makedirs(args.base)

        for r in [r for r in Read.reads if r.passed]:
            sub_i = min(len(r.subreads), 5)
            index = r.name[0:2]
            path = args.base + '/' + index
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
            if (r.subread_count < Read.MIN_COVERAGE or r.subread_count > Read.MAX_COVERAGE) and r.passed:
                r.passed = False
                r.fail_step = 'Minimum Coverage'

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

class Hit:
    '''A Hit with nhmmer'''

    def __init__(self, id, hsps):
        self.id = id
        self.hsps = hsps

def parse_reads(file):
    return {v.id:v for v in SeqIO.parse(file, 'fastq')}

def parse_hmm(file):

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
    parse_hmm(args.hmm)
    #dir_name = args.base + '-subreads'
    dir_name = args.base

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

    # build summary

    # print summary
