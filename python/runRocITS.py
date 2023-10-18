'''Try to run the entire Roc-ITS pipeline on a dataset.'''
import sys
import argparse
import os
import subprocess
import time
import shutil
import glob
import re

import multiprocessing as mp

if sys.version_info[0] < 3:
    raise Exception('Python 3 or a more recent version required.')

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--reads', required=True, help='Fastq file containing Nanopore reads')
parser.add_argument('--project', required=True, help='Project name used for output')
parser.add_argument('--joint-hmm', required=True, help='HMM for the joint sequence')
parser.add_argument('--max-cores', type=int, default=4, help='Maximum number of cores to use')
args = parser.parse_args()

base = f'{args.project}/{args.project}'

#
# Worker Class
#
class Worker(mp.Process):
    def __init__(self, task_queue: mp.Queue):
        super(Worker,self).__init__()
        self.task_queue = task_queue

    def run(self):
        for (function, *args) in iter(self.task_queue.get, None):
            print(f'Running: {function.__name__}({*args,})')
            
            function(*args)
            self.task_queue.task_done()

#
# Functions to run in the queue
#

def run_vsearch(task_queue: mp.Queue) -> None:
    '''Convert the Nanopore Fastq file to a 0-width fasta file using vsearch'''
    bf = f'{base}-vsearch.benchmark'
    if not os.path.exists(bf):
        start_time = time.time()
        shutil.copy2(args.reads, f'{base}.fq')
        fout = f'{base}.fa'
        subprocess.run(['vsearch', '--fastq_filter', args.reads, '--fastq_qmax',
                        '90', '--fasta_width', '0', '--fastaout', fout], check=True)
        write_benchmark(start_time, time.time(), bf)
    task_queue.put((run_nhmmer, task_queue))

def run_nhmmer(task_queue: mp.Queue) -> None:
    '''Identify joint sequences in reads using nhmmer'''
    bf = f'{base}-nhmmer.benchmark'
    if not os.path.exists(bf):
        start_time = time.time()
        f = open(f'{base}-nhmm.out', 'w')
        subprocess.run(['nhmmer', '--cpu', str(args.max_cores), '--tblout', 
                        f'{base}-tbl.out', '--noali', args.joint_hmm,
                        f'{base}.fa'], stdout=f, check=True)
        write_benchmark(start_time, time.time(), bf)
    task_queue.put((run_extract_subseqs, task_queue))

def run_extract_subseqs(task_queue: mp.Queue) -> None:
    '''Split reads into subseqs'''
    subprocess.run(['python', 'extract-subseqs.py', '--project', args.project, '--pipeline'],
                   check=True)
    for x in glob.glob(f'{base}/*/*-subreads.fa'):        
        subread = x.split('/')[-1]
        subread = re.sub('-subreads.fa', '', subread)
        task_queue.put((run_consent, subread, task_queue))

def run_consent(subread: str, task_queue: mp.Queue) -> None:
    '''Use CONSENT-correct to polish subreads'''
    base2 = f'{base}/{subread[:2]}/{subread}'
    bf = f'{base2}-CONSENT-correct.benchmark'
    if not os.path.exists(bf):
        start_time = time.time()
        fout = f'{base2}-corrected.fa'
        subprocess.run(['CONSENT-correct', '--in', f'{base2}-subreads.fq', '--out',
                        fout, '--type', 'ONT'], check=True)
        write_benchmark(start_time, time.time(), bf)
    task_queue.put((run_prank, subread, task_queue))

def run_prank(subread: str, task_queue: mp.Queue) -> None:
    '''Run prank to align the subreads for this read'''
    base2 = f'{base}/{subread[:2]}/{subread}'
    bf = f'{base2}-prank.benchmark'
    if not os.path.exists(bf):
        start_time = time.time()
        fout = f'{base2}-corrected'
        subprocess.run(['prank', '-F', f'-o={fout}', f'-d={base2}-corrected.fa'],
                       stdout=subprocess.DEVNULL, check=True)
        write_benchmark(start_time, time.time(), bf)
    task_queue.put((run_smush_msa, subread, task_queue))

def run_smush_msa(subread: str, task_queue: mp.Queue) -> None:
    '''Remove gaps and other issues from MSA'''
    base2 = f'{base}/{subread[:2]}/{subread}'
    subprocess.run(['python', 'smush-msa.py', '--pipeline', 
                    '--base', base2], stdout=subprocess.DEVNULL, check=True)
    task_queue.put((run_improve_msa, subread, task_queue))

def run_improve_msa(subread: str, task_queue: mp.Queue) -> None:
    '''Improve the MSA'''
    base2 = f'{base}/{subread[:2]}/{subread}'
    subprocess.run(['python', 'improve-msa.py', '--pipeline', '--base', base2],
                   stdout=subprocess.DEVNULL, check=True)
    task_queue.put((run_make_consensus, subread, task_queue))

def run_make_consensus(subread: str, task_queue: mp.Queue) -> None:
    '''Generate a consensus from the improved MSA'''
    base2 = f'{base}/{subread[:2]}/{subread}'
    subprocess.run(['python', 'make-consensus.py', '--pipeline', '--base', base2,
                    '--read', subread], stdout=subprocess.DEVNULL, check=True)

#
# helper functions
#
def make_project_dir() -> None:
    '''Check for existing project dir, create one otherwise'''
    if os.path.exists(args.project):
        if not os.path.isdir(args.project):
            raise Exception(f'{args.project} exists and is not a directory')
    else:
        os.mkdir(args.project)

def write_benchmark(start_time: float, end_time: float, file: str) -> None:
    '''Write a benchmark file'''
    with open(file, 'w') as bf:
        bf.write(f'{start_time}\n')
        bf.write(f'{end_time}\n')

if __name__ == '__main__':

    workers = []
    mgr = mp.Manager()
    queue = mgr.Queue()

    for i in range(args.max_cores):
        worker = Worker(queue)
        workers.append(worker)
        worker.start()

    make_project_dir()

    # start us off with vsearch
    queue.put((run_vsearch, queue))

    # Block until all items in queue are done
    queue.join()

    for _ in workers:
        queue.put(None)

    # Block until workers complete and join main process
    for worker in workers:
        worker.join()

    # Concatenate consensi
    with open(f'{args.project}/{args.project}-consensi.fa', 'w') as outfile:
        for x in glob.glob(f'{args.project}/{args.project}/*/*-consensus.fa'):
            with open(x) as infile:
                outfile.write(infile.read())

    print('Program Completed.')
