# Generating Consensus Reads 

For these examples, we will be working with the fastq output of the guppy basecaller - in this case __barcode06.fq__.

## Prerequisites

### Python Libraries

- biopython : https://github.com/biopython/biopython
- seaborn : https://seaborn.pydata.org/

### External programs 

 - MAFFT : https://mafft.cbrc.jp/alignment/software/ (no extensions are necessary)
 - CONSENT-correct : https://github.com/morispi/CONSENT/releases/tag/v2.2.2
 - PRANK : http://wasabiapp.org/software/prank/
 - PROBCONS : http://probcons.stanford.edu/

First, we must alter the fastq file so that the sequence and quality strings are each on a single line. Guppy generates quality scores outside of the range expected by vsearch by default, so we increase the maximum quality score to 90.

```bash
vsearch \
  --fastq_filter barcode06.fq \
  --fastq_qmax 90 \
  --fasta_width 0 \
  --fastaout barcode06.fa
```
Next we use hmmsearch to split the nanopore reads on the library joint sequence, effectively separating the full reads into subreads.

```bash
hmmsearch \
 --cpu 15 \
 --domtblout barcode06-dom.out \
 joint.hmm \
 barcode06.fa \
 barcode06-hmm.out
```
The following step performs quality control on the generated subreads. Subreads longer than --max-insert or shorter than --min-insert are excluded from further analysis. Then for each subread, the 5' and 3' barcodes from the adjoining joint sequences are extracted when processed and aligned via MAFFT. This alignment is then analyzed to identify and remove reads containing mixed subreads. 

Ths script generates a report listing the fate of each read and subread in the dataset, as well as up to 5 directories containing hashed details for all reads with one or more remaining subread. We recommend using only reads with 5 or more subreads (the reads5 directory), but reads with fewer subreads are reported and can be processed.

```bash
python skim-subseqs.py \
  --base barcode06 \
  --hmm barcode06-hmm.out \
  --reads barcode06.fq \
  --min-insert 1500 \
  --max-insert 3500 \
  --barcode-consensus-cutoff .65
```
This will result in a directory with the following structure:
```bash
filtered-subread-length-hist.png
read-list.tsv
reads1
reads2
reads3
reads4
reads5
subread-length-hist.png
subread-list.tsv
```
The png files show histograms describing the length of the sub-reads. The subread-list.tsv is a tab-delimited file describing each of the sub-reads with the following columns:
1) unique ID of nanopore read
2) begin coordinate of the sub-read
3) end coordinate of the sub-read
4) distance between start and end of sub-read
5) length of joint found
6) p-value of joint-match
7) 5' barcode (None if not found)
8) 3' barcode (None if not found)

Here is some example data:
```bash
ec64aa44-08ca-4858-b35f-acb8e06a9051    53321   55848   2194    335     3.2e-167        GTCGC   GGGTT
ec64aa44-08ca-4858-b35f-acb8e06a9051    55848   58383   2202    332     8.2e-168        None    AAGGT
ec64aa44-08ca-4858-b35f-acb8e06a9051    58383   60908   2188    335     3.1e-161        GTCGC   AAGGT
b7ab8f72-a2c3-4f42-b57f-159717d38c7e    0       1160    825     335     4.9e-169        CCAGC   GCAAT
b7ab8f72-a2c3-4f42-b57f-159717d38c7e    1160    3621    2126    335     4.5e-160        CCAAC   GCAAT
```

The read-list.tsv is a tab-delimited file with the following columns:
1) unique ID of nanopore read
2) something
3) something
4) average length of sub-reads
5) something
6) flag
7) something
8) 5' barcade consensus
9) something
10) 3' barcode consensus

Here is an example:
```bash
f548e84f-d909-4250-8c90-b958b4631b42    20      0       2178.8947368421054      13.457318276858818      True            19      tacag   19      gacaa
8872dcbe-68c5-4bea-bb34-67900b977b97    20      0       2099.5789473684213      18.524205754374716      True            18      acctg   18      gctga
d2f19bb4-1627-4323-9921-080e5e67f28a    0       52      None    None    False   Minimum Barcode Coverage        0       None    0       None
```

The reads[12345] directories contain the subreads with 1, 2, 3, 4, or 5+ subreads in them. Within each of the directories, the script will create a fasta and fastq file for the subreads in the following format (where the first two letters of the Nanopore read_id are used to create a series of sub-directories):

- 01/0152ca78-23d6-4944-a0a0-ca4fda370bd8-subreads.fa
- 01/0152ca78-23d6-4944-a0a0-ca4fda370bd8-subreads.fq

Each of these sub-reads needs to be processed to be aligned to generate a consensus sequence that corresponds to a RoC-ITS read. This is done through a series of processing steps (described in detail below) but that can be done in bulk and in parallel with the generateConsensi.sh bash script.

```bash
bash generateConsensi barcode06/reads5 10 output
```
The parameters are the directory containing the sub-read information (in this case the barcode06 sub-reads where 5+ subreads were found), the number of processors to use, and the output file name. This script is restartable in case it is terminated prematurely.

Within the script the following commands are used to fine-tune, align, and generate a consensus sequence. We use CONSENT-correct to use the overlapping nature of the subreads to correct each subread.

```bash
CONSENT-correct --in 0152ca78-23d6-4944-a0a0-ca4fda370bd8-subreads.fq\
  --type ONT \
  --out 0152ca78-23d6-4944-a0a0-ca4fda370bd8-corrected.fa
```

Then PRANK is used to further improve the subreads

```bash
prank -F -o=0152ca78-23d6-4944-a0a0-ca4fda370bd8-corrected \
  -d=0152ca78-23d6-4944-a0a0-ca4fda370bd8-subreads-corrected.fa > /dev/null
```
Next a python script is used to remove several types of non-informative gaps from the PRANK alignment, speeding up subsequent steps.

```bash
python smush-msa.py --msa 0152ca78-23d6-4944-a0a0-ca4fda370bd8-corrected-smushed-best.fas \
  > 0152ca78-23d6-4944-a0a0-ca4fda370bd8-corrected-smushed.best.fas
```
This step identifies highly conserved anchors of --word-size in the script and then uses PROBCONS to improve the alignment of regions between anchors.

```bash
python improve-msa.py --msa 0152ca78-23d6-4944-a0a0-ca4fda370bd8-corrected-smushed.best.fas \
  --max-cores 2 \
  --word-size 12 \
  > 0152ca78-23d6-4944-a0a0-ca4fda370bd8-smushed-improved.fa

Finally, the corrected subreads are used to generate a putative consensus sequence. All subreads are compared to this consensus, and any that fail to match the consensus at certain threshold are removed and a new consensus sequence is generated and the final consensus sequence is reported.

```bash
python make-consensus.py --msa 0152ca78-23d6-4944-a0a0-ca4fda370bd8-smushed-improved.fa \
  --read 0152ca78-23d6-4944-a0a0-ca4fda370bd8-smushed-improved.fa \
  > 0152ca78-23d6-4944-a0a0-ca4fda370bd8-consensus.fa
```
