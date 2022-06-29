# Generating Consensus Reads 

For these examples, we will be working with the fastq output of the guppy basecaller - in this case __barcode06.fq__.

## Prerequisits

### Python Libraries

- biopython : https://github.com/biopython/biopython
- seaborn : https://seaborn.pydata.org/

### External programs 

 - MAFFT : https://mafft.cbrc.jp/alignment/software/
 - CONSENT-correct : https://github.com/morispi/CONSENT
 - PRANK : https://github.com/morispi/CONSENT
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
The output will create a fasta and fastq file for the subreads in the following format:

- 01/0152ca78-23d6-4944-a0a0-ca4fda370bd8-subreads.fa
- 01/0152ca78-23d6-4944-a0a0-ca4fda370bd8-subreads.fq

Next we use CONSENT-correct to use the overlapping nature of the subreads to correct each subread.

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

Below is a bash subroutine that can be used to distribute these steps for all reads across N processors on a single machine that can be used in lieu of a more robust workflow/pipeline system. 

```bash
task(){
  BASE=`basename $1 -subreads.fq`
  DIR=`dirname $1`

  CONSENT-correct --in ${DIR}/${BASE}-subreads.fq \
    --out ${DIR}/${BASE}-corrected.fa --type ONT

  prank -F -o=${DIR}/${BASE}-corrected \
    -d=${DIR}/${BASE}-corrected.fa > /dev/null

  python smush-msa.py --msa ${DIR}/${BASE}-corrected.best.fas \
    > ${DIR}/${BASE}-corrected-smushed.best.fas

  python improve-msa.py --msa ${DIR}/${BASE}-corrected-smushed.best.fas \
    --max-cores 2 --word-size 12 > ${DIR}/${BASE}-smushed-improved.fa

  echo "working on ${DIR}/${BASE}-subreads"
  python make-consensus.py --msa ${DIR}/${BASE}-smushed-improved.fa \
    --read ${BASE} > ${DIR}/${BASE}-consensus.fa
}

N=25

for FILE in reads5/*/*-subreads.fq
  do
    ((i=i%N));
    ((i++==0)) && wait
    task $FILE&
  done
  cat reads5/*/*-consensus.fa > consensi-5plus-subreads.fa
done

```
