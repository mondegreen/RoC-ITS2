# RoC-ITS2

## Prerequisites

### Perl packages
(be sure to update your PERL5LIB appropriately):<br>
* Carp<br>
* FileHandle<br>
* Getopt::Long<br>
* LWP::UserAgent<br>
* Pod::Usage<br>
### Python Libraries

- biopython : https://github.com/biopython/biopython
- seaborn : https://seaborn.pydata.org/

### External programs

Be sure to add these programs to your PATH

 - MAFFT : https://mafft.cbrc.jp/alignment/software/ (no extensions are necessary)
 - CONSENT-correct : https://github.com/morispi/CONSENT/releases/tag/v2.2.2
 - PRANK : http://wasabiapp.org/software/prank/
 - PROBCONS : http://probcons.stanford.edu/
 - Muscle : https://github.com/rcedgar/muscle/
 - Hmmer : http://hmmer.org/
 - BLAST : https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download
 - CD-Hit : https://github.com/weizhongli/cdhit
 - VSearch : https://github.com/torognes/vsearch
 - QIIME2 : https://qiime2.org/

### Processing

There is a multi-stage processing approach that takes raw RoC-ITS reads and ultimately generates clusters of distinct 16S-ITS operons associated with a particular genera:

1) Run python processes to generate RoC-ITS sequences from RoC-ITS reads
2) Run perl to organize RoC-ITS sequences roughly into genera
3) Run perl to further cluster individual genera groupings into distinct RoC-ITS operons

