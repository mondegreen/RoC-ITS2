## Clustering RoC-ITS sequences

The initial clustering step takes the output of the python pipeline (in this case the file provided "test_RoCITS.fa" but you could also use the "output" file from the python validation run), extracts the 16S portion of each RoC-ITS sequence and uses the 16S portion to cluster the full length RoC-ITS sequences roughly at the genera-level. The following command will take a fasta file of RoC-ITS sequences and cluster them. Here is an example command:

```bash
perl generateGeneraClusters.pl -input test_RoCITS.fa -outDir rocItsClustersDir -clusterPid 0.95 -basename test -minsubreads 4 -rocitsPath . -minclstrSize 10
```

Note that the generateGeneraClusters.pl script is designed to, but not guarranted to, cluster by genera. Clustering occurs based on the 16S gene alone. The inclusiveness of the clustering is largely controlled by the clusterPid parameter which defaults to clustering at 95%. We strongly recommend the -minsubreads to be 5 for real data but for this example we have used 4 to provide a quick example case for testing/validation purposes.

This script will create a directory named by the -outDir option. Within that directory will be a series of fasta files with the format:<br>

test.v1.num=v2.fa<br>

where v1 is the cluster UID and v2 is the number of RoC-ITS sequences that are in it. Given the test data, this will produce 1 cluster with 162 sequences in it.

Internally the script performs several steps:

1. Filter the RoC-ITS sequences based on the number of sub-reads used during construction (default = 10; controlled by the -minsubreads parameter).
2. Identifying 16S gene within RoC-ITS sequences using regular expressions. This is controlled by two parameters:
  - -16s_fiveprime (defaults to "CTGAGCCAGGATCAAACTCT")
  - -16s_threeprime (defaults to "TACGG.TACCTTGTTACGACTT")
3. Blast step to identify 16S genes that could not be detected by the regular expression but that were nonetheless present
4. Extraction of the 16S genes
5. Clustering of the 16S genes using cd-hit-est. This clustering step is largely controlled by the -clusterPid paramter (default = 0.95)
6. Construction of the output data in a directory format specified by the -outDir parameter. Note that the entire RoC-ITS read is present in the output clusters even though only the 16S gene is used to determine cluster assignment.


## Clustering to identify distinct 16S-ITS operons

Further processing and clustering of the genera clusters now takes place using the full-length RoC-ITS sequences. This will further remove sequences with high numbers of errors and/or chimeric sequences.

```bash
perl processSampleSet.pl -fasta cluster.v1.num=v2.fa -targetDir targetV1 -pid 990 -identifier test -rocitsPath $PWD -log 1 -chimeraDb reference.fa -minSize 7
```

Many of these parameters will be cluster and/or organism.<br>
For example, we recommend that the minSize be set at roughly 1/10th to 1/20th of the total number of sequences. Most organisms have 10 or fewer ribosomal operons and this should filter out erroneous clusters. However, if multiple organisms are present in a sample, there may be no ideal minimum size that will not result in the loss of clusters that represent real organisms.
The pid parameter will dictate how many clusters are generated and how much variation is tolerated. 16S operons can be identical or show large amounts of variation between them. We have found that rough clustering with a pid parameter of 900 is good if you want to remove many of the chimera automatically. Manual inspection of the resulting multiple sequence alignments (MSA) may be necessary to refine this parameter if cluster corresponding to distinct operons are desired.
The chimeraDb is an optional parameter where the user can provide trusted reference sequences to identify chimeric RoC-ITS sequences.

Internally the processSampleSet.pl script performs several steps:

1. (optional) If a set of trusted 16S-ITS sequences is provided with the -chimeraDb option, vsearch is used to identify and eliminate chimera
2. Multiple sequence alignment (MSA) is performed using muscle
3. The MSA is improved. Columns that occur in 3 or fewer sequences will be ignored
4. The MSA is corrected. In this case, errors in columns will be corrected if they occur rarely (in less than 5% of the sequences or 3 sequences, which ever is the larger number). The minimum number of sequences can be overridden by the -overrideMinCnt parameter which is useful if a cluster is thought to contain a number of different organisms (for example, the Pseudomonas cluster in the publication).
5. The ends of the MSA are trimmed to eliminate noise
6. The variable columns are extracted (along with neighboring columns) from the MSA and used to construct a informative subsequence for further clustering. This makes the ratio of informative to uninformative positions uniform so that clustering will not behave wildly depending on the number of variant positions or the length of the 16S-ITS sequence.
7. cd-hit clustering of the extracted variant positions occurs
8. The full length RoC-ITS sequences are then clustered into a directory structure similar to the input based on the clusters found in step 7.
9. For each sub cluster an MSA is generated

Depending on the identity of the clustering (controlled by the -pid paramter) the clusters will reflect the number of ribosomal operons found in the organism.
