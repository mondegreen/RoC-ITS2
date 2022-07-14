## Clustering RoC-ITS sequences

The initial clustering step takes the output of the python pipeline (in this case the file test_RoCITS.fa but you could also use the "output" file from the python validation run), extracts the 16S portion of each RoC-ITS sequence and uses the 16S portion to cluster the full length RoC-ITS sequences roughly at the genera-level. The following command will take a fasta file of RoC-ITS sequences and cluster them. Here is an example command:

```bash
perl generateGeneraClusters.pl -input test_RoCITS.fa -outDir rocItsClustersDir -basename test -minsubreads 4 -rocitsPath . -minclstrSize 10
```

designed to, but not guarranted to, cluster by genera. We strongly recommend the -minsubreads to be 5 for real data but for this example we have used 4 to provide a quick example case for testing/validation purposes.

This script will create a directory named by the -outDir option. Within that directory will be a series of fasta files with the format:<br>

test.v1.num=v2.fa<br>

where v1 is the cluster UID and v2 is the number of RoC-ITS sequences that are in it. Given the test data, this will produce 1 cluster with 162 sequences in it.

## Clustering to identify distinct 16S-ITS operons

Further clustering of the genera clusters, but now using the full-length RoC-ITS sequences, will further remove sequences with high numbers of errors and/or chimeric sequences.

```bash
perl software4release/processSampleSet.pl -fasta cluster.v1.num=v2.fa -targetDir targetV1 -pid 990 -identifier test -rocitsPath $PWD -log 1 -chimeraDb reference.fa -minSize 7
```

Many of these parameters will be cluster and/or organism.<br>
For example, we recommend that the minSize be set at roughly 1/10th to 1/20th of the total number of sequences. Most organisms have 10 or fewer ribosomal operons and this should filter out erroneous clusters. However, if multiple organisms are present in a sample, there may be no ideal minimum size that will not result in the loss of clusters that represent real organisms.
The pid parameter will dictate how many clusters are generated and how much variation is tolerated. 16S operons can be identical or show large amounts of variation between them. We have found that rough clustering with a pid parameter of 900 is good if you want to remove many of the chimera automatically. Manual inspection of the resulting multiple sequence alignments (MSA) may be necessary to refine this parameter if cluster corresponding to distinct operons are desired.
The chimeraDb is an optional parameter where the user can provide trusted reference sequences to identify chimeric RoC-ITS sequences.

{describe outcome of processSampleSet}
