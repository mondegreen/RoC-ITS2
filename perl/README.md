#Clustering RoC-ITS sequences

The initial clustering step takes the output of the python pipeline. The following command will take a fasta file of RoC-ITS sequences and cluster them. Here is an example command:

```bash
perl generateGeneraClusters.pl -input test_RoCITS.fa -outDir rocItsClustersDir -basename test -minsubreads 4 -rocitsPath . -minclstrSize 10
```

designed to, but not guarranted to, cluster by genera.

This script will create a directory named by the -outDir option. Within that directory will be a series of fasta files with the format:<br>

test.v1.num=v2.fa<br>

where v1 is the cluster UID and v2 is the number of RoC-ITS sequences that are in it.


