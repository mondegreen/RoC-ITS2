#Clustering RoC-ITS sequences

The initial clustering step takes the output of the python pipeline. The following command will take a fasta file of RoC-ITS sequences and cluster them. Here is an example command:

```bash
perl generateGeneraClusters.pl -input test16.fa -outDir rocItsClustersDir -basename test -minsubreads 4 -rocitsPath . -minclstrSize 10
```

designed to, but not guarranted to, cluster by genera.

