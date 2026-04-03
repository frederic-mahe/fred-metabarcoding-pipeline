# fred-metabarcoding-pipeline

Fred's metabarcoding pipeline: from raw fastq files to an occurrence
table, at light-speed

This is a work in progress. An older pipeline is also [described
here](https://github.com/frederic-mahe/swarm/wiki/Fred's-metabarcoding-pipeline).

The pipeline described here is a swarm-based metabarcoding
pipeline. It relies on free and open-source tools such as
[cutadapt](https://cutadapt.readthedocs.io/en/stable/),
[swarm](https://github.com/torognes/swarm), and
[vsearch](https://github.com/torognes/vsearch).

This pipeline aims to produce single-nucleotide resolution clustering
results, equivalent to Amplicon Sequence Variants (ASVs). Clustering
results, chimera detection results, and taxonomic assignment results
are grouped in a ready-to-use occurrence table.

```sh
((BASH_VERSINFO[0] <= 4)) && echo "error: the shell you are using is too old"
```
