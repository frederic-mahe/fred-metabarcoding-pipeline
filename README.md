# fred-metabarcoding-pipeline

Fred's metabarcoding pipeline: from raw fastq files to an occurrence
table at light-speed

This is a work in progress. A less up-to-date pipeline is [described
here](https://github.com/frederic-mahe/swarm/wiki/Fred's-metabarcoding-pipeline).

The goal is to identify sequence variants, possibly down to a
single-nucleotide difference when necessary, and to produce a clean,
ready-to-use occurrence table (number of observations of each sequence
variant in each sample).


```sh
((BASH_VERSINFO[0] <= 4)) && echo "error: the shell you are using is too old"
```
