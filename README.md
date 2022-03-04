# submerse :aquarius::dna::aquarius:
###### Calculate genome depth and randomly sub-sample reads


_By Tom Stanton_ (he/him) :scientist: 
[email me!](mailto:tomdstanton@gmail.com?subject=[submerse]) \
[![alt text][1.1]][1] [![alt text][6.1]][6]

---
### About
This tool combines two excellent programs I was previously using to calculate
genome coverage and subsample reads: 
[fastq-info](https://github.com/raymondkiu/fastq-info) and 
[seqtk sample](https://github.com/lh3/seqtk).

Additional functionality is the automatic random sub-sampling of reads at 
chosen depths relative to the calculated depth of the genome.

<p align="center">
    <img src="https://render.githubusercontent.com/render/math?math=C = LN / G">
</p>

The depth calculation is based the Lander/Waterman equation, where coverage (C) based on read length (L), 
number of reads (N), and genome size (G) [1](#1).

### Usage
```commandline
usage: submerse.py [-h] [--version] -a ASSEMBLIES [ASSEMBLIES ...] -r READS [READS ...] [-s SEED] [-d DEPTHS [DEPTHS ...]] [-o OUT] [-t THREADS]

submerse

optional arguments:
  -h, --help            show this help message and exit
  --version             show version number and exit
  -a ASSEMBLIES [ASSEMBLIES ...], --assemblies ASSEMBLIES [ASSEMBLIES ...]
                        assembly fasta files (default: None)
  -r READS [READS ...], --reads READS [READS ...]
                        read fastq files, can be gzipped (default: None)
  -s SEED, --seed SEED  seed for random sub-sampling (default: None)
  -d DEPTHS [DEPTHS ...], --depths DEPTHS [DEPTHS ...]
                        depth(s) to subsample at (default: None)
  -o OUT, --out OUT     output directory for sub-sampled reads (default: None)
  -t THREADS, --threads THREADS
                        number of threads (default: 8)
```

### Examples
Calculate genome depth for paired-end reads and pipe to a tab-separated file.
```commandline
./submerse.py -a assemblies/*.fasta -r reads/*_{1,2}.fastq.gz >> genome_depths.txt
```
Same as above, but calculate N reads for ```10x```, ```20x```, and ```30x``` depth.
```commandline
./submerse.py -a assemblies/*.fasta -r reads/*_{1,2}.fastq.gz -d 10 20 30 >> genome_depths.txt
```
Same as above, but randomly subsample the reads to ```subsampled_reads/``` using a random seed of ```100```.
```commandline
./submerse.py -a assemblies/*.fasta -r reads/*_{1,2}.fastq.gz -d 10 20 30 -o subsampled_reads -s 100
```
This will generate a directory tree like this, where the sub-sampled reads have the **same file
name** as the original read file:

```
ls subsampled_reads/
subsampled_reads/submerse_10X_subsampled_reads
subsampled_reads/submerse_20X_subsampled_reads
subsampled_reads/submerse_30X_subsampled_reads
```
N.B. The output files will be **uncompressed** due to slow performance in  Python.

Future updates will improve speed and safety, but I will eventually port this to Rust.

---

### References
<a id="1">[1]</a>
Eric S. Lander, Michael S. Waterman (1988).
Genomic mapping by fingerprinting random clones: A mathematical analysis.
Genomics, 2, 3, 4 1988

[1]: http://twitter.com/tomstantonmicro
[1.1]: http://i.imgur.com/tXSoThF.png (twitter icon with padding)
[6]: http://www.github.com/tomdstanton
[6.1]: http://i.imgur.com/0o48UoR.png (github icon with padding)
