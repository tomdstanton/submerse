# submerse :aquarius::dna::aquarius:
###### Calculate genome depth and randomly sub-sample reads


_By Tom Stanton_ (he/him)


[![LindIn](https://img.shields.io/badge/LinkedIn-0077B5?style=for-the-badge&logo=linkedin&logoColor=white)](https://uk.linkedin.com/in/tom-stanton-676556100)
[![Holtlab](https://img.shields.io/badge/-Holt%20Lab-black?style=for-the-badge&logo=square&logoColor=white)](https://holtlab.net)
[![Twitter](https://img.shields.io/badge/Twitter-1DA1F2?style=for-the-badge&logo=twitter&logoColor=white)](https://twitter.com/tomstantonmicro)

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)


---
### About
This tool combines two excellent programs I was previously using to calculate
genome coverage and subsample reads: 
[fastq-info](https://github.com/raymondkiu/fastq-info) and 
[seqtk sample](https://github.com/lh3/seqtk).

Additional functionality is the automatic random sub-sampling of reads at 
chosen depths relative to the calculated depth of the genome.

The depth calculation is based the Lander/Waterman equation, where coverage (C) based on read length (L), 
number of reads (N), and genome size (G) [1](#1). 

<p align="center">
    <img src="https://render.githubusercontent.com/render/math?math=C = LN / G">
</p>

Read lengths (insert sizes) are determined by calculating the most frequent read length, and
averaged if there are paired read files. Unfortunately this requires iterating over all
the reads which increases run time. While this can be sped up with external libraries,
this program aims to operate from a single script with no dependencies.

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
                        number of threads (default: #cpus)
```
- [x] Automatic pairing of read and assembly files
- [x] Relative, random subsampling
- [x] Automatic insert size calculation
- [x] Pure Python, no dependencies, single script
- [x] No mapping required
- [ ] Speed optimisation (suggestions welcome!)
- [ ] Writing gzipped files (too slow)


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

Command 2 will produce a table like this:

**Sample**|**Assembly\_file**|**Genome\_size**|**Total\_reads**|**Total\_bases**|**Average\_insert\_size**|**Depth**|**Read\_file**|**Insert\_size**|**N\_reads**|**N\_bases**|**N\_reads\_at\_10X\_depth**|**N\_reads\_at\_20X\_depth**|**N\_reads\_at\_30X\_depth**
:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:
genome|genome.fasta|5284985|3048096|448774957|151|87|genome\_1.fastq.gz|151|1524048|224428859|175178|350356|525534
genome|genome.fasta|5284985|3048096|448774957|151|87|genome\_2.fastq.gz|151|1524048|224346098|175178|350356|525534
Future updates will improve speed and safety, but I will eventually port this to Rust.

---

### References
<a id="1">[1]</a>
Eric S. Lander, Michael S. Waterman (1988).
Genomic mapping by fingerprinting random clones: A mathematical analysis.
Genomics, 2, 3, 4 1988

[2]: http://twitter.com/tomstantonmicro
[2.1]: http://i.imgur.com/tXSoThF.png (twitter icon with padding)
[3]: http://www.github.com/tomdstanton
[3.1]: http://i.imgur.com/0o48UoR.png (github icon with padding)
[4]: mailto:tomdstanton@gmail.com?subject=[submerse]
[4.1]: https://i.imgur.com/vltiL8c.png
