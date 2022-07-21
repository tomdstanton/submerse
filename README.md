# submerse :aquarius::dna:
###### Calculate genome depth and randomly sub-sample reads


_By Tom Stanton_ (he/him)


[![LindIn](https://img.shields.io/badge/LinkedIn-0077B5?style=for-the-badge&logo=linkedin&logoColor=white)](https://uk.linkedin.com/in/tom-stanton-676556100)
[![Holtlab](https://img.shields.io/badge/-Holt%20Lab-black?style=for-the-badge&logo=square&logoColor=white)](https://holtlab.net)
[![Twitter](https://img.shields.io/badge/Twitter-1DA1F2?style=for-the-badge&logo=twitter&logoColor=white)](https://twitter.com/tomstantonmicro)

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)


---
### About
This tool is basically an amalgamation two excellent programs I was previously using to calculate
genome coverage and subsample reads: 
[fastq-info](https://github.com/raymondkiu/fastq-info) and 
[seqtk sample](https://github.com/lh3/seqtk).

Additional functionality is the automatic random sub-sampling of reads at 
chosen depths relative to the calculated depth of the genome.

Coverage (C) is calculated based on either one of three equations:
1. C = N bases / Genome size (from assembly)
2. C = (N reads * Mean read length) / Genome size (from assembly)
3. C = (N reads * Mode read length) / Genome size (from assembly)

### Dependencies and installation
* Python >= 3.9.
* No installation required, just clone the script and run!

### Usage
• Give `submerse` an assembly fasta and read fastq to calculate depth of each read file. If there is more than one
read file for the assembly (paired-end or hybrid) it will average these values. 

• You then have the option to randomly subsample at various depths based on the relative average 
depth `submerse` calculates.

• `submerse` will try to match reads and assemblies based on the _name_ of the file 
(minus the path and extensions). You can override this behavour and match assembly files to read files in the **order** in
which they were given to the arguments using the `--force` option.
```commandline
usage: python submerse.py -a assembly.fasta -1 sr_1.fastq.gz -2 sr_2.fastq.gz -s lr.fastq.gz -d 10 20 30 -out subsampled_reads -r 678

submerse

optional arguments:
  -h, --help            show this help message and exit
  --version             print version and exit
  -a ASSEMBLY [ASSEMBLY ...], --assembly ASSEMBLY [ASSEMBLY ...]
                        assembly fasta (default: None)
  -1 FORWARD [FORWARD ...], --forward FORWARD [FORWARD ...]
                        forward fastq(.gz) (default: None)
  -2 REVERSE [REVERSE ...], --reverse REVERSE [REVERSE ...]
                        reverse fastq(.gz) (default: None)
  -s SINGLE [SINGLE ...], --single SINGLE [SINGLE ...]
                        non-paired fastq(.gz); use for long-reads, combine with -1 -2 for hybrid assembly (default: None)
  -r RANDOM_SEED, --random-seed RANDOM_SEED
                        random integer seed for sub-sampling, e.g. 10 (default: None)
  -d DEPTHS [DEPTHS ...], --depths DEPTHS [DEPTHS ...]
                        depth(s) to subsample at (default: None)
  -e {1,2,3}, --equation {1,2,3}
                        equation used to calculate coverage (default: 1)
  -o OUT, --out OUT     output directory for sub-sampled reads (default: None)
  -f, --force           force given file order (default: False)
  -v, --verbose         verbose statements to stderr (default: None)
  -vv, --extra_verbose  debug statements to stderr (default: 30)
```
- [x] Automatic pairing of read and assembly files
- [x] Relative, random subsampling
- [x] Base Python, no dependencies, single script
- [x] No mapping/alignment required
- [ ] TODO Speed optimisation + threading (suggestions welcome!)
- [ ] TODO Writing gzipped sub-sampled reads (slow performance)

### Output table
Regardless of whether you are sub-sambling reads or just calculating genome covergae stats, the stdout will 
print a table like this, where each line corresponds to each read file provided:

|Sample            |Assembly_file|Genome_size|Total_reads|Total_bases|Ave_mean_read_length|Ave_mode_read_length|Total_coverage|Read_file        |Mean_read_length|Mode_read_length|Read_coverage|Max_read_length|Min_read_length|N_reads|N_bases  |
|------------------|-------------|-----------|-----------|-----------|--------------------|--------------------|--------------|-----------------|----------------|----------------|-------------|---------------|---------------|-------|---------|
|genome            |genome.fasta |5284985    |3048096    |448774957  |147                 |151                 |85            |genome_1.fastq.gz|147             |151             |42           |151            |35             |1524048|224428859|
|genome            |genome.fasta |5284985    |3048096    |448774957  |147                 |151                 |85            |genome_2.fastq.gz|147             |151             |42           |151            |35             |1524048|224346098|

### Examples
Calculate genome depth for paired-end reads and pipe to a tab-separated file.
```commandline
./submerse.py -a assembly.fasta -1 sr_1.fastq.gz -2 sr_2.fastq.gz >> genome_depths.txt
```
Same as above, but calculate N reads for ```10x```, ```20x```, and ```30x``` depth.
```commandline
./submerse.py -a assembly.fasta -1 sr_1.fastq.gz -2 sr_2.fastq.gz -d 10 20 30 >> genome_depths.txt
```
Same as above, but randomly subsample the reads to ```subsampled_reads/``` using a random seed of ```100```.
```commandline
./submerse.py -a assembly.fasta -1 sr_1.fastq.gz -2 sr_2.fastq.gz -d 10 20 30 -o . -r 100 >> genome_depths.txt
```
This will generate a directory tree like this, where the sub-sampled reads have the **same file
name** as the original read file:

```
submerse_10X_subsampled_reads/read_1.fastq
submerse_10X_subsampled_reads/read_2.fastq
submerse_20X_subsampled_reads/read_1.fastq
submerse_20X_subsampled_reads/read_2.fastq
submerse_30X_subsampled_reads/read_1.fastq
submerse_30X_subsampled_reads/read_2.fastq
```


### Future plans
Future updates will improve speed and safety, but I will eventually port this to Rust.

---

### References
<a id="1">[1]</a>
Eric S. Lander, Michael S. Waterman.
Genomic mapping by fingerprinting random clones: A mathematical analysis.
Genomics, 2-4, 1988

[2]: http://twitter.com/tomstantonmicro
[2.1]: http://i.imgur.com/tXSoThF.png (twitter icon with padding)
[3]: http://www.github.com/tomdstanton
[3.1]: http://i.imgur.com/0o48UoR.png (github icon with padding)
[4]: mailto:tomdstanton@gmail.com?subject=[submerse]
[4.1]: https://i.imgur.com/vltiL8c.png
