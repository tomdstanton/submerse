#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
submerse
-----------------
Genome depth calculation and relative random sub-sampling of read files.
-----------------
Requires:
- Python >=3.9
-----------------
Tom Stanton, 2022
"""

# ............... Imports ............... #
import os
from random import seed, choices
from argparse import ArgumentParser, RawTextHelpFormatter
import gzip
import sys
import logging
from os import cpu_count, path, makedirs
from re import compile, IGNORECASE, search, findall
from itertools import islice
from time import time
from statistics import mode

# ............... Attributes ............... #
__version__ = '0.0.4'
__author__ = "Tom Stanton"
__maintainer__ = "Tom Stanton"
__email__ = "tomdstanton@gmail.com"
__status__ = "Development"
__title__ = 'submerse'
__author_email__ = 'tomdstanton@gmail.com'
__description__ = 'Genome depth calculation and relative random sub-sampling of read files'
__license__ = 'gpl-3.0'

# ............... Globals ............... #
READ_REGEX = compile('|'.join([
    '_R[12]\.(f(?:ast)?q(?:\.gz)?)$',
    '_R[12]_[0-9]+?\.(f(?:ast)?q(?:\.gz)?)$',
    '_R[12].[0-9]+?\.(f(?:ast)?q(?:\.gz)?)$',
    '_[12]\.(f(?:ast)?q(?:\.gz)?)$',
    '_[12]_[0-9]+?\.(f(?:ast)?q(?:\.gz)?)$',
    '_[12].[0-9]+?\.(f(?:ast)?q(?:\.gz)?)$',
    '\.(f(?:ast)?q(?:\.gz)?)$']), IGNORECASE)

LOGGER = logging.getLogger(__name__)  # https://stackoverflow.com/a/45065655/10771025

ASCII = fr'''
           _                                  
 ___ _   _| |__  _ __ ___   ___ _ __ ___  ___ 
/ __| | | | '_ \| '_ ` _ \ / _ \ '__/ __|/ _ \
\__ \ |_| | |_) | | | | | |  __/ |  \__ \  __/
|___/\__,_|_.__/|_| |_| |_|\___|_|  |___/\___|

==============================================
Estimate genome depth and randomly sub-sample reads
github.com/tomdstanton/submerse
'''


# ............... Functions ............... #
def main():
    init_time = time()
    args = parse_args(sys.argv[1:])

    # Check arguments ##############################################################
    if not any([args.forward, args.reverse, args.single]):
        quit_with_error('No fastq files provided')

    if not any([args.assembly, args.genome_size]):
        quit_with_error('No genomes/genome-sizes provided')

    if args.force and len(args.assembly) not in [len(i) for i in [args.forward, args.reverse, args.single] if i]:
        quit_with_error(
            'The --force argument requires the number of reads in at least one of the read sets to be the '
            'same as the number of assemblies')

    if args.genome_size:
        args.force = True
        if len(args.genome_size) not in [len(i) for i in [args.forward, args.reverse, args.single] if i]:
            quit_with_error(
                'The --genome_size argument requires the number of reads in at least one of the read-sets to be the '
                'same as the number of genome sizes')

    if args.out and not args.depths:
        quit_with_error('Can only output sub-sampled reads at selected depths with the -d option')

    if args.out and not args.random_seed:
        quit_with_error('Can only randomly sub-sample reads with a random seed, try -r 100')

    # Set up logger ##############################################################
    seed(args.random_seed)  # Set random seed using the args
    line_width = 70  # Nice general terminal width
    logging.basicConfig(format="[%(asctime)s %(levelname)s %(funcName)s] %(message)s",
                        datefmt='%H:%M:%S', level=args.loglevel)
    LOGGER.info(format_str(f' {path.basename(__file__)} {__version__} '.center(line_width, '='), '36'))
    for i in ASCII.split('\n'):
        LOGGER.info(format_str(i, '36'))
    LOGGER.info(f'Platform: {sys.platform}')
    if sys.version_info[0] < 3:
        quit_with_error(f"Python version needs to be 3, not {sys.version_info[0]}")
    elif sys.version_info[1] < 9:
        quit_with_error(f"Python version needs to be 3.9, not {sys.version_info[0]}.{sys.version_info[1]}")
    LOGGER.info(f'Python: {sys.version_info[0]}.{sys.version_info[1]}')

    # Process input files ##############################################################
    samples = process_input(args.forward, args.reverse, args.single, args.assembly, args.genome_size, args.force)

    header = f'Sample\tAssembly_file\tGenome_size\tTotal_reads\tTotal_bases\tAve_mean_read_length\t' \
             f'Ave_mode_read_length\tTotal_coverage\tRead_file\tMean_read_length\tMode_read_length\t' \
             f'Read_coverage\tMax_read_length\tMin_read_length\tN_reads\tN_bases'
    if args.depths:
        header += '\t' + "\t".join([f"N_reads_at_{depth}X_depth" for depth in args.depths])
    print(header)

    # Iterate over samples ##############################################################
    for sample in samples:
        if not sample.genome_size:
            sample.genome_size = get_genome_size(sample.assembly)
        for read in sample.reads:
            read.genome_size = sample.genome_size
            read.count_bases(args.equation)

        sample.total_reads += sum([i.n_reads for i in sample.reads])
        sample.total_bases += sum([i.n_bases for i in sample.reads])
        sample.ave_mean_read_length += sum([i.mean_read_length for i in sample.reads]) / len(sample.reads)
        sample.ave_mode_read_length += sum([i.mode_read_length for i in sample.reads]) / len(sample.reads)

        if args.equation == 1:
            sample.coverage += sample.total_bases / sample.genome_size
        elif args.equation == 2:
            sample.coverage += (sample.total_reads * sample.ave_mean_read_length) / sample.genome_size
        elif args.equation == 3:
            sample.coverage += (sample.total_reads * sample.ave_mode_read_length) / sample.genome_size

        for read in sample.reads:
            if args.depths:
                for depth in args.depths:
                    read.subsamples.append(Subsample(sample.coverage, depth, read, args.out))
            print(f'{sample.get_output_string()}\t{read.get_output_string()}')

    LOGGER.info(format_str(f' Completed in {time() - init_time:.2f}s '.center(line_width, '='), '92'))
    sys.exit(0)


def parse_args(arguments):
    parser = ArgumentParser(formatter_class=RawTextHelpFormatter, description=format_str(ASCII, '36'),
                            usage='python submerse.py -a assembly.fasta -1 sr_1.fastq.gz -2 sr_2.fastq.gz '
                                  '-s lr.fastq.gz -d 10 20 30 -out subsampled_reads -r 678', add_help=False)

    input_parser = parser.add_argument_group("Input options")
    input_parser.add_argument('-a', '--assembly', nargs='+', type=str, metavar='.fasta',
                              required=False, help='assembly fasta, not needed if --genome-size set')
    input_parser.add_argument('-1', '--forward', nargs='+', type=str, metavar='fastq(.gz)',
                              required=False, help='forward reads')
    input_parser.add_argument('-2', '--reverse', nargs='+', type=str, metavar='fastq(.gz)',
                              required=False, help='reverse reads')
    input_parser.add_argument('-s', '--single', nargs='+', type=str, required=False, metavar='fastq(.gz)',
                              help='non-paired reads; use for long-reads, combine with -1 -2 for hybrid assembly')
    input_parser.add_argument('-f', '--force', action='store_true', default=False,
                              help='force assembly-read pairing by given file order rather than by filename')

    calculation_parser = parser.add_argument_group("Calculation options")
    calculation_parser.add_argument('-e', '--equation', type=int, default=1, choices=[1, 2, 3],
                                    help='equation used to calculate coverage:\n'
                                         '1. C = N bases / Genome size (from assembly)\n'
                                         '2. C = (N reads * Mean read length) / Genome size (from assembly)\n'
                                         '3. C = (N reads * Mode read length) / Genome size (from assembly)')
    calculation_parser.add_argument('-g', '--genome_size', type=int, nargs='+', metavar='bp',
                                    help='manually set genome size(s) to be paired with the read-set inputs; '
                                         'when used --assembly not required')

    subsample_parser = parser.add_argument_group("Sub-sampling options")
    subsample_parser.add_argument('-o', '--out', type=str, metavar='subsampled_reads/',
                                  help='trigger sub-sampling and set output directory for reads; REQUIRES: -r')
    subsample_parser.add_argument('-r', '--random_seed', type=int, metavar='100',
                                  help='random integer seed for sub-sampling')
    subsample_parser.add_argument('-d', '--depths', nargs='+', type=int, metavar='10 20 30',
                                  help='relative depth(s) to sub-sample at, eg 10x, 20x and 30x')
    # subsample_parser.add_argument('-n', '--nreads', nargs='+', type=int,
    #                               help='use a specified number of reads to sub-sample instead of relative depth; REQUIRES: -r -o')
    # subsample_parser.add_argument('-x', '--subsample_only', action='store_true', default=False,
    #                               help='skip any calculation and just randomly sub-sample reads; REQUIRES: -r -n -o')

    other_parser = parser.add_argument_group("Other options")
    other_parser.add_argument('--version', action='version', version=f'submerse v{__version__}',
                              help="print version and exit")
    other_parser.add_argument('-v', '--verbose', help="verbose statements to stderr",
                              action="store_const", dest="loglevel", const=logging.INFO)
    other_parser.add_argument('-vv', '--extra_verbose', help="debug statements to stderr",
                              action="store_const", dest="loglevel", const=logging.DEBUG, default=logging.WARNING)
    other_parser.add_argument('-h', '--help', action='help', help='show this help message and exit')
    return parser.parse_args(arguments)


def get_read_sizes(read_file: str, gzipped: bool):  # -> list[int]:  ## This might generate extra overhead
    # This function opens the fastq-formatted file to count the number of reads and gets the length of every read
    # -1 removes newline char, islice iterates over 2nd line in steps of 4
    # Maybe try: https://groverj3.github.io/articles/2019-08-22_just-write-your-own-python-parsers-for-fastq-files.html
    start = time()
    if gzipped:
        with gzip.open(read_file, 'rb') as f:
            sizes = [len(line) - 1 for line in islice(f, 1, None, 4)]
    else:
        with open(read_file, 'rb') as f:
            sizes = [len(line) - 1 for line in islice(f, 1, None, 4)]
    LOGGER.debug(f"Counted reads and lengths for {read_file} {format_str(f'[{time() - start:.2f}s]', '94')}")
    return sizes


def subsample_reads(read_file, out_dir, n_reads: int, gzipped: bool) -> str:
    start = time()
    makedirs(out_dir, exist_ok=True)
    out_file = f'{out_dir}/{path.basename(read_file)}'.removesuffix('.gz')
    openfile = gzip.open(read_file, 'rb') if gzipped else open(read_file, 'rb')  # Splitting by '@' doesn't work
    # Writing data once is faster than writing in pieces,
    # so generate sub-sampled reads in list comprehension then write the resulting joined chunk
    open(out_file, 'wb').write(
        b'\n'.join(read for read in choices(findall(b"\n".join([b"[^\n]+"] * 4), openfile.read()), k=n_reads)))
    # This does assume that the sequence lines are not wrapped and there are 4 lines per read
    openfile.close()
    LOGGER.debug(f"Written sub-sampled reads to {out_file} {format_str(f'[{time() - start:.2f}s]', '94')}")
    return out_file


def process_reads(sample, files, read_type, force: bool, order: int):
    """Assigns read files to an assembly and add to the Sample object"""
    for index, file in enumerate(files):
        extension_match = search(READ_REGEX, file)
        if extension_match:
            sample_name = path.basename(file).replace(extension_match.group(0), '')
            if force and order == index:
                sample.reads.append(ReadFile(file, extension_match.group(0), read_type))
                if sample.sample_name == 'unknown':  # Set the sample name with read name if --genome_size used
                    sample.sample_name = sample_name
            elif sample_name == sample.sample_name:
                sample.reads.append(ReadFile(file, extension_match.group(0), read_type))


def process_input(forward_files, reverse_files, single_files, assembly_files, genome_sizes, force: bool):
    """Processes input files and stores in Sample object"""
    for files in [i for i in [assembly_files, forward_files, reverse_files, single_files] if i]:
        for file in files:
            if not path.isfile(file):
                quit_with_error(f'{file} does not exist')

    sample_list = []  # List of sample objects to return
    genomes = assembly_files if assembly_files else genome_sizes

    for index, file in enumerate(genomes):  # Loop over assemblies/genome sizes and add reads
        # We need enumerate in case the --force option is used to add reads based on order
        # Init the sample class, if genome sizes used, set the name to unknown and populate later with read name
        genome_name = path.basename(file).rsplit(".", 1)[0] if assembly_files else "unknown"
        sample = Sample(genome_name, file, genome_sizes)
        # Add read files to the sample objects
        if forward_files:
            process_reads(sample, forward_files, 'forward', force, index)
        if reverse_files:
            process_reads(sample, reverse_files, 'reverse', force, index)
        if single_files:
            process_reads(sample, single_files, 'single', force, index)
        # Check sample objects
        if not sample.reads:
            LOGGER.warning(format_str(f"WARNING: No read files found for {sample.sample_name}", '93'))
        elif len(sample.reads) > 3:
            LOGGER.warning(format_str(f"WARNING: More than 3 files for {sample.sample_name}: "
                                      f"{' '.join(sample.reads)}", '93'))
        else:
            sample_list.append(sample)

    if not sample_list:
        quit_with_error('No files to analyse')
    return sample_list


def quit_with_error(message: str):
    """Displays the given message and ends the program's execution."""
    LOGGER.error(format_str(f"ERROR: {message}", '91'))
    sys.exit(1)


def format_str(string: str, format_number: str) -> str:
    """Add colour to string using a format string number"""
    return f"\033[{format_number}m{string}\033[0m"


def get_genome_size(assembly: str) -> int:
    """Quickly determine genome size from fasta file"""
    start = time()
    with open(assembly, 'rt') as f:  # Can we open in binary mode?
        size = sum([sum([len(x) for x in i.split('\n')[1:]]) for i in f.read().split('>')])
    LOGGER.debug(f"Calculated size of {assembly} ({size}bp) {format_str(f'[{time() - start:.2f}s]', '94')}")
    return size


# ............... Classes ............... #
class Sample(object):
    def __init__(self, sample_name, assembly_file, genome_sizes):
        self.sample_name = sample_name
        self.assembly = "None" if genome_sizes else assembly_file
        self.reads = []
        self.total_reads, self.total_bases, self.ave_mean_read_length, \
        self.ave_mode_read_length, self.coverage = [0] * 5
        self.genome_size = assembly_file if genome_sizes else 0

    def get_output_string(self) -> str:
        return f'{self.sample_name}\t{path.basename(self.assembly)}\t' \
               f'{self.genome_size}\t{self.total_reads}\t{self.total_bases}\t' \
               f'{round(self.ave_mean_read_length)}\t{round(self.ave_mode_read_length)}\t{round(self.coverage)}'


class ReadFile(object):
    def __init__(self, filepath, extension, read_type):
        self.path = filepath
        self.extension = extension
        self.read_type = read_type
        self.gzipped = True if extension.endswith('.gz') else False
        self.n_reads, self.mean_read_length, self.mode_read_length, self.n_bases, \
        self.min_read_length, self.coverage, self.max_read_length, self.genome_size = [0] * 8
        self.subsamples = []

    def get_output_string(self):
        output_string = f'{path.basename(self.path)}\t{round(self.mean_read_length)}\t{round(self.mode_read_length)}\t{round(self.coverage)}\t' \
                        f'{round(self.max_read_length)}\t{round(self.min_read_length)}\t{round(self.n_reads)}\t{round(self.n_bases)}'
        for s in self.subsamples:
            output_string += f'\t{s.n_subsampled_reads}'
        return output_string

    def count_bases(self, equation: int):
        read_sizes = get_read_sizes(self.path, self.gzipped)
        # Start timing here
        start = time()
        self.n_reads = len(read_sizes)
        self.n_bases = sum(read_sizes)
        self.mean_read_length = self.n_bases / self.n_reads
        self.mode_read_length = mode(read_sizes)
        self.max_read_length = max(read_sizes)
        self.min_read_length = min(read_sizes)
        if equation == 1:
            self.coverage = self.n_bases / self.genome_size
        elif equation == 2:
            self.coverage = (self.n_reads * self.mean_read_length) / self.genome_size
        elif equation == 3:
            self.coverage = (self.n_reads * self.mode_read_length) / self.genome_size

        LOGGER.debug(f"Performed read size stats for {self.path} {format_str(f'[{time() - start:.2f}s]', '94')}")


class Subsample(object):
    def __init__(self, coverage, depth, read, out_dir):
        if depth:
            self.depth = depth
            self.depth_string = f'submerse_{depth}X_subsampled_reads'
            self.n_subsampled_reads = round((depth / coverage) * read.n_reads)
        self.path = '-'
        if out_dir:
            if read.n_reads <= self.n_subsampled_reads:
                LOGGER.warning(format_str(
                    f'{read.path} has {read.n_reads} reads and sub-sampling at {depth}X depth'
                    f' would require {self.n_subsampled_reads} reads', '93'))
            else:
                self.path = subsample_reads(read.path, f'{out_dir}/{self.depth_string}',
                                            self.n_subsampled_reads, read.gzipped)


# ............... Main Program ............... #
if __name__ == "__main__":
    main()
