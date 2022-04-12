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
from random import seed, choices
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import gzip
import sys
import logging
from os import cpu_count, path, makedirs
from re import compile, IGNORECASE, search, findall
from itertools import islice
from concurrent import futures

# ............... Attributes ............... #
__version__ = '0.0.2'
__author__ = "Tom Stanton"
__maintainer__ = "Tom Stanton"
__email__ = "tomdstanton@gmail.com"
__status__ = "Development"
__title__ = 'submerse'
__author_email__ = 'tomdstanton@gmail.com'
__description__ = 'Genome depth calculation and relative random sub-sampling of read files'
__license__ = 'gpl-3.0'

# ............... Globals ............... #
read_regex = compile('|'.join([
    '_R[12]\.(f(?:ast)?q(?:\.gz)?)$',
    '_R[12]_[0-9]+?\.(f(?:ast)?q(?:\.gz)?)$',
    '_R[12].[0-9]+?\.(f(?:ast)?q(?:\.gz)?)$',
    '_[12]\.(f(?:ast)?q(?:\.gz)?)$',
    '_[12]_[0-9]+?\.(f(?:ast)?q(?:\.gz)?)$',
    '_[12].[0-9]+?\.(f(?:ast)?q(?:\.gz)?)$',
    '\.(f(?:ast)?q(?:\.gz)?)$']), IGNORECASE)

logger = logging.getLogger(__name__)  # https://stackoverflow.com/a/45065655/10771025


# ............... Functions ............... #
def main():
    args = parse_args(sys.argv[1:])
    # args = parse_args("-a /Users/tom/Bioinformatics/submerse/R-AJ101A_S62_L002.fasta -1 /Users/tom/Bioinformatics/submerse/R-AJ101A_S62_L002_1.fastq.gz -2 /Users/tom/Bioinformatics/submerse/R-AJ101A_S62_L002_2.fastq.gz -f".split())

    if not any([args.forward, args.reverse, args.single]):
        quit_with_error('No fastq files provided')

    if args.force and len(args.assembly) not in [len(i) for i in [args.forward, args.reverse, args.single] if i]:
        quit_with_error('The --force argument requires the number of reads in at least one of the read sets to be the '
                        'same as the number of assemblies')

    if args.out and not args.depths:
        quit_with_error('Can only output sub-sampled reads at selected depths with the --depths option')

    if args.out and not args.seed:
        quit_with_error('Can only randomly sub-sample reads with a random seed, try --seed 100')

    args.threads = min(args.threads, cpu_count())
    seed(args.random_seed)

    line_width = 70  # Nice general terminal width
    logging.basicConfig(format="[%(asctime)s::%(levelname)s::%(funcName)s] %(message)s",
                        datefmt='%H:%M:%S', level=args.loglevel)
    logger.info(f' {path.basename(__file__)} {__version__} '.center(line_width, '='))
    logo()  # Logger output logo if verbose
    logger.info(f'Platform: {sys.platform}')
    logger.info(f'Python: {".".join([str(i) for i in sys.version_info[:3]])}')
    logger.info(f'Using {args.threads} threads')

    samples = process_input(args.forward, args.reverse, args.single, args.assembly, args.force)

    header = f'Sample\tAssembly_file\tGenome_size\tTotal_reads\tTotal_bases\tAverage_insert_size\tDepth\tRead_file\t' \
             f'Insert_size\tN_reads\tN_bases'
    if args.depths:
        header += '\t' + "\t".join([f"N_reads_at_{depth}X_depth" for depth in args.depths])
    print(header)
    executor = futures.ProcessPoolExecutor(args.threads)
    ftures = [executor.submit(s.calculate_coverage(args.depths, args.out)) for s in samples]
    futures.wait(ftures)
    sys.exit(0)


def parse_args(arguments):
    parser = ArgumentParser(description='submerse', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', action='version', version=f'submerse v{__version__}',
                        help="print version and exit")
    parser.add_argument('-a', '--assembly', nargs='+', type=str, required=True, help='assembly fasta')
    parser.add_argument('-1', '--forward', nargs='+', type=str, required=False, help='forward fastq(.gz)')
    parser.add_argument('-2', '--reverse', nargs='+', type=str, required=False, help='reverse fastq(.gz)')
    parser.add_argument('-s', '--single', nargs='+', type=str, required=False,
                        help='non-paired fastq(.gz); use for long-reads, combine with -1 -2 for hybrid assembly')
    parser.add_argument('-r', '--random-seed', type=int, help='random integer seed for sub-sampling, e.g. 10')
    parser.add_argument('-d', '--depths', nargs='+', type=int, help='depth(s) to subsample at')
    parser.add_argument('-o', '--out', type=str, help='output directory for sub-sampled reads')
    parser.add_argument('-t', '--threads', type=int, default=8, help='number of threads')
    parser.add_argument('-f', '--force', action='store_true', default=False, help='force given file order')
    parser.add_argument('-v', '--verbose', help="verbose statements to stderr",
                        action="store_const", dest="loglevel", const=logging.INFO)
    parser.add_argument('-vv', '--extra_verbose', help="debug statements to stderr",
                        action="store_const", dest="loglevel", const=logging.DEBUG, default=logging.WARNING)
    return parser.parse_args(arguments)


def get_reads_and_insert(read_file, gzipped):
    # This function opens the fastq-formatted file to count the number of reads and gets the length of every
    # read to calculate the most frequently occurring insert size.
    # It's probably the slowest part of the program, and pure-python optimisations are encouraged.
    openfile = gzip.open(read_file, 'rb') if gzipped else open(read_file, 'rb')
    # Size is the most common insert size, however this is slower as requires collection of all seq lengths
    # islice iterates over every 2nd line in steps of 4
    # -1 removes newline char count
    sizes = [len(line) - 1 for line in islice(openfile, 1, None, 4)]
    openfile.close()
    return len(sizes), max(set(sizes), key=sizes.count), sum(sizes)


def subsample_reads(read_file, out_dir, n_reads, gzipped):
    makedirs(out_dir, exist_ok=True)
    out_file = f'{out_dir}/{path.basename(read_file)}'.removesuffix('.gz')
    openfile = gzip.open(read_file, 'rb') if gzipped else open(read_file, 'rb')  # Splitting by '@' doesn't work
    # Writing data once is faster than writing in pieces,
    # so generate sub-sampled reads in list comprehension then write the resulting joined chunk
    open(out_file, 'wb').write(
        b'\n'.join(read for read in choices(findall(b"\n".join([b"[^\n]+"] * 4), openfile.read()), k=n_reads)))
    # This does assume that the sequence lines are not wrapped and there are 4 lines per read
    openfile.close()
    return out_file


def process_reads(sample, files, read_type, force, order):
    """Assigns read files to an assembly and add to the Sample object"""

    for index, file in enumerate(files):
        extension_match = search(read_regex, file)
        if extension_match:
            if force and order == index:
                sample.reads.append(ReadFile(file, extension_match.group(0), read_type))
            elif path.basename(file).replace(extension_match.group(0), '') == sample.sample_name:
                sample.reads.append(ReadFile(file, extension_match.group(0), read_type))


def process_input(forward_files, reverse_files, single_files, assembly_files, force):
    """Processes input files and stores in Sample object"""

    for files in [assembly_files, forward_files, reverse_files, single_files]:
        if files:
            for file in files:
                if not path.isfile(file):
                    quit_with_error(f'{file} does not exist')
    sample_list = []

    for index, file in enumerate(assembly_files):
        sample = Sample(path.basename(file).rsplit(".", 1)[0], file)
        if forward_files:
            process_reads(sample, forward_files, 'forward', force, index)
        if reverse_files:
            process_reads(sample, reverse_files, 'reverse', force, index)
        if single_files:
            process_reads(sample, single_files, 'single', force, index)

        if not sample.reads:
            logger.warning(format_str(f"WARNING: No read files found for {sample}", '93'))
        elif len(sample.reads) > 3:
            logger.warning(format_str(f"WARNING: More than 3 files for {sample}: {' '.join(sample.reads)}", '93'))
        else:
            sample_list.append(sample)

    if not sample_list:
        quit_with_error('No files to analyse')
    return sample_list


def quit_with_error(message):
    """Displays the given message and ends the program's execution."""
    logger.error(format_str(f"ERROR: {message}", '91'))
    sys.exit(1)


def format_str(string, format_number):
    return f"\033[{format_number}m{string}\033[0m"


def logo():
    logger.info(format_str('''███████╗██╗   ██╗██████╗ ███╗   ███╗███████╗██████╗ ███████╗███████╗''', '36'))
    logger.info(format_str('''██╔════╝██║   ██║██╔══██╗████╗ ████║██╔════╝██╔══██╗██╔════╝██╔════╝''', '36'))
    logger.info(format_str('''███████╗██║   ██║██████╔╝██╔████╔██║█████╗  ██████╔╝███████╗█████╗  ''', '36'))
    logger.info(format_str('''╚════██║██║   ██║██╔══██╗██║╚██╔╝██║██╔══╝  ██╔══██╗╚════██║██╔══╝  ''', '36'))
    logger.info(format_str('''███████║╚██████╔╝██████╔╝██║ ╚═╝ ██║███████╗██║  ██║███████║███████╗''', '36'))
    logger.info(format_str('''╚══════╝ ╚═════╝ ╚═════╝ ╚═╝     ╚═╝╚══════╝╚═╝  ╚═╝╚══════╝╚══════╝''', '36'))


def get_genome_size(assembly):
    with open(assembly, 'rt') as f:  # Can we open in binary mode?
        return sum([sum([len(x) for x in i.split('\n')[1:]]) for i in f.read().split('>')])


# ............... Classes ............... #
class Sample(object):
    def __init__(self, sample_name, assembly_file):
        self.sample_name = sample_name
        self.assembly = assembly_file
        self.reads = []
        self.genome_size = 0
        self.total_reads = 0
        self.total_bases = 0
        self.ave_insert_size = 0
        self.coverage = 0
        self.assembly_method = ''

    def get_output_string(self):
        return f'{self.sample_name}\t{path.basename(self.assembly)}\t' \
               f'{self.genome_size}\t{self.total_reads}\t{self.total_bases}\t' \
               f'{self.ave_insert_size}\t{self.coverage}'

    # We perform all the I/O operations here so we can do it concurrently
    def calculate_coverage(self, subsample_depths, out_dir):
        self.genome_size = get_genome_size(self.assembly)
        for read in self.reads:
            read.n_reads, read.insert_size, read.n_bases = get_reads_and_insert(read.path, read.gzipped)

        self.total_reads += sum([i.n_reads for i in self.reads])
        self.total_bases += sum([i.n_bases for i in self.reads])
        self.ave_insert_size += round(sum([i.insert_size for i in self.reads]) / len(self.reads))
        self.coverage += round((self.total_reads * self.ave_insert_size) / self.genome_size)

        for read in self.reads:
            if subsample_depths:
                for depth in subsample_depths:
                    read.subsamples.append(Subsample(self.coverage, depth, read, out_dir))
            print(f'{self.get_output_string()}\t{read.get_output_string()}')


class ReadFile(object):
    def __init__(self, filepathpath, extension, read_type):
        self.path = filepathpath
        self.extension = extension
        self.read_type = read_type
        self.gzipped = True if extension.endswith('.gz') else False
        self.n_reads, self.insert_size, self.n_bases = 0, 0, 0
        self.subsamples = []

    def get_output_string(self):
        output_string = f'{path.basename(self.path)}\t{self.insert_size}\t{self.n_reads}\t{self.n_bases}'
        for s in self.subsamples:
            output_string += f'\t{s.n_subsampled_reads}'
        return output_string


class Subsample(object):
    def __init__(self, coverage, depth, read, out_dir):
        self.depth = depth
        self.depth_string = f'submerse_{depth}X_subsampled_reads'
        self.n_subsampled_reads = round((depth / coverage) * read.n_reads)
        self.path = '-'
        if out_dir:
            if read.n_reads <= self.n_subsampled_reads:
                logger.warning(format_str(
                    f'{read.path} has {read.n_reads} reads and sub-sampling at {depth}X depth'
                    f' would require {self.n_subsampled_reads} reads', '93'))
            else:
                self.path = subsample_reads(read.path, f'{out_dir}/{self.depth_string}',
                                            self.n_subsampled_reads, read.gzipped)


# ............... Main Program ............... #
if __name__ == "__main__":
    main()

