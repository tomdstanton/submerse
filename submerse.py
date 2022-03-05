#!/usr/bin/env python3

from random import seed, choices
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import gzip
import sys
from os import cpu_count, path, makedirs
from re import compile, IGNORECASE, search, findall
from itertools import islice
from concurrent import futures

__version__ = '0.0.1'


def parse_args(arguments):
    parser = ArgumentParser(description='submerse', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', action='version', version=f'submerse v{__version__}', help="print version and exit")
    parser.add_argument('-a', '--assemblies', nargs='+', type=str, required=True, help='assembly fasta files')
    parser.add_argument('-r', '--reads', nargs='+', type=str, required=True, help='read fastq files, can be gzipped')
    parser.add_argument('-s', '--seed', type=int, help='seed for random sub-sampling')
    parser.add_argument('-d', '--depths', nargs='+', type=int, help='depth(s) to subsample at')
    parser.add_argument('-o', '--out', type=str, help='output directory for sub-sampled reads')
    parser.add_argument('-t', '--threads', type=int, default=cpu_count(), help='number of threads')
    return parser.parse_args(arguments)


def get_reads_and_insert(read_file, gzipped):
    openfile = gzip.open(read_file, 'rb') if gzipped else open(read_file, 'rb')
    sizes = [len(line.strip()) for line in islice(openfile, 1, None, 4)]  # Is collecting all seq lengths slower?
    openfile.close()
    return len(sizes), max(set(sizes), key=sizes.count), sum(sizes)  # Size is the most common insert size


def subsample_reads(read_file, out_dir, n_reads, gzipped):
    makedirs(out_dir, exist_ok=True)
    out_file = f'{out_dir}/{path.basename(read_file)}'.removesuffix('.gz')
    openfile = gzip.open(read_file, 'rt') if gzipped else open(read_file, 'rt')  # Splitting by '@' so open as rt?
    # Writing data once is faster than writing in pieces,
    # so generate sub-sampled reads in list comprehension then write the resulting joined chunk
    open(out_file, 'w').write(
        '\n'.join(read for read in choices(findall("\n".join(["[^\n]+"] * 4), openfile.read()), k=n_reads)))
    # This does assume that the sequence lines are not wrapped and there are 4 lines per read
    openfile.close()
    return out_file


def process_reads(read_files, assembly_files):
    for file in assembly_files + read_files:
        if not path.isfile(file):
            quit_with_error(f'{file} does not exist')

    samples = []
    read_regex = compile('|'.join([
        '_R[12]\.(f(?:ast)?q(?:\.gz)?)$',
        '_R[12]_[0-9]+?\.(f(?:ast)?q(?:\.gz)?)$',
        '_R[12].[0-9]+?\.(f(?:ast)?q(?:\.gz)?)$',
        '_[12]\.(f(?:ast)?q(?:\.gz)?)$',
        '_[12]_[0-9]+?\.(f(?:ast)?q(?:\.gz)?)$',
        '_[12].[0-9]+?\.(f(?:ast)?q(?:\.gz)?)$',
        '\.(f(?:ast)?q(?:\.gz)?)$']), IGNORECASE)

    for file in assembly_files:
        sample = Sample(path.basename(file).rsplit(".", 1)[0], file)
        for read_file in read_files:
            extension_match = search(read_regex, read_file)
            if extension_match:
                if path.basename(read_file).replace(extension_match.group(0), '') == sample.sample_name:
                    sample.reads.append(ReadFile(read_file, extension_match.group(0)))

        if not sample.reads:
            print(f'WARNING: No read files found for {sample}', file=sys.stderr)
        elif len(sample.reads) > 2:
            print(f'WARNING: More than 2 files for {sample}: {" ".join(sample.reads)}', file=sys.stderr)
        else:
            samples.append(sample)

    if not samples:
        quit_with_error('No files to analyse')
    return samples


def quit_with_error(message):
    print(f'ERROR: {message}', file=sys.stderr)
    sys.exit(1)


def get_genome_size(assembly):
    with open(assembly, 'rt') as f:  # Can we open in binary mode?
        return sum([sum([len(x) for x in i.split('\n')[1:]]) for i in f.read().split('>')])


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
    def __init__(self, path, extention):
        self.path = path
        self.extention = extention
        self.gzipped = True if extention.endswith('.gz') else False
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
                print(f'WARNING: {read.path} has {read.n_reads} reads and sub-sampling at '
                      f'{depth}X depth would require {self.n_subsampled_reads} reads', file=sys.stderr)
            else:
                self.path = subsample_reads(read.path, f'{out_dir}/{self.depth_string}',
                                            self.n_subsampled_reads, read.gzipped)


if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    if args.out and not args.depths:
        quit_with_error('Can only output sub-sampled reads at selected depths with the --depths option')

    if args.out and not args.seed:
        quit_with_error('Can only randomly sub-sample reads with a random seed, try --seed 100')

    args.threads = min(args.threads, cpu_count())
    seed(args.seed)
    samples = process_reads(args.reads, args.assemblies)
    header = f'Sample\tAssembly_file\tGenome_size\tTotal_reads\tTotal_bases\tAverage_insert_size\tDepth\tRead_file\t' \
             f'Insert_size\tN_reads\tN_bases'
    if args.depths:
        header += '\t' + "\t".join([f"N_reads_at_{depth}X_depth" for depth in args.depths])
    print(header)

    executor = futures.ProcessPoolExecutor(args.threads)
    ftures = [executor.submit(sample.calculate_coverage(args.depths, args.out)) for sample in samples]
    futures.wait(ftures)
    sys.exit(0)
