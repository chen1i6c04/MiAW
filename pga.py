import os
import argparse
from tempfile import TemporaryDirectory
from src.cmd import *

CURRENT_DIR = os.path.dirname(os.path.realpath(__file__))
DB = os.path.join(CURRENT_DIR, 'database')
BIN_DIR = os.path.join(CURRENT_DIR, 'bin')


def workflow(forward_reads, reverse_reads, outdir, threads):
    cline = FastqcCommandline(forward_reads, reverse_reads, outdir=outdir, threads=2)
    cline()

    tmp = TemporaryDirectory(prefix='kmerfinder_', dir=outdir)
    tmp_dir = tmp.name
    cline = KmerFinderCommandline(infile=" ".join((forward_reads, reverse_reads)), outdir=tmp_dir,
                                  database=os.path.join(DB, 'kmerfinder_db', 'bacteria.ATG'),
                                  taxonomy_file=os.path.join(DB, 'kmerfinder_db', 'bacteria.tax'),
                                  extented_output=True)
    cline()
    os.renames(os.path.join(tmp_dir, 'results.txt'), os.path.join(outdir, 'kmerfinder.txt'))

    tmp = TemporaryDirectory(prefix='spades_', dir=outdir)
    tmp_dir = tmp.name
    cline = SpadesCommandline(reads_1=forward_reads, reads_2=reverse_reads, output_dir=tmp_dir, tmp='/dev/shm',
                              threads=threads, careful=True)
    cline()
    seqfile = os.path.join(outdir, 'spades.fa')
    os.renames(os.path.join(tmp_dir, 'contigs.fasta'), seqfile)

    tmp = TemporaryDirectory(prefix='plasmidfinder_', dir=outdir)
    tmp_dir = tmp.name
    cline = PlasmidFinderCommandline(infile=seqfile, outdir=tmp_dir,
                                     database=os.path.join(DB, 'plasmidfinder_db'),
                                     tmp=tmp_dir, extented_output=True)
    cline()
    os.renames(os.path.join(tmp_dir, 'results_tab.tsv'), os.path.join(outdir, 'plasmidfinder.txt'))

    cline = AMRFinderCommandline(cmd=os.path.join(BIN_DIR, 'amr', 'amrfinder'),
                                 nuc_fasta=seqfile, output_file=os.path.join(outdir, 'amrfinder.txt'), threads=threads)
    cline()


def main():
    parser = argparse.ArgumentParser(usage='This program is for analyze MiSeq output.')
    parser.add_argument("-r1", required=True, help="file with forward paired-end reads")
    parser.add_argument("-r2", required=True, help="file with reverse paired-end reads")
    parser.add_argument("-o", "--outdir", required=True, help="directory to store all the resulting files")
    parser.add_argument("-t", "--threads", default=8, type=int, required=False, help="number of threads [Default 8]")
    args = parser.parse_args()

    forward_reads = args.r1
    reverse_reads = args.r2
    outdir = args.outdir
    threads = args.threads

    os.makedirs(outdir, exist_ok=True)
    workflow(forward_reads, reverse_reads, outdir, threads)


if __name__ == '__main__':
    main()