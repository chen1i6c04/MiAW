import os
import argparse
import subprocess
from tempfile import TemporaryDirectory
from src.logs import create_logger
from src.cmd import FastqcCommandline, KmerFinderCommandline, AMRFinderCommandline, PlasmidFinderCommandline,\
    ShovillCommandline

CURRENT_DIR = os.path.dirname(os.path.realpath(__file__))
DB = os.path.join(CURRENT_DIR, 'db')
BIN_DIR = os.path.join(CURRENT_DIR, 'bin')

# trimmomatic
ADAPTERS = os.path.join(DB, 'trimmomatic.fa')
MIN_BQ = 3
TRIM_OPT = f"ILLUMINACLIP:{ADAPTERS}:2:30:10 LEADING:{MIN_BQ} TRAILING:{MIN_BQ} SLIDINGWINDOW:4:20 MINLEN:36 TOPHRED33"


def trimming(forward_reads, reverse_reads, outdir, threads):
    forward_paired = os.path.join(outdir, 'paired_R1.fq.gz')
    reverse_paired = os.path.join(outdir, 'paired_R2.fq.gz')
    cmd = f"trimmomatic PE -threads {threads} {forward_reads} {reverse_reads} {forward_paired} /dev/null" \
          f" {reverse_paired} /dev/null {TRIM_OPT}"
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    stdout_str, stderr_str = p.communicate()
    return_code = p.returncode
    if return_code:
        raise subprocess.CalledProcessError(return_code, cmd, stdout_str, stderr_str)
    return forward_paired, reverse_paired


def workflow(forward_reads, reverse_reads, outdir, threads):
    logger = create_logger(os.path.join(outdir, 'pga.log'))

    logger.info("Quality control checks.")
    cline = FastqcCommandline(forward_reads, reverse_reads, outdir=outdir, threads=2)
    cline()

    logger.info("Species identification.")
    tmp = TemporaryDirectory(prefix='kmerfinder_', dir=outdir)
    tmp_dir = tmp.name
    cline = KmerFinderCommandline(infile=" ".join((forward_reads, reverse_reads)), outdir=tmp_dir,
                                  database=os.path.join(DB, 'kmerfinder_db', 'bacteria.ATG'),
                                  taxonomy_file=os.path.join(DB, 'kmerfinder_db', 'bacteria.tax'),
                                  extented_output=True)
    cline()
    os.renames(os.path.join(tmp_dir, 'results.txt'), os.path.join(outdir, 'kmerfinder.txt'))

    logger.info("Trimming reads.")
    tmp = TemporaryDirectory(prefix='trimmomatic_', dir=outdir)
    tmp_dir = tmp.name
    forward_paired, reverse_paired = trimming(forward_reads, reverse_reads, tmp_dir, threads)

    logger.info("Assembling reads with SPAdes.")
    cline = ShovillCommandline(r1=forward_paired, r2=reverse_paired, outdir=outdir, depth=0, tmpdir='/dev/shm',
                               cpus=threads, ram=32, nostitch=True, force=True)
    cline()

    contigs_file = os.path.join(outdir, 'contigs.fa')

    logger.info("Plasmid type prediction.")
    tmp = TemporaryDirectory(prefix='plasmidfinder_', dir=outdir)
    tmp_dir = tmp.name
    cline = PlasmidFinderCommandline(infile=contigs_file, outdir=tmp_dir,
                                     database=os.path.join(DB, 'plasmidfinder_db'),
                                     tmp=tmp_dir, extented_output=True)
    cline()
    os.renames(os.path.join(tmp_dir, 'results_tab.tsv'), os.path.join(outdir, 'plasmidfinder.txt'))

    logger.info("Identifies AMR genes.")
    cline = AMRFinderCommandline(cmd=os.path.join(BIN_DIR, 'amr', 'amrfinder'),
                                 nuc_fasta=contigs_file, output_file=os.path.join(outdir, 'amrfinder.txt'), threads=threads)
    cline()
    logger.info("Done.")


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
