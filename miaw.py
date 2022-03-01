import os
import shutil
import logging
import argparse
import subprocess
from tempfile import TemporaryDirectory
from Bio import SeqIO
from Bio.SeqUtils.lcc import lcc_simp
from src.assembly import assembly_pipeline
from src.utils import estimate_genome_size, count_bases


class LoggerFactory:

    FMT = "(PID:%(process)d)%(asctime)-20s[%(levelname)s] %(message)s"
    DATEFMT = "%Y-%m-%d %H:%M:%S"

    def __init__(self):
        self._logger = logging.getLogger(__name__)
        self._logger.setLevel(logging.INFO)

    def addLogBoxHandler(self, logbox):
        log_handler = logging.Handler(logbox)
        formatter = logging.Formatter(self.FMT, self.DATEFMT)
        log_handler.setFormatter(formatter)
        self._logger.addHandler(log_handler)

    def addFileHandler(self, logfile):
        file_handler = logging.FileHandler(logfile, mode="a", encoding=None, delay=False)
        formatter = logging.Formatter(self.FMT, self.DATEFMT)
        file_handler.setFormatter(formatter)
        self._logger.addHandler(file_handler)

    def addConsoleHandler(self):
        stream_handler = logging.StreamHandler()
        stream_handler.setLevel(logging.INFO)
        formatter = logging.Formatter(self.FMT, self.DATEFMT)
        stream_handler.setFormatter(formatter)
        self._logger.addHandler(stream_handler)

    def redircet_stdout(self):
        pass

    def create(self):
        return self._logger


def trimming(read_1, read_2, outdir, threads):
    paired_1 = os.path.join(outdir, 'R1.fastq.gz')
    paired_2 = os.path.join(outdir, 'R2.fastq.gz')
    json_report = os.path.join(outdir, 'fastp.json')
    html_report = os.path.join(outdir, 'fastp.html')
    cmd = [
        'fastp',
        '-i', read_1,
        '-I', read_2,
        '-o', paired_1,
        '-O', paired_2,
        '--length_required', '36',
        '--cut_front', '3',
        '--thread', str(threads),
        '--detect_adapter_for_pe',
        '-j', json_report,
        '-h', html_report,
        '-t', '1',
    ]
    subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return paired_1, paired_2


def species_identify(input_files, database, outfile):
    idx_file = database + '.ATG'
    tax_file = database + '.tax'
    with TemporaryDirectory() as tmpdir:
        subprocess.run(
            ['kmerfinder.py', '-i', input_files, '-o', tmpdir, '-db', idx_file, '-tax', tax_file, '-x']
        )
        shutil.copyfile(os.path.join(tmpdir, 'results.txt'), outfile)


def remove_low_complexity_sequence(source, destination):
    with open(destination, 'w') as out_f:
        for record in SeqIO.parse(source, 'fasta'):
            if lcc_simp(record.seq) > 0:
                SeqIO.write(record, out_f, 'fasta')


def main():
    parser = argparse.ArgumentParser(prog="MiAW", description='Central lab MiSeq Analysis Workflow')
    parser.add_argument("-1", "--short_1", required=True, help="file with forward paired-end reads")
    parser.add_argument("-2", "--short_2", required=True, help="file with reverse paired-end reads")
    parser.add_argument("-o", "--outdir", required=True, help="directory to store all the resulting files")
    parser.add_argument("-t", "--threads", default=8, type=int, required=False, help="number of threads [Default 8]")
    parser.add_argument("--tax_db", default='', required=False, help="Path of kmerfinder database")
    parser.add_argument("--qc", action="store_true", help="Quality check with FastQC [Default OFF]")

    args = parser.parse_args()
    outdir = args.outdir
    threads = args.threads

    os.makedirs(outdir, exist_ok=True)
    logfile = os.path.join(outdir, 'maiw.log')
    logger_factory = LoggerFactory()
    logger_factory.addFileHandler(logfile)
    logger_factory.addConsoleHandler()
    logger = logger_factory.create()

    assembly_dirname = os.path.join(outdir, 'spades')
    kmerfinder_filename = os.path.join(outdir, 'kmerfinder.txt')
    os.makedirs(assembly_dirname, exist_ok=True)

    if args.qc:
        logger.info("Quality control checks")
        subprocess.run(['fastqc', '-o', outdir, '-t', '2', args.short_1, args.short_2])

    if args.tax_db:
        logger.info("Identify species")
        species_identify(f'{args.short_1} {args.short_2}', args.tax_db, kmerfinder_filename)

    logger.info("Trim raw-reads.")
    r1, r2 = trimming(args.short_1, args.short_2, args.outdir, threads)
    total_bases = count_bases(r1) + count_bases(r2)
    gsize = estimate_genome_size(r1, threads)
    logger.info(f"Estimated genome size was {gsize}bp.")
    origin_depth = int(total_bases / gsize)
    logger.info(f"Estimated sequencing depth: {origin_depth}x.")

    depth = 100
    if origin_depth > depth:
        fraction = depth / origin_depth
        logger.info(f"Subsampling reads by factor {fraction:.3f} to get from {origin_depth}x to {depth}x")
        sub_r1 = os.path.join(assembly_dirname, 'R1.sub.fastq')
        with open(sub_r1, 'w+b') as handle:
            subprocess.run(['seqtk', 'sample', r1, str(fraction)], stdout=handle)
        sub_r2 = os.path.join(assembly_dirname, 'R2.sub.fastq')
        with open(sub_r2, 'w+b') as handle:
            subprocess.run(['seqtk', 'sample', r2, str(fraction)], stdout=handle)
        _r1, _r2 = sub_r1, sub_r2
    else:
        _r1, _r2 = r1, r2

    logger.info("Assembly with SPAdes")
    result_dirname = assembly_pipeline(_r1, _r2, assembly_dirname, threads)
    remove_low_complexity_sequence(
        os.path.join(result_dirname, 'pilon.fasta'),
        os.path.join(outdir, 'assembly.fasta'),
    )
    shutil.rmtree(assembly_dirname)
    # os.remove(r1)
    # os.remove(r2)
    logger.info("Done")


if __name__ == '__main__':
    main()
