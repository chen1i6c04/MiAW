import os
import sys
import shutil
import logging
import argparse
import subprocess
from tempfile import TemporaryDirectory
from Bio import SeqIO
from Bio.SeqUtils.lcc import lcc_simp
from src.assembly import assembly_pipeline



class LoggerFactory:

    FMT = "%(asctime)-20s[%(levelname)s] %(message)s"
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

    def create(self):
        return self._logger


def estimate_genome_size(filename, threads):
    with TemporaryDirectory() as tmpdir:
        output = subprocess.getoutput(
            f"kmc -sm -t{threads} -k21 -ci10 {filename} {tmpdir}/kmc {tmpdir} | "
            f"grep 'No. of unique counted k-mers' | "
            f"awk '{{print $NF}}'",
        )
    return int(output.strip().split()[-1])


def count_bases(filename):
    output = subprocess.getoutput(
        f"seqtk fqchk {filename} | grep ALL | awk '{{print $2}}'",
    )
    return int(output.strip())


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
    subprocess.run(cmd)
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


def taxonomic_classification(r1, r2, database, kraken2_report, bracken_report, threads):
    kraken2_cmd = f"kraken2 --db {database} --output /dev/null --threads {threads} " \
                  f"--report {kraken2_report} --memory-mapping --paired {r1} {r2}"
    bracken_cmd = f"bracken -i {kraken2_report} -d {database} -w /dev/null -r 300 -o {bracken_report}"
    cmd = kraken2_cmd + " && " + bracken_cmd
    subprocess.run(cmd, shell=True)


def check_dependency():
    version = {
        "FastQC" : "fastqc -v",
        "fastp"  : "fastp -v 2>&1",
        "Kraken2": "kraken2 -v | grep version",
        "Bracken": "grep VERSION $(which bracken)",
        "SPAdes" : "spades.py -v",
        "KMC"    : "kmc | grep K-Mer",
        "seqtk"  : "seqtk 2>&1 | grep Version",
        "Pilon"  : "pilon | grep version",
    }
    logger = logging.getLogger(__name__)
    for program_name, cmd in version.items():
        child_process = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE)
        if child_process.returncode:
            logger.error(msg=f"Could not determine version of {program_name}")
            sys.exit("Abort")
        else:
            version = child_process.stdout.decode().strip()
            logger.info(msg=f"Using {program_name:8} | {version}")


def main():
    parser = argparse.ArgumentParser(prog="MiAW", description='Central lab MiSeq Analysis Workflow')
    parser.add_argument("-1", "--short_1", required=True, help="file with forward paired-end reads")
    parser.add_argument("-2", "--short_2", required=True, help="file with reverse paired-end reads")
    parser.add_argument("-o", "--outdir", required=True, help="directory to store all the resulting files")
    parser.add_argument("-t", "--threads", default=8, type=int, required=False, help="number of threads [Default 8]")
    parser.add_argument("--tax_db", default='', required=False, help="Path of kraken2/bracken database")
    parser.add_argument(
        "--no-quality-check", action="store_true", help="Disable quality check [default: OFF]"
    )
    parser.add_argument(
        "--no-assembly", action="store_true", help="Disable denovo assembly [default: OFF]"
    )

    args = parser.parse_args()
    outdir = args.outdir
    threads = args.threads

    os.makedirs(outdir, exist_ok=True)
    logfile = os.path.join(outdir, 'maiw.log')
    logger_factory = LoggerFactory()
    logger_factory.addFileHandler(logfile)
    logger_factory.addConsoleHandler()
    logger = logger_factory.create()
    
    check_dependency()
    
    logger.info(f"""Illumina sequences
        {args.short_1}
        {args.short_2}""")

    assembly_dirname = os.path.join(outdir, 'spades')
    # kmerfinder_filename = os.path.join(outdir, 'kmerfinder.txt')
    kraken2_filename = os.path.join(outdir, 'kraken2.txt')
    bracken_filename = os.path.join(outdir, 'bracken.txt')

    if args.no_quality_check is False:
        logger.info("Quality control checks")
        subprocess.run(['fastqc', '-o', outdir, '-t', '2', args.short_1, args.short_2])

    if args.tax_db:
        # logger.info("Identify species")
        # species_identify(f'{args.short_1} {args.short_2}', args.tax_db, kmerfinder_filename)
        logger.info("Classify")
        taxonomic_classification(args.short_1, args.short_2, args.tax_db, kraken2_filename, bracken_filename, threads)

    logger.info("Trim raw-reads.")
    trim_r1, trim_r2 = trimming(args.short_1, args.short_2, args.outdir, threads)
    total_bases = count_bases(trim_r1) + count_bases(trim_r2)
    gsize = estimate_genome_size(trim_r1, threads)
    logger.info(f"Estimated genome size was {gsize}bp.")
    origin_depth = int(total_bases / gsize)
    logger.info(f"Estimated sequencing depth: {origin_depth}x.")

    if args.no_assembly is False:
        os.makedirs(assembly_dirname, exist_ok=True)
        depth = 100
        if origin_depth > depth:
            fraction = depth / origin_depth
            logger.info(f"Subsampling reads by factor {fraction:.3f} to get from {origin_depth}x to {depth}x")
            sub_r1 = os.path.join(assembly_dirname, 'R1.sub.fastq')
            subprocess.run(f"seqtk sample {trim_r1} {fraction} > {sub_r1}", shell=True)
            sub_r2 = os.path.join(assembly_dirname, 'R2.sub.fastq')
            subprocess.run(f"seqtk sample {trim_r2} {fraction} > {sub_r2}", shell=True)
            paired_1, paired_2 = sub_r1, sub_r2
        else:
            paired_1, paired_2 = trim_r1, trim_r2
        logger.info("Assembly with SPAdes")
        result_dirname = assembly_pipeline(paired_1, paired_2, assembly_dirname, threads)
        remove_low_complexity_sequence(
            os.path.join(result_dirname, 'pilon.fasta'),
            os.path.join(outdir, 'assembly.fasta'),
        )
        shutil.rmtree(assembly_dirname)
    os.remove(r1)
    os.remove(r2)
    logger.info("Done")


if __name__ == '__main__':
    main()
