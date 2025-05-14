__author__ = 'B.H Chen'
__email__ = 'chen1i6c04@gmail.com'

import os
import re
import sys
import gzip
import shutil
import argparse
import subprocess
from tempfile import TemporaryDirectory, gettempdir
from loguru import logger
from Bio import SeqIO
from Bio.SeqUtils.lcc import lcc_simp
from src.database import check_database_installation

current_location = os.path.dirname(os.path.abspath(__file__))


def estimate_genome_size(filename, threads):
    with TemporaryDirectory() as tmpdir:
        cmd = (f"kmc -sm -t{threads} -k21 -ci10 {filename} {tmpdir}/kmc {tmpdir} | "
               f"grep -oP 'No. of unique counted k-mers\s+:\s+\K\d+'")
        child_process = syscall(cmd, stdout=True)
    return int(child_process.stdout)


def parse_genome_size(string):
    unit_map = {'K': 1e3, 'M': 1e6, 'G': 1e9, 'k': 1e3, 'm': 1e6, 'g': 1e9}
    result = re.fullmatch(r'^([\d.]+)([KMGkmg])', string)
    if result is None:
        logger.error(f"Couldn't parse {string}")
    else:
        value, unit = result.groups()
        return int(float(value) * unit_map[unit])


def count_bases(filename):
    output = subprocess.getoutput(
        f"seqtk fqchk {filename} | grep ALL | cut -f 2",
    )
    return int(output.strip())


def validate_fastq(file):
    """Checks the input fastq is really a fastq
    """
    def is_gzip(f):
        with open(f, 'rb') as f:
            signature = f.read(2)
        return signature == b'\x1f\x8b'
    # to get extension
    open_function = gzip.open if is_gzip(file) else open
    with open_function(file, 'rt') as handle:
        if any(SeqIO.parse(handle, 'fastq')):
            logger.info(f"FASTQ {file} checked")
        else:
            logger.error(f"Input file {file} is not in the FASTQ format.")
            sys.exit('Abort')


def syscall(cmd, stdout=False, stderr=False):
    if stdout:
        stdout_str = subprocess.PIPE
    else:
        stdout_str = None
    if stderr:
        stderr_str = subprocess.PIPE
    else:
        stderr_str = None
    shell = True if isinstance(cmd, str) else False
    child_process = subprocess.run(
        cmd, shell=shell, stdout=stdout_str, stderr=stderr_str, universal_newlines=True,
    )
    if child_process.returncode:
        logger.error(f"Command '{cmd}' failed")
        sys.exit('Abort')
    return child_process


def run_fastp(input_1, input_2, outdir, threads):
    output_1 = os.path.join(outdir, 'R1.trim.fastq.gz')
    output_2 = os.path.join(outdir, 'R2.trim.fastq.gz')
    json_report = os.path.join(outdir, 'fastp.json')
    html_report = os.path.join(outdir, 'fastp.html')
    cmd = (f"fastp -i {input_1} -I {input_2} -o {output_1} -O {output_2} -l 50 -w {16 if threads > 16 else threads} "
           f"--detect_adapter_for_pe -t 1 -z 6 -5 -3 -j {json_report} -h {html_report}")
    logger.info(f"Running: {cmd}")
    syscall(cmd)
    return output_1, output_2


def run_spades(paired_1, paired_2, outdir, threads):
    assembly = os.path.join(outdir, 'contigs.fasta')
    cmd = (f"spades.py -1 {paired_1} -2 {paired_2} -o {outdir} -t {threads} --tmp-dir {gettempdir()} --cov-cutoff auto "
           f"--only-assembler --isolate")
    logger.info(f"Running: {cmd}")
    syscall(cmd)
    if not os.path.exists(assembly):
        logger.error("Assembly fail!")
    return assembly


def rename_and_filter_contigs(src, dst):
    with open(dst, 'w') as handle:
        for record in SeqIO.parse(src, 'fasta'):
            node, length, coverage = re.search('NODE_(\d+)_length_(\d+)_cov_(\d*\.*\d+)', record.id).groups()
            node, length, coverage = int(node), int(length), float(coverage)
            if length < 200:
                logger.info(f"Removing short contig (< 200 bp) : {record.id}")
            elif coverage < 2:
                logger.info(f"Removing low coverage contig (< 2 x) : {record.id}")
            elif lcc_simp(record.seq) == 0:
                logger.info(f"Removing homopolymer contig : {record.id}")
            else:
                origname = record.id
                record.id = f"contig{node:04}"
                record.description = f"length={length} cov={coverage:.2f} origname={origname}"
                SeqIO.write(record, handle, 'fasta')


def extract_unmapped_sequences(input_1, input_2, output_1, output_2, target, threads):
    cmd = (f"bwa-mem2 mem -t {threads} {target} {input_1} {input_2} | "
           f"samtools view -@ {threads} -b -h -f 12 - | "
           f"samtools fastq -@ {threads} -c 6 -1 {output_1} -2 {output_2} -0 {os.devnull} -s {os.devnull} -n -")
    syscall(cmd)


def run_kraken2_and_bracken(short_one, short_two, database, outdir, threads):
    kraken2_report = os.path.join(outdir, 'kraken2.txt')
    bracken_report = os.path.join(outdir, 'bracken.txt')
    target = os.path.join(current_location, 'database', 'ncbi_plasmids')

    with TemporaryDirectory(dir=outdir) as tmpdir:
        unmapped_one = os.path.join(tmpdir, 'unmapped_1.fastq.gz')
        unmapped_two = os.path.join(tmpdir, 'unmapped_2.fastq.gz')
        extract_unmapped_sequences(short_one, short_two, unmapped_one, unmapped_two, target, threads)
        kraken2_cmd = (f"kraken2 --db {database} --output {os.devnull} --threads {threads} --confidence 0.6 "
                       f"--report {kraken2_report} --memory-mapping --paired {unmapped_one} {unmapped_two}")
        logger.info(f"Running: {kraken2_cmd}")
        syscall(kraken2_cmd)
    bracken_cmd = f"bracken -i {kraken2_report} -d {database} -w /dev/null -o {bracken_report} -l G"
    logger.info(f"Running: {bracken_cmd}")
    syscall(bracken_cmd)


def run_busco(input_file, output_dir, database, num_threads=1):
    specific_pattern = re.compile("^short_summary.specific.*.busco.json$")
    generic_pattern = re.compile("^short_summary.generic.*.busco.json$")
    with TemporaryDirectory() as tmpdir:
        cmd = f"busco -o busco --out_path {tmpdir} -m geno --offline --download_path {database} " \
              f"-i {input_file} -c {num_threads} --auto-lineage-prok"
        logger.info(f"Running: {cmd}")
        syscall(cmd)
        out_path = os.path.join(tmpdir, "busco")
        for f in os.listdir(out_path):
            match = specific_pattern.fullmatch(f)
            if match:
                shutil.copyfile(
                    os.path.join(out_path, f),
                    os.path.join(output_dir, 'busco_specific_summary.json')
                )
            match = generic_pattern.fullmatch(f)
            if match:
                shutil.copyfile(
                    os.path.join(out_path, f),
                    os.path.join(output_dir, 'busco_generic_summary.json')
                )


def run_skani(input_file, output_dir, database, num_threads=1):
    cmd = f"skani search -q {input_file} -d {database} -o {output_dir} -n 5 -t {num_threads}"
    syscall(cmd)


def check_dependency():
    version = {
        "FastQC": "fastqc -v",
        "fastp": "fastp -v 2>&1",
        "Kraken2": "kraken2 -v | grep version",
        "Bracken": "grep VERSION $(which bracken)",
        "SPAdes": "spades.py -v",
        "KMC": "kmc | grep K-Mer",
        "seqtk": "seqtk 2>&1 | grep Version",
        "BUSCO": "busco -v",
        "bwa-mem2": "bwa-mem2 version",
        "samtools": "samtools version | head -n 1",
        "skANI": "skani -V",
        "rasusa": "rasusa -V",
        "pypolca": "pypolca -V"
    }
    for program_name, cmd in version.items():
        child_process = syscall(cmd, stdout=True)
        if child_process.returncode:
            logger.error(f"Could not determine version of {program_name}")
            sys.exit("Abort")
        else:
            version = child_process.stdout.strip()
            logger.info(f"Using {program_name:8} | {version}")


def main():
    parser = argparse.ArgumentParser(prog="MiAW", description='Central lab MiSeq Analysis Workflow')
    parser.add_argument("-1", "--short_1", required=True, help="file with forward paired-end reads")
    parser.add_argument("-2", "--short_2", required=True, help="file with reverse paired-end reads")
    parser.add_argument("-o", "--outdir", required=True, help="directory to store all the resulting files")
    parser.add_argument("-t", "--threads", default=8, type=int, required=False, help="number of threads [Default 8]")
    parser.add_argument('-g', '--genome-size', default=None, 
                        help='Estimated genome size eg. 3.2M or 5m <blank=AUTO> (default: "")')
    parser.add_argument('-x', '--depth', default=100, type=int,
                        help='Sub-sample reads to this depth. Disable with --depth 0 (default: 50)')
    parser.add_argument(
        "--kraken2_db",
        default='',
        required=False,
        help="Path of kraken2/bracken database"
    )
    parser.add_argument(
        "--busco_db",
        default='',
        required=False,
        help="Path of busco database"
    )
    parser.add_argument(
        "--skani_db",
        default='',
        required=False,
        help="Path of skANI database"
    )
    parser.add_argument(
        "--disable-fastqc", action="store_true", help="Disable FastQC [default: OFF]"
    )

    args = parser.parse_args()
    outdir = args.outdir
    threads = args.threads

    os.makedirs(outdir, exist_ok=True)
    logfile = os.path.join(outdir, 'maiw.log')
    fmt = "{time:YYYY-MM-DD HH:mm:ss} [{level}] {message}"
    logger.add(logfile, format=fmt, level='INFO')
    logger.add(sys.stderr, format=fmt, level='ERROR')
    logger.add(lambda _: sys.exit(1), level="ERROR")

    logger.info("Checking dependencies")
    check_dependency()
    check_database_installation(os.path.join(current_location, 'database'))
    
    validate_fastq(args.short_1)
    validate_fastq(args.short_2)

    spades_output = os.path.join(outdir, 'spades')
    final_assembly = os.path.join(outdir, 'assembly.fasta')

    if args.disable_fastqc is False:
        logger.info("Quality control checks")
        syscall(f'fastqc -o {outdir} -t 2 {args.short_1} {args.short_2}')

    logger.info("Trim raw-reads.")
    trim_1, trim_2 = run_fastp(args.short_1, args.short_2, args.outdir, threads)
    if args.kraken2_db:
        logger.info("Running Kraken2 & Bracken")
        logger.info(f"Kraken2 database is {args.kraken2_db}")
        run_kraken2_and_bracken(trim_1, trim_2, args.kraken2_db, outdir, threads)
    total_bases = count_bases(trim_1) + count_bases(trim_2)

    if args.genome_size:
        gsize = parse_genome_size(args.genome_size)
        logger.info(f"Using genome size was {gsize}bp.")
    else:
        gsize = estimate_genome_size(trim_1, threads)
        logger.info(f"Estimated genome size was {gsize}bp.")
    origin_depth = int(total_bases / gsize)
    logger.info(f"Estimated sequencing depth: {origin_depth}x.")

    os.makedirs(spades_output, exist_ok=True)
    if args.depth:
        if origin_depth > args.depth:
            fraction = args.depth / origin_depth
            logger.info(f"Subsampling reads by factor {fraction:.3f} to get from {origin_depth}x to {args.depth}x")
            sub_1 = os.path.join(spades_output, 'R1.sub.fastq.gz')
            sub_2 = os.path.join(spades_output, 'R2.sub.fastq.gz')
            syscall(f"rasusa reads -f {fraction} -O g -o {sub_1} -o {sub_2} {trim_1} {trim_2}")
            paired_1, paired_2 = sub_1, sub_2
        else:
            logger.info("No read depth reduction requested or necessary.")
            paired_1, paired_2 = trim_1, trim_2
    else:
        paired_1, paired_2 = trim_1, trim_2
    logger.info("Running SPAdes")
    spades_assembly = run_spades(paired_1, paired_2, spades_output, threads)
    shutil.copyfile(spades_assembly, os.path.join(outdir, 'spades.fasta'))
    pypolca_output = os.path.join(outdir, 'pypolca')
    logger.info("Running pypolca")
    syscall(f"pypolca run -a {spades_assembly} -1 {paired_1} -2 {paired_2} -t {threads} -o {pypolca_output} --careful")
    rename_and_filter_contigs(
        os.path.join(pypolca_output, 'pypolca_corrected.fasta'),
        final_assembly,
    )
    shutil.rmtree(spades_output)

    if args.busco_db:
        logger.info("Running busco")
        run_busco(final_assembly, outdir, args.busco_db, threads)

    if args.skani_db:
        logger.info("Running skANI")
        run_skani(final_assembly, os.path.join(outdir, 'skani.txt'), args.skani_db, threads)

    os.remove(trim_1)
    os.remove(trim_2)
    logger.info("Done")


if __name__ == '__main__':
    main()

