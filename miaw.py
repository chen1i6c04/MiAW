__author__ = 'B.H Chen'
__email__ = 'chen1i6c04@gmail.com'

import os
import re
import sys
import json
import shutil
import logging
import argparse
import subprocess
from tempfile import TemporaryDirectory
from loguru import logger
from Bio import SeqIO
from Bio.SeqUtils.lcc import lcc_simp


LOCATION = os.path.dirname(os.path.abspath(__file__))
DATABASE = os.path.join(LOCATION, 'database')


def estimate_genome_size(filename, threads):
    with TemporaryDirectory() as tmpdir:
        output = subprocess.getoutput(
            f"kmc -sm -t{threads} -k21 -ci10 {filename} {tmpdir}/kmc {tmpdir} | "
            f"grep 'No. of unique counted k-mers' | "
            f"awk '{{print $NF}}'",
        )
    return int(output.strip().split()[-1])


def parse_genome_size(pattern):
    unit_map = {'K': 1e3, 'M': 1e6, 'G': 1e9}
    result = re.fullmatch(r'^([\d.]+)([KMG])', pattern)
    if result is None:
        sys.exit(f"Couldn't parse {pattern}")
    else:
        value, unit = result.groups()
        return int(float(value) * unit_map[unit])


def count_bases(filename):
    output = subprocess.getoutput(
        f"seqtk fqchk {filename} | grep ALL | awk '{{print $2}}'",
    )
    return int(output.strip())


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
    cmd = ['fastp', '-i', input_1, '-I', input_2, '-o', output_1, '-O', output_2, '--length_required', '50',
           '--thread', str(threads), '--detect_adapter_for_pe', '-t', '1', '-z', '6', '-5', '-3',
           '-j', json_report, '-h', html_report]
    logger.info(f"fastp command: {' '.join(cmd)}")
    syscall(cmd)
    return output_1, output_2


def run_spades(paired_1, paired_2, outdir, threads):
    assembly = os.path.join(outdir, 'contigs.fasta')
    cmd = ['spades.py', '-1', paired_1, '-2', paired_2, '-o', outdir, '-t', str(threads), '--tmp-dir', '/tmp',
           '--cov-cutoff', 'auto', '--only-assembler', '--disable-gzip-output', '--isolate']
    syscall(cmd)
    return assembly


def remove_low_complexity_and_shorter_sequence(source, destination):
    with open(destination, 'w') as handle:
        for record in SeqIO.parse(source, 'fasta'):
            if lcc_simp(record.seq) > 0 and len(record.seq) >= 200:
                SeqIO.write(record, handle, 'fasta')
            else:
                logger.info(f"Removing contig {record.id}")


def remove_target_sequence(input_1, input_2, output_1, output_2, target, threads):
    cmd = f"bwa-mem2 mem -t {threads} {target} {input_1} {input_2} | " \
          f"samtools sort -@ {threads} -n -O BAM - | " \
          f"samtools view -@ {threads} -f 12 -O BAM - | " \
          f"samtools fastq -@ {threads} -c 6 -1 {output_1} -2 {output_2} -0 /dev/null -s /dev/null -n -"
    syscall(cmd)


def run_kraken2_and_bracken(r1, r2, database, outdir, threads):
    kraken2_report = os.path.join(outdir, 'kraken2.txt')
    bracken_report = os.path.join(outdir, 'bracken.txt')
    with TemporaryDirectory(dir=outdir) as tmpdir:
        query_1 = os.path.join(tmpdir, 'R1.fq.gz')
        query_2 = os.path.join(tmpdir, 'R2.fq.gz')
        index = os.path.join(DATABASE, 'ncbi_plasmids')
        remove_target_sequence(r1, r2, query_1, query_2, index, threads)
        kraken2_cmd = f"kraken2 --db {database} --output /dev/null --threads {threads} " \
                      f"--report {kraken2_report} --memory-mapping --paired {query_1} {query_2}"
        bracken_cmd = f"bracken -i {kraken2_report} -d {database} -w /dev/null -o {bracken_report} -l G"
        cmd = kraken2_cmd + " && " + bracken_cmd
        syscall(cmd)


def run_busco(input_file, output_dir, database, num_threads=1):
    specific_pattern = re.compile("^short_summary.specific.*.busco.json$")
    generic_pattern = re.compile("^short_summary.generic.*.busco.json$")
    with TemporaryDirectory() as tmpdir:
        cmd = f"busco -o busco --out_path {tmpdir} -m geno --offline --download_path {database} " \
              f"-i {input_file} -c {num_threads} --auto-lineage-prok"
        logger.info(f"busco command: {cmd}")
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
        "bwa-mem2": "bwa-mem2 version"
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
                        help='Estimated genome size eg. 3.2M <blank=AUTO> (default: "")')
    parser.add_argument('-x', '--depth', default=50, type=int,
                        help='Sub-sample reads to this depth. Disable with --depth 0 (default: 50)')
    parser.add_argument(
        "--kraken2_db",
        default='',
        required=False,
        help="Path of kraken2/bracken database. If not provide, won't run kraken2"
    )
    parser.add_argument(
        "--busco_db",
        default='',
        required=False,
        help="Path of busco database. If not provide, won't run busco"
    )
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
    fmt = "{time:YYYY-MM-DD HH:mm:ss} [{level}] {message}"
    logger.add(logfile, format=fmt, level='INFO')
    logger.add(sys.stderr, format=fmt, level='ERROR')
    
    check_dependency()
    
    logger.info(f"Read 1: {args.short_1}")
    logger.info(f"Read 2: {args.short_2}")

    spades_output = os.path.join(outdir, 'spades')
    final_assembly = os.path.join(outdir, 'assembly.fasta')

    if args.no_quality_check is False:
        logger.info("Quality control checks")
        syscall(f'fastqc -o {outdir} -t 2 {args.short_1} {args.short_2}')

    logger.info("Trim raw-reads.")
    trim_1, trim_2 = run_fastp(args.short_1, args.short_2, args.outdir, threads)
    if args.kraken2_db:
        logger.info("Running Kraken2/Bracken")
        logger.info(f"Kraken2/Bracken database is {args.kraken2_db}")
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

    if not args.no_assembly:
        os.makedirs(spades_output, exist_ok=True)
        if args.depth:
            if origin_depth > args.depth:
                fraction = args.depth / origin_depth
                logger.info(f"Subsampling reads by factor {fraction:.3f} to get from {origin_depth}x to {args.depth}x")
                sub_1 = os.path.join(spades_output, 'R1.sub.fastq')
                syscall(f"seqtk sample {trim_1} {fraction} > {sub_1}")
                sub_2 = os.path.join(spades_output, 'R2.sub.fastq')
                syscall(f"seqtk sample {trim_2} {fraction} > {sub_2}")
                paired_1, paired_2 = sub_1, sub_2
            else:
                logger.info("No read depth reduction requested or necessary.")
                paired_1, paired_2 = trim_1, trim_2
        else:
            paired_1, paired_2 = trim_1, trim_2
        logger.info("Assembly with SPAdes")
        spades_assembly = run_spades(paired_1, paired_2, spades_output, threads)
        remove_low_complexity_and_shorter_sequence(
            spades_assembly,
            final_assembly,
        )
        shutil.rmtree(spades_output)
        if args.busco_db:
            logger.info("Running busco")
            run_busco(final_assembly, outdir, args.busco_db, threads)
    os.remove(trim_1)
    os.remove(trim_2)
    logger.info("Done")


if __name__ == '__main__':
    main()
