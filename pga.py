import os
import re
import yaml
import shutil
import argparse
from tempfile import TemporaryDirectory
from src.logs import create_logger
from src.cmd import run_cmd


TMPDIR = str()
CURRENT_DIR = os.path.dirname(os.path.realpath(__file__))
DB = os.path.join(CURRENT_DIR, 'db')
BIN_DIR = os.path.join(CURRENT_DIR, 'bin')

# trimmomatic
ADAPTERS = os.path.join(DB, 'trimmomatic.fa')
MIN_BQ = 3
TRIM_OPT = f"ILLUMINACLIP:{ADAPTERS}:2:30:10 HEADCROP:10 LEADING:{MIN_BQ} TRAILING:{MIN_BQ} SLIDINGWINDOW:4:20 " \
           f"MINLEN:36 TOPHRED33"


def trimming(forward_reads, reverse_reads, outdir, threads, compress=False):
    extension = {True: '.fq.gz', False: '.fq'}
    forward_paired = os.path.join(outdir, 'paired_R1' + extension[compress])
    reverse_paired = os.path.join(outdir, 'paired_R2' + extension[compress])
    cmd = f"trimmomatic PE -threads {threads} {forward_reads} {reverse_reads} {forward_paired} /dev/null" \
          f" {reverse_paired} /dev/null {TRIM_OPT}"
    run_cmd(cmd)
    return forward_paired, reverse_paired


def estimate_genome_size(infile, threads):
    with TemporaryDirectory(dir=TMPDIR) as tmp:
        cmd = f'kmc -sm -t{threads} -k21 -ci10 {infile} {tmp}/kmc {tmp}'
        stdout, stderr = run_cmd(cmd)
    pattern = r'No. of unique counted k-mers\s+\W\s+(\d+)'
    result = re.search(pattern, stdout.decode())
    return int(result.group(1))


def count_bases(infile):
    stdout, stderr = run_cmd(f'seqtk fqchk {infile} | grep ALL')
    return int(stdout.decode().split()[1])


def species_identify(input_files, outdir):
    database = os.path.join(DB, 'kmerfinder_db', 'bacteria.ATG')
    taxonomy_file = os.path.join(DB, 'kmerfinder_db', 'bacteria.tax')
    with TemporaryDirectory(dir=TMPDIR) as tmp:
        run_cmd(f'kmerfinder.py -i {input_files} -o {tmp} -db {database} -tax {taxonomy_file} -x')
        shutil.copy(os.path.join(tmp, 'results.txt'), os.path.join(outdir, 'kmerfinder.txt'))


def correct_reads(forward_reads, reverse_reads, outdir, threads):
    run_cmd(f'spades.py -1 {forward_reads} -2 {reverse_reads} -o {outdir} -t {threads} --tmp-dir {TMPDIR} '
            f'--only-error-correction')
    yaml_file = os.path.join(outdir, 'corrected', 'corrected.yaml')
    with open(yaml_file) as handle:
        data = yaml.safe_load(handle)
    return data[0]['left reads'][0], data[0]['right reads'][0]


def assembly(forward_reads, reverse_reads, outdir, threads):
    forward_trim_reads, reverse_trim_reads = trimming(forward_reads, reverse_reads, outdir, threads)
    assembler = 'spades.py'
    run_cmd(f'{assembler} -1 {forward_trim_reads} -2 {reverse_trim_reads} -o {outdir} --tmp-dir {TMPDIR} -t {threads} '
            f'--careful --disable-gzip-output')


def plasmid_identify(seqfile, outdir):
    database = os.path.join(DB, 'plasmidfinder_db')
    with TemporaryDirectory(dir=TMPDIR) as tmp:
        run_cmd(f"plasmidfinder.py -i {seqfile} -o {tmp} -p {database} -tmp {tmp} -x")
        shutil.copy(os.path.join(tmp, 'results_tab.tsv'), os.path.join(outdir, 'plasmidfinder.txt'))


def resistance_gene_detection(seqfile, output_file, threads):
    database = os.path.join(DB, 'amr', 'latest')
    with TemporaryDirectory(dir=TMPDIR) as tmp:
        tmp_file = os.path.join(tmp, 'locus.faa')
        run_cmd(f'prodigal -i {seqfile} -a {tmp_file}')
        run_cmd(f'amrfinder -p {tmp_file} -o {output_file} -d {database} --threads {threads}')


def main():
    parser = argparse.ArgumentParser(usage='This program is for analyze MiSeq output.')
    parser.add_argument("-1", required=True, help="file with forward paired-end reads")
    parser.add_argument("-2", required=True, help="file with reverse paired-end reads")
    parser.add_argument("-o", "--outdir", required=True, help="directory to store all the resulting files")
    parser.add_argument("-t", "--threads", default=8, type=int, required=False, help="number of threads [Default 8]")
    parser.add_argument("--tmp-dir", default='/tmp', required=False,
                        help="directory for temporary files. [Default /tmp]")
    args = parser.parse_args().__dict__
    forward_reads, reverse_reads = args['1'], args['2']
    outdir = args['outdir']
    threads = args['threads']
    tmp_dir = args['tmp_dir']

    global TMPDIR
    TMPDIR = tmp_dir

    os.makedirs(outdir, exist_ok=True)
    logger = create_logger(os.path.join(outdir, 'pga.log'))

    logger.info("Quality control checks.")
    run_cmd(f'fastqc -o {outdir} -t 2 {forward_reads} {reverse_reads}')

    logger.info("Species identification.")
    species_identify(f'{forward_reads} {reverse_reads}', outdir)

    total_bases = count_bases(forward_reads) + count_bases(reverse_reads)
    gsize = estimate_genome_size(forward_reads, threads)
    depth = int(total_bases/gsize)

    with TemporaryDirectory(dir=outdir) as spades_dir:
        if depth > 100:
            factor = '%.3f' % (100 / depth)
            logger.info(f"Subsampling reads by factor {factor} to get from {depth}x to 100x")
            r1 = os.path.join(spades_dir, 'R1.sub.fq.gz')
            run_cmd(f"seqtk sample {forward_reads} {factor} | pigz --fast -c -p {threads} > {r1}")
            r2 = os.path.join(spades_dir, 'R2.sub.fq.gz')
            run_cmd(f"seqtk sample {reverse_reads} {factor} | pigz --fast -c -p {threads} > {r2}")
        else:
            r1, r2 = forward_reads, reverse_reads
        logger.info("Assembling reads with SPAdes.")
        assembly(r1, r2, spades_dir, threads)
        shutil.copy(os.path.join(spades_dir, 'contigs.fasta'), outdir)
    seqfile = os.path.join(outdir, 'contigs.fasta')

    logger.info("Plasmid type prediction.")
    plasmid_identify(seqfile, outdir)

    logger.info("Identifies AMR genes.")
    output_file = os.path.join(outdir, 'amrfinder.txt')
    resistance_gene_detection(seqfile, output_file, threads)
    logger.info("Done.")


if __name__ == '__main__':
    main()
