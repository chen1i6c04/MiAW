import os
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
            f'--careful')


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
    parser.add_argument("-r1", required=True, help="file with forward paired-end reads")
    parser.add_argument("-r2", required=True, help="file with reverse paired-end reads")
    parser.add_argument("-o", "--outdir", required=True, help="directory to store all the resulting files")
    parser.add_argument("-t", "--threads", default=8, type=int, required=False, help="number of threads [Default 8]")
    parser.add_argument("--tmp-dir", default='/tmp', required=False,
                        help="directory for temporary files. [Default /tmp]")
    args = parser.parse_args()

    global TMPDIR
    TMPDIR = args.tmp_dir

    os.makedirs(args.outdir, exist_ok=True)
    logger = create_logger(os.path.join(args.outdir, 'pga.log'))

    logger.info("Quality control checks.")
    run_cmd(f'fastqc -o {args.outdir} -t {args.threads} {args.r1} {args.r2}')

    logger.info("Species identification.")
    species_identify(f'{args.r1} {args.r2}', args.outdir)

    logger.info("Assembling reads with SPAdes.")
    with TemporaryDirectory(dir=args.outdir) as outdir:
        assembly(args.r1, args.r2, outdir, args.threads)
        shutil.copy(os.path.join(outdir, 'contigs.fasta'), args.outdir)
    seqfile = os.path.join(args.outdir, 'contigs.fa')

    logger.info("Plasmid type prediction.")
    plasmid_identify(seqfile, args.outdir)

    logger.info("Identifies AMR genes.")
    output_file = os.path.join(args.outdir, 'amrfinder.txt')
    resistance_gene_detection(seqfile, output_file, args.threads)
    logger.info("Done.")


if __name__ == '__main__':
    main()
