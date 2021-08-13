import os
import shutil
import logging
import argparse
from tempfile import TemporaryDirectory
from src.utils import run_cmd, estimate_genome_size, count_bases
from src.polish import pilon_polish


TMPDIR = str()
CURRENT_DIR = os.path.dirname(os.path.realpath(__file__))
DB = os.path.join(CURRENT_DIR, 'db')
BIN_DIR = os.path.join(CURRENT_DIR, 'bin')

# trimmomatic
ADAPTERS = os.path.join(DB, 'trimmomatic.fa')
MIN_BQ = 3


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


def crop_position(filepath, window_size=3, gap=10):
    p = run_cmd(f"seqtk fqchk {filepath}")
    lines = p.stdout.decode().strip().split('\n')[3:]
    per_base_content = (line.split()[2:6] for line in lines)
    content_gaps = []
    for line in per_base_content:
        a, c, g, t = [float(value) for value in line]
        content_gaps.append(max(abs(a - t), abs(c - g)))
    for start in range(len(content_gaps), 0, -1):
        end = start - window_size
        window = content_gaps[end:start]
        if max(window) < gap:
            crop = start
            break
    else:
        crop = len(content_gaps)
    return crop


def trimming(short_1, short_2, outdir, threads):
    opt = f"CROP:{crop_position(short_1)} ILLUMINACLIP:{ADAPTERS}:2:30:10 LEADING:{MIN_BQ} TRAILING:{MIN_BQ} " \
          f"SLIDINGWINDOW:4:20 MINLEN:36 TOPHRED33"
    r1 = os.path.join(outdir, 'R1.fq')
    r2 = os.path.join(outdir, 'R2.fq')
    cmd = f"trimmomatic PE -threads {threads} {short_1} {short_2} {r1} /dev/null {r2} /dev/null {opt}"
    run_cmd(cmd)
    return r1, r2


def species_identify(input_files, outfile):
    database = os.path.join(DB, 'kmerfinder_db', 'bacteria.ATG')
    taxonomy_file = os.path.join(DB, 'kmerfinder_db', 'bacteria.tax')
    with TemporaryDirectory(dir=TMPDIR) as tmpdir:
        run_cmd(f'kmerfinder.py -i {input_files} -o {tmpdir} -db {database} -tax {taxonomy_file} -x')
        shutil.copyfile(os.path.join(tmpdir, 'results.txt'), outfile)


def spadesCommandline(prog='spades.py', **kwargs):
    parameters = {
        'short_1': '-1',
        'short_2': '-2',
        'outdir': '-o',
        'threads': '--threads',
        'memory': '--memory',
        'tmpdir': '--tmp-dir',
        'disable_gzip_output': '--disable-gzip-output',
        'isolate': '--isolate',
        'careful': '--careful',
        'cov_cutoff': '--cov-cutoff',
        'only_assembler': '--only-assembler'
    }
    cmd = [prog]
    for key, value in kwargs.items():
        parameter = parameters.get(key)
        if parameter:
            if isinstance(value, bool) and value:
                cmd += [parameter]
            elif isinstance(value, bool) and value is False:
                pass
            else:
                cmd += [parameter, value]
        else:
            raise ValueError('Parameter %s is invalided variable' % key)
    cmd = [str(term) for term in cmd]
    cmd = ' '.join(cmd)
    return cmd


def main():
    parser = argparse.ArgumentParser(prog="MiWA", description='This program is for analyze MiSeq output.')
    parser.add_argument("-1", "--short_1", required=True, help="file with forward paired-end reads")
    parser.add_argument("-2", "--short_2", required=True, help="file with reverse paired-end reads")
    parser.add_argument("-o", "--outdir", required=True, help="directory to store all the resulting files")
    parser.add_argument("-t", "--threads", default=8, type=int, required=False, help="number of threads [Default 8]")
    parser.add_argument("--qc", action="store_true", help="Quality check with FastQC [Default OFF]")
    parser.add_argument("--id", action="store_true", help="Predicted species with Kmerfinder [Default OFF]")
    parser.add_argument("--tmp-dir", default='/tmp', required=False,
                        help="directory for temporary files. [Default /tmp]")
    args = parser.parse_args()
    outdir = args.outdir
    threads = args.threads

    global TMPDIR
    TMPDIR = args.tmp_dir

    os.makedirs(outdir, exist_ok=True)
    logfile = os.path.join(outdir, 'pga.log')
    logger_factory = LoggerFactory()
    logger_factory.addFileHandler(logfile)
    logger_factory.addConsoleHandler()
    logger = logger_factory.create()

    spades_dirname = os.path.join(outdir, 'spades')
    pilon_dirname = os.path.join(spades_dirname, 'pilon')
    kmerfinder_filename = os.path.join(outdir, 'kmerfinder.txt')
    spades_filename = os.path.join(spades_dirname, 'contigs.fasta')
    pilon_filename = os.path.join(pilon_dirname, 'pilon.fasta')
    os.makedirs(spades_dirname, exist_ok=True)

    if args.qc:
        logger.info("Quality control checks")
        run_cmd(f'fastqc -o {outdir} -t 2 {args.short_1} {args.short_2}')

    if args.id:
        logger.info("Identify species")
        species_identify(f'{args.short_1} {args.short_2}', kmerfinder_filename)

    logger.info("Trim raw-reads.")
    r1, r2 = trimming(args.short_1, args.short_2, spades_dirname, threads)
    total_bases = count_bases(r1) + count_bases(r2)
    gsize = estimate_genome_size(r1, threads, TMPDIR)
    logger.info(f"Estimated genome size was {gsize}bp.")
    origin_depth = int(total_bases / gsize)
    logger.info(f"Estimated sequencing depth: {origin_depth}x.")

    depth = 80
    if origin_depth > depth:
        fraction = depth / origin_depth
        logger.info(f"Subsampling reads by factor {fraction:.3f} to get from {origin_depth}x to {depth}x")
        sub_r1 = os.path.join(spades_dirname, 'R1.sub.fastq')
        run_cmd(f"seqtk sample {r1} {fraction} > {sub_r1}")
        sub_r2 = os.path.join(spades_dirname, 'R2.sub.fastq')
        run_cmd(f"seqtk sample {r2} {fraction} > {sub_r2}")
        _r1, _r2 = sub_r1, sub_r2
    else:
        _r1, _r2 = r1, r2

    logger.info("SPAdes assemblies")
    spades_cmd = spadesCommandline(short_1=_r1, short_2=_r2, outdir=spades_dirname, tmpdir=TMPDIR,
                                   threads=threads, disable_gzip_output=True, cov_cutoff='auto', only_assembler=True)
    run_cmd(spades_cmd)

    logger.info("Polishing assembly with Pilon.")
    change_count = pilon_polish(_r1, _r2, spades_filename, pilon_dirname, threads)
    logger.info(f"Total number of changes: {change_count}")
    asm_filename = os.path.join(outdir, 'assembly.fasta')
    shutil.copyfile(pilon_filename, asm_filename)
    shutil.rmtree(spades_dirname)
    logger.info("Done")


if __name__ == '__main__':
    main()
