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
    # for start in range(len(content_gaps)):
    #     end = start + window_size
    #     window = content_gaps[start: end]
    #     if max(window) < gap:
    #         head_crop = start
    #         break
    # else:
    #     head_crop = 0
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


def species_identify(input_files, outdir):
    database = os.path.join(DB, 'kmerfinder_db', 'bacteria.ATG')
    taxonomy_file = os.path.join(DB, 'kmerfinder_db', 'bacteria.tax')
    with TemporaryDirectory(dir=TMPDIR) as tmpdir:
        run_cmd(f'kmerfinder.py -i {input_files} -o {tmpdir} -db {database} -tax {taxonomy_file} -x')
        shutil.copy(os.path.join(tmpdir, 'results.txt'), os.path.join(outdir, 'kmerfinder.txt'))


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
    parser = argparse.ArgumentParser(prog="Pathogen Analyze", description='This program is for analyze MiSeq output.')
    parser.add_argument("-1", "--short_1", required=True, help="file with forward paired-end reads")
    parser.add_argument("-2", "--short_2", required=True, help="file with reverse paired-end reads")
    parser.add_argument("-o", "--outdir", required=True, help="directory to store all the resulting files")
    parser.add_argument("-t", "--threads", default=8, type=int, required=False, help="number of threads [Default 8]")
    parser.add_argument("--qc", action="store_true", help="Quality check with FastQC")
    parser.add_argument("--id", action="store_true", help="Predicted species with Kmerfinder")
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

    if args.qc:
        logger.info("Quality control checks")
        run_cmd(f'fastqc -o {outdir} -t 2 {args.short_1} {args.short_2}')

    if args.id:
        logger.info("Identify species")
        species_identify(f'{args.short_1} {args.short_2}', outdir)

    spades_dir = os.path.join(outdir, 'spades_assembly')
    os.makedirs(spades_dir)

    logger.info("Trim raw-reads.")

    trim_short_1, trim_short_2 = trimming(args.short_1, args.short_2, spades_dir, threads)
    total_bases = count_bases(trim_short_1) + count_bases(trim_short_2)
    gsize = estimate_genome_size(trim_short_1, threads, TMPDIR)
    logger.info(f"Estimated genome size was {gsize}bp.")
    origin_depth = int(total_bases / gsize)
    logger.info(f"Estimated sequencing depth: {origin_depth}x.")

    depth = 80
    if origin_depth > depth:
        fraction = depth / origin_depth
        logger.info(f"Subsampling reads by factor {fraction:.3f} to get from {origin_depth}x to {depth}x")
        sub_short_1 = os.path.join(spades_dir, 'sub_R1.fastq')
        run_cmd(f"seqtk sample {trim_short_1} {fraction} > {sub_short_1}")
        sub_short_2 = os.path.join(spades_dir, 'sub_R2.fastq')
        run_cmd(f"seqtk sample {trim_short_2} {fraction} > {sub_short_2}")
        _short_1, _short_2 = sub_short_1, sub_short_2
    else:
        _short_1, _short_2 = trim_short_1, trim_short_2

    logger.info("SPAdes assemblies")
    spades_cmd = spadesCommandline(short_1=_short_1, short_2=_short_2, outdir=spades_dir, tmpdir=TMPDIR,
                                   threads=threads, disable_gzip_output=True, cov_cutoff='auto', only_assembler=True)
    run_cmd(spades_cmd)
    spades_output = os.path.join(spades_dir, 'contigs.fasta')

    logger.info("Polishing assembly with Pilon.")
    pilon_dir = os.path.join(spades_dir, 'pilon_polish')
    change_count = pilon_polish(_short_1, _short_2, spades_output, pilon_dir, threads)
    logger.info(f"Total number of changes: {change_count}")
    pilon_output = os.path.join(pilon_dir, 'pilon.fasta')
    finall_assembly = os.path.join(outdir, 'assembly.fasta')
    shutil.copyfile(pilon_output, finall_assembly)
    shutil.copy(os.path.join(spades_dir, 'spades.log'), outdir)
    shutil.rmtree(spades_dir)
    logger.info("Done")


if __name__ == '__main__':
    main()
