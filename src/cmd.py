from .application import Option, Switch, AbstractCommandline


class SpadesCommandline(AbstractCommandline):
    """
    parameters
    reads_1           : file with forward paired-end reads
    reads_2           : file with reverse paired-end reads
    single_reads      : file with unpaired reads
    output_dir        : directory to store all the resulting files (required)
    tmp               : directory for temporary files
    threads           : number of threads
    memory            : RAM limit for SPAdes in Gb (terminates if exceeded) [default: 250]
    trusted_contigs   : file with trusted contigs
    untrusted_contigs : file with untrusted contigs
    careful           : tries to reduce number of mismatches and short indels [default: False]
    """
    def __init__(self, cmd='spades.py', **kwargs):
        self.parameters = [
            Option(['-1', 'reads_1']),
            Option(['-2', 'reads_2']),
            Option(['-s', 'single_reads']),
            Option(['-o', 'output_dir']),
            Option(['--tmp-dir', 'tmp']),
            Option(['-t', 'threads'], lambda x: isinstance(x, int)),
            Option(['-m', 'memory'], lambda x: isinstance(x, int)),
            Option(['--trusted-contigs', 'trusted_contigs']),
            Option(['--untrusted-contigs', 'untrusted_contigs']),
            Switch(['--careful', 'careful'])
        ]
        AbstractCommandline.__init__(self, cmd, **kwargs)


class FastqcCommandline(AbstractCommandline):
    """
    outdir  : Create all output files in the specified output directory.
              Please note that this directory must exist as the program
              will not create it.  If this option is not set then the
              output file for each sequence file is created in the same
              directory as the sequence file which was processed.
    threads : Specifies the number of files which can be processed simultaneously.
    quiet   : Supress all progress messages on stdout and only report errors.
    """
    def __init__(self, *args, **kwargs):
        self.parameters = [
            Option(['--outdir', 'outdir']),
            Option(['--threads', 'threads']),
            Switch(['--noextract', 'noextract']),
            Switch(['--extract', 'extract']),
        ]
        AbstractCommandline.__init__(self, 'fastqc', *args, **kwargs)


class KmerFinderCommandline(AbstractCommandline):
    def __init__(self, cmd='kmerfinder.py', **kwargs):
        self.parameters = [
            Option(['-i', 'infile']),
            Option(['-o', 'outdir']),
            Option(['-db', 'database']),
            Option(['-tax', 'taxonomy_file']),
            Option(['-kp', 'kma_path']),
            Switch(['-x', 'extented_output']),
            Switch(['-q', 'quiet'])
        ]
        AbstractCommandline.__init__(self, cmd, **kwargs)


class PlasmidFinderCommandline(AbstractCommandline):
    def __init__(self, cmd='plasmidfinder.py', **kwargs):
        self.parameters = [
            Option(['-i', 'infile']),
            Option(['-o', 'outdir']),
            Option(['-p', 'database']),
            Option(['-tmp', 'tmp']),
            Option(['-l', 'mincov']),
            Option(['-t', 'threshold']),
            Switch(['-x', 'extented_output']),
            Switch(['-q', 'quiet'])
        ]
        AbstractCommandline.__init__(self, cmd, **kwargs)


class AMRFinderCommandline(AbstractCommandline):
    def __init__(self, cmd='amrfinder', **kwargs):
        self.parameters = [
            Option(['-n', 'nuc_fasta']),
            Option(['-o', 'output_file']),
            Option(['-d', 'database']),
            Option(['--threads', 'threads']),
        ]
        AbstractCommandline.__init__(self, cmd, **kwargs)


class ShovillCommandline(AbstractCommandline):
    def __init__(self, cmd='shovill', **kwargs):
        self.parameters = [
            Option(['--R1', 'r1']),
            Option(['--R2', 'r2']),
            Option(['--outdir', 'outdir']),
            Option(['--depth', 'depth']),
            Option(['--gsize', 'gsize']),
            Option(['--tmpdir', 'tmpdir']),
            Option(['--cpus', 'cpus']),
            Option(['--ram', 'ram']),
            Switch(['--trim', 'trim']),
            Switch(['--force', 'force']),
            Switch(['--noreadcorr', 'noreadcorr']),
            Switch(['--nostitch', 'nostitch']),
        ]
        AbstractCommandline.__init__(self, cmd, **kwargs)