import os
import subprocess


def run_spades(paired_1, paired_2, outdir, threads):
    cmd = [
        'spades.py',
        '-1', paired_1,
        '-2', paired_2,
        '-o', outdir,
        '--threads', str(threads),
        '--tmp-dir', '/tmp',
        '--cov-cutoff', 'auto',
        '--only-assembler',
        '--disable-gzip-output',
        '--isolate'
    ]
    subprocess.run(cmd)


def reads_alignment(paired_1, paired_2, assembly, outfile, threads):
    subprocess.run(['bwa-mem2', 'index', assembly])
    bwa_process = subprocess.Popen(
        [
            'bwa-mem2', 'mem',
            '-v', '3',
            '-x', 'intractg',
            '-t', str(threads),
            assembly, paired_1, paired_2,
        ], stdout=subprocess.PIPE
    )
    subprocess.run(
        [
            'samtools', 'sort',
            '--threads', str(threads),
            '-m', '500m',
            '--reference', assembly,
            '-T', '/tmp',
            '-o', outfile,
        ],
        stdin=bwa_process.stdout)
    subprocess.run(['samtools', 'index', outfile])


def run_pilon(assembly, alignments, outdir, threads):
    cmd = [
        'pilon',
        '--genome', assembly,
        '--frags', alignments,
        '--minmq', '60',
        '--minqual', '3',
        '--fix', 'bases',
        '--output', 'pilon',
        '--outdir', outdir,
        '--threads', str(threads),
        '--mindepth', '0.25',
        '--changes',
    ]
    subprocess.run(cmd)


def assembly_pipeline(paired_1, paired_2, outdir, threads):
    pilon_dirname = os.path.join(outdir, 'pilon')
    asm_filename = os.path.join(outdir, 'contigs.fasta')
    bam_filename = os.path.join(pilon_dirname, 'alignments.sort.bam')
    os.makedirs(pilon_dirname, exist_ok=True)
    run_spades(paired_1, paired_2, outdir, threads)
    reads_alignment(paired_1, paired_2, asm_filename, bam_filename, threads)
    run_pilon(asm_filename, bam_filename, pilon_dirname, threads)
    return pilon_dirname
