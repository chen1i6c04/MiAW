import os
import subprocess


def reads_alignment(paired_1, paired_2, assembly, outfile, threads):
    subprocess.run(f'bwa-mem2 index {assembly}', shell=True, check=True)
    subprocess.run(
        f'bwa-mem2 mem -v 3 -x intractg -t {threads} {assembly} {paired_1} {paired_2} | '
        f'samtools sort --threads {threads} -m 500m --reference {assembly} -T /tmp -o {outfile}',
        shell=True, check=True
    )
    subprocess.run(f'samtools index {outfile}', shell=True, check=True)


def run_polca(assembly, paired_1, paired_2, output_dir, num_threads, logfile):
    assembly = os.path.abspath(assembly)
    os.chdir(output_dir)
    cmd = f"polca.sh -a {assembly} -r '{paired_1} {paired_2}' -t {num_threads} 2>&1 >> {logfile}"
    subprocess.run(cmd, shell=True, check=True)


def assembly_pipeline(paired_1, paired_2, outdir, threads, logfile):
    polca_dirname = os.path.join(outdir, 'polca')
    spades_assembly = os.path.join(outdir, 'contigs.fasta')
    polca_assembly = os.path.join(polca_dirname, 'contigs.fasta.PolcaCorrected.fa')
    os.makedirs(polca_dirname)
    subprocess.run(
        f"spades.py -1 {paired_1} -2 {paired_2} -o {outdir} -t {threads} --tmp-dir /tmp --cov-cutoff auto "
        f"--only-assembler --disable-gzip-output --isolate",
        shell=True, check=True,
    )
    run_polca(spades_assembly, paired_1, paired_2, polca_dirname, threads, logfile)
    return polca_assembly
