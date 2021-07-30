import os
from .utils import run_cmd


def pilon_polish(fq_1, fq_2, reference, outdir, threads=4):
    os.makedirs(outdir)
    run_cmd(f"bwa index {reference}")
    run_cmd(f"samtools faidx  {reference}")
    bam = os.path.join(outdir, 'alignments.sort.bam')
    cmd = f"bwa mem -v 3 -x intractg -t {threads} {reference} {fq_1} {fq_2} | " \
          f"samclip --ref {reference}.fai | " \
          f"samtools sort --threads {threads} -m 500m --reference {reference} -T /tmp/ -o {bam}"
    run_cmd(cmd)
    run_cmd(f"samtools index  {bam}")
    cmd = f"pilon --genome {reference} --frags {bam} --minmq 60 --minqual 3 --fix bases --output pilon" \
          f" --outdir {outdir} --threads {threads} --changes --mindepth 0.25"
    run_cmd(cmd)

    pilon_changes_file = os.path.join(outdir, 'pilon.changes')
    change_count = 0
    with open(pilon_changes_file) as pilon_change:
        for line in pilon_change:
            change_count += 1
    return change_count
