import subprocess
from tempfile import TemporaryDirectory


def estimate_genome_size(filename, threads):
    with TemporaryDirectory() as tmpdir:
        kmc_process = subprocess.Popen(
            ['kmc', '-sm', f'-t{threads}', '-k21', '-ci10', filename, f'{tmpdir}/kmc', tmpdir],
            stdout=subprocess.PIPE,
        )
        grep_process = subprocess.Popen(
            ['grep', 'No. of unique counted k-mers'],
            stdin=kmc_process.stdout,
            stdout=subprocess.PIPE
        )
        output = subprocess.check_output(
            ['awk', '{print $NF}'],
            stdin=grep_process.stdout
        )
    return int(output.decode().strip().split()[-1])


def count_bases(filename):
    seqtk_process = subprocess.Popen(
        ['seqtk', 'fqchk', filename],
        stdout=subprocess.PIPE
    )
    grep_process = subprocess.Popen(
        ['grep', 'ALL'],
        stdin=seqtk_process.stdout,
        stdout=subprocess.PIPE
    )
    output = subprocess.check_output(
        ['awk', '{print $2}'],
        stdin=grep_process.stdout
    )
    return int(output.decode().strip())
