import subprocess
from tempfile import TemporaryDirectory


def run_cmd(cmd):
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, check=True)
    return p


def estimate_genome_size(filepath, threads, tmpdir='/tmp'):
    with TemporaryDirectory(dir=tmpdir) as tmp:
        cmd = f"kmc -sm -t{threads} -k21 -ci10 {filepath} {tmp}/kmc {tmp} | grep 'No. of unique counted k-mers' | " \
              f"awk '{{print $NF}}'"
        p = run_cmd(cmd)
    return int(p.stdout.decode().strip())


def count_bases(filepath):
    p = run_cmd(f"seqtk fqchk {filepath} | grep ALL | awk '{{print $2}}'")
    return int(p.stdout.decode().strip())
