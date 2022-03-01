import subprocess
from tempfile import TemporaryDirectory


def estimate_genome_size(filename, threads):
    with TemporaryDirectory() as tmpdir:
        output = subprocess.check_output(
            f"kmc -sm -t{threads} -k21 -ci10 {filename} {tmpdir}/kmc {tmpdir} | "
            f"grep 'No. of unique counted k-mers' | "
            f"awk '{{print $NF}}'",
            shell=True
        )
    return int(output.decode().strip().split()[-1])


def count_bases(filename):
    output = subprocess.check_output(
        f"seqtk fqchk {filename} | grep ALL | awk '{{print $2}}'",
        shell=True
    )
    return int(output.decode().strip())
